% Channel-specific expected-reward Cox analysis for IT, RT, and BR


clear; clc; close all;

inputFolderName = 'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';
modelingFolderName = 'D:\Nill\data\BART\0_0_new_IED\context_modeling\param_recovery_1_modeling\';
outputFolderName = 'D:\Nill\code\BART\IED\0_0_new_IED\IED7_Cox_expected_reward\';
if ~exist(outputFolderName,'dir'), mkdir(outputFolderName); end

settings.postIEDWindowMillisecondsRT = 500;
settings.postIEDWindowMillisecondsIT = 1000;
settings.postIEDWindowMillisecondsBR = 1000;
settings.maximumRTSeconds = 20;
settings.defaultSamplingFrequencyHz = 1000;
settings.combineLeftAndRight = true;
settings.minimumIEDsPerChannel = 3;
settings.minimumEventsInsideWindow = 5;
settings.minimumEventsOutsideWindow = 5;
settings.channelFDRLevel = 0.05; 
settings.numberOfPermutations = 500;
settings.randomSeed = 20260713;

rng(settings.randomSeed,'twister');

analysisTypes = ["IT","RT","BR"];
allResults = table();

for analysisType = analysisTypes
    fprintf('\n============================================================\n');
    fprintf('Starting %s channel analysis\n',analysisType);
    fprintf('============================================================\n');
    R = runChannelAnalysis(inputFolderName,modelingFolderName,analysisType,settings);
    writetable(R,fullfile(outputFolderName,analysisType + "_all_channels.csv"));
    writetable(R(R.significantChannel,:),fullfile(outputFolderName, ...
        analysisType + "_significant_channels_only.csv"));
    allResults = [allResults; R]; 
end

writetable(allResults,fullfile(outputFolderName,'IT_RT_BR_all_channels.csv'));
writetable(allResults(allResults.significantChannel,:),fullfile(outputFolderName, ...
    'IT_RT_BR_significant_channels_only.csv'));

areaSummary = makeAreaSummary(allResults);
writetable(areaSummary,fullfile(outputFolderName, ...
    'IT_RT_BR_significant_channel_area_summary.csv'));

plotAreaSummary(areaSummary,outputFolderName);
save(fullfile(outputFolderName,'IT_RT_BR_channel_results.mat'), ...
    'allResults','areaSummary','settings');

fprintf('\nFinished. Significant channels: %d of %d fitted channels.\n', ...
    sum(allResults.significantChannel),sum(allResults.status=="fitted"));

function results = runChannelAnalysis(inputFolder,modelingFolder,analysisType,S)
files = dir(fullfile(inputFolder,'*.LFPIED.mat'));
if isempty(files), error('No .LFPIED.mat files found.'); end
results = initializeResults();

for ff = 1:numel(files)
    parts = split(string(files(ff).name),'.');
    patientID = parts(1);
    fprintf('\n[%s] Participant %d/%d: %s\n',analysisType,ff,numel(files),patientID);
    D = load(fullfile(inputFolder,files(ff).name));
    if ~isfield(D,'LFPIED'), fprintf('  skipped: no LFPIED\n'); continue; end
    L = D.LFPIED;

    modelFiles = dir(fullfile(modelingFolder,patientID + "*TDdataParamRecovery.mat"));
    if isempty(modelFiles)
        fprintf('  skipped: expected-reward file not found\n'); continue
    elseif numel(modelFiles)>1
        fprintf('  skipped: multiple expected-reward files matched patient ID\n'); continue
    end
    M = load(fullfile(modelingFolder,modelFiles(1).name));
    if ~isfield(M,'TDdataParamRecovery') || ...
            ~all(isfield(M.TDdataParamRecovery,{'bestApIdx','bestAnIdx','expectedReward'}))
        fprintf('  skipped: expected-reward fields missing\n'); continue
    end
    TD=M.TDdataParamRecovery;
    expectedReward=double(squeeze(TD.expectedReward(TD.bestApIdx,TD.bestAnIdx,:)));
    expectedReward=expectedReward(:);

    common = {'selectedChans','anatomicalLocs','RTs','isControl','balloonType'};
    if analysisType=="RT", needed=[common,{'IED_occurance_RT'}];
    elseif analysisType=="IT", needed=[common,{'ITs','IED_occurance_IT'}];
    else, needed=[common,{'ITs','BankedTrials','IED_occurance_IT'}]; end
    if any(~isfield(L,needed)), fprintf('  skipped: missing required field\n'); continue; end

    selectedChans = round(double(L.selectedChans(:)));
    if isempty(selectedChans), continue; end
    labels = selectedLabels(L.anatomicalLocs,selectedChans);
    labels = cleanLabels(labels,S.combineLeftAndRight);
    excluded = excludedLabels(labels);

    RT = double(L.RTs(:)); control=double(L.isControl(:));
    balloon=mapColor(double(L.balloonType(:)));
    if analysisType=="RT"
        duration=RT; endpoint=true(size(RT)); IED=L.IED_occurance_RT;
        n=min([numel(RT),numel(control),numel(balloon),numel(expectedReward)]);
    elseif analysisType=="IT"
        IT=double(L.ITs(:)); duration=IT; endpoint=true(size(IT));
        IED=L.IED_occurance_IT;
        n=min([numel(RT),numel(IT),numel(control),numel(balloon),numel(expectedReward)]);
    else
        IT=double(L.ITs(:)); bank=double(L.BankedTrials(:));
        duration=IT; endpoint=bank==1; IED=L.IED_occurance_IT;
        n=min([numel(RT),numel(IT),numel(bank),numel(control),numel(balloon),numel(expectedReward)]);
    end
    RT=RT(1:n); duration=duration(1:n); endpoint=endpoint(1:n);
    control=control(1:n); balloon=balloon(1:n);
    expectedReward=expectedReward(1:n);
    valid=control==0 & isfinite(RT) & RT>0 & RT<=S.maximumRTSeconds & ...
        isfinite(duration) & duration>0 & ismember(balloon,[1 2 3]) & ...
        isfinite(expectedReward);
    if analysisType=="BR", valid=valid & ismember(bank(1:n),[0 1]); end
    Fs=getFs(L,S.defaultSamplingFrequencyHz);
    window=getWindow(analysisType,S)/1000;

    for cc=1:numel(selectedChans)
        fprintf('  channel %d/%d\n',cc,numel(selectedChans));
        if excluded(cc), continue; end
        channelIED=IED(isfinite(IED(:,2)) & round(IED(:,2))==cc,:);
        validIED=validIEDRows(channelIED,valid);
        nIED=sum(validIED);
        status="fitted"; beta=NaN; se=NaN; z=NaN; p=NaN;
        HR=NaN; lo=NaN; hi=NaN; nIn=0; nOut=0; nRows=0;
        if nIED<S.minimumIEDsPerChannel
            status="skipped: too few IEDs";
        else
            trials=find(valid);
            observedRows=cell(numel(trials),1);
            for tt=1:numel(trials)
                tr=trials(tt);
                observedRows{tt}=makeRows(tr,duration(tr),endpoint(tr), ...
                    balloon(tr),expectedReward(tr),channelIED,Fs,window);
            end
            CP=vertcat(observedRows{:});
            nRows=height(CP);
            nIn=sum(CP.eventAtStop & CP.postIED==1);
            nOut=sum(CP.eventAtStop & CP.postIED==0);
            if nIn<S.minimumEventsInsideWindow
                status="skipped: too few events inside window";
            elseif nOut<S.minimumEventsOutsideWindow
                status="skipped: too few events outside window";
            else
                try
                    [beta,se,z]=fitInteractionCox( ...
                        CP,S.minimumEventsInsideWindow);
                    fprintf('    running %d permutations',S.numberOfPermutations);
                    permBeta=NaN(S.numberOfPermutations,1);
                    for permutationIndex=1:S.numberOfPermutations
                        permutedRows=cell(numel(trials),1);
                        for trialIndex=1:numel(trials)
                            tr=trials(trialIndex);
                            shiftedIED=circularShiftIEDsWithinTrial( ...
                                channelIED,tr,duration(tr),Fs);
                            permutedRows{trialIndex}=makeRows(tr,duration(tr), ...
                                endpoint(tr),balloon(tr),expectedReward(tr), ...
                                shiftedIED,Fs,window);
                        end
                        permutedCP=vertcat(permutedRows{:});
                        try
                            permBeta(permutationIndex)=fitInteractionCox( ...
                                permutedCP,S.minimumEventsInsideWindow);
                        catch
                            permBeta(permutationIndex)=NaN;
                        end
                        if mod(permutationIndex,50)==0
                            fprintf(' ... %d',permutationIndex);
                        end
                    end
                    fprintf('\n');
                    validPermBeta=permBeta(isfinite(permBeta));
                    if numel(validPermBeta)<max(100,ceil(0.80*S.numberOfPermutations))
                        error('Too few successful permutation fits (%d of %d).', ...
                            numel(validPermBeta),S.numberOfPermutations);
                    end
                    p=(1+sum(abs(validPermBeta)>=abs(beta)))/(numel(validPermBeta)+1);
                    HR=exp(beta); lo=exp(beta-1.96*se); hi=exp(beta+1.96*se);
                    if ~all(isfinite([beta,se,p,HR,lo,hi])), status="failed: nonfinite estimate"; end
                catch ME
                    status="failed: "+string(ME.message);
                end
            end
        end
        row=table(analysisType,patientID,selectedChans(cc),cc,labels(cc),status, ...
            sum(valid),nIED,nRows,nIn,nOut,beta,se,z,p,NaN,HR,lo,hi,false, ...
            'VariableNames',results.Properties.VariableNames);
        results=[results;row]; 
    end
end

% Correct within each participant and outcome, not across the entire cohort.
patients=unique(results.patientID);
for pp=1:numel(patients)
    idx=results.patientID==patients(pp) & results.status=="fitted" & ...
        isfinite(results.pValue_IEDxExpectedReward);
    if any(idx)
        results.pValueFDR_IEDxExpectedReward(idx)= ...
            BH(results.pValue_IEDxExpectedReward(idx));
        results.significantChannel(idx)= ...
            results.pValueFDR_IEDxExpectedReward(idx)<S.channelFDRLevel;
    end
end
results=sortrows(results,{'significantChannel','pValueFDR_IEDxExpectedReward'}, ...
    {'descend','ascend'});
end

function T=initializeResults()
T=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
    strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),false(0,1),'VariableNames', ...
    {'analysisType','patientID','originalChannelNumber','localSelectedChannelIndex', ...
    'anatomicalArea','status','nValidTrials','nIEDs','nCountingRows', ...
    'nEventsInsideWindow','nEventsOutsideWindow','beta_IEDxExpectedReward', ...
    'standardError_IEDxExpectedReward','zStatistic_IEDxExpectedReward', ...
    'pValue_IEDxExpectedReward','pValueFDR_IEDxExpectedReward', ...
    'interactionHazardRatio','interactionHazardRatioCILow', ...
    'interactionHazardRatioCIHigh','significantChannel'});
end

function S=makeAreaSummary(R)
R=R(R.significantChannel,:);
if isempty(R)
    S=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
        'VariableNames',{'analysisType','anatomicalArea','nSignificantChannels', ...
        'nParticipants','medianInteractionHazardRatio'}); return
end
[G,type,area]=findgroups(R.analysisType,R.anatomicalArea);
S=table(type,area,splitapply(@numel,R.pValue_IEDxExpectedReward,G), ...
    splitapply(@(x)numel(unique(x)),R.patientID,G), ...
    splitapply(@median,R.interactionHazardRatio,G),'VariableNames', ...
    {'analysisType','anatomicalArea','nSignificantChannels','nParticipants', ...
    'medianInteractionHazardRatio'});
S=sortrows(S,{'analysisType','nSignificantChannels'},{'ascend','descend'});
end

function plotAreaSummary(S,out)
if isempty(S), return; end
types=["IT","RT","BR"];
f=figure('Color','w','Position',[100 100 1500 700]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i=1:3
    nexttile; A=S(S.analysisType==types(i),:);
    if isempty(A), title(types(i)+": no significant channels"); axis off; continue; end
    A=sortrows(A,'nSignificantChannels','ascend');
    barh(categorical(A.anatomicalArea,A.anatomicalArea),A.nSignificantChannels);
    xlabel('Number of significant channels'); title(types(i)); 
    set(gca,'FontSize',11,'FontWeight','bold');
end
exportgraphics(f,fullfile(out,'IT_RT_BR_significant_channels_by_area.pdf'),'ContentType','vector');
close(f);
end

function labels=selectedLabels(raw,selected)
if iscell(raw), labels=string(raw); elseif ischar(raw), labels=string(cellstr(raw)); else, labels=string(raw); end
labels=labels(:);
if numel(labels)>=max(selected), labels=labels(selected);
elseif numel(labels)~=numel(selected), error('Cannot map anatomical labels to selected channels.'); end
end

function x=cleanLabels(x,combine)
x=strip(string(x)); x(ismissing(x)|strlength(x)==0)="Unknown";
if combine
    x=regexprep(x,'^\s*\(?\s*(Left|Right|LH|RH|L|R)\s*\)?[_\-\s]+','','ignorecase');
    x=regexprep(x,'[_\-\s]+\(?\s*(Left|Right|LH|RH|L|R)\s*\)?\s*$','','ignorecase');
end
x=strip(regexprep(x,'[_\-]+',' '));
end

function tf=excludedLabels(x)
n=lower(regexprep(strip(string(x)),'[_\-]+',' ')); c=regexprep(n,'[^a-z0-9]','');
tf=contains(n,'nac')|contains(n,'accumbens')|contains(n,'white matter')| ...
    startsWith(c,'wm')|contains(c,'lateralventricle')|contains(c,'ventraldc');
end

function c=mapColor(x)
c=NaN(size(x)); ok=isfinite(x)&ismember(round(x),[1 2 3 11 12 13]);
c(ok)=mod(round(x(ok))-1,10)+1;
end

function fs=getFs(L,d)
if isfield(L,'Fs')&&isscalar(L.Fs)&&isfinite(L.Fs)&&L.Fs>0, fs=double(L.Fs); else, fs=d; end
end

function w=getWindow(type,S)
if type=="RT", w=S.postIEDWindowMillisecondsRT;
elseif type=="IT", w=S.postIEDWindowMillisecondsIT; else, w=S.postIEDWindowMillisecondsBR; end
end

function v=validIEDRows(IED,valid)
if isempty(IED)||size(IED,2)<3, v=false(0,1); return; end
tr=round(IED(:,1)); v=isfinite(tr)&tr>=1&tr<=numel(valid)&isfinite(IED(:,3))&IED(:,3)>=1;
i=find(v); v(i)=valid(tr(i));
end

function T=emptyCP()
T=table(zeros(0,1),zeros(0,1),false(0,1),false(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'tStart','tStop','censored','eventAtStop','postIED','balloonColorCode','expectedReward'});
end

function R=makeRows(tr,duration,event,color,reward,IED,fs,window)
if isempty(IED), times=[]; else
    ok=isfinite(IED(:,1))&round(IED(:,1))==tr&isfinite(IED(:,3))&IED(:,3)>=1;
    times=double(IED(ok,3))/fs; times=sort(times(times>0&times<duration));
end
breaks=unique([0;times(:);min(times(:)+window,duration);duration],'sorted');
start=breaks(1:end-1); stop=breaks(2:end); ok=stop>start; start=start(ok); stop=stop(ok);
mid=(start+stop)/2; post=false(numel(start),1);
for k=1:numel(start), post(k)=any(mid(k)>=times & mid(k)<=times+window); end
cens=true(numel(start),1); ev=false(numel(start),1);
if event, cens(end)=false; ev(end)=true; end
R=table(start,stop,cens,ev,double(post),repmat(color,numel(start),1), ...
    repmat(reward,numel(start),1),'VariableNames', ...
    {'tStart','tStop','censored','eventAtStop','postIED','balloonColorCode','expectedReward'});
end

function [beta,se,z]=fitInteractionCox(CP,minimumEventsPerGroup)
postIED=double(CP.postIED);
rawV=double(CP.expectedReward);
eventRows=CP.eventAtStop==1;

nEventsInside=sum(eventRows & postIED==1);
nEventsOutside=sum(eventRows & postIED==0);
if nEventsInside<minimumEventsPerGroup || ...
        nEventsOutside<minimumEventsPerGroup
    error('Insufficient outcome events inside or outside the post-IED window.');
end

insideEventV=rawV(eventRows & postIED==1);
outsideEventV=rawV(eventRows & postIED==0);
if numel(unique(insideEventV))<3 || numel(unique(outsideEventV))<3 || ...
        std(insideEventV)<=sqrt(eps)*max(1,abs(mean(insideEventV))) || ...
        std(outsideEventV)<=sqrt(eps)*max(1,abs(mean(outsideEventV)))
    error(['Expected reward has insufficient event-level variation ' ...
        'inside or outside the post-IED window.']);
end

V=rawV-mean(rawV,'omitnan');
if ~all(isfinite(V))
    error('Expected reward contains a nonfinite value.');
end

observedColors=unique(CP.balloonColorCode);
colorEventCounts=zeros(numel(observedColors),1);
for colorIndex=1:numel(observedColors)
    colorEventCounts(colorIndex)=sum(eventRows & ...
        CP.balloonColorCode==observedColors(colorIndex));
end
if any(colorEventCounts<minimumEventsPerGroup)
    error('At least one balloon-color level has insufficient outcome events.');
end

[~,referenceIndex]=max(colorEventCounts);
referenceColor=observedColors(referenceIndex);
X=[postIED,V];
for color=reshape(observedColors(observedColors~=referenceColor),1,[])
    X(:,end+1)=double(CP.balloonColorCode==color);
end

X(:,end+1)=postIED.*V;
interactionIndex=size(X,2);
if ~any(postIED==0) || ~any(postIED==1)
    error('Post-IED predictor does not vary.');
end

varying=any(X~=X(1,:),1);
if ~varying(interactionIndex)
    error('The post-IED x expected-reward interaction has no variation.');
end
X=X(:,varying);
interactionIndex=sum(varying(1:interactionIndex));
if rank(X)<size(X,2)
    error('Cox design matrix is rank deficient.');
end

informationScale=X'*X;
if rcond(informationScale)<1e-12
    error('Cox design matrix is numerically ill-conditioned.');
end

T=[CP.tStart CP.tStop];
cens=logical(CP.censored);
opts=statset('coxphfit');
opts.Display='off';
opts.MaxIter=500;

singularState=warning('error','MATLAB:singularMatrix');
nearSingularState=warning('error','MATLAB:nearlySingularMatrix');
warningCleanup=onCleanup(@()restoreWarningStates( ...
    singularState,nearSingularState)); 

[b,~,~,stats]=coxphfit(X,T,'Censoring',cens, ...
    'Ties','efron','Baseline',0,'Options',opts);
beta=b(interactionIndex);
se=stats.se(interactionIndex);
z=stats.z(interactionIndex);
if ~all(isfinite([beta,se,z])) || se<=0
    error('Cox fit returned an unstable interaction estimate.');
end
end

function restoreWarningStates(singularState,nearSingularState)
warning(singularState);
warning(nearSingularState);
end

function shiftedIED=circularShiftIEDsWithinTrial(channelIED,tr,duration,fs)
shiftedIED=channelIED;
rows=isfinite(channelIED(:,1)) & round(channelIED(:,1))==tr & ...
    isfinite(channelIED(:,3)) & channelIED(:,3)>=1;
if ~any(rows) || ~isfinite(duration) || duration<=0
    return;
end
times=double(channelIED(rows,3))/fs;
offset=rand*duration;
shiftedTimes=mod(times+offset,duration);
shiftedTimes(shiftedTimes<=0)=eps(duration);
shiftedIED(rows,3)=shiftedTimes*fs;
end

function q=BH(p)
p=p(:); m=numel(p); [s,o]=sort(p); adj=s.*m./(1:m)';
adj=flipud(cummin(flipud(adj))); adj=min(adj,1); q=NaN(m,1); q(o)=adj;
end