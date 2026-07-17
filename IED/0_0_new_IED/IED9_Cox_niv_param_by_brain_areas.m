% Channel-specific Niv-asymmetry Cox analysis for IT, RT, and BR


clear; clc; close all;

inputFolderName = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';
nivParameterFile = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\context_modeling\param_recovery_4_param_recovery\alpha_comparison.csv';
outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED9_Cox_niv_param_by_brain_areas\';
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
settings.minimumParticipantsPerArea = 5;
settings.areaFDRLevel = 0.10;
settings.numberOfPermutations = 1000;
settings.randomSeed = 20260713;

rng(settings.randomSeed,'twister');

[nivTable,nivSummary] = readNivParameterTable(nivParameterFile);
writetable(nivTable,fullfile(outputFolderName,'niv_parameters_used.csv'));

analysisTypes = ["IT","RT","BR"];
allChannelResults = table();
allAreaResults = table();

for analysisType = analysisTypes
    fprintf('\n============================================================\n');
    fprintf('Starting %s channel analysis\n',analysisType);
    fprintf('============================================================\n');

    channelResults = runChannelAnalysis( ...
        inputFolderName,analysisType,settings,nivTable);
    areaResults = makeAreaNivResults( ...
        channelResults,analysisType,settings);

    writetable(channelResults,fullfile(outputFolderName, ...
        analysisType + "_all_channels.csv"));
    writetable(areaResults,fullfile(outputFolderName, ...
        analysisType + "_all_anatomical_areas.csv"));
    writetable(areaResults(areaResults.significantArea,:), ...
        fullfile(outputFolderName, ...
        analysisType + "_significant_anatomical_areas_only.csv"));

    allChannelResults = [allChannelResults; channelResults]; 
    allAreaResults = [allAreaResults; areaResults]; 
end

writetable(allChannelResults,fullfile(outputFolderName, ...
    'IT_RT_BR_all_channels.csv'));
writetable(allAreaResults,fullfile(outputFolderName, ...
    'IT_RT_BR_all_anatomical_areas.csv'));
writetable(allAreaResults(allAreaResults.significantArea,:), ...
    fullfile(outputFolderName, ...
    'IT_RT_BR_significant_anatomical_areas_only.csv'));

plotAreaSummary(allAreaResults,outputFolderName);
save(fullfile(outputFolderName,'IT_RT_BR_Niv_channel_results.mat'), ...
    'allChannelResults','allAreaResults','nivTable','nivSummary','settings');

fprintf('\nFinished. Significant anatomical areas: %d of %d fitted areas.\n', ...
    sum(allAreaResults.significantArea),sum(allAreaResults.status=="fitted"));

function results = runChannelAnalysis(inputFolder,analysisType,S,nivTable)
files = dir(fullfile(inputFolder,'*.LFPIED.mat'));
if isempty(files), error('No .LFPIED.mat files found.'); end
results = initializeChannelResults();

for ff = 1:numel(files)
    parts = split(string(files(ff).name),'.');
    patientID = parts(1);
    fprintf('\n[%s] Participant %d/%d: %s\n', ...
        analysisType,ff,numel(files),patientID);

    [nivFound,alphaPlus,alphaMinus,nivRaw,nivZ] = ...
        getNivParametersForPatient(nivTable,patientID);
    if ~nivFound
        fprintf('  skipped: Niv parameters not found\n');
        continue;
    end

    D = load(fullfile(inputFolder,files(ff).name));
    if ~isfield(D,'LFPIED')
        fprintf('  skipped: no LFPIED\n');
        continue;
    end
    L = D.LFPIED;

    common = {'selectedChans','anatomicalLocs','RTs', ...
        'isControl','balloonType'};
    if analysisType=="RT"
        needed = [common,{'IED_occurance_RT'}];
    elseif analysisType=="IT"
        needed = [common,{'ITs','IED_occurance_IT'}];
    else
        needed = [common,{'ITs','BankedTrials','IED_occurance_IT'}];
    end
    if any(~isfield(L,needed))
        fprintf('  skipped: missing required field\n');
        continue;
    end

    selectedChans = round(double(L.selectedChans(:)));
    if isempty(selectedChans), continue; end
    labels = selectedLabels(L.anatomicalLocs,selectedChans);
    labels = cleanLabels(labels,S.combineLeftAndRight);
    excluded = excludedLabels(labels);

    RT = double(L.RTs(:));
    control = double(L.isControl(:));
    balloon = mapColor(double(L.balloonType(:)));

    if analysisType=="RT"
        duration = RT;
        endpoint = true(size(RT));
        IED = L.IED_occurance_RT;
        n = min([numel(RT),numel(control),numel(balloon)]);
    elseif analysisType=="IT"
        IT = double(L.ITs(:));
        duration = IT;
        endpoint = true(size(IT));
        IED = L.IED_occurance_IT;
        n = min([numel(RT),numel(IT),numel(control),numel(balloon)]);
    else
        IT = double(L.ITs(:));
        bank = double(L.BankedTrials(:));
        duration = IT;
        endpoint = bank==1;
        IED = L.IED_occurance_IT;
        n = min([numel(RT),numel(IT),numel(bank), ...
            numel(control),numel(balloon)]);
    end

    RT = RT(1:n);
    duration = duration(1:n);
    endpoint = endpoint(1:n);
    control = control(1:n);
    balloon = balloon(1:n);
    valid = control==0 & isfinite(RT) & RT>0 & ...
        RT<=S.maximumRTSeconds & isfinite(duration) & duration>0 & ...
        ismember(balloon,[1 2 3]);
    if analysisType=="BR"
        bank = bank(1:n);
        valid = valid & ismember(bank,[0 1]);
    end

    Fs = getFs(L,S.defaultSamplingFrequencyHz);
    window = getWindow(analysisType,S)/1000;

    for cc = 1:numel(selectedChans)
        fprintf('  channel %d/%d\n',cc,numel(selectedChans));
        if excluded(cc), continue; end

        channelIED = IED(isfinite(IED(:,2)) & ...
            round(IED(:,2))==cc,:);
        validIED = validIEDRows(channelIED,valid);
        nIED = sum(validIED);

        status = "fitted";
        beta = NaN; se = NaN; z = NaN;
        HR = NaN; lo = NaN; hi = NaN;
        nIn = 0; nOut = 0; nRows = 0;

        if nIED<S.minimumIEDsPerChannel
            status = "skipped: too few IEDs";
        else
            trials = find(valid);
            observedRows = cell(numel(trials),1);
            for tt = 1:numel(trials)
                tr = trials(tt);
                observedRows{tt} = makeRows( ...
                    tr,duration(tr),endpoint(tr),balloon(tr), ...
                    channelIED,Fs,window);
            end
            CP = vertcat(observedRows{:});
            nRows = height(CP);
            nIn = sum(CP.eventAtStop & CP.postIED==1);
            nOut = sum(CP.eventAtStop & CP.postIED==0);

            if nIn<S.minimumEventsInsideWindow
                status = "skipped: too few events inside window";
            elseif nOut<S.minimumEventsOutsideWindow
                status = "skipped: too few events outside window";
            else
                try
                    [beta,se,z] = fitChannelCox( ...
                        CP,S.minimumEventsInsideWindow);
                    HR = exp(beta);
                    lo = exp(beta-1.96*se);
                    hi = exp(beta+1.96*se);
                    if ~all(isfinite([beta,se,HR,lo,hi]))
                        status = "failed: nonfinite estimate";
                    end
                catch ME
                    status = "failed: "+string(ME.message);
                end
            end
        end

        row = table(analysisType,patientID,selectedChans(cc),cc, ...
            labels(cc),status,sum(valid),nIED,nRows,nIn,nOut, ...
            alphaPlus,alphaMinus,nivRaw,nivZ,beta,se,z,HR,lo,hi, ...
            'VariableNames',results.Properties.VariableNames);
        results = [results; row]; 
    end
end

results = sortrows(results, ...
    {'status','patientID','anatomicalArea','localSelectedChannelIndex'}, ...
    {'descend','ascend','ascend','ascend'});
end

function T = initializeChannelResults()
T = table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
    strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1), ...
    'VariableNames',{'analysisType','patientID','originalChannelNumber', ...
    'localSelectedChannelIndex','anatomicalArea','status', ...
    'nValidTrials','nIEDs','nCountingRows','nEventsInsideWindow', ...
    'nEventsOutsideWindow','alphaPlus','alphaMinus','nivAsymmetryRaw', ...
    'nivAsymmetryZ','beta_postIED','standardError_postIED', ...
    'zStatistic_postIED','hazardRatio_postIED', ...
    'hazardRatioCILow','hazardRatioCIHigh'});
end

function A = makeAreaNivResults(R,analysisType,S)
fitted = R.status=="fitted" & isfinite(R.beta_postIED) & ...
    isfinite(R.nivAsymmetryZ);
P = aggregateParticipantAreaEffects(R(fitted,:));
areas = unique(P.anatomicalArea);
n = numel(areas);

A = table(repmat(analysisType,n,1),areas,zeros(n,1),zeros(n,1), ...
    NaN(n,1),NaN(n,1),NaN(n,1),NaN(n,1),NaN(n,1),NaN(n,1), ...
    NaN(n,1),false(n,1),strings(n,1), ...
    'VariableNames',{'analysisType','anatomicalArea','nParticipants', ...
    'nChannels','beta_IEDxNiv','standardError_IEDxNiv', ...
    'betaCILow','betaCIHigh','correlation','pValue_IEDxNiv', ...
    'pValueFDR_IEDxNiv','significantArea','status'});

rng(S.randomSeed+getAnalysisSeedOffset(analysisType),'twister');
for aa = 1:n
    D = P(P.anatomicalArea==areas(aa),:);
    x = D.nivAsymmetryZ;
    y = D.medianChannelBeta;
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);

    A.nParticipants(aa) = numel(x);
    A.nChannels(aa) = sum(D.nChannels(ok));

    if numel(x)<S.minimumParticipantsPerArea
        A.status(aa) = "skipped: too few participants";
        continue;
    elseif std(x)<=sqrt(eps)*max(1,abs(mean(x))) || ...
            std(y)<=sqrt(eps)*max(1,abs(mean(y)))
        A.status(aa) = "skipped: insufficient variation";
        continue;
    end

    [slope,se,r] = simpleSlope(x,y);
    fprintf('  [%s] area %d/%d: %s; running %d permutations', ...
        analysisType,aa,n,areas(aa),S.numberOfPermutations);
    nullSlope = NaN(S.numberOfPermutations,1);
    for permutationIndex = 1:S.numberOfPermutations
        nullSlope(permutationIndex) = ...
            simpleSlope(x(randperm(numel(x))),y);
        if mod(permutationIndex,100)==0
            fprintf(' ... %d',permutationIndex);
        end
    end
    fprintf('\n');

    p = (1+sum(abs(nullSlope)>=abs(slope))) / ...
        (S.numberOfPermutations+1);
    A.beta_IEDxNiv(aa) = slope;
    A.standardError_IEDxNiv(aa) = se;
    A.betaCILow(aa) = slope-1.96*se;
    A.betaCIHigh(aa) = slope+1.96*se;
    A.correlation(aa) = r;
    A.pValue_IEDxNiv(aa) = p;
    A.status(aa) = "fitted";
end

fitted = A.status=="fitted" & isfinite(A.pValue_IEDxNiv);
if any(fitted)
    A.pValueFDR_IEDxNiv(fitted) = BH(A.pValue_IEDxNiv(fitted));
    A.significantArea(fitted) = ...
        A.pValueFDR_IEDxNiv(fitted)<S.areaFDRLevel;
end
A = sortrows(A,{'significantArea','pValueFDR_IEDxNiv'}, ...
    {'descend','ascend'});
end

function P = aggregateParticipantAreaEffects(R)
if isempty(R)
    P = table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
        zeros(0,1),zeros(0,1),zeros(0,1), ...
        'VariableNames',{'patientID','anatomicalArea', ...
        'medianChannelBeta','medianChannelSE','nChannels', ...
        'nivAsymmetryRaw','nivAsymmetryZ'});
    return;
end

[G,patient,area] = findgroups(R.patientID,R.anatomicalArea);
P = table(patient,area, ...
    splitapply(@median,R.beta_postIED,G), ...
    splitapply(@median,R.standardError_postIED,G), ...
    splitapply(@numel,R.beta_postIED,G), ...
    splitapply(@(x)x(1),R.nivAsymmetryRaw,G), ...
    splitapply(@(x)x(1),R.nivAsymmetryZ,G), ...
    'VariableNames',{'patientID','anatomicalArea', ...
    'medianChannelBeta','medianChannelSE','nChannels', ...
    'nivAsymmetryRaw','nivAsymmetryZ'});
end

function [slope,se,r] = simpleSlope(x,y)
x = x(:); y = y(:);
xc = x-mean(x); yc = y-mean(y);
denominator = sum(xc.^2);
slope = sum(xc.*yc)/denominator;
intercept = mean(y)-slope*mean(x);
residual = y-(intercept+slope*x);
se = sqrt(sum(residual.^2)/max(numel(x)-2,1)/denominator);
r = sum(xc.*yc)/sqrt(sum(xc.^2)*sum(yc.^2));
end

function plotAreaSummary(S,out)
f = figure('Color','w','Position',[100 100 1600 900]);
layout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
types = ["IT","RT","BR"];
colors = {[0.847 0.333 0.153],[0.204 0.459 0.702], ...
    [0.250 0.600 0.250]};

for ii = 1:3
    ax = nexttile;
    A = S(S.analysisType==types(ii) & S.status=="fitted",:);
    if isempty(A)
        title(ax,types(ii)+": no fitted areas");
        axis(ax,'off');
        continue;
    end
    A = sortrows(A,'beta_IEDxNiv','ascend');
    y = (1:height(A))';
    hold(ax,'on');
    xline(ax,0,'--','Color',[0.35 0.35 0.35],'LineWidth',1.2);
    for rr = 1:height(A)
        plot(ax,[A.betaCILow(rr),A.betaCIHigh(rr)],[y(rr),y(rr)], ...
            '-','Color',colors{ii},'LineWidth',1.5);
    end
    scatter(ax,A.beta_IEDxNiv,y,55,colors{ii},'filled');
    significantRows = find(A.significantArea);
    scatter(ax,A.beta_IEDxNiv(significantRows),y(significantRows), ...
        95,colors{ii},'filled','MarkerEdgeColor','k','LineWidth',1.2);
    for rr = reshape(significantRows,1,[])
        text(ax,A.betaCIHigh(rr),y(rr),"  *", ...
            'FontSize',13,'FontWeight','bold','VerticalAlignment','middle');
    end
    ax.YTick = y;
    ax.YTickLabel = A.anatomicalArea;
    ax.YDir = 'reverse';
    ax.FontSize = 11;
    ax.FontWeight = 'bold';
    box(ax,'off');
    xlabel(ax,'Change in post-IED beta per 1-SD Niv');
    title(ax,types(ii),'FontSize',14,'FontWeight','bold');
end
title(layout,'Anatomical localization of Niv moderation', ...
    'FontSize',16,'FontWeight','bold');
exportgraphics(f,fullfile(out, ...
    'IT_RT_BR_Niv_moderation_by_anatomical_area.pdf'), ...
    'ContentType','vector');
close(f);
end

function labels = selectedLabels(raw,selected)
if iscell(raw)
    labels = string(raw);
elseif ischar(raw)
    labels = string(cellstr(raw));
else
    labels = string(raw);
end
labels = labels(:);
if numel(labels)>=max(selected)
    labels = labels(selected);
elseif numel(labels)~=numel(selected)
    error('Cannot map anatomical labels to selected channels.');
end
end

function x = cleanLabels(x,combine)
x = strip(string(x));
x(ismissing(x)|strlength(x)==0) = "Unknown";
if combine
    x = regexprep(x, ...
        '^\s*\(?\s*(Left|Right|LH|RH|L|R)\s*\)?[_\-\s]+','', ...
        'ignorecase');
    x = regexprep(x, ...
        '[_\-\s]+\(?\s*(Left|Right|LH|RH|L|R)\s*\)?\s*$','', ...
        'ignorecase');
end
x = strip(regexprep(x,'[_\-]+',' '));
end

function tf = excludedLabels(x)
n = lower(regexprep(strip(string(x)),'[_\-]+',' '));
c = regexprep(n,'[^a-z0-9]','');
tf = contains(n,'nac') | contains(n,'accumbens') | ...
    contains(n,'white matter') | startsWith(c,'wm') | ...
    contains(c,'lateralventricle') | contains(c,'ventraldc') | ...
    n=="unknown";
end

function c = mapColor(x)
c = NaN(size(x));
ok = isfinite(x) & ismember(round(x),[1 2 3 11 12 13]);
c(ok) = mod(round(x(ok))-1,10)+1;
end

function fs = getFs(L,d)
if isfield(L,'Fs') && isscalar(L.Fs) && isfinite(L.Fs) && L.Fs>0
    fs = double(L.Fs);
else
    fs = d;
end
end

function w = getWindow(type,S)
if type=="RT"
    w = S.postIEDWindowMillisecondsRT;
elseif type=="IT"
    w = S.postIEDWindowMillisecondsIT;
else
    w = S.postIEDWindowMillisecondsBR;
end
end

function v = validIEDRows(IED,valid)
if isempty(IED) || size(IED,2)<3
    v = false(0,1);
    return;
end
tr = round(IED(:,1));
v = isfinite(tr) & tr>=1 & tr<=numel(valid) & ...
    isfinite(IED(:,3)) & IED(:,3)>=1;
idx = find(v);
v(idx) = valid(tr(idx));
end

function R = makeRows(tr,duration,event,color,IED,fs,window)
if isempty(IED)
    times = [];
else
    ok = isfinite(IED(:,1)) & round(IED(:,1))==tr & ...
        isfinite(IED(:,3)) & IED(:,3)>=1;
    times = double(IED(ok,3))/fs;
    times = sort(times(times>0 & times<duration));
end

breaks = unique([0;times(:);min(times(:)+window,duration);duration], ...
    'sorted');
start = breaks(1:end-1);
stop = breaks(2:end);
ok = stop>start;
start = start(ok); stop = stop(ok);
mid = (start+stop)/2;
post = false(numel(start),1);
for kk = 1:numel(start)
    post(kk) = any(mid(kk)>=times & mid(kk)<=times+window);
end
cens = true(numel(start),1);
ev = false(numel(start),1);
if event
    cens(end) = false;
    ev(end) = true;
end
R = table(start,stop,cens,ev,double(post), ...
    repmat(color,numel(start),1), ...
    'VariableNames',{'tStart','tStop','censored','eventAtStop', ...
    'postIED','balloonColorCode'});
end

function [beta,se,z] = fitChannelCox(CP,minimumEventsPerGroup)
postIED = double(CP.postIED);
eventRows = CP.eventAtStop==1;
if sum(eventRows & postIED==1)<minimumEventsPerGroup || ...
        sum(eventRows & postIED==0)<minimumEventsPerGroup
    error('Insufficient outcome events inside or outside the post-IED window.');
end

observedColors = unique(CP.balloonColorCode);
colorEventCounts = zeros(numel(observedColors),1);
for colorIndex = 1:numel(observedColors)
    colorEventCounts(colorIndex) = sum(eventRows & ...
        CP.balloonColorCode==observedColors(colorIndex));
end
if any(colorEventCounts<minimumEventsPerGroup)
    error('At least one balloon-color level has insufficient outcome events.');
end

[~,referenceIndex] = max(colorEventCounts);
referenceColor = observedColors(referenceIndex);
X = postIED;
for color = reshape(observedColors(observedColors~=referenceColor),1,[])
    X(:,end+1) = double(CP.balloonColorCode==color); 
end

if ~any(postIED==0) || ~any(postIED==1)
    error('Post-IED predictor does not vary.');
end
if rank(X)<size(X,2)
    error('Cox design matrix is rank deficient.');
end
if rcond(X'*X)<1e-12
    error('Cox design matrix is numerically ill-conditioned.');
end

opts = statset('coxphfit');
opts.Display = 'off';
opts.MaxIter = 500;
singularState = warning('error','MATLAB:singularMatrix');
nearSingularState = warning('error','MATLAB:nearlySingularMatrix');
warningCleanup = onCleanup(@()restoreWarningStates( ...
    singularState,nearSingularState)); 

[b,~,~,stats] = coxphfit(X,[CP.tStart CP.tStop], ...
    'Censoring',logical(CP.censored),'Ties','efron', ...
    'Baseline',0,'Options',opts);
beta = b(1);
se = stats.se(1);
z = stats.z(1);
if ~all(isfinite([beta,se,z])) || se<=0
    error('Cox fit returned an unstable post-IED estimate.');
end
end

function restoreWarningStates(singularState,nearSingularState)
warning(singularState);
warning(nearSingularState);
end

function [nivTable,nivSummary] = readNivParameterTable(nivParameterFile)
if ~exist(nivParameterFile,'file')
    error('Niv parameter CSV was not found: %s',nivParameterFile);
end

fileID = fopen(nivParameterFile,'r');
if fileID==-1
    error('Could not open Niv parameter file: %s',nivParameterFile);
end
fileCleanup = onCleanup(@()fclose(fileID));
headerLine = fgetl(fileID);
if ~ischar(headerLine) || isempty(strtrim(headerLine))
    error('The Niv parameter file has no readable header row.');
end

if contains(headerLine,sprintf('\t'))
    delimiter = sprintf('\t');
elseif contains(headerLine,',')
    delimiter = ',';
elseif contains(headerLine,';')
    delimiter = ';';
else
    error('Could not identify the delimiter in the Niv parameter CSV.');
end

headerCells = strsplit(headerLine,delimiter);
headerNames = lower(strip(string(headerCells(:))));
headerNames = erase(headerNames,'"');
headerNames = regexprep(headerNames,'[^a-z0-9]','');
idColumn = find(headerNames=="ptid",1);
alphaPlusColumn = find(headerNames=="fitalphaplus",1);
alphaMinusColumn = find(headerNames=="fitalphaminus",1);
if isempty(idColumn) || isempty(alphaPlusColumn) || isempty(alphaMinusColumn)
    error(['Niv CSV must contain ptID, fit_alpha_plus, ' ...
        'and fit_alpha_minus.']);
end

formatSpecification = repmat('%q',1,numel(headerCells));
rawColumns = textscan(fileID,formatSpecification, ...
    'Delimiter',delimiter,'ReturnOnError',false);
clear fileCleanup;

patientID = strip(erase(erase(string(rawColumns{idColumn}), ...
    '"'),char(39)));
alphaPlus = str2double(string(rawColumns{alphaPlusColumn}));
alphaMinus = str2double(string(rawColumns{alphaMinusColumn}));
denominator = alphaMinus+alphaPlus;
nivRaw = (alphaMinus-alphaPlus)./denominator;

valid = ~ismissing(patientID) & strlength(patientID)>0 & ...
    isfinite(alphaPlus) & isfinite(alphaMinus) & denominator>0 & ...
    isfinite(nivRaw);
patientID = patientID(valid);
alphaPlus = alphaPlus(valid);
alphaMinus = alphaMinus(valid);
nivRaw = nivRaw(valid);
normalizedPatientID = normalizeParticipantID(patientID);

[~,uniqueRows] = unique(normalizedPatientID,'stable');
if numel(uniqueRows)<numel(normalizedPatientID)
    warning(['Duplicate participant IDs were found in the Niv CSV. ' ...
        'Only the first row for each participant will be used.']);
end
patientID = patientID(uniqueRows);
normalizedPatientID = normalizedPatientID(uniqueRows);
alphaPlus = alphaPlus(uniqueRows);
alphaMinus = alphaMinus(uniqueRows);
nivRaw = nivRaw(uniqueRows);

nivMean = mean(nivRaw,'omitnan');
nivSD = std(nivRaw,0,'omitnan');
if ~isfinite(nivSD) || nivSD<=eps
    error('Niv asymmetry has zero or invalid variance.');
end
nivZ = (nivRaw-nivMean)./nivSD;

nivTable = table(patientID,normalizedPatientID,alphaPlus,alphaMinus, ...
    nivRaw,nivZ,'VariableNames',{'patientID','normalizedPatientID', ...
    'alphaPlus','alphaMinus','nivAsymmetryRaw','nivAsymmetryZ'});
nivSummary = struct('mean',nivMean,'sd',nivSD, ...
    'nParticipants',height(nivTable));
end

function [found,alphaPlus,alphaMinus,nivRaw,nivZ] = ...
    getNivParametersForPatient(nivTable,patientID)
matchingRows = nivTable.normalizedPatientID== ...
    normalizeParticipantID(patientID);
found = any(matchingRows);
alphaPlus = NaN; alphaMinus = NaN; nivRaw = NaN; nivZ = NaN;
if ~found, return; end
row = find(matchingRows,1,'first');
alphaPlus = nivTable.alphaPlus(row);
alphaMinus = nivTable.alphaMinus(row);
nivRaw = nivTable.nivAsymmetryRaw(row);
nivZ = nivTable.nivAsymmetryZ(row);
end

function normalizedID = normalizeParticipantID(patientID)
normalizedID = lower(strip(string(patientID)));
normalizedID = regexprep(normalizedID,'[^a-z0-9]','');
end

function offset = getAnalysisSeedOffset(analysisType)
if analysisType=="IT"
    offset = 1000;
elseif analysisType=="RT"
    offset = 2000;
else
    offset = 3000;
end
end

function q = BH(p)
p = p(:);
m = numel(p);
[s,o] = sort(p);
adj = s.*m./(1:m)';
adj = flipud(cummin(flipud(adj)));
adj = min(adj,1);
q = NaN(m,1);
q(o) = adj;
end
