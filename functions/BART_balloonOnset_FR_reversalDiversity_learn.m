function [BARTstats] = BART_balloonOnset_FR_reversalDiversity_learn(ptID)   % BARTstats
% BART_BALLOONONSET_FR_reversalDiversity_learn analyzes neuronal data for the BART task and looks for evidence of asymmetric learning as predicted
% in Distributional RL.
%
%   [BARTstats] = BART_balloonOnset_FR_reversalDiversity_learn(ptID) analyzes single unit data for
%   the patient specified in the string ptID. Output is raw an z-scored
%   firing rate.
%

% author: TAP20210721



%% TODO:::
% "to test for asymmetric learning, we modeled neuron responses w/ classic
% and distributional RL models to test which as a better fit to the data"
% -- mueller et al. "in all cases the model is used for each neuron, to
% predict the firing rate on each trial."
% fit linear models where y is FR and x is RPE using fitlm or other(max ...?) . 
% ?? where exactly are cue and outcome aligned differently?
% 

% ptID = '202118'; % This line for debugging only

%% Initialize Variables
rewardProbPvals = [];
outProbPvals = [];

%% loading neural and beahvioral data.
if ~isstring(ptID), ptID = num2str(ptID);,end;

nevList = dir(sprintf('//155.100.91.44/d/Data/preProcessed/BART_units/%s/Data/*.nev',ptID));
% nevList2 = dir(sprintf('X:/Data/preProcessed/BART_units/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
    nevFile(strfind(nevFile,'\'))='/'; % had to flip slashes so openNEV would load the correct file.
end
% [trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;

% loading behavioral matFile
matFile = sprintf('//155.100.91.44/d/Data/preProcessed/BART_units/%s/Data/%s.bartBHV.mat',ptID,ptID);
load(matFile)
pointsEarned = [data.points];
results = arrayfun(@(x) x.result, data, 'UniformOutput',false);

% standard 3-D [chan, unit, timestamp (seconds)] matrix.
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];

% channel deets.
inclChans = unique(ChanUnitTimestamp(:,1));
[myLabels,~,microLabels,microLocsMNI,generalLabels] = microLabelsBART(ptID);
inclChans(inclChans-96>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.\
nChans = length(inclChans);

%% task parameters in chronological order..
% There aren't any trigs that == 4
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);

% task identifiers
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
isCTRL = balloonIDs>10;

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
[~,sortedOutcomeIdcs] = sort(outcomeType);

% overall numbers of outcome types
nBanks = sum(outcomeType==1);
nPops = sum(outcomeType==2);

% adjusting for trial numbers
% only including complete trials; generally => excluding the last trial.
nTrials = min([length(outcomeType) length(balloonIDs)]);
inflateTimes = inflateTimes(1:nTrials);
balloonTimes = balloonTimes(1:nTrials);
balloonIDs = balloonIDs(1:nTrials);
isCTRL = isCTRL(1:nTrials);

% timing parameters.
pre = 2;
post = 3;


% firign rate timing parameters
binWidth = 50;
Fspikes = 1000;

%% color map
% reward probability colormap
cMap(1,:) = [0.5 0.5 0.5];  % gray
cMap(2,:) = [1 0 0];        % red
cMap(3,:) = [1 0.5 0];      % orange
cMap(4,:) = [1 0.9 0];      % yellow
cMap(5,:) = [1 0 1];        % controls

% risk and reward colormap
rcMap(1,:) = [0.5 0.5 0.5];
rcMap(2,:) = [1 0 0];
rcMap(3,:) = [1 0.5 0];
rcMap(4,:) = [1 0.9 0];
rcMap(5,:) = rgb('lightcoral');
rcMap(6,:) = rgb('rosybrown');
rcMap(7,:) = rgb('violet');

% colormap per trial. 
balloonEdgeColorMap = ones(length(balloonIDs),3)*.5;
riskRewardEdgeColorMap = zeros(length(balloonIDs),3);

% populating balloon color map
for x = 1:3
	balloonEdgeColorMap(balloonIDs==x,:) = repmat(cMap(x+1,:),sum(balloonIDs==x),1);
end

%% setting up Reward Variables
% regressors for linear models.
% [1 2 3 4 11 12 13 14] = [Y O R G Yc Oc Rc Gc]

% reward probability variable.
rewardProbability = zeros(nTrials,1);
% gray balloons
rewardProbability(balloonIDs==14) = 0;
% colored balloons
rewardProbability(balloonIDs==1) = cumsum(outcomeType(balloonIDs==1)==1)./(1:sum(balloonIDs==1))';
rewardProbability(balloonIDs==2) = cumsum(outcomeType(balloonIDs==2)==1)./(1:sum(balloonIDs==2))';
rewardProbability(balloonIDs==3) = cumsum(outcomeType(balloonIDs==3)==1)./(1:sum(balloonIDs==3))';
% collored controls
rewardProbability(balloonIDs==11 | balloonIDs==12 | balloonIDs==13) = 1;

% reward probability categories: gray, red, orange, yellow, controls (purple)
% increasing probability of reward (0 < red < orange < yellow < 1)
rewardProbCats = zeros(nTrials,1);
rewardProbCats(balloonIDs==14) = 1; % gray
rewardProbCats(balloonIDs==1) = 4; % yellow
rewardProbCats(balloonIDs==2) = 3; % orange
rewardProbCats(balloonIDs==3) = 2; % red
rewardProbCats(balloonIDs==11 | balloonIDs==12 | balloonIDs==13) = 5; % pink
categorical(rewardProbCats);
rs = unique(rewardProbCats);

%% RPEs
% According to Mueller paper, Rewarded trials have an r of 1 and unrewarded
% trials have an r of 0. RPEs are calculated by taking the result of the
% trial (1 if rewarded (banked), 0 if unrewarded (popped)) and subtracting
% the probability of reward of that trial.

% RPE = NaN(length(rewardProbability),1); % Initialize RPE vector
r = double(outcomeType); % convert outcomeType to double so it can be used in calculation w/ rewardProbability
r(r==2)=0; % set popped trials to a value of 0 instead of an index of 2.
RPE = r-rewardProbability; % this seems to check out. there are positive, negative and RPEs w/ a value of zero.
% [length(find(RPE>0)),length(find(RPE==0)),length(find(RPE<0))] % quick check

% alternatively, RPE r could be the number of points.
% RPE = data.points-rewardProbability;


unitcount = 1;

%% Loops to generate firing rates -- need to do modeling/B+ & B- in the loop for every unit.
% looping over Channels
for ch = nChans:-1:1
    % looping over number of units in the AP data. This just looks for
    % the number of unique units from channel ch.
    nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
    for un = 1:nUnits
        fprintf('\nprocessing and plotting for channel %d unit %d of %d',ch,un,nUnits)

        % getting unit times for the current channel and unit. getting the
        % timestamp from the correct channel/unit combo.
        unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
        if ~isempty(unitTimes) || length(unitTimes>5) 
            unitcount=unitcount+1;
            
            %% cue aligned spikes. 
            % loooping over trials
            for tt = 1:nTrials
                % putting the data in a structure
                spikes.channel(ch).unit(un).trial(tt).cueTimes = unitTimes(unitTimes>balloonTimes(tt)-pre & unitTimes<balloonTimes(tt)+post) - repmat(balloonTimes(tt)-pre,length(unitTimes(unitTimes>balloonTimes(tt)-pre & unitTimes<balloonTimes(tt)+post)),1);
                % generates FR in each bin.
                Rstc(tt,:) = psthBins(round(spikes.channel(ch).unit(un).trial(tt).cueTimes*Fspikes), binWidth, Fspikes, 1, (pre+post)*Fspikes);
            end % looping over trials
            FR(ch,un).Cue = Rstc;
            rawCueFR = Rstc;
% keyboard            

            % timing
            tSecCue = linspace(-pre, post, size(Rstc,2));

            % baseline correcting cue-aligned firing rates. when to do this?
            baseline = true;
            if baseline
                % define baseline time window (in seconds)
                bWin = [-1.2 -0.2];
                Rstc = Rstc-repmat(mean(Rstc(:,tSecCue>bWin(1) & tSecCue<bWin(2)),2),1,size(Rstc,2));
            else
                bWin = [-1.2 -0.2];
            end

            % setting up variables for statistics -- should I use this for
            % stuff besides ANOVAs??
            statTimeWin = [0.25 1.5]; % in seconds.
            %statTimeWin = [0.2 3.2]; % another option in seconds.

            % Get columns for analysis!
            win = [find(tSecCue == statTimeWin(1)), find(tSecCue==statTimeWin(2))];

            % z-score in window of analysis (?) -- will this need to be demeaned(?)
            % I think I'm going to need to do this manually as the matlab
            % function z scores everything instead of w/ respect to the
            % baseline...
%             zCueFR{ch,un} = zscore(Rstc(:,tSecCue>-pre+1 & tSecCue<post-1));
%             zCueFR{ch,un} = zscore(Rstc(:,tSecCue>=statTimeWin(1) & tSecCue<=statTimeWin(2)));
            bMean = mean(mean(Rstc(:,tSecCue>=bWin(1) & tSecCue<=bWin(2))));
            bSTD = std2(Rstc(:,tSecCue>=bWin(1) & tSecCue<=bWin(2))); % not confident this is the best way to do this...
            zCueFR{ch,un} = (Rstc-bMean)/bSTD;
            FR(ch,un).zCue = (Rstc-bMean)/bSTD;
            


            % Asymmetric Scaling (RPEs at Cue)
            % these? or more aptly put, use the z-scored data?
%             banksCue = zCueFR{ch,un}(outcomeType == 1,:);
%             banksCueFRTrials = mean(banksCue(:,tSecCue>=statTimeWin(1) & tSecCue <=statTimeWin(2)),2);
%             banksCueFRTime = mean(banksCue(:,tSecCue>=statTimeWin(1) & tSecCue <=statTimeWin(2)),1);
%             popsCue = zCueFR{ch,un}(outcomeType==2,:);
%             popsCueFRTrials = mean(popsCue(:,tSecCue>=statTimeWin(1) & tSecCue <=statTimeWin(2)),2);
%             popsCueFRTime = mean(popsCue(:,tSecCue>=statTimeWin(1) & tSecCue <=statTimeWin(2)),1);
            
            % B+ and B-. Regress the chosen cue probability against the firing rate at feedback
            % separately for rewarded and unrewarded trials. exclude FR = zero if in window -- averaged across window of interest (statTimeWin)
%             mdlBanksCue{ch,un} = fitlm(banksCueFRTrials, rewardProbability(outcomeType==1),'Exclude',banksCueFRTrials == min(banksCueFRTrials));
%             mdlPopsCue{ch,un} = fitlm(popsCueFRTrials, rewardProbability(outcomeType==2),'Exclude',popsCueFRTrials == min(popsCueFRTrials));
%             betaPlusCue{ch,un} = abs(mdlBanksCue{ch,un}.Coefficients{2,1}); %
%             betaMinusCue{ch,un} = abs(mdlPopsCue{ch,un}.Coefficients{2,1});
            % betas = [betaPlus,betaMinus]
            % tau's are for outputs.
%             tauCue{ch,un} = betaPlusCue{ch,un}/(betaPlusCue{ch,un}+betaMinusCue{ch,un}); % Tau is a measure of optimism.
%             tau{ch,un} = tauCue; % probably need to change to not be a struct so we can append them.


            % Cue aligned plots (Mueller Fig. 1b) ?? should I smooth this?
            % Mueller only plots these for individual units as an example
            % of optimistic/neutral and pessimistic neurons.
%             grayCueFR = Rstc(rewardProbCats==1,:);
%             meanGrayCueFR = mean(grayCueFR);
%             errGrayCueFR = std(grayCueFR)./sqrt(sum(rewardProbCats==1));
%             yellowCueFR = Rstc(rewardProbCats==4,:);
%             meanYellowCueFR = mean(yellowCueFR);
%             errYellowCueFR = std(yellowCueFR)./sqrt(sum(rewardProbCats==2));
%             orangeCueFR = Rstc(rewardProbCats==3,:);
%             meanOrangeCueFR = mean(orangeCueFR);
%             errOrangeCueFR = std(orangeCueFR)./sqrt(sum(rewardProbCats==3));
%             redCueFR = Rstc(rewardProbCats==2,:);
%             meanRedCueFR = mean(redCueFR);
%             errRedCueFR = std(redCueFR)./sqrt(sum(rewardProbCats==4));
%             controlCueFR = Rstc(rewardProbCats==5,:);
%             meanControlCueFR = mean(controlCueFR);
%             errControlCueFR = std(controlCueFR)./sqrt(sum(rewardProbCats==5));
% 
%             figure(101)
%             plot(tSecCue,meanGrayCueFR,'color',cMap(1,:))
%             hold on
%             patch([tSecCue fliplr(tSecCue)], [meanGrayCueFR-errGrayCueFR fliplr(meanGrayCueFR+errGrayCueFR)], cMap(1,:),'facealpha',0.5,'edgecolor','none')
%             plot(tSecCue,meanRedCueFR,'color',cMap(2,:))
%             patch([tSecCue fliplr(tSecCue)],[meanRedCueFR-errRedCueFR fliplr(meanRedCueFR + errRedCueFR)],cMap(2,:),'facealpha',0.5,'edgecolor','none')
%             plot(tSecCue,meanOrangeCueFR,'color',cMap(3,:))
%             patch([tSecCue fliplr(tSecCue)],[meanOrangeCueFR-errOrangeCueFR fliplr(meanOrangeCueFR+errOrangeCueFR)],cMap(3,:),'facealpha',0.5,'edgecolor','none')
%             plot(tSecCue,meanYellowCueFR,'color',cMap(4,:))
%             patch([tSecCue fliplr(tSecCue)],[meanYellowCueFR-errYellowCueFR fliplr(meanYellowCueFR+errYellowCueFR)],cMap(4,:),'facealpha',0.5,'edgecolor','none')
%             plot(tSecCue,meanControlCueFR,'color',cMap(5,:))
%             patch([tSecCue fliplr(tSecCue)],[meanControlCueFR-errControlCueFR fliplr(meanControlCueFR+errControlCueFR)],cMap(5,:),'facealpha',0.5,'edgecolor','none')
%             title('Firing Rate @ Cue')
%             xlabel('time relative to balloon appearance (s)')
%             xticks([-2 -1 0 1 2 3])
%             ylabel('Firing Rate (Hz)')
%             text(-1,-4,sprintf('Channel %i, unit %i',ch,un))

            % Mueller figure 2. Plots normalized firing rate with plot.
            % they look at normalized firing rate at choice and
            % banked/popped trials separately.

            %% outcome aligned spikes.
            % loooping over trials
            for tt2 = 1:nTrials
                % putting the data in a structure
                spikes.channel(ch).unit(un).trial(tt2).outTimes = unitTimes(unitTimes>outcomeTimes(tt2)-pre & unitTimes<outcomeTimes(tt2)+post) - repmat(outcomeTimes(tt2)-pre,length(unitTimes(unitTimes>outcomeTimes(tt2)-pre & unitTimes<outcomeTimes(tt2)+post)),1);

                Rsto(tt2,:) = psthBins(round(spikes.channel(ch).unit(un).trial(tt2).outTimes*Fspikes), binWidth, Fspikes, 1, (pre+post)*Fspikes);
            end % looping over trials
            FR(ch,un).Out= Rsto;
            rawOutFR = Rsto;
            

            % Outcome aligned plots (Mueller Fig. 1b) ?? should I smooth this?
            % probably should z score as well.
%             grayOutFR = Rsto(rewardProbCats==1,:);
%             meanGrayOutFR = mean(grayOutFR);
%             errGrayOutFR = std(grayOutFR)./sqrt(sum(rewardProbCats==1));
%             yellowOutFR = Rsto(rewardProbCats==4,:);
%             meanYellowOutFR = mean(yellowOutFR);
%             errYellowOutFR = std(yellowOutFR)./sqrt(sum(rewardProbCats==2));
%             orangeOutFR = Rsto(rewardProbCats==3,:);
%             meanOrangeOutFR = mean(orangeOutFR);
%             errOrangeOutFR = std(orangeOutFR)./sqrt(sum(rewardProbCats==3));
%             redOutFR = Rsto(rewardProbCats==2,:);
%             meanRedOutFR = mean(redOutFR);
%             errRedOutFR = std(redOutFR)./sqrt(sum(rewardProbCats==4));
%             controlOutFR = Rsto(rewardProbCats==5,:);
%             meanControlOutFR = mean(controlOutFR);
%             errControlOutFR = std(controlOutFR)./sqrt(sum(rewardProbCats==5));

%             figure(102)
%             plot(tSecCue,meanGrayOutFR,'color',cMap(1,:))
%             hold on
%             patch([tSecCue fliplr(tSecCue)], [meanGrayOutFR-errGrayOutFR fliplr(meanGrayOutFR+errGrayOutFR)], cMap(1,:),'facealpha',0.5,'edgecolor','none')
%             plot(tSecCue,meanRedOutFR,'color',cMap(2,:))
%             patch([tSecCue fliplr(tSecCue)],[meanRedOutFR-errRedOutFR fliplr(meanRedOutFR + errRedOutFR)],cMap(2,:),'facealpha',0.5,'edgecolor','none')
%             plot(tSecCue,meanOrangeOutFR,'color',cMap(3,:))
%             patch([tSecCue fliplr(tSecCue)],[meanOrangeOutFR-errOrangeOutFR fliplr(meanOrangeOutFR+errOrangeOutFR)],cMap(3,:),'facealpha',0.5,'edgecolor','none')
%             plot(tSecCue,meanYellowOutFR,'color',cMap(4,:))
%             patch([tSecCue fliplr(tSecCue)],[meanYellowOutFR-errYellowOutFR fliplr(meanYellowOutFR+errYellowOutFR)],cMap(4,:),'facealpha',0.5,'edgecolor','none')
%             plot(tSecCue,meanControlOutFR,'color',cMap(5,:))
%             patch([tSecCue fliplr(tSecCue)],[meanControlOutFR-errControlOutFR fliplr(meanControlOutFR+errControlOutFR)],cMap(5,:),'facealpha',0.5,'edgecolor','none')
%             title('Firing Rate @ outcome')
%             xlabel('time relative to outcome (s)')
%             xticks([-2 -1 0 1 2 3])
%             ylabel('Firing Rate (Hz)')

            %baseline
            if baseline
                % baseline correcting outcome-aligned data with the same window as the cue.
                Rsto = Rsto-repmat(mean(Rsto(:,tSecCue>bWin(1) & tSecCue<bWin(2)),2),1,size(Rsto,2));
            end
            
            
            

            % timing
            tSecOut = linspace(-pre, post, size(Rsto,2));

            % z-score -- will need to be demeaned(?) -- also, this like
            % doesn't set the baseline window to a score of 0, right?
%             zOutFR{ch,un} = zscore(Rsto(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)));
            % make own z score that whose baseline is pre-cue, right?
            zOutFR{ch,un} = (Rsto-bMean)/bSTD;
            FR(ch,un).zOut = (Rsto-bMean)/bSTD;

%             keyboard
           

            % find average z-scored FR across different reward prob. cats.
            % should these be timeStatWins?
            for z = length(unique(rewardProbCats)):-1:1
                rProbFRbarC{ch,un}(z) = mean(zCueFR{ch,un}(rewardProbCats==rs(z)));
                rProbFRerrC{ch,un}(z) = std(zCueFR{ch,un}(rewardProbCats==rs(z)))./sum(rewardProbCats==rs(z));
                rProbFRbarO{ch,un}(z) = mean(zOutFR{ch,un}(rewardProbCats==rs(z)));
                rProbFRerrO{ch,un}(z) = std(zOutFR{ch,un}(rewardProbCats==rs(z)))./sum(rewardProbCats==rs(z));
            end

            % attempt to demean -- from Elliot but didn't use it if
            % z-scored... Also, when I demeaned the z-scored data the mean
            % was 0
            %             rProbFRbarC{ch,un} = rProbFRbarC{ch,un}-repmat(mean(rProbFRbarC{ch,un}),1,length(rProbFRbarC{ch,un}));
            %             rProbFRbarO{ch,un} = rProbFRbarO{ch,un}-repmat(mean(rProbFRbarO{ch,un}),1,length(rProbFRbarO{ch,un}));

            % We're going to focus the initial analysis on unit w/
            % monotonic RPE/z-scored data. for 202118 only unit 3 seems to
            % do anything...

            %%%%% REVERSAL POINT CALC %%%%%%%%%%%%%%
            % 		RPy{ch,un} = mean(rProbFRbar{ch,un})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



            % Asymmetric Scaling (RPEs at analysis) -- should I z-score
            % these? or more aptly put, use the z-scored data?
%             banksOutcome = zOutFR{ch,un}(outcomeType == 1,:);
%             banksOutcomeFRTrials = mean(banksOutcome(:,tSecOut>=statTimeWin(1) & tSecOut <=statTimeWin(2)),2);
%             banksOutcomeFRTime = mean(banksOutcome(:,tSecOut>=statTimeWin(1) & tSecOut <=statTimeWin(2)),1);
%             popsOutcome = zOutFR{ch,un}(outcomeType==2,:);
%             popsOutcomeFRTrials = mean(popsOutcome(:,tSecOut>=statTimeWin(1) & tSecOut <=statTimeWin(2)),2);
%             popsOutcomeFRTime = mean(popsOutcome(:,tSecOut>=statTimeWin(1) & tSecOut <=statTimeWin(2)),1);
            % B+ and B-. Regress the chosen cue probability against the firing rate at feedback
            % separately for rewarded and unrewarded trials. exclude FR = zero if in window -- averaged across window of interest (statTimeWin)

            % Outcome models. Elliot used rewardProbCats not rewardProbability
%             mdlBanksOutcome{ch,un} = fitlm(banksOutcomeFRTrials, rewardProbCats(outcomeType==1),'Exclude',banksOutcomeFRTrials == min(banksOutcomeFRTrials));
%             mdlPopsOutcome{ch,un} = fitlm(popsOutcomeFRTrials, rewardProbCats(outcomeType==2),'Exclude',popsOutcomeFRTrials == min(popsOutcomeFRTrials));
%             betaPlusOutcome{ch,un} = abs(mdlBanksOutcome{ch,un}.Coefficients{2,1}); %
%             betaMinusOutcome{ch,un} = abs(mdlPopsOutcome{ch,un}.Coefficients{2,1});
            % betas = [betaPlus,betaMinus]
%             tauOutcome{ch,un} = betaPlusOutcome{ch,un}/(betaPlusOutcome{ch,un}+betaMinusOutcome{ch,un}); % Tau is a measure of optimism.
%             tau{ch,un} = tauOutcome;


            tbl = table(mean(FR(ch,un).zCue(:,tSecOut>=statTimeWin(1) & tSecOut <=statTimeWin(2)),2),rewardProbCats,'VariableNames',{'ZFR','pRewardCats'});
            tblOut = table(mean(FR(ch,un).zOut(:,tSecOut>=statTimeWin(1) & tSecOut <=statTimeWin(2)),2),rewardProbCats,'VariableNames',{'ZFROut','pRewardCats'});

            basic = true;
            if basic
                try
                    LM = fitglme(tbl,'ZFR ~ pRewardCats^2','Exclude',mean(FR(ch,un).zCue(:,tSecCue>=statTimeWin(1) & tSecCue <=statTimeWin(2)),2) == min(mean(FR(ch,un).zCue(:,tSecCue>=statTimeWin(1) & tSecCue <=statTimeWin(2)),2)));
                catch
                    LM = fitglme(tbl,'ZFR ~ pRewardCats^2');
                end
                %then do outcome aligned. check if betas invert!
                try
                    LMOut = fitglme(tblOut, 'ZFROut ~ pRewardCats^2','Exclude',mean(FR(ch,un).zOut(:,tSecOut>=statTimeWin(1) & tSecOut <=statTimeWin(2)),2) == min(mean(FR(ch,un).zOut(:,tSecOut>=statTimeWin(1) & tSecOut <=statTimeWin(2)),2)));
                catch
                    LMOut = fitglme(tblOut,'ZFROut ~ pRewardCats^2');
                end
            end % EHS has all kinds of options if basic isn't true. Could be added!
            BARTstats(unitcount).patientID = ptID;
            BARTstats(unitcount).ChanUnit = [ch un];
            BARTstats(unitcount).prePostTime = [-pre post];
            BARTstats(unitcount).unitlocation.generalLabels = generalLabels{floor((inclChans(ch)-96)./9)+1};
            BARTstats(unitcount).Cue = rawCueFR;
            BARTstats(unitcount).zCue = (Rstc-bMean)/bSTD;
            BARTstats(unitcount).Out = rawOutFR;
             BARTstats(unitcount).zOut = (Rsto-bMean)/bSTD;
             BARTstats(unitcount).tSecCue = tSecCue;
            BARTstats(unitcount).tSecOut = tSecOut;
            BARTstats(unitcount).statTimeWin = statTimeWin;
            BARTstats(unitcount).rewardProbLM = LM;
            BARTstats(unitcount).rewardProbANOVA = anova(LM);
            BARTstats(unitcount).rewardProbRsquared = LM.Rsquared;
            rewardProbPvals = cat(1,rewardProbPvals,BARTstats(unitcount).rewardProbANOVA.pValue(1));
            BARTstats(unitcount).rewardProbPvals = rewardProbPvals;
            BARTstats(unitcount).OutLM = LMOut;
            BARTstats(unitcount).OutLMANOVA = anova(LMOut);
            BARTstats(unitcount).outProbRSquared = LMOut.Rsquared;
            outProbPvals = cat(1,outProbPvals,BARTstats(unitcount).OutLMANOVA.pValue(1));
            BARTstats(unitcount).OutProbPvals = outProbPvals;
            BARTstats(unitcount).Points = [data.points];
            BARTstats(unitcount).Results =results;
            BARTstats(unitcount).Score = [data.score];
            BARTstats(unitcount).BalloonIDs = balloonIDs;


        else
            fprintf('\n... no spikes for this unit')
            BARTstats(1).patientID = ptID;
            BARTstats(1).prePostTime = [-pre post];
            BARTstats(1).ChanUnit = [0 0];
            BARTstats(1).unitlocation.NMM = {''};
            BARTstats(1).unitlocation.myLabels = {''};
            BARTstats(1).Cue = [];
            BARTstats(1).zCue = [];
            BARTstats(1).Out = [];
             BARTstats(1).zOut = [];
             BARTstats(1).tSecCue = [];
            BARTstats(1).tSecOut = [];
            BARTstats(1).statTimeWin = [];
            BARTstats(1).rewardProbLM =[];
            BARTstats(1).rewardProbANOVA = [];
            BARTstats(1).rewardProbRsquared =[];
%             rewardProbPvals = cat(1,rewardProbPvals,BARTstats(unitcount).rewardProbANOVA.pValue(1));
            BARTstats(1).rewardProbPvals = [];
            BARTstats(1).OutLM = [];
            BARTstats(1).OutLMANOVA = [];
            BARTstats(1).outProbRSquared = [];
%             outProbPvals = cat(1,outProbPvals,BARTstats(unitcount).outLMANOVA.pValue(1));
            BARTstats(1).OutProbPvals = [];
            BARTstats(1).Points = [];
            BARTstats(1).Results = [];
            BARTstats(1).Score = [];
            BARTstats(1).BalloonIDs = [];

                
        end % end ~isempty(unitTimes)
    end % end units loop
end % end channels loop

BARTstats(1) = []; % should get rid of the first 'line' of the struct for units w/ no spikes!

%% Met with Tim Muller 20231122. Helped me realize that the model fitting for asymmetric learning was fitting behavioral data.
% S, alpha+ and alpha- must be fit simultaneously. the TDlearn function already does a lot of this. Need to fit variables to actual data.
% not sure the best way to do so. Like I need to compare the predicted FR
% to the actual FR for RPE neurons! Somewhere I need to do a grid search of
% S and the two alpha values to see which is ideal!

% so reading the updated version, I think we need to fit on model to neural
% data and the other to the behavioral data. maybe?

% neural data is FR @ outcome (their window was 200-600 ms).

% pull out firing rates for every neuron on every trial % may need to Z-score...
if exist('FR','var')
    for k = 1:size(FR,1)
        for q = 1:size(FR,2)
            idxFR(k,q) = ~isempty(FR(k,q).Cue);
        end
    end
    unitsFR = FR(idxFR); % pull out the ch/un combinations that had spikes

    [c,u] = find(idxFR==1); % this just builds an index of chans/units that
    z = 1:nChans;

    % calculate delta (delta/RPE are the same thing, right) which is FR.Out - FR.Cue
    for k = 1:length(unitsFR)
        deltaFR{k} = unitsFR(k).Out - unitsFR(k).Cue;
        deltaFRz{k} = unitsFR(k).zOut - unitsFR(k).zCue;
    end
end

end % end function

%% Below this line __________________________________________________________________ Needs to go into across units analysis!
% models={'SSSL','SSAL','ASSL','ASAL'}; % okay, I think we have to do value prediction first... get RPE and feed that into scaling!
% clear predFR predFRRPE V RewardPE
% for mod = 1:length(models)
%     % Initalize variables.
%     R(1) = data(1).points;
%     X(1) = 0.5;								% starting with a coin flip for risk
%     V = zeros(1,nTrials);
%     Vrisk = zeros(1,nTrials);
%     a = .15; % alpha for SL model(s)
%     models{mod}
%     % start w/ learning!
%     for k = 1:numel(u)
%         if strcmp(models{mod}(3:4),'SL') % START Symmetric Learning
%             for t = 2:nTrials
%                 % reward varriable.
%                 if strcmp(data(t).result,'banked')
%                     R(t) = data(t).points;			% outcome on current trial.
%                 else
%                     R(t) = 0;
%                 end
% 
%                 % risk variable.
%                 if isCTRL(t)
%                     X(t) = 0;
%                 else
%                     X(t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
%                     % current trial risk is defined as P(pop) on previously observed balloons.
%                 end
% 
%                 % updating risk and reward PE
%                 RewardPE(t) = R(t) - V(t-1,1);
%                 XPE(t) = X(t) - Vrisk(t-1);
% 
%                 % updating value
%                 Vrisk(t) = Vrisk(t-1) + a*XPE(t-1);
%                 V(t,1) = V(t-1,1) + a*RewardPE(t-1);
%                 V(t,2) = V(t-1,2) + a*RPE(t-1);
%                 % if t==nTrials; keyboard; end;
%             end
% 
%             % Plot Symmetric Value Prediction!
%             figure % I don't know what to compare this too! Need to learn about 'ground truth'
%             hold on
%             plot(V(:,1))
%             plot(V(:,2))
%             legend('RewardPE','RPE')
%             title(sprintf('Value Estimates, Model: %s',models{mod}))
%             %                 subtitle(sprintf('S = %.3f, \alpha^+ = %.3f, \alpha^- = %.3f',S(q),alphaPlus(w),alphaMinus(w)))
%             if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}), 'dir')
%                 mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}))
%             end
%             halfMaximize(gcf,'page')
%             saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s/%s_%d_%d_Value_Model%s.pdf',ptID,models{mod},ptID,z(c(k)),z(u(k)),models{mod}))
%             close(gcf)
% 
% 
%             vanillaMDL{1,k} = fitlm(RewardPE,mean(unitsFRz.Out{k,1}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2));
%             beta0(1,k)=vanillaMDL{1,k}.Coefficients{1,1};
%             beta1(1,k)= vanillaMDL{1,k}.Coefficients{2,1};
%             vanillaMDL{2,k} = fitlm(RPE,mean(unitsFRz.Out{k,1}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2));
%             beta0(2,k) = vanillaMDL{2,k}.Coefficients{1,1};
%             beta1(2,k) = vanillaMDL{2,k}.Coefficients{2,1};
% 
%             posidx(:,1) = RewardPE>0; posidx(:,2) = RPE >0;
%             negidx(:,1) = RewardPE<=0; negidx(:,2) = RPE<=0;
%             if strcmp(models{mod}(1:2),'SS')
%                 predFR = beta0(1,k) + (beta1(1,k).*RewardPE);
%                 predFRRPE = beta0(2,k) + (beta1(2,k).*RPE);
% 
%                 % Plot Symmetric Scaling FR Prediction
%                 figure
%                 hold on
%                 plot(mean(unitsFRz.Out{k}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2))
%                 plot(predFR)
%                 plot(predFRRPE)
%                 legend('Actual Z-Scored','RewardPRE','RPE')
%                 title(sprintf('Firing Rate Model: %s',models{mod}))
%                 if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}), 'dir')
%                     mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}))
%                 end
%                 halfMaximize(gcf,'page')
%                 saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s/%s_%d_%d_FR_Model_%s.pdf',ptID,models{mod},ptID,z(c(k)),z(u(k)),models{mod}))
%                 close(gcf)
% 
%             elseif strcmp(models{mod(1:2),'AS'})
%                 S = 0:.025:1;
%                 for q = 1:length(S)
%                     predFR(posidx(:,1)) = beta0(1,k) + (beta1(1,k).*RewardPE(posidx(:,1))*S(q));
%                     predFR(negidx(:,1)) = beta0(1,k) + (beta1(1,k).*RewardPE(negidx(:,1))*(1-S(q)));
%                     predFRRPE(posidx(:,2)) = beta0(2,k) + (beta1(2,k).*RPE(posidx(:,2))*S(q));
%                     predFRRPE(negidx(:,2)) = beta0(2,k) + (beta1(2,k).*RPE(negidx(:,2))*(1-S(q)));
% 
%                     % Plot Asymmetric Scaling FR Prediction
%                     figure
%                     hold on
%                     plot(mean(unitsFRz.Out{k}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2))
%                     plot(predFR)
%                     plot(predFRRPE)
%                     legend('Actual Z-Scored','RewardPRE','RPE')
%                     title(sprintf('Firing Rate Model: %s (S = %.3f, 1-S = %.3f)',models{mod},S(q),1-S(q)))
%                     if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}), 'dir')
%                         mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}))
%                     end
%                     halfMaximize(gcf,'page')
%                     saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s/%s_%d_%d_FR_Model_%s_%.3fS.pdf',ptID,models{mod},ptID,z(c(k)),z(u(k)),models{mod},S(q)))
%                     close(gcf)
% 
%                 end
%             end % end scaling for SL condition
% 
%         elseif strcmp(models{mod}(3:4),'AL')
%             alphaPlus = 0:0.025:1;
%             alphaMinus = 0:0.025:1;
%             for w = 1:length(alphaPlus)
%                 for m = 1:length(alphaMinus)
%                     for t = 2:nTrials
%                         % reward varriable.
%                         if strcmp(data(t).result,'banked')
%                             R(t) = data(t).points;			% outcome on current trial.
%                         else
%                             R(t) = 0;
%                         end
% 
%                         % risk variable.
%                         if isCTRL(t)
%                             X(t) = 0;
%                         else
%                             X(t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
%                             % current trial risk is defined as P(pop) on previously observed balloons.
%                         end
% 
%                         % updating risk and reward PE
%                         RewardPE(t) = R(t) - V(t-1,1);
%                         XPE(t) = X(t) - Vrisk(t-1);
% 
%                         % updating value
%                         Vrisk(t) = Vrisk(t-1) + a*XPE(t-1);
%                         %                         V(t) = V(t-1) + a*RewardPE(t-1);
%                         % need to identify if trial has +/- RPE.
%                         if RewardPE(t) >0
%                             V(t,1) = V(t-1,1)+alphaPlus(w)*RewardPE(t-1);
%                         else
%                             V(t,1) = V(t-1,1)+alphaMinus(m)*RewardPE(t-1);
%                         end
%                         if RPE(t)>0
%                             V(t,2) = V(t-1,2)+alphaPlus(w)*RPE(t-1);
%                         else
%                             V(t,2) = V(t-1,2)+alphaMinus(m)*RPE(t-1);
%                         end
%                     end
% 
%                     % Plot Asymmetric Value Prediction!
%                     figure % I don't know what to compare this too!
%                     hold on
%                     plot(V(:,1))
%                     plot(V(:,2))
%                     legend('RewardPE','RPE')
%                     title(sprintf('Value Estimates, Model: %s',models{mod}))
%                     subtitle(sprintf('{\alpha}^+ = %.3f, {\alpha}^- = %.3f',alphaPlus(w),alphaMinus(m)))
%                     if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}), 'dir')
%                         mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}))
%                     end
%                     halfMaximize(gcf,'page')
%                     saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s/%s_%d_%d_Value_Model%s_alphaplus%.3f_and_%.3fminus.pdf',ptID,models{mod},ptID,z(c(k)),z(u(k)),models{mod},alphaPlus(w),alphaMinus(m)))
%                     close(gcf)
% 
% 
%                     % Use RPEs from AL model to get Betas for linear models
%                     vanillaMDL{1,k} = fitlm(RewardPE,mean(unitsFRz.Out{k,1}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2));
%                     beta0(1,k)=vanillaMDL{1,k}.Coefficients{1,1};
%                     beta1(1,k)= vanillaMDL{1,k}.Coefficients{2,1};
%                     vanillaMDL{2,k} = fitlm(RPE,mean(unitsFRz.Out{k,1}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2));
%                     beta0(2,k) = vanillaMDL{2,k}.Coefficients{1,1};
%                     beta1(2,k) = vanillaMDL{2,k}.Coefficients{2,1};
% 
%                     posidx(:,1) = RewardPE>0; posidx(:,2) = RPE >0;
%                     negidx(:,1) = RewardPE<=0; negidx(:,2) = RPE<=0;
% 
%                     if strcmp(models{mod}(1:2),'SS')
%                         predFR = beta0(1,k) + (beta1(1,k).*RewardPE);
%                         predFRRPE = beta0(2,k) + (beta1(2,k).*RPE);
% 
%                         % Plot Symmetric Scaling FR Prediction
%                         figure
%                         hold on
%                         plot(mean(unitsFRz.Out{k}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2))
%                         plot(predFR)
%                         plot(predFRRPE)
%                         legend('Actual Z-Scored','RewardPRE','RPE')
%                         title(sprintf('Firing Rate Model: %s',models{mod}))
%                         subtitle(sprintf('{\alpha}^+ = %.3f, {\alpha}^- = %.3f',alphaPlus(w),alphaMinus(m)))
%                         if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}), 'dir')
%                             mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}))
%                         end
%                         halfMaximize(gcf,'page')
%                         saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s/%s_%d_%d_FR_Model%s_alphaplus%.3f_and_%.3fminus.pdf',ptID,models{mod},ptID,z(c(k)),z(u(k)),models{mod},alphaPlus(w),alphaMinus(m)))
%                         close(gcf)
% 
% 
%                     elseif strcmp(models{mod(1:2),'AS'})
%                         S = 0:.025:1;
%                         for q = 1:length(S)
%                             predFR(posidx(:,1)) = beta0(1,k) + (beta1(1,k).*RewardPE(posidx(:,1))*S(q));
%                             predFR(negidx(:,1)) = beta0(1,k) + (beta1(1,k).*RewardPE(negidx(:,1))*(1-S(q)));
%                             predFRRPE(posidx(:,2)) = beta0(2,k) + (beta1(2,k).*RPE(posidx(:,2))*S(q));
%                             predFRRPE(negidx(:,2)) = beta0(2,k) + (beta1(2,k).*RPE(negidx(:,2))*(1-S(q)));
% 
%                             % Plot Asymmetric Scaling FR Prediction
%                             figure % I don't know what to compare this too!
%                             hold on
%                             plot(mean(unitsFRz.Out{k}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2))
%                             plot(predFR)
%                             plot(predFRRPE)
%                             legend('Actual Z-Scored','RewardPRE','RPE')
%                             title(sprintf('Firing Rate, Model: %s',models{mod}))
%                             subtitle(sprintf('S = %.3f, {\alpha}^+ = %.3f, {\alpha}^- = %.3f',S(q),alphaPlus(w),alphaMinus(m)))
%                             if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}), 'dir')
%                                 mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s',ptID,models{mod}))
%                             end
%                             halfMaximize(gcf,'page')
%                             saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s/%s_%d_%d_FR_Model%s_S%.3f_alphas%.3fand%.3fplusminus.pdf',ptID,models{mod},ptID,z(c(k)),z(u(k)),models{mod},S(q),alphaPlus(w),alphaMinus(m)))
%                             close(gcf)
%                         end % end S
%                     end % end scaling for AL condition
%                 end % end alpha minus
%             end %end alpha plus
%         end % end Learning
%     end % end looping through units
% end % end models
% %%
% end % end function
 

 %% Below is backwards!
% 
%         % this has to go inside each Learning Model
%     for k = 1:length(deltaFR)
%         vanillaMDL{1,k} = fitlm(RewardPE,mean(unitsFRz.Out{k,1}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2));
%         beta0(1,k)=vanillaMDL{1,k}.Coefficients{1,1};
%         beta1(1,k)= vanillaMDL{1,k}.Coefficients{2,1};
%         vanillaMDL{2,k} = fitlm(RPE,mean(unitsFRz.Out{k,1}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2));
%         beta0(2,k) = vanillaMDL{2,k}.Coefficients{1,1};
%         beta1(2,k) = vanillaMDL{2,k}.Coefficients{2,1};
%     end
% 
% 
%     %     for k = 1:length(deltaAsymmetricLearning)
%     %         posidx{k} = squeeze(mean(deltaAsymmetricLearning{k}(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2))>0;
%     %         negidx{k} = squeeze(mean(deltaAsymmetricLearning{k}(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2))<=0;
%     %     end
% 
%     posidx(:,1) = RewardPE>0; posidx(:,2) = RPE >0;
%     negidx(:,1) = RewardPE<=0; negidx(:,2) = RPE<=0;
% 
%     a=.5; % chosen alpha value for symmetric learning
% 
%     for k = 1:length(beta1)
%         if strcmp(models{mod}(1:2),'SS')
%             predFR = beta0(1,k) + (beta1(1,k).*RewardPE);
%             predFRRPE = beta0(2,k) + (beta1(2,k).*RPE);
% 
%             if strcmp(models{mod}(3:4),'SL')
%                 for t = 2:nTrials
%                                         
%                 end
% 
%                 figure % I don't know what to compare this too!
%                 hold on
%                 plot(V(:,1))
%                 plot(V(:,2))
%                 legend('RewardPE','RPE')
%                 title(sprintf('Value Estimates, Model: %s',models{mod}))
%                 %                 subtitle(sprintf('S = %.3f, \alpha^+ = %.3f, \alpha^- = %.3f',S(q),alphaPlus(w),alphaMinus(w)))
%                 if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s',ptID), 'dir')
%                     mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s',ptID))
%                 end
%                 halfMaximize(gcf,'page')
%                 saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s_%d_%d_Value_Model%s.pdf',ptID,ptID,z(c(k)),z(u(k)),models{mod}))
%                 close(gcf)
% 
%                 
% 
%             elseif strcmp(models{mod}(3:4),'AL')
%                 alphaPlus = 0:0.025:1;
%                 alphaMinus = 0:0.025:1;
%                 for w = 1:length(alphaPlus)
%                     for m = 1:length(alphaMinus)
%                         for t = 2:nTrials
%                             % reward varriable.
%                             %                     if strcmp(data(t).result,'banked')
%                             %                         R(t) = data(t).points;			% outcome on current trial.
%                             %                     else
%                             %                         R(t) = 0;
%                             %                     end
% 
%                             % risk variable.
%                             if isCTRL(t)
%                                 X(t) = 0;
%                             else
%                                 X(t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
%                                 % current trial risk is defined as P(pop) on previously observed balloons.
%                             end
% 
%                             % updating risk and reward PE
%                             %                     RewardPE(t) = R(t) - V(t-1);
%                             XPE(t) = X(t) - Vrisk(t-1);
% 
%                             % updating value
%                             Vrisk(t) = Vrisk(t-1) + a*XPE(t-1);
%                             %                             V(t) = V(t-1) + a*RewardPE(t-1);
%                             % need to identify if trial has +/- RPE.
%                             if RewardPE(t) >0
%                                 V(t,1) = V(t-1,1)+alphaPlus(w)*predFR(t-1);
%                             else
%                                 V(t,1) = V(t-1,1)+alphaMinus(m)*predFR(t-1);
%                             end
%                             if RPE(t)>0
%                                 V(t,2) = V(t-1,2)+alphaPlus(w)*predFRRPE(t-1);
%                             else
%                                 V(t,2) = V(t-1,2)+alphaMinus(m)*predFRRPE(t-1);
%                             end
% 
%                             % if t==nTrials; keyboard; end;
%                         end
%                         
% 
%                         figure
%                         hold on
%                         plot(mean(unitsFRz.Out{k}(:,tSecOut>=statTimeWin(1) & tSecOut<=statTimeWin(2)),2))
%                         plot(predFR)
%                         plot(predFRRPE)
%                         legend('Actual Z-Scored','RewardPRE','RPE')
%                         title(sprintf('Firing Rate Model: %s',models{mod}))
%                         if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s',ptID), 'dir')
%                             mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s',ptID))
%                         end
%                         halfMaximize(gcf,'page')
%                         saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s_%d_%d_FR_Model_%s.pdf',ptID,ptID,z(c(k)),z(u(k)),models{mod}))
%                         close(gcf)
%                     end % end alphaMinus
%                 end % end alphaPlus
%             end
% 
%         elseif strcmp(models{mod}(1:2),'AS')
%             S = 0:.025:1;
%             for q = 1:length(S)
%                 predFR(posidx(:,1)) = beta0(1,k) + (beta1(1,k).*RewardPE(posidx(:,1))*S(q));
%                 predFR(negidx(:,1)) = beta0(1,k) + (beta1(1,k).*RewardPE(negidx(:,1))*(1-S(q)));
%                 predFRRPE(posidx(:,2)) = beta0(2,k) + (beta1(2,k).*RPE(posidx(:,2))*S(q));
%                 predFRRPE(negidx(:,2)) = beta0(2,k) + (beta1(2,k).*RPE(negidx(:,2))*(1-S(q)));
% 
%                 if strcmp(models{mod}(3:4),'SL')
%                     for t = 2:nTrials
%                         % reward varriable.
%                         %                         if strcmp(data(t).result,'banked')
%                         %                             R(t) = data(t).points;			% outcome on current trial.
%                         %                         else
%                         %                             R(t) = 0;
%                         %                         end
% 
%                         % risk variable.
%                         if isCTRL(t)
%                             X(t) = 0;
%                         else
%                             X(t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
%                             % current trial risk is defined as P(pop) on previously observed balloons.
%                         end
% 
%                         % updating risk and reward PE
%                         %                         RewardPE(t) = R(t) - V(t-1);
%                         XPE(t) = X(t) - Vrisk(t-1);
% 
%                         % updating value
%                         Vrisk(t) = Vrisk(t-1) + a*XPE(t-1);
%                         %         V(t) = V(t-1) + a*RewardPE(t-1);
%                         V(t,1) = V(t-1,1) + a*RewardPE(t-1);
%                         V(t,2) = V(t-1,2) + a*RPE(t-1);
% 
% 
%                         % if t==nTrials; keyboard; end;
%                     end
%                 elseif strcmp(models{mod}(3:4),'AL')
%                     alphaPlus = 0:0.025:1;
%                     alphaMinus = 0:0.025:1;
%                     for w = 1:length(alphaPlus)
%                         for m = 1:length(alphaMinus)
%                             for t = 2:nTrials
%                                 % reward varriable.
%                                 %                         if strcmp(data(t).result,'banked')
%                                 %                             R(t) = data(t).points;			% outcome on current trial.
%                                 %                         else
%                                 %                             R(t) = 0;
%                                 %                         end
% 
%                                 % risk variable.
%                                 if isCTRL(t)
%                                     X(t) = 0;
%                                 else
%                                     X(t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
%                                     % current trial risk is defined as P(pop) on previously observed balloons.
%                                 end
% 
%                                 % updating risk and reward PE
%                                 %                         RewardPE(t) = R(t) - V(t-1);
%                                 XPE(t) = X(t) - Vrisk(t-1);
% 
%                                 % updating value
%                                 Vrisk(t) = Vrisk(t-1) + a*XPE(t-1);
%                                 %         V(t) = V(t-1) + a*RewardPE(t-1);
%                                 if RewardPE(t) >0
%                                     V(t,1) = V(t-1,1)+alphaPlus(w)*predFR(t-1);
%                                 else
%                                     V(t,1) = V(t-1,1)+alphaMinuw(m)*predFR(t-1);
%                                 end
%                                 if RPE(t)>0
%                                     V(t,2) = V(t-1,2)+alphaPlus(w)*predFRRPE(t-1);
%                                 else
%                                     V(t,2) = V(t-1,2)+alphaMinus(m)*predFRRPE(t-1);
%                                 end
% 
% 
%                                 % if t==nTrials; keyboard; end;
%                             end % end trials of value estimation
%                             figure % I don't know what to compare this too!
%                             hold on
%                             plot(V(:,1))
%                             plot(V(:,2))
%                             legend('RewardPE','RPE')
%                             title(sprintf('Value Estimates, Model: %s',models{mod}))
%                             subtitle(sprintf('S = %.3f, \alpha^+ = %.3f, \alpha^- = %.3f',S(q),alphaPlus(w),alphaMinus(m)))
%                             if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s',ptID), 'dir')
%                                 mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s',ptID))
%                             end
%                             halfMaximize(gcf,'page')
%                             saveas(gcf,sprintf('D:/Data/Alex/BART/MuellerAnalysis/AsymmetricLearning/%s/%s_%d_%d_Value_Model%s_S%.3f_alphas%.3fand%.3fplusminus.pdf',ptID,ptID,z(c(k)),z(u(k)),models{mod},S(q),alphaPlus(w),alphaMinus(m)))
%                             close(gcf)
%                         end % end iterating through alpha minus
%                     end % end alpha plus
%                 end % end selection of learning model
%                 
%                 % TODO::: look at R^2
%             end % end S
%         end % end selection of scaling model
%     end % end units
% end % end models (4 types)







% Believe below is junk...
%% Asymmetric Learning!
% Should I just use TDlearn func?
% V<- V + alpha*delta; delta = r-V, alpha = learning rate
% for DistRL: V<- V + alpha+*delta, delta>0 & V<- V + alpha-*delta. delta
% <=0
% First, Asymmetric Scaling
% to fit neural data to RL models they predicted the firing rate at outcome
% from RPE. (FR = B0+B1*delta) For dist: FR = B0+B+*delta, delta >0 and FR
% = B0+B-*delta, delta <= 0;

% Replaced asymmetric scaling equations.
% FR = B0+B1*delta*S, delta >0
% FR = B0+B1*delta(1-S), delta <=0
% S is bound between 0 and 1. if S is near 1, ositive RPE 
% B0 and B1 are the same in both equation and were fir in the same
% regression model. GET BETA0 AND BETA1. Inputs are going to be delta (1
% per trial, so average across time)


% Estimating the parameters:
% regressors are generated by passing through the model the option chosen
% and the reward observed on each trial of the training set (would that mean
% balloon for BART?)

%% after talking w/ Elliot we are going to regress the FR = B0+B1*delta*S
% WE NEED S TO BE A SLOPE NOT REALLY A SCALAR.... means a fitLM is probs
% the way to go.
% 
% 
% x= -10:.1:10;
% xpos = 0:.1:15;
% xneg = -15:.1:0;
% for q = 1:length(S)
%     for k = 1:length(deltaAsymmetricLearning)
%         y{k} = polyval(P(k,:),x);
%         figure(k+50)
%         scatter(squeeze(mean(deltaAsymmetricLearning{k}(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)),squeeze(mean(unitsFR(k).Out(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)))
%         hold on
%         plot(x,y{k})
% 
%         figure((k+60)*q)
%         scatter(squeeze(mean(deltaAsymmetricLearning{k}(posidx{k},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)),squeeze(mean(unitsFR(k).Out(posidx{k},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)))
%         hold on
%         ypos{q,k} = polyval([(beta1(k)*S(q)) beta0(k)],xpos);
%         title('Pos')
%         plot(ypos{q,k})
% 
%         figure((k+70)*q)
%         scatter(squeeze(mean(deltaAsymmetricLearning{k}(negidx{k},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)),squeeze(mean(unitsFR(k).Out(negidx{k},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)))
%         hold on
%         yneg{q,k} = polyval([(beta1(k)*(1-S(q))) beta0(k)],xneg);
%         title('Neg')
%     end
% end
% 
% 
% % SO I think my issue was that I was visualizing a line w/ coeefs that spanned a new x range, but was calculating R^2 for the data that got scaled by S instead...
% S = 0:.025:1;
% 
% from there we get beta0 and beta1. then we can grid search s!
% for k = 1:length(deltaAsymmetricLearning)
%     mdlTrad{k}=fitlm(squeeze(mean(deltaAsymmetricLearning{k}(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)),squeeze(mean(unitsFR(k).Out(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
%     beta0(k) = mdlTrad{k}.Coefficients{1,1};
%     beta1(k) = mdlTrad{k}.Coefficients{2,1};
% end
% 
% now that we have betas, we can fit values of s -- need to add an outer loop that goes through every ch/un
% for q = 1:length(deltaAsymmetricLearning)
%     keyboard
%     close all
%     for k = 1:length(S)
%         feedpp{q,k} = beta1(q)*S(k);
%         feedpn{q,k} = beta1(q)*(1-S(k));
%         checkppx{q,k} = (0:.1:14)*beta1(q);
%         checkpnx{q,k} = (-16:.1:0)*beta1(q);
%         pp{q,k} = polyval([feedpp{q,k},beta0(q)], mean(deltaAsymmetricLearning{q}(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
%         pn{q,k} = polyval([feedpn{q,k},beta0(q)], mean(deltaAsymmetricLearning{q}(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
%         checkpp{q,k} = polyval([S(k),beta0(q)], checkppx{q,k});
%         checkpn{q,k} = polyval([(1-S(k)),beta0(q)], checkpnx{q,k});
%         recheckpp{q,k} = polyval([S(k)*beta1(q), beta0(q)],0:.1:14);
%         recheckpn{q,k} = polyval([(1-S(k))*beta1(q),beta0(q)], -16:.1:0);
% 
%         feedmdlpp{q,k} = beta0(q)+(beta1(q)*mean(deltaAsymmetricLearning{q}(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
%         feedmdlpn{q,k} = beta0(q)+(beta1(q)*mean(deltaAsymmetricLearning{q}(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
%         mdlSp{q,k} = fitlm(mean(deltaAsymmetricLearning{q}(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)*S(k),feedmdlpp{q,k});
%         mdlSn{q,k} = fitlm(mean(deltaAsymmetricLearning{q}(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)*(1-S(k)),feedmdlpn{q,k});
%         R2pnmdl(q,k) = mdlSn{q,k}.Rsquared.Ordinary ;
%         R2ppmdl(q,k) = mdlSp{q,k}.Rsquared.Ordinary;
% 
%         figure(k+500)
%         subplot(3,1,1)
%         scatter(mean(deltaAsymmetricLearning{q}(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2),mean(unitsFR(q).Out(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
%         xlabel('\delta')
%         ylabel('Mean FR (Outcome)')
%         
%         subplot(3,1,2)
%         plot(mean(deltaAsymmetricLearning{q}(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2),pp{q,k},'LineWidth',2.5)
%         hold on
%         scatter(mean(deltaAsymmetricLearning{q}(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2),mean(unitsFR(q).Out(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
%         plot(checkppx{q,k},checkpp{q,k},'--k')
%         plot(0:.1:14,recheckpp{q,k},'--g')
%         title('PP FR \delta >0')
%         text(1,1,sprintf('S = %.03f',S(k)))
%         legend('pp','','\beta1*S','\beta1')
%         xlabel('\delta')
%         ylabel('Mean FR (Outcome)')
% 
%         subplot(3,1,3)
%         plot(mean(deltaAsymmetricLearning{q}(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2),pn{q,k},'LineWidth',2.5)
%         hold on
%         scatter(mean(deltaAsymmetricLearning{q}(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2), mean(unitsFR(q).Out(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2))
%         plot(checkpnx{q,k},checkpn{q,k},'--k')
%         plot(-16:.1:0,recheckpn{q,k},'--g')
%         title('PN FR \delta <=0')
%         legend('pn','','\beta1*(1-S)','\beta1')
%         xlabel('\delta')
%         ylabel('Mean FR (Outcome)')
% 
%         halfMaximize(k+500,'page')
%         if ~exist(sprintf('D:/Data/Alex/BART/MuellerAnalysis/S/%s',ptID), 'dir')
%             mkdir(sprintf('D:/Data/Alex/BART/MuellerAnalysis/S/%s',ptID))
%         end
%         [c,u] = find(idxFR==1);
%         z = 1:nChans;
%         saveas(k+500,sprintf('D:/Data/Alex/BART/MuellerAnalysis/S/%s/%s_chan%d_unit%d_%.3f_SPlots.pdf',ptID,ptID,z(c(q)),z(u(q)),S(k)))
% 
%             figure(600+k)
%             mdlSp{k}.plot;
%             figure(700+k)
%             mdlSn{k}.plot;
% 
%         calculate R2
%         calculate residuals sum((actual-predicted).^2)
%         ppRes(q,k) = sum((mean(unitsFR(q).Out(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)-pp{q,k}).^2);
%         pnRes(q,k) = sum((mean(unitsFR(q).Out(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)-pn{q,k}).^2);
% 
%         calculate sum((y-Ybar).^2)
%         ppDifMean = sum((mean(unitsFR(q).Out(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)-mean(mean(unitsFR(q).Out(posidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2))).^2);
%         pnDifMean = sum((mean(unitsFR(q).Out(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)-mean(mean(unitsFR(q).Out(negidx{q},tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2))).^2);
% 
%         do final R2 calculation
%         R2pp(q,k) = 1 - (ppRes(q,k)/ppDifMean);
%         R2pn(q,k) = 1 - (pnRes(q,k)/pnDifMean);
%         
%     end
%     BestS(q,:) =[S(find(R2pp(q,:)==max(R2pp(q,:)))), S(find(R2pn(q,:)==max(R2pn(q,:))))]
%     BestSmdl(q,:) = [S(find(R2ppmdl(q,:)==max(R2ppmdl(q,:)),1)), S(find(R2pnmdl(q,:)==max(R2pnmdl(q,:)),1))]
%     BARTstats(q).BestS = BestS(q,:);
%     BARTstats(q).BestSmdl = BestSmdl(q,:);
% end
% 
% 
% try to find best S value for both pos/neg RPEs. Make a mesh of S x (1-S)
% w/ a geometry of R squared values.
% 
% %


