function [reversals,taus,STATS, halfReversals, halfTaus, testTaus, halftestTaus] = BART_balloonOnset_FR_reversalDiversity_riskModel(ptID,plotFlag)   % BARTstats
% BART_BALLOONONSET_FR_
% author: EHS20220711
%  input args
% ptID = '202202';
% plotFlag = true; 
% Output:
% reversals cell structure
                % reversals{:,1} = RPc (this is an attempt at quantifying the reversal points for units w/ more than one crossing of zero when cue aligned.
                % reversals{:,2} = RPo (this is an attempt at quantifying the reversal points for units w/ more than one crossing of zero when outcome aligned.
                        % RPc(1) = max/min
                        % RPc(n) = linear model of first time when RPE crosses 0
                        % RPc(n_end) = final linear model of first time when RPE crosses 0
                        % RPc(end) = mean of RP(n)s to compare against RPc(1)
                % reversals{:,3} = RP (this is an attempt at estimating the reversal point when neuron has only 1 crossing of zero when cue aligned.
                % reversals{:,4} = RPOut (this is an attempt at estimating the reversal point when neuron has only 1 crossing of zero when outcome aligned.
                % reversals{:,5} = logical indices aiming to identify 'monotonic' units when cue aligned
                % reversals{:,6} = logical indices aiming to identify 'monotonic' units when outcome aligned
                % reversals{:,7} = patient ID
                % reversals{:,8} = chan ID
                % reversals{:,9} = unit ID
                % reversals{:,10} = CuePolyOrder
                % reversals{:,11} = OutPolyOrder
               
%halfReversals Cell structure
%                 halfReversals{:,1} = RPcg1
%                 halfReversals{:,2} = RPcg2
%                 halfReversals{:,3} = RPog1
%                 halfReversals{:,4} = RPog2
%                 halfReversals{:,5} = linear encoding for g1 cue
%                 halfReversals{:,6} = linear encoding for g2 cue
%                 halfReversals{:,7} = linear encoding for g1 out
%                 halfReversals{:,8} = linear encoding for g2 out
%                 halfReversals{:,9} = RPg1
%                 halfReversals{:,10} = RPg2
%                 halfReversals{:,11} =RPOutg1
%                 halfReversals{:,12} = RPOutg2
%                 halfReversals{:,13} = channel
%                 halfReversals{:,14} = unit
%                 halfReversals{:,15} = ptID
%                 halfReversals{:,16} = BARTstats.bestCuePolyOrderg1;
%                 halfReversals{:,17} = BARTstats.bestOutPolyOrderg1;
%                 halfReversals{:,18} = BARTstats.bestCuePolyOrderg2;
%                 halfReversals{:,19} = BARTstats.bestOutPolyOrderg2;
                


%halfTaus Structure
%                 halfTau(:,1) = taug1
%                 halfTau(:,2) = taug2
%                 halfTau(:,3) = betaPlus
%                 halfTau(:,4) = betaMinus
%                 halfTau(:,5) = betaPlusg1
%                 halfTau(:,6) = betaMinusg1
%                 halfTau(:,7) = betaPlusg1
%                 halfTau(:,8) = betaminusg2


psthType = 'bins'; %'chronux' or 'bins'

% initializing reversal point var
taus = [];
reversals = {};
STATS = struct();
rewardProbPvals = [];
halfReversals = {};
halfTaus = [];
testTaus = [];
halftestTaus=[];

% reversals cell structure
                % reversals{:,1} = RPc (this is an attempt at quantifying the reversal points for units w/ more than one crossing of zero when cue aligned.
                % reversals{:,2} = RPo (this is an attempt at quantifying the reversal points for units w/ more than one crossing of zero when outcome aligned.
                        % RPc(1) = max/min
                        % RPc(n) = linear model of first time when RPE crosses 0
                        % RPc(n_end) = final linear model of first time when RPE crosses 0
                        % RPc(end) = mean of RP(n)s to compare against RPc(1)
                % reversals{:,3} = RP (this is an attempt at estimating the reversal point when neuron has only 1 crossing of zero when cue aligned.
                % reversals{:,4} = RPOut (this is an attempt at estimating the reversal point when neuron has only 1 crossing of zero when outcome aligned.
                % reversals{:,5} = logical indices aiming to identify 'monotonic' units when cue aligned
                % reversals{:,6} = logical indices aiming to identify 'monotonic' units when outcome aligned
                % reversals{:,7} = patient ID
                % reversals{:,8} = chan ID
                % reversals{:,9} = unit ID
                % reversals{:,10} = CuePolyOrder
                % reversals{:,11} = OutPolyOrder
               
%halfReversals Cell structure
%                 halfReversals{:,1} = RPcg1
%                 halfReversals{:,2} = RPcg2
%                 halfReversals{:,3} = RPog1
%                 halfReversals{:,4} = RPog2
%                 halfReversals{:,5} = linear encoding for g1 cue
%                 halfReversals{:,6} = linear encoding for g2 cue
%                 halfReversals{:,7} = linear encoding for g1 out
%                 halfReversals{:,8} = linear encoding for g2 out
%                 halfReversals{:,9} = RPg1
%                 halfReversals{:,10} = RPg2
%                 halfReversals{:,11} =RPOutg1
%                 halfReversals{:,12} = RPOutg2
%                 halfReversals{:,13} = channel
%                 halfReversals{:,14} = unit
%                 halfReversals{:,15} = ptID
%                 halfReversals{:,16} = BARTstats.bestCuePolyOrderg1;
%                 halfReversals{:,17} = BARTstats.bestOutPolyOrderg1;
%                 halfReversals{:,18} = BARTstats.bestCuePolyOrderg2;
%                 halfReversals{:,19} = BARTstats.bestOutPolyOrderg2;
                


%halfTaus Structure
%                 halfTau(:,1) = taug1
%                 halfTau(:,2) = taug2
%                 halfTau(:,3) = betaPlus
%                 halfTau(:,4) = betaMinus
%                 halfTau(:,5) = betaPlusg1
%                 halfTau(:,6) = betaMinusg1
%                 halfTau(:,7) = betaPlusg1
%                 halfTau(:,8) = betaminusg2

% default plotFlag to true
if nargin<2
    plotFlag = true;
end


%% TODO:: exception for 202117
% TODO::: change when to calculate an RP. previously this was done by
% looking at crossings of y-axis. Now need to rely on best fit polynomial!


%% loading neural and beahvioral data.
% nevList = dir(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/*.nev',ptID));
nevList = dir(sprintf('//155.100.91.44/d/Data/preProcessed/BART_units/%s/Data/*.nev',ptID));
% nevList = dir(sprintf('X:/Data/preProcessed/BART_units/%s/Data/*.nev',ptID));
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
% matFile = sprintf('~/data/BART/BART_EMU/%s/Data/%s.bartBHV.mat',ptID,ptID);
matFile = sprintf('//155.100.91.44/d/Data/preProcessed/BART_units/%s/Data/%s.bartBHV.mat',ptID,ptID);
load(matFile)
pointsEarned = [data.points];

% standard 3-D [chan, unit, timestamp (seconds)] matrix.
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];

% channel deets.
inclChans = unique(ChanUnitTimestamp(:,1));
[myLabels,~,microLabels,~,generalLabels] = microLabelsBART(ptID); % these are actually atlas labels. for neurology labels, use the first input.
inclChans(inclChans-96>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.\
nChans = length(inclChans);


%% task parameters
% % task parameters (out of order)
% alignName = 'banksandpops';
% bankTimes = trigTimes(trigs==25);
% popTimes = trigTimes(trigs==26);
% nBanks = length(bankTimes);
% nPops = length(popTimes);
% outcomeType = [ones(1,nBanks) 2*ones(1,nPops)]; % 1 = bank, 2 = pop
% nTrials = length(outcomeType);
% %1: trial start ::      [1 2 3 4 11 12 13 14] = [Y O R G Yc Oc Rc Gc]
ballonIdcs = trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14;
balloonTimes = trigTimes(ballonIdcs);
% inflateTimes = trigTimes(trigs==23);
% respTimes = trigTimes(trigs==26 | trigs==25);

% task identifiers
balloonIDs = trigs(ballonIdcs);
% isCTRL = balloonIDs>10;


%% now for a a brief trial identity validation
% cross referencing the balloonIDs from the NSPdata and the ballon IDs from
% the task matlab data
validation = false;
if validation
    A = balloonIDs;
    B = [data.trial_type]+([data.is_control]*10);
    figure
    hold on
    plot(1:length(A),A)
    plot(1:length(B),B)
end


% task parameters in chronological order..
% respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
%% TODO::: DOUBLE CHECK 
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop 

% adjusting for trial numbers
% only including complete trials; generally => excluding the last trial.
nTrials = min([length(outcomeType) length(balloonIDs)]);
% inflateTimes = inflateTimes(1:nTrials);
balloonTimes = balloonTimes(1:nTrials);
balloonIDs = balloonIDs(1:nTrials);
% isCTRL = isCTRL(1:nTrials);
pointsEarned = pointsEarned(1:nTrials);
% pointsEarned = pointsEarned(1:nTrials);
% poppedTrials = logical(trigs(trigs==25 | trigs==26)-25); % 0 = bank; 1 = pop;
% [~,sortedTrialIdcs] = sort(balloonIDs);
% 
% % duration of balloon inflation
% inflateDurations = outcomeTimes-inflateTimes;

% timing parameters.
pre = 3;
post = 5;


%% color map
% reward probability colormap
cMap(1,:) = [0.5 0.5 0.5];  % gray
cMap(2,:) = [1 0 0];        % red
cMap(3,:) = [1 0.5 0];      % orange
cMap(4,:) = [1 0.9 0];      % yellow
cMap(5,:) = [1 0 1];        % controls


%% setting up regressors for linear models.
% [1 2 3 4 11 12 13ba 14] = [Y O R Yc Oc Rc Gc]

% cumulative estimate of reward probability variable.
rewardProbability = zeros(nTrials,1);
% gray balloons
rewardProbability(balloonIDs==14) = 0;
% colored balloons
rewardProbability(balloonIDs==1) = cumsum(outcomeType(balloonIDs==1)==1)./(1:sum(balloonIDs==1))';
rewardProbability(balloonIDs==2) = cumsum(outcomeType(balloonIDs==2)==1)./(1:sum(balloonIDs==2))';
rewardProbability(balloonIDs==3) = cumsum(outcomeType(balloonIDs==3)==1)./(1:sum(balloonIDs==3))';
% colored controls
rewardProbability(balloonIDs==11 | balloonIDs==12 | balloonIDs==13) = 1;

% regressor for risk
% cumulativeRiskEstimate = rewardProbability.^2;

% reward probability categories: gray, red, orange, yellow, controls (purple)
% increasing probability of reward (0 < red < orange < yellow < 1)
rewardProbCats = zeros(nTrials,1);
rewardProbCats(balloonIDs==14) = 1;
rewardProbCats(balloonIDs==1) = 4;
rewardProbCats(balloonIDs==2) = 3;
rewardProbCats(balloonIDs==3) = 2;
rewardProbCats(balloonIDs==11 | balloonIDs==12 | balloonIDs==13) = 5;


% also setting up a regressor for risk.
riskCats = zeros(nTrials,1);
riskCats(balloonIDs>10) = 1;
riskCats(balloonIDs==1) = 2;
riskCats(balloonIDs==2) = 3;
riskCats(balloonIDs==3) = 4;

% reward magnitude regressor - running average
rewardMagnitude = zeros(nTrials,1);
rewardMagnitude(balloonIDs==14) = 0;
rewardMagnitude(balloonIDs==1) = (cumsum(pointsEarned(balloonIDs==1))./(1:sum(balloonIDs==1)))';
rewardMagnitude(balloonIDs==2) = (cumsum(pointsEarned(balloonIDs==2))./(1:sum(balloonIDs==2)))';
rewardMagnitude(balloonIDs==3) = (cumsum(pointsEarned(balloonIDs==3))./(1:sum(balloonIDs==3)))';
rewardMagnitude(balloonIDs==11) = (cumsum(pointsEarned(balloonIDs==11))./(1:sum(balloonIDs==11)))';
rewardMagnitude(balloonIDs==12) = (cumsum(pointsEarned(balloonIDs==12))./(1:sum(balloonIDs==12)))';
rewardMagnitude(balloonIDs==13) = (cumsum(pointsEarned(balloonIDs==13))./(1:sum(balloonIDs==13)))';

% reward magnitude regressor - categorical
rewardMagCats = zeros(nTrials,1);
rewardMagCats(balloonIDs==14) = 1;
rewardMagCats(balloonIDs==3 | balloonIDs==13) = 2;
rewardMagCats(balloonIDs==2 | balloonIDs==12) = 3;
rewardMagCats(balloonIDs==1 | balloonIDs==11) = 4;

% setting up variables for statistics
statTimeWin = [0.2 2.2]; % in seconds.
% [20220603]:: [0.2 1.7] is a good window.
% RTs = [data.rt];
% preRTwins = [zeros(length(RTs),1) RTs'];
% cueWin = [0.2 1.5*median(RTs)];
% statTimeWin = cueWin;


% for removing trial categories in order to examine a subset of trial
% categories (e.g., no gray balloons...)
removeIdcs = false;
if removeIdcs
    rmvIdcs = rewardProbCats==5 | rewardProbCats==1; % removing positive controls
    rewardProbability(rmvIdcs) = [];
    rewardProbCats(rmvIdcs) = [];
    outcomeType(rmvIdcs) = [];
    outcomeTimes(rmvIdcs) = [];
    balloonIDs(rmvIdcs) = [];
    balloonTimes(rmvIdcs) = [];
    isCTRL(rmvIdcs) = [];
    inflateTimes(rmvIdcs) = [];
    inflateDurations(rmvIdcs) = [];
    riskCats(rmvIdcs) = [];
    rewardMagnitude(rmvIdcs) = [];
    rewardMagCats(rmvIdcs) = [];
end

% making categorical variables categorical, (with order)
categorical(rewardProbCats);
categorical(riskCats);
categorical(rewardMagCats);

%redefining nTrials as the length of the new regressor vectors.
nTrials = length(rewardProbCats);

% Index to concatenate reversal outputs!
unitInd =1;
% looping over Channels
for ch = nChans:-1:1
    % looping over number of units in the AP data
    nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
    for un = nUnits:-1:1
        fprintf('\nprocessing and plotting for unit %d of %d on channel %d in patient %s....',un,nUnits,ch,ptID)
        
        % kernel width for firing rates.
        if strcmp(psthType,'chronux')
            kernelWidth = 100 ./1000;
        elseif strcmp(psthType,'bins')
            binWidth = 50; % in ms now???
        end
        
        % getting unit times for the current channel and unit.
        unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
        if ~isempty(unitTimes) || length(unitTimes>5)
            
            %% cue aligned spikes.
            for tt = nTrials:-1:1
                % putting the data in a structure
                cueSpikes.channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>balloonTimes(tt)-pre & unitTimes<balloonTimes(tt)+post) - repmat(balloonTimes(tt)-pre,length(unitTimes(unitTimes>balloonTimes(tt)-pre & unitTimes<balloonTimes(tt)+post)),1);
                cueSpikesCell{tt,1} = cueSpikes.channel(ch).unit(un).trial(tt).times';
            end % looping over trials
            
            % calculating psths
            if strcmp(psthType,'chronux')
                [R,t,~] = psth(cueSpikes.channel(ch).unit(un).trial, kernelWidth, 'n', [0 pre+post]);
            elseif strcmp(psthType,'bins')
                [R,~,t] = psthBins(cueSpikes.channel(ch).unit(un).trial(1).times, binWidth, 1, 1, pre+post); % CHECK THIS?!
            end
            


            % single trial cue-aligned firing rates.
            for t2 = 1:nTrials
                if isempty(cueSpikes.channel(ch).unit(un).trial(t2).times)
                    RC(t2,:) = zeros(1,length(t));
                    %               fprintf('\n     zero spikes recorded in trial %d! \n', t2)
                else
                    if strcmp(psthType,'chronux')
                        [RC(t2,:)] = psth(cueSpikes.channel(ch).unit(un).trial(t2), kernelWidth, 'n', [0 pre+post]);
                    elseif strcmp(psthType,'bins')
                        [rSmo(t2,:),RC(t2,:),t] = psthBins(cueSpikes.channel(ch).unit(un).trial(t2).times, binWidth, 1, 1, pre+post);
                    end
                end
            end 

            % looping over trials
            
            % timing
            % tSecCue = t-repmat(pre,1,length(t));
            tSecCue = linspace(-pre,post,length(t));
            
            % baseline correcting cue-aligned firing rates.
            baseline = true;
            if baseline
                % define baseline time window (in seconds)
                bWin = [-1.2 -0.2];
                RC = RC-repmat(mean(RC(:,tSecCue>bWin(1) & tSecCue<bWin(2)),2),1,size(RC,2));
            end

            % split trials into 2 groups randomly. placed here so whether baseline is true/false it gets incorporated into group 1/2
            % was getting an error because sometimes one of the groups only had 1 pop outcome in it and so when the linear model wasattempted to be fitted, there was insufficient inputs.
            verify = false;
            while verify == false
                randTrials = randperm(nTrials);
                group1Trials = randTrials(1:floor(end/2));
                group2Trials = randTrials(floor(end/2)+1:end);
                if sum(outcomeType(group1Trials)==2) >1 & sum(outcomeType(group2Trials)==2) >1
                    verify = true;
                end
            end
            RC1 = RC(group1Trials,:);
            RC2 = RC(group2Trials,:);
                      
            
            %% [20230320] risk realted firing rates. 
            
            %%%%%%%% CUE-ALIGNED FIRING VARIABLES
            %%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% added lines for 2 groups
            nonNeg = false;
            if nonNeg
                % nonegative normalization
                tmpFR = RC-repmat(mean(RC),size(RC,1),1);
                tmptsec = linspace(-pre+1,post-1,size(tmpFR,2));
                cueNoNeg = squeeze(mean(tmpFR(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),2));
                ZcueNoNeg = squeeze(mean(tmpFR(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),2)./sqrt(std(tmpFR(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),[],2)));
            	
                FRmeasure = squeeze(mean(RC(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),2));  

                tmpFRg1 = RC1-repmat(mean(RC1),size(RC1,1),1);
                tmptsecg1 = linspace(-pre+1,post-1,size(tmpFRg1,2));
                cueNoNeg = squeeze(mean(tmpFRg1(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),2));
                ZcueNoNegg1 = squeeze(mean(tmpFRg1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),2)./sqrt(std(tmpFRg1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),[],2)));

                FRmeasureg1 = squeeze(mean(RC1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),2));

                tmpFRg2 = RC2-repmat(mean(RC2),size(RC2,1),1);
                tmptsecg2 = linspace(-pre+1, post-1, size(tmpFRg2,2));
                cueNoNeg = squeeze(mean(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),2));
                ZcueNoNegg2 = squeeze(mean(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),2)./sqrt(std(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),[],2)));

                FRmeasureg2 = squeeze(mean(RC2(:,tmptsecg2>statTimeWin(2) & tmptsecg2<statTimeWin(2)),2));            

            else
                % mean firing rates in the stat window.
                cueVar = squeeze(mean(RC(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
                % z-scored mean firing rates in the stat window.
                ZcueVar = zscore(squeeze(mean(RC(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                % softmax normalization
                ScueVar = softmax(squeeze(mean(RC(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                
                
                % means of z-scored firing rates...
                tmpFR = (RC(:,tSecCue>-pre+1 & tSecCue<post-1) - mean(RC(:,tSecCue>-pre+1 & tSecCue<post-1)))./(std(RC(:,tSecCue>-pre+1 & tSecCue<post-1)));
                tmptsec = linspace(-pre+1,post-1,size(tmpFR,2));
                cueVarZ = squeeze(mean(tmpFR(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),2));
                
                % change variable names
                FRmeasure = cueVarZ;

                cueVarg1 = squeeze(mean(RC1(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
                ZcueVarg1 = zscore(squeeze(mean(RC1(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                ScueVarg1 = softmax(squeeze(mean(RC1(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                tmpFRg1 = (RC1(:,tSecCue>-pre+1 & tSecCue<post-1) - mean(RC1(:,tSecCue>-pre+1 & tSecCue<post-1)))./(std(RC1(:,tSecCue>-pre+1 & tSecCue<post-1)));
                tmptsecg1 = linspace(-pre+1,post-1,size(tmpFRg1,2));
                cueVarZg1 = squeeze(mean(tmpFRg1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),2));
                FRmeasureg1 = cueVarZg1;

                cueVarg2 = squeeze(mean(RC2(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
                ZcueVarg2 = zscore(squeeze(mean(RC2(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                ScueVarg2 = softmax(squeeze(mean(RC2(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                tmpFRg2 = (RC2(:,tSecCue>-pre+1 & tSecCue<post-1) - mean(RC2(:,tSecCue>-pre+1 & tSecCue<post-1)))./(std(RC2(:,tSecCue>-pre+1 & tSecCue<post-1)));
                tmptsecg2 = linspace(-pre+1,post-1,size(tmpFRg2,2));
                cueVarZg2 = squeeze(mean(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),2));
                FRmeasureg2 = cueVarZg2; 
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % anova from linear model below. NOTE:: can use either of the
            % pevious two variables.
            tbl = table(FRmeasure,rewardProbability,rewardProbCats,riskCats,rewardMagnitude,rewardMagCats...
                ,'VariableNames',{'normFiring','pReward','pRewardCats','riskCats','rewardMagnitude','rewardMagCats'});

            tblg1 = table(FRmeasureg1,rewardProbability(group1Trials),rewardProbCats(group1Trials),riskCats(group1Trials),rewardMagnitude(group1Trials),rewardMagCats(group1Trials)...
                ,'VariableNames',{'normFiring','pReward','pRewardCats','riskCats','rewardMagnitude','rewardMagCats'});

            tblg2 = table(FRmeasureg2,rewardProbability(group2Trials),rewardProbCats(group2Trials),riskCats(group2Trials),rewardMagnitude(group2Trials),rewardMagCats(group2Trials)...
                ,'VariableNames',{'normFiring','pReward','pRewardCats','riskCats','rewardMagnitude','rewardMagCats'});
            % excluding minimum value for trials that have many zero-FR
            % firing rates.s
            
            
            basic = true;
            if basic
                try
                    LM = fitglme(tbl,'normFiring ~ pRewardCats^2','Exclude',FRmeasure==min(FRmeasure))
                catch
                    LM = fitglme(tbl,'normFiring ~ pRewardCats^2')
                end
                try
                    LMg1 = fitglme(tblg1,'normFiring ~ pRewardCats^2','Exclude',FRmeasureg1 == min(FRmeasureg1));
                catch
                    LMg1 = fitglme(tblg1,'normFiring ~ pRewardCats^2');    
                end
                try
                    LMg2 = fitglme(tblg2,'normFiring ~ pRewardCats^2','Exclude',FRmeasureg2 == min(FRmeasureg2));
                catch
                    LMg2 = fitglme(tblg2,'normFiring ~ pRewardCats^2');
                end
                %                 keyboard
            else % will mess with the rest if we want to turn basic off.
                %% linear regression to determine probability selectivity.
                robusto = false;
                inclIntercepts = true;
                if ~robusto
                    if ~inclIntercepts
                        % 2) if you just want to exclude intercepts...
                        if sum(FRmeasure==min(FRmeasure))>=2
                            try
                                LM = fitlm(tbl,'normFiring ~ pRewardCats','Exclude',FRmeasure==min(FRmeasure),'intercept',false);
                            catch
                                LM = fitlm(tbl,'normFiring ~ pRewardCats','intercept',false);
                            end
                        else
                            LM = fitlm(tbl,'normFiring ~ pRewardCats','intercept',false);
                        end
                        
                    else
                        % 3) if you want to include intercepts.
                        if sum(FRmeasure==min(FRmeasure))>=2
                            try
                                LM = fitglme(tbl,'normFiring ~ 1 + riskCats)','Exclude',FRmeasure==min(FRmeasure));
                            catch
                                LM = fitglme(tbl,'normFiring ~ 1 + riskCats)');
                            end
                        else
                            LM = fitglme(tbl,'normFiring ~ pRewardCats*pReward + (1 | riskCats)');
                        end
                    end
                else
                    if ~inclIntercepts
                        % 2) if you just want to exclude intercepts...
                        if sum(FRmeasure==min(FRmeasure))>=2
                            try
                                LM = fitlm(tbl,'normFiring ~ pRewardCats','Exclude',FRmeasure==min(FRmeasure),'intercept',false,'RobustOpts','on');
                            catch
                                LM = fitlm(tbl,'normFiring ~ pRewardCats','intercept',false,'RobustOpts','on');
                            end
                        else
                            LM = fitlm(tbl,'normFiring ~ pRewardCats','intercept',false,'RobustOpts','on');
                        end
                        
                    else
                        % 3) if you want to include intercepts.
                        if sum(FRmeasure==min(FRmeasure))>=2
                            try
                                LM = fitlm(tbl,'normFiring ~ pRewardCats','Exclude',FRmeasure==min(FRmeasure),'RobustOpts','on');
                            catch
                                LM = fitlm(tbl,'normFiring ~ pRewardCats','RobustOpts','on');
                            end
                        else
                            LM = fitlm(tbl,'normFiring ~ pRewardCats','RobustOpts','on');
                        end
                    end
                end
            end
            
            % Saving regression data to output structure.
            rewardProbPvalsg1 = [];
            rewardProbPvalsg2 = [];

            BARTstats.rewardProbLM = LM;
            BARTstats.rewardProbANOVA = anova(LM);
            BARTstats.rewardProbRsquared = LM.Rsquared;
            rewardProbPvals = cat(1,rewardProbPvals,BARTstats.rewardProbANOVA.pValue(1));

            BARTstats.rewardProbLMg1 = LMg1;
            BARTstats.rewardProbANOVAg1 = anova(LMg1);
            BARTstats.rewardProbRsquaredg1 = LMg1.Rsquared;
            rewardProbPvalsg1 = cat(1,rewardProbPvalsg1,BARTstats.rewardProbANOVAg1.pValue(1));

            BARTstats.rewardPRobLMg2 = LMg2;
            BARTstats.rewardProbANOVAg2 = anova(LMg2);
            BARTstats.rewardProbRsquaredg2 = LMg2.Rsquared;
            rewardprobPvalsg2 = cat(1,rewardProbPvalsg2,BARTstats.rewardProbANOVAg2.pValue(1));
            
            % doing more precise curve fitting.
            xs = 1:0.5:5;
            [p1,S1,mu1] = polyfit(tbl.pRewardCats,tbl.normFiring,1);
            [yf1,d1] = polyval(p1,xs,S1,mu1); % found the yf from eval!!!!!!!
            [p2,S2,mu2] = polyfit(tbl.pRewardCats,tbl.normFiring,2);
            [yf2,d2] = polyval(p2,xs,S2,mu2);
            [p3,S3,mu3] = polyfit(tbl.pRewardCats,tbl.normFiring,3);
            [yf3,d3] = polyval(p3,xs,S3,mu3);

            [p1g1,S1g1,mu1g1] = polyfit(tblg1.pRewardCats,tblg1.normFiring,1);
            [yf1g1,d1g1] = polyval(p1g1,xs,S1g1,mu1g1);
            [p2g1,S2g1,mu2g1] = polyfit(tblg1.pRewardCats,tblg1.normFiring,2);
            [yf2g1,d2g1] = polyval(p2g1,xs,S2g1,mu2g1);
            [p3g1,S3g1,mu3g1] = polyfit(tblg1.pRewardCats,tblg1.normFiring,3);
            [yf3g1,d3g1] = polyval(p3g1,xs,S3g1,mu3g1);

            [p1g2,S1g2,mu1g2] = polyfit(tblg2.pRewardCats,tblg2.normFiring,1);
            [yf1g2,d1g2] = polyval(p1g2,xs,S1g2,mu1g2);
            [p2g2,S2g2,mu2g2] = polyfit(tblg2.pRewardCats,tblg2.normFiring,2);
            [yf2g2,d2g2] = polyval(p2g2,xs,S2g2,mu2g2);
            [p3g2,S3g2,mu3g2] = polyfit(tblg2.pRewardCats,tblg2.normFiring,3);
            [yf3g2,d3g2] = polyval(p3g2,xs,S3g2,mu3g2);
            
            % finding minima of the norm residuals. 
            [~,bestCuePoly] = min([S1.normr S2.normr S3.normr]);
            [~,bestCuePolyg1] = min([S1g1.normr S2g1.normr S3g1.normr]);
            [~,bestCuePolyg2] = min([S1g2.normr S2g2.normr S3g2.normr]);
            			
			% estimating reversal point from the best fit polynomial
			% solve for X = 0...
% 			RPs = 

            if LM.Coefficients{2,6}<=LM.Coefficients{3,6}
                bestCuePoly = 1;
            elseif LM.Coefficients{2,6}>=LM.Coefficients{3,6}
                bestCuePoly = 2;
            end

            if LMg1.Coefficients{2,6}<=LMg1.Coefficients{3,6}
                bestCuePolyg1 = 1;
            elseif LMg1.Coefficients{2,6} >= LMg1.Coefficients{3,6}
                bestCuePolyg1 = 2;
            end

            if LMg2.Coefficients{2,6}<=LMg2.Coefficients{3,6}
                bestCuePolyg2 = 1;
            elseif LMg2.Coefficients{2,6} >= LMg2.Coefficients{3,6}
                bestCuePolyg2 = 2;
            end

			% saving the poly fit for cue-aligned firing rates.
			BARTstats.polyXvals = xs;
			BARTstats.bestCuePolyOrder = bestCuePoly;
			eval(['BARTstats.bestCuePoly = yf' num2str(bestCuePoly) ';']) % it's unclear to me what the yf does. is it like %d?
	%		BARTstats.polyRP = polyRP;
            
            BARTstats.bestCuePolyOrderg1 = bestCuePolyg1;
            BARTstats.bestCuePolyOrderg2 = bestCuePolyg2;
            
            plotPoly = false;
            if plotPoly
                figure(100*un)
                subplot(2,2,1)
                hold on
                plot(tbl.pRewardCats,tbl.normFiring,'k.')
                plot(xs,yf1,'r-')
                plot(xs,yf1+2*d1,'m--',xs,yf1-2*d1,'m--')
                hold off
                xlim([0 6])
                axis square
                
                subplot(2,2,2)
                hold on
                plot(tbl.pRewardCats,tbl.normFiring,'k.')
                plot(xs,yf2,'r-')
                plot(xs,yf2+2*d2,'m--',xs,yf2-2*d2,'m--')
                hold off
                xlim([0 6])
                axis square
                
                subplot(2,2,3)
                hold on
                plot(tbl.pRewardCats,tbl.normFiring,'k.')
                plot(xs,yf3,'r-')
                plot(xs,yf3+2*d3,'m--',xs,yf3-2*d3,'m--')
                hold off
                xlim([0 6])
                axis square
                
%                 saveas()
            end

            % probability selective neurons from Tim's paper:
            %   - p<0.05 in linear regression between probability level and
            %   mean firing rate on each trial in a 0.25-1.25 s window
            %   post-cue.
            
            % just getting number of reward probability categories.
            rs = unique(rewardProbCats);
            
            
            %% plotting
            if plotFlag
                if ishandle(un); close(un); end
                figure(un)
                
                
                %% cue aligned spikes [raster and rates].
                % TODO:: color code the raster.
                subplot(6,3,1)
                %             try
                plotSpikeRasterCats(cueSpikesCell,rewardProbCats,cMap,'PlotType','vertline');
                xlim([1 pre+post-1])
                set(gca,'XTickLabel',{})
                axis xy square
                %             catch
                %                 text(0,0,'not enough trials or spikes to plot raster')
                %                 axis off
                %             end
                
                subplot(6,3,4)
                hold on
                for rr = 1:length(rs)
                    patch([tSecCue fliplr(tSecCue)],[((mean(RC(rewardProbCats==rs(rr),:)))+((std(RC(rewardProbCats==rs(rr),:)))./sqrt(sum(rewardProbCats==rs(rr)))))...
                        fliplr((mean(RC(rewardProbCats==rs(rr),:)))-((std(RC(rewardProbCats==rs(rr),:)))./sqrt(sum(rewardProbCats==rs(rr)))))]...
                        ,cMap(rs(rr),:),'facealpha',0.3,'edgecolor','none')
                    plot(tSecCue,(mean(RC(rewardProbCats==rs(rr),:)))','color',cMap(rs(rr),:))
                end
                hold off
                
                % deets
                hold off
                axis tight square
                xlim([-pre+1 post-1])
                xlabel('time relative to balloon onset (s)')
                ylabel('firing rates')
                title('cue')
            end
            
            clear R t E rSmo
            %% outcome aligned spikes.
            % loooping over trials
            for tt2 = 1:nTrials
                % putting the data in a structure
                posSpikes.channel(ch).unit(un).trial(tt2).times = unitTimes(unitTimes>outcomeTimes(tt2)-pre & unitTimes<outcomeTimes(tt2)+post) - repmat(outcomeTimes(tt2)-pre,length(unitTimes(unitTimes>outcomeTimes(tt2)-pre & unitTimes<outcomeTimes(tt2)+post)),1);
                posSpikesCell{tt2,1} = posSpikes.channel(ch).unit(un).trial(tt2).times';
            end % looping over trials
            
            % calculating psth
            if strcmp(psthType,'chronux')
                [R,t,E] = psth(posSpikes.channel(ch).unit(un).trial, kernelWidth, 'n', [0 pre+post]);
            elseif strcmp(psthType,'bins')
                [R(t2,:),rSmo(t2,:),t] = psthBins(posSpikes.channel(ch).unit(un).trial(t2).times, binWidth, 1, 1, pre+post);
            end
            
            for t2 = nTrials:-1:1
                % single trial firing rates.
                if isempty(posSpikes.channel(ch).unit(un).trial(t2).times)
                    Rp(t2,:) = zeros(1,length(t));
                    %               fprintf('\n     zero spikes recorded in trial %d! \n', t2)
                else
                    if strcmp(psthType,'chronux')
                        [Rp(t2,:)] = psth(posSpikes.channel(ch).unit(un).trial(t2), kernelWidth, 'n', [0 pre+post]);
                    elseif strcmp(psthType,'bins')
                        [rSmo(t2,:),Rp(t2,:),t] = psthBins(posSpikes.channel(ch).unit(un).trial(t2).times, binWidth, 1, 1, pre+post);
                    end
                end
            end % looping over trials
            

            
            % timing
            %             tSecCue = t-repmat(pre,1,length(t));
            tSecCue = linspace(-pre,post,length(t));
            
            if baseline
                % baseline correcting outcome-alogned data with the same window as the cue.
                Rp = Rp-repmat(mean(Rp(:,tSecCue>bWin(1) & tSecCue<bWin(2)),2),1,size(Rp,2));
            end
            
            % Splitting the outcome aligned Firing rates into 2 groups. Not going to rerandomize so that the groups stay consistent.
            Rpg1 = Rp(group1Trials,:);
            Rpg2 = Rp(group2Trials,:);
            
            %%%%%%%% OUTCOME-ALIGNED FIRING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % largely for plotting
            if nonNeg
                % nonegative normalization
                tmpFR = Rp-repmat(mean(Rp),size(Rp,1),1);
                tmptsec = linspace(-pre+1,post-1,size(tmpFR,2));
                cueNoNeg = squeeze(mean(tmpFR(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),2));
                ZcueNoNeg = squeeze(mean(tmpFR(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),2)./sqrt(std(tmpFR(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),[],2)));
                
                postFRmeasure = ZcueNoNeg;

                tmpFRg1 = Rpg1-repmat(mean(Rpg1),size(Rpg1,1),1);
                tmptsecg1 = linspace(-pre+1,post-1,size(tmpFRg1,2));
                cueNoNegg1 = squeeze(mean(tmpFRg1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),2)./sqrt(std(tmpFRg1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),[],2)));
                ZcueNoNegg1 = squeeze(mean(tmpFRg1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),2)./sqrt(std(tmpFRg1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),[],2)));
                postFRmeasureg1 = ZcueNoNegg1

                tmpFRg2 = Rpg2-repmat(mean(Rpg2),size(Rpg2,1),1);
                tmptsecg2 = linspace(-pre+1,post-1,size(tmpFRg2,2));
                cueNoNegg2 = squeeze(mean(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),2)./sqrt(std(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),[],2)));
                ZcueNoNegg2 = squeeze(mean(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),2)./sqrt(std(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),[],2)));
                postFRmeasureg2 = ZcueNoNegg2
            else
                % mean firing rates in the stat window.
                outVar = squeeze(mean(Rp(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
                % z-scored mean firing rates in the stat window.
                ZoutVar = zscore(squeeze(mean(Rp(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                % softmax normalization
                SoutVar = softmax(squeeze(mean(Rp(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));

                outVarg1 = squeeze(mean(Rpg1(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
                ZoutVarg1 = zscore(squeeze(mean(Rpg1(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                SoutVarg1 = softmax(squeeze(mean(Rpg1(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                
                outVarg2 = squeeze(mean(Rpg2(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
                ZoutVarg2 = zscore(squeeze(mean(Rpg2(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                SoutVarg2 = softmax(squeeze(mean(Rpg2(:,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2)));
                
                % means of z-scored firing rates...
                tmpFR = (Rp(:,tSecCue>-pre+1 & tSecCue<post-1) - mean(Rp(:,tSecCue>-pre+1 & tSecCue<post-1)))./(std(Rp(:,tSecCue>-pre+1 & tSecCue<post-1)));
                tmptsec = linspace(-pre+1,post-1,size(tmpFR,2));
                outVarZ = squeeze(mean(tmpFR(:,tmptsec>statTimeWin(1) & tmptsec<statTimeWin(2)),2));
                
                % change variable names
                postFRmeasure = outVarZ;

                tmpFRg1 = (Rpg1(:,tSecCue>-pre+1 & tSecCue<post-1) - mean(Rpg1(:,tSecCue>-pre+1 & tSecCue<post-1)))./(std(Rpg1(:,tSecCue>-pre+1 & tSecCue<post-1)));
                tmptsecg1 = linspace(-pre+1,post-1,size(tmpFRg1,2));
                outVarZg1 = squeeze(mean(tmpFRg1(:,tmptsecg1>statTimeWin(1) & tmptsecg1<statTimeWin(2)),2));
                postFRmeasureg1 = outVarZg1;

                tmpFRg2 = (Rpg2(:,tSecCue>-pre+1 & tSecCue<post-1) - mean(Rpg2(:,tSecCue>-pre+1 & tSecCue<post-1)))./(std(Rpg2(:,tSecCue>-pre+1 & tSecCue<post-1)));
                tmptsecg2 = linspace(-pre+1,post-1,size(tmpFRg2,2));
                outVarZg2 = squeeze(mean(tmpFRg2(:,tmptsecg2>statTimeWin(1) & tmptsecg2<statTimeWin(2)),2));
                postFRmeasureg2 = outVarZg2;
                
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            
            
            %% outcome model.
            outVar = squeeze(mean(postFRmeasure,2)); outVarg1 = squeeze(mean(postFRmeasureg1,2)); outVarg2 = squeeze(mean(postFRmeasureg2,2));
            outTbl = table(outVar,rewardProbCats,'VariableNames',{'normFiring','pRewardCats'});
            outTblg1 = table(outVarg1,rewardProbCats(group1Trials),'VariableNames',{'normFiring','pRewardCats'});
            outTblg2 = table(outVarg2,rewardProbCats(group2Trials),'VariableNames',{'normFiring','pRewardCats'});
            try
                outLM = fitglme(outTbl,'normFiring ~ pRewardCats^2','Exclude',outVar==min(outVar))

            catch
                outLM = fitglme(outTbl,'normFiring ~ pRewardCats^2')
            end
            try
                 outLMg1 = fitglme(outTblg1,'normFiring~pRewardCats^2','Exclude',outVarg1 == min(outVarg1));
            catch
                outLMg1 = fitglme(outTblg1,'normFiring ~ pRewardCats^2');
            end
            try
                outLMg2 = fitglme(outTblg2,'normFiring~pRewardCats^2','Exclude',outVarg2 == min(outVarg2));
            catch
                outLMg2 = fitglme(outTblg2,'normFiring ~ pRewardCats^2');
            end
            BARTstats.outLM = outLM;
            BARTstats.outLMg1 = outLMg1;
            BARTstats.outLMg2 = outLMg2;
            
            
            %% anova for reward
            rewardVar = squeeze(mean(Rp(outcomeType==1,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
            rewardVarg1 = squeeze(mean(Rpg1(outcomeType(group1Trials)==1,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
            rewardVarg2 = squeeze(mean(Rpg2(outcomeType(group2Trials)==1,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));

            [BARTstats.pReward,BARTstats.tblReward,BARTstats.statsReward] =...
                anova1(rewardVar,rewardProbCats(outcomeType==1),'off');
            rewardTbl = table(rewardVar,rewardProbCats(outcomeType==1),'VariableNames',{'FR','probCats'});
            BARTstats.rewardLM = fitlm(rewardTbl,'FR ~ 1 + probCats');

            [BARTstats.pRewardg1, BARTstats.tblRewardg1, BARTstats.statsRewargdg1] = ...
                anova1(rewardVarg1,rewardProbCats(outcomeType(group1Trials)==1),'off');
            rewardTblg1 = table(rewardVarg1,rewardProbCats(outcomeType(group1Trials)==1),'VariableNames',{'FR','probCats'});
            BARTstats.rewardLMg1 = fitlm(rewardTblg1,'FR~1+probCats');

            [BARTstats.pRewardg2, BARTstats.tblRewardg2, BARTstats.statsRewargdg2] = ...
                anova1(rewardVarg2,rewardProbCats(outcomeType(group2Trials)==1),'off');
            rewardTblg2 = table(rewardVarg2,rewardProbCats(outcomeType(group2Trials)==1),'VariableNames',{'FR','probCats'});
            BARTstats.rewardLMg2 = fitlm(rewardTblg2,'FR~1+probCats');
            
			% fiitting curves for outcome-aligned firing rates
			[p1,S1,mu1] = polyfit(tbl.pRewardCats,postFRmeasure,1);
            [yf1,d1] = polyval(p1,xs,S1,mu1);
            [p2,S2,mu2] = polyfit(tbl.pRewardCats,postFRmeasure,2);
            [yf2,d2] = polyval(p2,xs,S2,mu2);
            [p3,S3,mu3] = polyfit(tbl.pRewardCats,postFRmeasure,3);
            [yf3,d3] = polyval(p3,xs,S3,mu3);

            [p1g1,S1g1,mug1] = polyfit(tblg1.pRewardCats,postFRmeasureg1,1);
            [yf1g1,d1g1] = polyval(p1g1,xs,S1g1,mu1g1);
            [p2g1,S2g1,mu2g1] = polyfit(tblg1.pRewardCats,postFRmeasureg1,2);
            [yf2g1,d2g1] = polyval(p2g1,xs,S2g1,mu2g1);
            [p3g1,S3g1,mu3g1] = polyfit(tblg1.pRewardCats,postFRmeasureg1,3);
            [yf3g1,d3g1] = polyval(p3g1,xs,S3g1,mu3g1);

            [p1g2,S1g2,mug2] = polyfit(tblg2.pRewardCats,postFRmeasureg2,1);
            [yf1g2,d1g2] = polyval(p1g2,xs,S1g2,mu1g2);
            [p2g2,S2g2,mu2g2] = polyfit(tblg2.pRewardCats,postFRmeasureg2,2);
            [yf2g2,d2g2] = polyval(p2g2,xs,S2g2,mu2g2);
            [p3g2,S3g2,mu3g2] = polyfit(tblg2.pRewardCats,postFRmeasureg2,3);
            [yf3g2,d3g2] = polyval(p3g2,xs,S3g2,mu3g2);
            
            if outLM.Coefficients{2,6}<=outLM.Coefficients{3,6}
                bestOutPoly = 1;
            elseif outLM.Coefficients{2,6}>=outLM.Coefficients{3,6}
                bestOutPoly = 2;
            end

            if outLMg1.Coefficients{2,6}<=outLMg1.Coefficients{3,6}
                bestOutPolyg1 = 1;
            elseif outLMg1.Coefficients{2,6} >= outLMg1.Coefficients{3,6}
                bestOutPolyg1 = 2;
            end
            
            if outLMg2.Coefficients{2,6}<=outLMg2.Coefficients{3,6}
                bestOutPolyg2 = 1;
            elseif outLMg2.Coefficients{2,6}>=outLMg2.Coefficients{3,6}
                bestOutPolyg2 = 2;
            end
			
			% saving the poly fit for cue-aligned firing rates.
			BARTstats.bestOutPolyOrder = bestOutPoly;
			eval(['BARTstats.bestOutPoly = yf' num2str(bestOutPoly) ';'])
            BARTstats.bestOutPolyOrderg1 = bestOutPolyg1;
            BARTstats.bestOutPolyOrderg2 = bestOutPolyg2;
            
            
            
            %% plotting reward aligned raster and rates for positive outcomes.
            if plotFlag
                figure(un)
                subplot(6,3,2)
                try
                    plotSpikeRasterCats(posSpikesCell(outcomeType==1),rewardProbCats(outcomeType==1),cMap,'PlotType','vertline');
                    xlim([1 pre+post-1])
                    set(gca,'XTickLabel',{})
                    axis xy square
                catch
                    text(0,0,'not enough trials or spikes to plot raster')
                    axis off
                    
                end
                
                subplot(6,3,5)
                hold on
                for rr = 1:length(rs)
                    if ismember(rr,unique(rewardProbCats(outcomeType==1)))
						% double-check & vs | [20230317]
                        patch([tSecCue fliplr(tSecCue)],[((mean(Rp(outcomeType==1 & rewardProbCats==rs(rr),:)))+((std(Rp(outcomeType==1 & rewardProbCats==rs(rr),:)))./sqrt(sum(outcomeType==1 & rewardProbCats==rs(rr)))))...
                            fliplr((mean(Rp(outcomeType==1 & rewardProbCats==rs(rr),:)))-((std(Rp(outcomeType==1 & rewardProbCats==rs(rr),:)))./sqrt(sum(outcomeType==1 & rewardProbCats==rs(rr)))))]...
                            ,cMap(rs(rr),:),'facealpha',0.3,'edgecolor','none')
                        plot(tSecCue,(mean(Rp(outcomeType==1 & rewardProbCats==rs(rr),:)))','color',cMap(rs(rr),:))
                    end
                end
                hold off
                
                % deets
                hold off
                axis tight square
                xlim([-pre+1 post-1])
                xlabel('time relative to positive outcome (s)')
                ylabel('firing rates')
                title('positive outcomes')
            end
            
            %% anova for NO reward
            noRewardVar = squeeze(mean(Rp(outcomeType==2,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
            [BARTstats.pNoReward,BARTstats.tblNoReward,BARTstats.statsNoReward] =...
                anova1(noRewardVar,rewardProbCats(outcomeType==2),'off');
            noRewardTbl = table(noRewardVar,rewardProbCats(outcomeType==2),'VariableNames',{'FR','probCats'});
            BARTstats.noRewardLM = fitlm(noRewardTbl,'FR ~ 1 + probCats');

            noRewardVarg1 = squeeze(mean(Rpg1(outcomeType(group1Trials)==2,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
            [BARTstats.pNoRewardg1,BARTstats.tblNoRewardg1,BARTstats.statsNoRewardg1] = ...
                anova1(noRewardVarg1, rewardProbCats(outcomeType(group1Trials)==2),'off');
            noRewardTblg1 = table(noRewardVarg1,rewardProbCats(outcomeType(group1Trials)==2),'VariableNames',{'FR','probCats'});
            BARTstats.noRewardLMg1 = fitlm(noRewardTblg1,'FR ~ 1 + probCats');

            noRewardVarg2 = squeeze(mean(Rpg2(outcomeType(group2Trials)==2,tSecCue>statTimeWin(1) & tSecCue<statTimeWin(2)),2));
            [BARTstats.pNoRewardg2,BARTstats.tblNoRewardg2,BARTstats.statsNoRewardg2] = ...
                anova1(noRewardVarg2, rewardProbCats(outcomeType(group2Trials)==2),'off');
            noRewardTblg2 = table(noRewardVarg2,rewardProbCats(outcomeType(group2Trials)==2),'VariableNames',{'FR','probCats'});
            BARTstats.noRewardLMg2 = fitlm(noRewardTblg2,'FR ~ 1 + probCats');
                       
            
             %%%%%%%%%%%%%%  TAU  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %% TAU calculation
             betaPlus = abs(BARTstats.rewardLM.Coefficients{2,1});     % beta for reward linear model.
             betaMinus = abs(BARTstats.noRewardLM.Coefficients{2,1});  % beta for unrewarded linear model
             tau = betaPlus/(betaPlus + betaMinus);                    % tau calculation
             taus = cat(1,taus,tau);									  % saving taus...

             testTau = (betaPlus - betaMinus)/(betaPlus + betaMinus);
             testTaus = cat(1,testTaus,testTau);

             betaPlusg1 = abs(BARTstats.rewardLMg1.Coefficients{2,1});
             betaMinusg1 = abs(BARTstats.noRewardLMg1.Coefficients{2,1});
             taug1 = betaPlusg1/(betaPlusg1 + betaMinusg1);
             testTaug1 = (betaPlusg1-betaMinusg1)/(betaPlusg1 + betaMinusg1);

             betaPlusg2 = abs(BARTstats.rewardLMg2.Coefficients{2,1});
             betaMinusg2 = abs(BARTstats.noRewardLMg2.Coefficients{2,1});
             taug2 = betaPlusg2/(betaPlusg2 + betaMinusg2);
             testTaug2 = (betaPlusg2 - betaMinusg2)/(betaPlusg2 + betaMinusg2);

             halfTau = [taug1, taug2 betaPlus betaMinus betaPlusg1 betaMinusg1 betaPlusg2 betaMinusg2];
             halfTaus = cat(1,halfTaus,halfTau);

             halftestTau = [testTaug1 testTaug2];
             halftestTaus= cat(1,halftestTaus,halftestTau);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %% plotting reward aligned raster and rates for negative outcomes.
            if plotFlag
                figure(un)
                subplot(6,3,3)
                try
                    plotSpikeRasterCats(posSpikesCell(outcomeType==2),rewardProbCats(outcomeType==2),cMap,'PlotType','vertline');
                    xlim([1 pre+post-1])
                    set(gca,'XTickLabel',{})
                    axis xy square
                catch
                    text(0,0,'unable to plot raster...')
                    axis off
                end
                
                subplot(6,3,6)
                hold on
                for rr = 1:length(rs)
                    if ismember(rr,unique(rewardProbCats(outcomeType==2)))
                        patch([tSecCue fliplr(tSecCue)],[((mean(Rp(outcomeType==2 & rewardProbCats==rs(rr),:),1))+((std(Rp(outcomeType==2 & rewardProbCats==rs(rr),:),[],1))./sqrt(sum(outcomeType==2 & rewardProbCats==rs(rr)))))...
                            fliplr((mean(Rp(outcomeType==2 & rewardProbCats==rs(rr),:),1))-((std(Rp(outcomeType==2 & rewardProbCats==rs(rr),:),[],1))./sqrt(sum(outcomeType==2 & rewardProbCats==rs(rr)))))]...
                            ,cMap(rs(rr),:),'facealpha',0.3,'edgecolor','none')
                        plot(tSecCue,(mean(Rp(outcomeType==2 & rewardProbCats==rs(rr),:),1))','color',cMap(rs(rr),:))
                    end
                end
                hold off
                
                % deets
                hold off
                axis tight square
                xlim([-pre+1 post-1])
                xlabel('time relative to negative outcome (s)')
                ylabel('firing rates')
                title('negative outcomes')
                
                %% second column of per-unit plot:
                %   1) text data summary
                %   2) reversal point plot
                %   3) scatter plot with linear regression
                
                % data summary
                subplot(3,2,[3 5])
                hold on
                % patient ID
                text(0,1,sprintf('patient: %s',ptID))
                % which unit and where it is, anatomically.
                if ch<=8
                    text(0,0.8,sprintf('chan: %d, unit: %d; %s (%s)',inclChans(ch),un,string(microLabels{1}),generalLabels{1}))
                    unitLoc = deblank(microLabels{1});
                    if strcmp(unitLoc,'Unknown') || strcmp(unitLoc,'Left Cerebral White Matter') || strcmp(unitLoc,'Right Cerebral White Matter')
                        unitLoc = deblank(myLabels{1});
                    end
                elseif ch>8 && ch<=16
                    text(0,0.8,sprintf('chan: %d, unit: %d; %s (%s)',inclChans(ch),un,string(microLabels{2}),generalLabels{2}))
                    unitLoc = deblank(microLabels{2});
                    if strcmp(unitLoc,'Unknown') || strcmp(unitLoc,'Left Cerebral White Matter') || strcmp(unitLoc,'Right Cerebral White Matter')
                        unitLoc = deblank(myLabels{2});
                    end
                elseif ch>16 && ch<=24
                    text(0,0.8,sprintf('chan: %d, unit: %d; %s (%s)',inclChans(ch),un,string(microLabels{3}),generalLabels{3}))
                    unitLoc = deblank(microLabels{3});
                    if strcmp(unitLoc,'Unknown') || strcmp(unitLoc,'Left Cerebral White Matter') || strcmp(unitLoc,'Right Cerebral White Matter')
                        unitLoc = deblank(myLabels{3});
                    end
                elseif ch>24 && ch<=32
                    text(0,0.8,sprintf('chan: %d, unit: %d; %s (%s)',inclChans(ch),un,string(microLabels{4}),generalLabels{4}))
                    unitLoc = deblank(microLabels{4});
                    if strcmp(unitLoc,'Unknown') || strcmp(unitLoc,'Left Cerebral White Matter') || strcmp(unitLoc,'Right Cerebral White Matter')
                        unitLoc = deblank(myLabels{4});
                    end
                else
                    error('trying to loop over more channels than recorded electrodes. something went wrong...')
                end
                
                % probability linear regression results.
                try
                    % looping over CUE model coefficients and listing them below the unit location
                    for cfs = 2:size(LM.Coefficients,1)
                        text(0,1-(0.2*cfs),sprintf('%s: t(%d) = %.2f, p = %.2f',LM.Coefficients{cfs,1},LM.Coefficients{cfs,5},LM.Coefficients{cfs,4},LM.Coefficients{cfs,6}));
                    end
                    
                    text(0,0,'~~~~ cue model above this line, outcome model below ~~~~')
                    
                    % looping over OUTCOME model coefficients and listing them below the unit location
                    for cf2 = 2:size(outLM.Coefficients,1)
                        text(0,0-(0.2*cf2),sprintf('%s: t(%d) = %.2f, p = %.2f',outLM.Coefficients{cf2,1},outLM.Coefficients{cf2,5},outLM.Coefficients{cf2,4},outLM.Coefficients{cf2,6}));
                    end
                    
                catch
                    text(0,0.6,'probability selectivity could not be determined...')
                end
                
                axis off
                hold off
                ylim([-1 1])
                
            else
                if ch<=8
                    unitLoc = deblank(microLabels{1});
                    if strcmp(unitLoc,'Unknown') || strcmp(unitLoc,'Left Cerebral White Matter') || strcmp(unitLoc,'Right Cerebral White Matter')
                        unitLoc = deblank(myLabels{1});
                    end
                elseif ch>8 && ch<=16
                    unitLoc = deblank(microLabels{2});
                    if strcmp(unitLoc,'Unknown') || strcmp(unitLoc,'Left Cerebral White Matter') || strcmp(unitLoc,'Right Cerebral White Matter')
                        unitLoc = deblank(myLabels{2});
                    end
                elseif ch>16 && ch<=24
                    unitLoc = deblank(microLabels{3});
                    if strcmp(unitLoc,'Unknown') || strcmp(unitLoc,'Left Cerebral White Matter') || strcmp(unitLoc,'Right Cerebral White Matter')
                        unitLoc = deblank(myLabels{3});
                    end
                elseif ch>24 && ch<=32
                    unitLoc = deblank(microLabels{4});
                    if strcmp(unitLoc,'Unknown') || strcmp(unitLoc,'Left Cerebral White Matter') || strcmp(unitLoc,'Right Cerebral White Matter')
                        unitLoc = deblank(myLabels{4});
                    end
                else
                    error('trying to loop over more channels than recorded electrodes. something went wrong...')
                end
                
                % probability linear regression results.
                try
                    % looping over model coefficients and listing them below the unit location
                    for cfs = 2:size(LM.Coefficients,1)
                        text(0,1-(0.2*cfs),sprintf('%s: t(%d) = %.2f, p = %.2f',LM.Coefficients{cfs,1},LM.Coefficients{cfs,5},LM.Coefficients{cfs,4},LM.Coefficients{cfs,6}));
                    end
                catch
                    text(0,0.6,'probability selectivity could not be determined...')
                end
            end

            
            %% calculating reversal points
            for z = length(unique(rewardProbCats)):-1:1
                % for Cue
                rProbFRbar(z) = mean(FRmeasure(rewardProbCats==rs(z)));
                rProbFRerr(z) = std(FRmeasure(rewardProbCats==rs(z)))./sum(rewardProbCats==rs(z));
                rProbFRbarg1(z) = mean(FRmeasureg1(rewardProbCats(group1Trials)==rs(z)));
                rProbFRerrg1(z) = std(FRmeasureg1(rewardProbCats(group1Trials)==rs(z)))./sum(rewardProbCats==rs(z));
                rProbFRbarg2(z) = mean(FRmeasureg2(rewardProbCats(group2Trials)==rs(z)));
                rProbFRerrg2(z) = std(FRmeasureg2(rewardProbCats(group2Trials)==rs(z)))./sum(rewardProbCats==rs(z));
                
                % for outcome
                rProbFRbarOut(z) = mean(postFRmeasure(rewardProbCats==rs(z)));
                rProbFRerrOut(z) = std(postFRmeasure(rewardProbCats==rs(z)))./sum(rewardProbCats==rs(z));
                rProbFRbarOutg1(z) = mean(postFRmeasureg1(rewardProbCats(group1Trials)==rs(z)));
                rProbFRerrOutg1(z) = std(postFRmeasureg1(rewardProbCats(group1Trials)==rs(z)))./sum(rewardProbCats==rs(z));
                rProbFRbarOutg2(z) = mean(postFRmeasureg2(rewardProbCats(group2Trials)==rs(z)));
                rProbFRerrOutg2(z) = std(postFRmeasureg2(rewardProbCats(group2Trials)==rs(z)))./sum(rewardProbCats==rs(z));
            end
            
            % de-meaning rProbFRbar for RP calculation
            if nonNeg
                rProbFRbar = rProbFRbar-repmat(mean(rProbFRbar),1,length(rProbFRbar));
                rProbFRbarOut = rProbFRbarOut-repmat(mean(rProbFRbarOut),1,length(rProbFRbarOut));

                rProbFRbarg1 = rProbFRbarg1-repmat(mean(rProbFRbarg1),1,length(rProbFRbarg1));
                rProbFRbarg2 = rProbFRbarg2-repmat(mean(rProbFRbarg2),1,length(rProbFRbarg2));

                rProbFRbarOutg1 = rProbFRbarOutg1-repmat(mean(rProbFRbarOutg1),1,length(rProbFRbarOutg1));
                rProbFRbarOutg2 = rProbFRbarOutg2-repmat(mean(rProbFRbarOutg2),1,length(rProbFRbarOutg2));
            end         
            
            % Cue
            sF = sign(rProbFRbar);
            dFsF = diff(sF); % this helps identify times when rProbFRbar crosses from negative to positive or vice versa
            dF = diff(rProbFRbar);  % I don't think these ever get used.    
            sFdF = sign(dF);

            sFg1 = sign(rProbFRbarg1);
            dFsFg1 = diff(sFg1);
            sFg2 = sign(rProbFRbarg2);
            dFsFg2 = diff(sFg2);

            % Out
            sFOut = sign(rProbFRbarOut);
            dFsFOut = diff(sFOut);
            dFOut = diff(rProbFRbarOut);
            sFdFout = sign(dFOut);

            sFOutg1 = sign(rProbFRbarOutg1);
            dFsFOutg1 = diff(sFOutg1);
            sFOutg2 = sign(rProbFRbarOutg2);
            dFsFOutg2 = diff(sFOutg2);
            
            clear RP RPOut RPc RPo RPcg1 RPcg2 RPog1 RPog2
            RP = -10; % leaving this here in case we don't catch a scenario!
            RPOut = -10; % leaving this here to see there are times when the if/ends don't catch something.
            RPg1 = -10;
            RPg2 = -10;
            RPOutg1 = -10;
            RPOutg2 = -10;
            RPc=-10;
            RPo=-10;
            reversals{unitInd,5} = false; % for some reason we weren't hitting some conditions where this was either false or true, so initialzing as false.
            reversals{unitInd,6} = false;
            %TODO::: Need to use bestCuePoly and bestOutPoly to decide which RPs to calcualte. The problem is sometimes we need both...
           
            % Cue aligned RPs
            if LM.Coefficients{3,6} <.05
%             if sum(dFsF==2 | dFsF == -2)>1 % I've been working on how to deal with quadratic encoding so we're starting there.
                reversals{unitInd,5} = false;
                % the thought is to have several measures of possible RP values.
                % RPc(1) = max/min
                % RPc(n) = linear model of first time when RPE crosses 0
                % RPc(n_end) = final linear model of first time when RPE crosses 0
                % RPc(end) = mean of RP(n)s to compare against RPc(1)

                % initialize RPc to make sure all of the options are correctly executed
                RPc = -10*ones(1,sum(dFsF==2 | dFsF == -2)+2);

                % Looking at functions that cross 0 more than once and the best model fit was quadratic
%                 if BARTstats.bestCuePolyOrder==2 
                    pCue = polyfit(1:5,rProbFRbar,2); % fit a polynomial to the data
                    xRP = linspace(1,5,100); % create an x vector
                    yRP = polyval(pCue,xRP); % execute the fit polynomial 
                    if pCue(1)>0 % look for respective max/min
                        RPc(1) = xRP(find(yRP==min(yRP)));
                    else
                        RPc(1) = xRP(find(yRP==max(yRP)));
                    end

%                 elseif BARTstats.bestCuePolyOrder == 1 % these cross 0 more than once but are best modeled by a line
%                     nposC = find(dFsF == 2 | dFsF ==-2); % find points when 0 is crossed.
%                     pc1 = polyfit([nposC(1) nposC(end)+1],[rProbFRbar(nposC(1)) rProbFRbar(nposC(end)+1)],1); % made the decision to fit across the first/last crossings of zero
%                     xc = linspace(1,5,100);
%                     yc1 = polyval(pc1,xc);
%                     [~, idxc] = min(abs(yc1));
%                     RPc(1) = xc(idxc); % this should be when rProbFRbar crosses 0.
%                 end
                % now to determine the RP everytime 0 is crossed!
                nRPs = 2:sum(dFsF==2 | dFsF == -2)+1; 
                RPposidx = find(dFsF ==2 | dFsF == -2); % find indicies when a crossing of 0 occurs
                for k = 1:length(RPposidx)
                    p1c = polyfit([RPposidx(k) RPposidx(k)+1],[rProbFRbar(RPposidx(k)) rProbFRbar(RPposidx(k)+1)],1); %fit a polynomial between crossings!
                    xc = linspace(1,5,100);
                    yc = polyval(p1c,xc); % evaluate the polynomial
                    [~,idxc] = min(abs(yc));
                    RPc(k+1) = xc(idxc); % record the RP
                end
                RPc(end) = mean(RPc(2:(end-1))); % take an average of the crossings. We're not including max/min in this calculation. Going to be used as comparison for max/min
            end
            if LM.Coefficients{2,6} < .05 % now we get to linear encoding!
                if length(nonzeros(dFsF))  == 1 % corresponds to a 'monotonic unit'
                    reversals{unitInd,5} = true; % report that this neuron only crosses 0 once
                    RPidx = find(dFsF==2 | dFsF == -2); % find when cross occurs.
                    p = polyfit([RPidx, RPidx+1], [rProbFRbar(RPidx) rProbFRbar(RPidx+1)], 1);
                    x1 = linspace(1,5,100);
                    y1 = polyval(p,x1);
                    [~,idx] = min(abs(y1));
                    RP = x1(idx);
                else % if unit isn't monotonic, we look at whole linear regression!
                    reversals{unitInd,5} = false;
                    tstRP = -10*ones(1,sum(dFsF==2|dFsF==-2));
%                     tst = fitlm(tbl,'normFiring~pRewardCats');           
                    % then calculate reversal points using entire linear model
                    % first fitting a line
%                     p = polyfit([1 5],[rProbFRbar(1) rProbFRbar(end)],1);
%                     p = [tst.Coefficients{2,1}, tst.Coefficients{1,1}];   % this uses the whole linear regression, instead of just making a model between the first and last points
                    % then evaluating for zero.
%                     x1 = linspace(1,5,100);
%                     y1 = polyval(p,x1);
%                     [~,idx] = min(abs(y1));
%                     RP = x1(idx);
                    RPposidx = find(dFsF ==2 | dFsF == -2);
                    for k = 1:length(RPposidx)
                    p1c = polyfit([RPposidx(k) RPposidx(k)+1],[rProbFRbar(RPposidx(k)) rProbFRbar(RPposidx(k)+1)],1); %fit a polynomial between crossings!
                    xc = linspace(1,5,100);
                    yc = polyval(p1c,xc); % evaluate the polynomial
                    [~,idxc] = min(abs(yc));
                    tstRP(k) = xc(idxc);
                    end
                    RP = mean(tstRP);
                end
           end
%             keyboard
            % Cue-aligned group 1!
            if LMg1.Coefficients{3,6} < .05
%             if sum(dFsFg1==2 | dFsFg1 == -2)>1 % I've been working on how to deal with quadratic encoding so we're starting there.
                halfReversals{unitInd,5} = false;
                % the thought is to have several measures of possible RP values.
                % RPc(1) = max/min
                % RPc(n) = linear model of first time when RPE crosses 0
                % RPc(n_end) = final linear model of first time when RPE crosses 0
                % RPc(end) = mean of RP(n)s to compare against RPc(1)

                % initialize RPc to make sure all of the options are correctly executed
                RPcg1 = -10*ones(1,sum(dFsFg1==2 | dFsFg1 == -2)+2);

                % Looking at functions that cross 0 more than once and the best model fit was quadratic
%                 if BARTstats.bestCuePolyOrderg1==2 
                    pCueg1 = polyfit(1:5,rProbFRbarg1,2); % fit a polynomial to the data
                    xRPg1 = linspace(1,5,100); % create an x vector
                    yRPg1 = polyval(pCueg1,xRPg1); % execute the fit polynomial 
                    if pCueg1(1)>0 % look for respective max/min
                        RPcg1(1) = xRPg1(find(yRPg1==min(yRPg1)));
                    else
                        RPcg1(1) = xRPg1(find(yRPg1==max(yRPg1)));
                    end

%                 elseif BARTstats.bestCuePolyOrderg1 == 1 % these cross 0 more than once but are best modeled by a line
%                     nposCg1 = find(dFsFg1 == 2 | dFsFg1 ==-2); % find points when 0 is crossed.
%                     pc1g1 = polyfit([nposCg1(1) nposCg1(end)+1],[rProbFRbarg1(nposCg1(1)) rProbFRbarg1(nposCg1(end)+1)],1); % made the decision to fit across the first/last crossings of zero
%                     xcg1 = linspace(1,5,100);
%                     yc1g1 = polyval(pc1g1,xcg1);
%                     [~, idxcg1] = min(abs(yc1g1));
%                     RPcg1(1) = xcg1(idxcg1); % this should be when rProbFRbar crosses 0.
%                 end
                % now to determine the RP everytime 0 is crossed!
                nRPsg1 = 2:sum(dFsFg1==2 | dFsFg1 == -2)+1; 
                RPposidxg1 = find(dFsFg1 ==2 | dFsFg1 == -2); % find indicies when a crossing of 0 occurs
                for k = 1:length(RPposidxg1)
                    p1cg1 = polyfit([RPposidxg1(k) RPposidxg1(k)+1],[rProbFRbarg1(RPposidxg1(k)) rProbFRbarg1(RPposidxg1(k)+1)],1); %fit a polynomial between crossings!
                    xcg1 = linspace(1,5,100);
                    ycg1 = polyval(p1cg1,xcg1); % evaluate the polynomial
                    [~,idxcg1] = min(abs(ycg1));
                    RPcg1(k+1) = xcg1(idxcg1); % record the RP
                end
                RPcg1(end) = mean(RPcg1(2:(end-1))); % take an average of the crossings. We're not including max/min in this calculation. Going to be used as comparison for max/min
            end
            if LMg1.Coefficients{2,6} < .05 % now we get to linear encoding!
                 % report that this neuron only crosses 0 once
                
                if length(nonzeros(dFsFg1)) == 1
                    halfReversals{unitInd,5} = true;
%                     tstg1.Coefficients{2,1} > 0  %was nonzeros(dFsFg1) == 2 % line w/ positive slope. I.e. from neg RPE to pos.
                    RPidxg1 = find(dFsFg1==2 | dFsFg1 == -2); % find when cross occurs.
                    % then calculate reversal points
                    % first fitting a line
                    pg1 = polyfit([RPidxg1 RPidxg1+1],[rProbFRbarg1(RPidxg1) rProbFRbarg1(RPidxg1+1)],1);
                    % then evaluating for zero.
                    x1g1 = linspace(1,5,100);
                    y1g1 = polyval(pg1,x1g1);
                    [~,idxg1] = min(abs(y1g1));
                    RPg1 = x1g1(idxg1);
                else % nonzeros(dFsFg1) == -2 % line w/ negative slope. I.e. from pos RPE to neg.
                    halfReversals{unitInd,5} = false;
%                     RPidxg1 = find(dFsFg1==-2);
                    tstg1 = fitlm(tblg1,'normFiring~pRewardCats');
                     pg1 = [tstg1.Coefficients{2,1}, tstg1.Coefficients{1,1}];%polyfit([RPidxg1 RPidxg1+1],[rProbFRbarg1(RPidxg1) rProbFRbarg1(RPidxg1+1)],1);
                    % then evaluating for zero.
                    x1g1 = linspace(1,5,100);
                    y1g1 = polyval(pg1,x1g1);
                    [~,idxg1] = min(abs(y1g1));
                    RPg1 = x1g1(idxg1);
                end
            end
            % Cue aligned group 2!
            if LMg2.Coefficients{3,6} < .05
%             if sum(dFsFg2==2 | dFsFg2 == -2)>1 % I've been working on how to deal with quadratic encoding so we're starting there.
                halfReversals{unitInd,6} = false;
                % the thought is to have several measures of possible RP values.
                % RPc(1) = max/min
                % RPc(n) = linear model of first time when RPE crosses 0
                % RPc(n_end) = final linear model of first time when RPE crosses 0
                % RPc(end) = mean of RP(n)s to compare against RPc(1)

                % initialize RPc to make sure all of the options are correctly executed
                RPcg2 = -10*ones(1,sum(dFsFg2==2 | dFsFg2 == -2)+2);

                % Looking at functions that cross 0 more than once and the best model fit was quadratic
%                 if BARTstats.bestCuePolyOrderg2==2 
                    pCueg2 = polyfit(1:5,rProbFRbarg2,2); % fit a polynomial to the data
                    xRPg2 = linspace(1,5,100); % create an x vector
                    yRPg2 = polyval(pCueg2,xRPg2); % execute the fit polynomial 
                    if pCueg2(1)>0 % look for respective max/min
                        RPcg2(1) = xRPg2(find(yRPg2==min(yRPg2)));
                    else
                        RPcg2(1) = xRPg2(find(yRPg2==max(yRPg2)));
                    end

%                 elseif BARTstats.bestCuePolyOrderg2 == 1 % these cross 0 more than once but are best modeled by a line
%                     nposCg2 = find(dFsFg2 == 2 | dFsFg2 ==-2); % find points when 0 is crossed.
%                     pc1g2 = polyfit([nposCg2(1) nposCg2(end)+1],[rProbFRbarg2(nposCg2(1)) rProbFRbarg2(nposCg2(end)+1)],1); % made the decision to fit across the first/last crossings of zero
%                     xcg2 = linspace(1,5,100);
%                     yc1g2 = polyval(pc1g2,xcg2);
%                     [~, idxcg2] = min(abs(yc1g2));
%                     RPcg2(1) = xcg2(idxcg2); % this should be when rProbFRbar crosses 0.
%                 end
                % now to determine the RP everytime 0 is crossed!
                nRPsg2 = 2:sum(dFsFg2==2 | dFsFg2 == -2)+1; 
                RPposidxg2 = find(dFsFg2 ==2 | dFsFg2 == -2); % find indicies when a crossing of 0 occurs
                for k = 1:length(RPposidxg2)
                    p1cg2 = polyfit([RPposidxg2(k) RPposidxg2(k)+1],[rProbFRbarg2(RPposidxg2(k)) rProbFRbarg2(RPposidxg2(k)+1)],1); %fit a polynomial between crossings!
                    xcg2 = linspace(1,5,100);
                    ycg2 = polyval(p1cg2,xcg2); % evaluate the polynomial
                    [~,idxcg2] = min(abs(ycg2));
                    RPcg2(k+1) = xcg2(idxcg2); % record the RP
                end
                RPcg2(end) = mean(RPcg2(2:(end-1))); % take an average of the crossings. We're not including max/min in this calculation. Going to be used as comparison for max/min
            end
            if LMg2.Coefficients{2,6} < .05 % now we get to linear encoding!
                if length(nonzeros(dFsFg2)) == 1 %monotonic unit
                    RPidxg2 = find(dFsFg2==2 | dFsFg2 == -2); % find when cross occurs.
                    % then calculate reversal points
                    % first fitting a line
                    pg2 = polyfit([RPidxg2 RPidxg2+1],[rProbFRbarg2(RPidxg2) rProbFRbarg2(RPidxg2+1)],1);
                    % then evaluating for zero.
                    x1g2 = linspace(1,5,100);
                    y1g2 = polyval(pg2,x1g2);
                    [~,idxg2] = min(abs(y1g2));
                    RPg2 = x1g2(idxg2);
                else %if nonzeros(dFsFg2) == -2 % line w/ negative slope. I.e. from pos RPE to neg.
                    %                     RPidxg2 = find(dFsFg2==-2);
                    tstg2 = fitlm(tblg2,'normFiring~pRewardCats');
                    pg2 = [tstg2.Coefficients{2,1}, tstg2.Coefficients{1,1}]; %polyfit([RPidxg2 RPidxg2+1],[rProbFRbarg2(RPidxg2) rProbFRbarg2(RPidxg2+1)],1);
                    % then evaluating for zero.
                    x1g2 = linspace(1,5,100);
                    y1g2 = polyval(pg2,x1g2);
                    [~,idxg2] = min(abs(y1g2));
                    RPg2 = x1g2(idxg2);
                end
            end


            % Outcome aligned RPs
            if outLM.Coefficients{3,6} < .05
%             if sum(dFsFOut==2 | dFsFOut == -2)>1 % I've been working on how to deal with quadratic encoding so we're starting there.
                reversals{unitInd,6} = false;
                % the thought is to have several measures of possible RP
                % values.
                % RPo(1) = max/min
                % RPo(n) = linear model of first time when RPE crosses 0
                % RPo(n_end) = final linear model of first time when RPE crosses 0
                % RPo(end) = mean of RP(n)s to compare against RPc(1)

                % initialize RPo to make sure all of the options are correctly executed
                RPo = -10*ones(1,sum(dFsFOut==2 | dFsFOut == -2)+2);

%                 if BARTstats.bestOutPolyOrder == 2
                    pOut = polyfit(1:5,rProbFRbarOut,2);
                    xRPOut = linspace(1,5,100);
                    yRPOut = polyval(pOut,xRPOut);
                    if pOut(1) < 0 % report max/min
                        RPo(1) = xRPOut(find(yRPOut==max(yRPOut)));
                    else
                        RPo(1) = xRPOut(find(yRPOut == min(yRPOut)));
                    end
%                 elseif BARTstats.bestOutPolyOrder == 1
%                     npos = find(dFsFOut == 2 | dFsFOut ==-2);
%                     po1 = polyfit([npos(1) npos(end)+1],[rProbFRbarOut(npos(1)) rProbFRbarOut(npos(end)+1)],1);
%                     xo = linspace(1,5,100);
%                     yo1 = polyval(po1,xo);
%                     [~, idxo] = min(abs(yo1));
%                     RPo(1) = xo(idxo);
%                 end
                nRPos = 2:sum(dFsFOut==2 | dFsFOut == -2)+1;
                RPoposidx = find(dFsFOut ==2 | dFsFOut == -2);
                for k = 1:length(RPoposidx)
                    p1o = polyfit([RPoposidx(k) RPoposidx(k)+1],[rProbFRbarOut(RPoposidx(k)) rProbFRbarOut(RPoposidx(k)+1)],1);
                    xo = linspace(1,5,100);
                    yo = polyval(p1o,xo);
                    [~,idxo] = min(abs(yo));
                    RPo(k+1) = xo(idxo);
                end
                RPo(end) = mean(RPo(2:(end-1)));
            end
            if outLM.Coefficients{2,6} < .05 % now we get to linear encoding!
                if length(nonzeros(dFsFOut)) == 1 % linear, monotonic neuron
                    reversals{unitInd,6} = true;
                    RPoidx = find(dFsFOut==2 | dFsFOut == -2);
                    % then calculate reversal points
                    % first fitting a line
                    po = polyfit([RPoidx RPoidx+1],[rProbFRbarOut(RPoidx) rProbFRbarOut(RPoidx+1)],1);
                    % then evaluating for zero.
                    xo = linspace(1,5,100);
                    yo = polyval(po,xo);
                    [~,idxo] = min(abs(yo));
                    RPOut = xo(idxo);
                else %if nonzeros(dFsFOut) == -2 % line w/ negative slope. I.e. from pos RPE to neg.
%                     RPoidx = find(dFsFOut==-2)
                    reversals{unitInd,6} = false;
%                     outtst = fitlm(outTbl,'normFiring~pRewardCats');
%                     po = [outtst.Coefficients{2,1}, outtst.Coefficients{1,1}];%polyfit([RPoidx RPoidx+1],[rProbFRbarOut(RPoidx) rProbFRbarOut(RPoidx+1)],1);
                    % then evaluating for zero.
%                     xo = linspace(1,5,100);
%                     yo = polyval(po,xo);
%                     [~,idxo] = min(abs(yo));
                    tstRPOut = -10*ones(1,sum(dFsFOut==2 | dFsFOut==-2));
                    RPposidxOut = find(dFsFOut == 2 | dFsFOut == -2);
                    for k = 1:length(RPposidxOut)
                        po = polyfit([RPposidxOut(k) RPposidxOut(k)+1],[rProbFRbarOut(RPposidxOut(k)) rProbFRbarOut(RPposidxOut(k)+1)],1);
                        xo = linspace(1,5,100);
                        yo = polyval(po,xo);
                        [~,idxo] = min(abs(yo));
                        tstRPOut(k) = xo(idxo);
                    end
                    RPOut = mean(tstRPOut) ;
                end
            end
            
            % Outcome aligned RPs -- group 1!
            if outLMg1.Coefficients{3,6} < .05
% if sum(dFsFOutg1==2 | dFsFOutg1 == -2)>1 % I've been working on how to deal with quadratic encoding so we're starting there.
                halfReversals{unitInd,7} = false;
                
                % initialize RPo to make sure all of the options are correctly executed
                RPog1 = -10*ones(1,sum(dFsFOutg1==2 | dFsFOutg1 == -2)+2);

%                 if BARTstats.bestOutPolyOrderg1 == 2
                    pOutg1 = polyfit(1:5,rProbFRbarOutg1,2);
                    xRPOutg1 = linspace(1,5,100);
                    yRPOutg1 = polyval(pOutg1,xRPOutg1);
                    if pOutg1(1) < 0 % report max/min
                        RPog1(1) = xRPOutg1(find(yRPOutg1==max(yRPOutg1)));
                    else
                        RPog1(1) = xRPOutg1(find(yRPOutg1 == min(yRPOutg1)));
                    end
%                 elseif BARTstats.bestOutPolyOrderg1 == 1
%                     nposg1 = find(dFsFOutg1 == 2 | dFsFOutg1 ==-2);
%                     po1g1 = polyfit([nposg1(1) nposg1(end)+1],[rProbFRbarOutg1(nposg1(1)) rProbFRbarOutg1(nposg1(end)+1)],1);
%                     xog1 = linspace(1,5,100);
%                     yo1g1 = polyval(po1g1,xog1);
%                     [~, idxog1] = min(abs(yo1g1));
%                     RPog1(1) = xog1(idxog1);
%                 end
                nRPosg1 = 2:sum(dFsFOutg1==2 | dFsFOutg1 == -2)+1;
                RPoposidxg1 = find(dFsFOutg1 ==2 | dFsFOutg1 == -2);
                for k = 1:length(RPoposidxg1)
                    p1og1 = polyfit([RPoposidxg1(k) RPoposidxg1(k)+1],[rProbFRbarOutg1(RPoposidxg1(k)) rProbFRbarOutg1(RPoposidxg1(k)+1)],1);
                    xog1 = linspace(1,5,100);
                    yog1 = polyval(p1og1,xog1);
                    [~,idxog1] = min(abs(yog1));
                    RPog1(k+1) = xog1(idxog1);
                end
                RPog1(end) = mean(RPog1(2:(end-1)));
            end
            if outLMg1.Coefficients{2,6} <.05 % now we get to linear encoding!
                if length(nonzeros(dFsFOutg1)) == 1 % linear monotonic
                    halfReversals{unitInd,7} = true;
                    RPoidxg1 = find(dFsFOutg1==2| dFsFOutg1 == -2);
                    % then calculate reversal points
                    % first fitting a line
                    pog1 = polyfit([RPoidxg1 RPoidxg1+1],[rProbFRbarOutg1(RPoidxg1) rProbFRbarOutg1(RPoidxg1+1)],1);
                    % then evaluating for zero.
                    xog1 = linspace(1,5,100);
                    yog1 = polyval(pog1,xog1);
                    [~,idxog1] = min(abs(yog1));
                    RPOutg1 = xog1(idxog1);
                else %if nonzeros(dFsFOutg1) == -2 % line w/ negative slope. I.e. from pos RPE to neg.
%                     RPoidxg1 = find(dFsFOutg1==-2)
                    halfReversals{unitInd,7} = false;
                    outtstg1 = fitlm(outTblg1,'normFiring~pRewardCats');
                    pog1 = [outtstg1.Coefficients{2,1}, outtstg1.Coefficients{1,1}];%polyfit([RPoidxg1 RPoidxg1+1],[rProbFRbarOutg1(RPoidxg1) rProbFRbarOutg1(RPoidxg1+1)],1);
                    % then evaluating for zero.
                    xog1 = linspace(1,5,100);
                    yog1 = polyval(pog1,xog1);
                    [~,idxog1] = min(abs(yog1));
                    RPOutg1 = xog1(idxog1);
                end
            end

            % Outcome aligned RPs -- group 2!
            if outLMg2.Coefficients{3,6} < .05
%                   if sum(dFsFOutg2==2 | dFsFOutg2 == -2)>1 % I've been working on how to deal with quadratic encoding so we're starting there.
                halfReversals{unitInd,8} = false;
                
                % initialize RPo to make sure all of the options are correctly executed
                RPog2 = -10*ones(1,sum(dFsFOutg2==2 | dFsFOutg2 == -2)+2);

%                 if BARTstats.bestOutPolyOrderg2 == 2
%                     pOutg2 = polyfit(1:5,rProbFRbarOutg2,2);
%                     xRPOutg2 = linspace(1,5,100);
%                     yRPOutg2 = polyval(pOutg2,xRPOutg2);
%                     if pOutg2(1) < 0 % report max/min
%                         RPog2(1) = xRPOutg2(find(yRPOutg2==max(yRPOutg2)));
%                     else
%                         RPog2(1) = xRPOutg2(find(yRPOutg2 == min(yRPOutg2)));
%                     end
%                 elseif BARTstats.bestOutPolyOrderg2 == 1
%                     nposg2 = find(dFsFOutg2 == 2 | dFsFOutg2 ==-2);
%                     po1g2 = polyfit([nposg2(1) nposg2(end)+1],[rProbFRbarOutg2(nposg2(1)) rProbFRbarOutg2(nposg2(end)+1)],1);
%                     xog2 = linspace(1,5,100);
%                     yo1g2 = polyval(po1g2,xog2);
%                     [~, idxog2] = min(abs(yo1g2));
%                     RPog2(1) = xog2(idxog2);
%                 end
                nRPosg2 = 2:sum(dFsFOutg2==2 | dFsFOutg2 == -2)+1;
                RPoposidxg2 = find(dFsFOutg2 ==2 | dFsFOutg2 == -2);
                for k = 1:length(RPoposidxg2)
                    p1og2 = polyfit([RPoposidxg2(k) RPoposidxg2(k)+1],[rProbFRbarOutg2(RPoposidxg2(k)) rProbFRbarOutg2(RPoposidxg2(k)+1)],1);
                    xog2 = linspace(1,5,100);
                    yog2 = polyval(p1og2,xog2);
                    [~,idxog2] = min(abs(yog2));
                    RPog2(k+1) = xog2(idxog2);
                end
                RPog2(end) = mean(RPog2(2:(end-1)));
            end
            if outLMg2.Coefficients{2,6} < .05 % now we get to linear encoding!
                if length(nonzeros(dFsFOutg2)) == 1 % linear monotonic
                    halfReversals{unitInd,8} = true;
                    RPoidxg2 = find(dFsFOutg2==2);
                    % then calculate reversal points
                    % first fitting a line
                    pog2 = polyfit([RPoidxg2 RPoidxg2+1],[rProbFRbarOutg2(RPoidxg2) rProbFRbarOutg2(RPoidxg2+1)],1);
                    % then evaluating for zero.
                    xog2 = linspace(1,5,100);
                    yog2 = polyval(pog2,xog2);
                    [~,idxog2] = min(abs(yog2));
                    RPOutg2 = xog2(idxog2);
                else %if nonzeros(dFsFOutg2) == -2 % line w/ negative slope. I.e. from pos RPE to neg.
%                     RPoidxg2 = find(dFsFOutg2==-2)
                    halfReversals{unitInd,8} = false;
                    outtstg2 = fitlm(outTblg2,'normFiring~pRewardCats');
                    pog2 = [outtstg2.Coefficients{2,1}, outtstg2.Coefficients{1,1}];%polyfit([RPoidxg2 RPoidxg2+1],[rProbFRbarOutg2(RPoidxg2) rProbFRbarOutg2(RPoidxg2+1)],1);
                    % then evaluating for zero.
                    xog2 = linspace(1,5,100);
                    yog2 = polyval(pog2,xog2);
                    [~,idxog2] = min(abs(yog2));
                    RPOutg2 = xog2(idxog2);
                end
            end
            %% saving statistics structure.
            loopVar = length(taus);
            BARTstats.patientID = ptID;
            BARTstats.timeWin = statTimeWin;
            BARTstats.unitLocation = unitLoc;
            BARTstats.group1Trials = group1Trials;
            BARTstats.group2Trials = group2Trials;
            STATS(loopVar).stats = BARTstats;


            %             %% specifying output variable
            RPs = [RP, RPOut];

            if exist("RPc",'var')
                reversals{unitInd,1} = RPc;
            end
            if exist("RPo",'var')
                reversals{unitInd,2} = RPo;
            end
            if exist("RP",'var')
                reversals{unitInd,3} = RP;
            end
            if exist("RPOut",'var')
                reversals{unitInd,4} = RPOut;
            end
            reversals{unitInd,8}= ch;
            reversals{unitInd,9}=un;
            reversals{unitInd,7} = ptID; %str2double(ptID);
            reversals{unitInd,10} = BARTstats.bestCuePolyOrder;
            reversals{unitInd,11} = BARTstats.bestOutPolyOrder;

            if exist("RPcg1",'var')
                halfReversals{unitInd,1} = RPcg1;
            end
            if exist("RPcg2",'var')
                halfReversals{unitInd,2} = RPcg2;
            end
            if exist("RPog1",'var')
                halfReversals{unitInd,3} = RPog1;
            end
            if exist("RPog2",'var')
                halfReversals{unitInd,4} = RPog2;
            end
            if exist("RPg1",'var')
                halfReversals{unitInd,9} = RPg1;
            end
            if exist("RPg2",'var')
                halfReversals{unitInd,10} = RPg2;
            end
            if exist("RPOutg1",'var')
                halfReversals{unitInd,11} = RPOutg1;
            end
            if exist("RPoutg2",'var')
                halfReversals{unitInd,12} = RPOutg2;
            end
            halfReversals{unitInd,13}= ch;
            halfReversals{unitInd,14}=un;
            halfReversals{unitInd,15} = str2double(ptID);
            halfReversals{unitInd,16} = BARTstats.bestCuePolyOrderg1;
            halfReversals{unitInd,17} = BARTstats.bestOutPolyOrderg1;
            halfReversals{unitInd,18} = BARTstats.bestCuePolyOrderg2;
            halfReversals{unitInd,19} = BARTstats.bestOutPolyOrderg2;

            % plotting the reversal point data
            if plotFlag
                markerSize = 15;
                subplot(3,2,4)
                hold on
%                 BOXPLOT = false;
%                 if BOXPLOT
%                     for rp = 1:length(rs)
%                         betterBoxplot(rs(rp),FRmeasure(rewardProbCats==rs(rp)),cMap(rs(rp),:),3,'.')
%                         betterBoxplot(rs(rp),postFRmeasure(rewardProbCats==rs(rp)),cMap(rs(rp),:),3,'.')
%                     end
%                     
%                     xlim([min(unique(rewardProbCats))-1  max(unique(rewardProbCats))+1])
%                     
%                     % plotting mean line (assuming data is de-meaned or z-scored.
%                     line([0 5],[0 0],'color',[0.3 0.3 0.3],'linestyle',':')
%                     line(rs,rProbFRbar,'color','k','linestyle','-')
%                     line(rs,rProbFRbarOut,'color','k','linestyle','--')
%                     
%                 else

                    % plotting polyfit from abov
                    plot(xs,BARTstats.bestCuePoly,'color',rgb('darkgray'),'linewidth',2)
                    plot(xs,BARTstats.bestOutPoly,'color',rgb('black'),'linewidth',2)
                    text(0,0.01,'cue','color',rgb('darkgray'),'fontweight','bold')
                    text(0,-0.01,'outcome','color',rgb('black'),'fontweight','bold')

                    % plotting errors
                    for rp = 1:length(rs)
                        line([rs(rp) rs(rp)],[rProbFRbar(rp)-rProbFRerr(rp) rProbFRbar(rp)+rProbFRerr(rp)],'color',cMap(rs(rp),:),'linewidth',6)
                        line([rs(rp) rs(rp)],[rProbFRbarOut(rp)-rProbFRerrOut(rp) rProbFRbarOut(rp)+rProbFRerrOut(rp)],'color',cMap(rs(rp),:),'linewidth',6)
                    end
                    % plotting mean line (assuming data is de-meaned or z-scored.
                    line([0 5],[0 0],'color',[0.3 0.3 0.3],'linestyle',':')
                    line(rs,rProbFRbar,'color','k','linestyle',':')
                    line(rs,rProbFRbarOut,'color','k','linestyle',':')
                    
                    
%                 end
                hold off
                
                xlabel('p(reward) categories')
                ylabel('normalized firing rates')
%                 title(sprintf('reversal point = %.2f',RP))
                % maybe replace reversal point with local max/min
                axis square tight
                xlim([0 6])
                
                
                halfMaximize(un,'page')
                
                %             suptitle('choice prob coding -------------- reward prob coding -------------- no reward prob coding')
                
                % saving figures
%                 saveDir = sprintf('~/Dropbox/BART_firing_RISK/');
                saveDir = sprintf('//155.100.91.44/d/Data/Alex/BART/Units/BART_firing/Risk');
%                 saveas(un,fullfile(saveDir,sprintf('pt%s_ch%d_un%d_firingRate_riskModel_%dcats.pdf',ptID,ch,un,length(rs))))
                saveas(un,fullfile(saveDir,sprintf('pt%s_ch%d_un%d_firingRate_riskModel_%dcats.pdf',ptID,ch,un,length(rs))))
                close(un)
                
            end





            unitInd = unitInd+1;
        end % if no spikes recorded for this unit.
    end % looping over units
end % looping over channels




























