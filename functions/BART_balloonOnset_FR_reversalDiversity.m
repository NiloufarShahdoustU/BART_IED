function [] = BART_balloonOnset_FR_reversalDiversity(ptID)   % BARTstats
% BART_BALLOONONSET_FR_RISK analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_balloonOnset_FR_risk(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20181005
% edited: TAP20210721

% input args
% ptID = '202118'; % This line for debugging only

%% loading neural and beahvioral data.
% nevList = dir(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/*.nev',ptID));
nevList = dir(sprintf('D:/Data/preProcessed/BART_units/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
% [trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;

% loading behavioral matFile
% matFile = sprintf('~/data/BART/BART_EMU/%s/Data/%s.bartBHV.mat',ptID,ptID);
matFile = sprintf('D:/Data/preProcessed/BART_units/%s/Data/%s.bartBHF.mat',ptID,ptID);
load(matFile)
pointsEarned = [data.points];

% standard 3-D [chan, unit, timestamp (seconds)] matrix.
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];

% channel deets.
inclChans = unique(ChanUnitTimestamp(:,1));
microLabels = microLabelsBART(ptID);
inclChans(inclChans-96>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.\
nChans = length(inclChans);


%% task parameters in chronological order..
% There aren't any trigs that == 4
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
respTimes = trigTimes(trigs==26 | trigs==25);

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
pointsEarned = pointsEarned(1:nTrials);
poppedTrials = logical(trigs(trigs==25 | trigs==26)-25); % 0 = bank; 1 = pop;
[~,sortedTrialIdcs] = sort(balloonIDs);

% duration of balloon inflation
inflateDurations = outcomeTimes-inflateTimes;

% timing parameters.
pre = 2;
post = 3;

% firign rate timing parameters
binWidth = 50;
Fspikes = 1000;


%% [20210326]:: reframing this as in Tim's study:
%  1) ANOVA for reward probability aligned on choice.
%  2) ANOVA for reward probability aligned on reward.
%  3) ANOVA for reward probability aligned on no reward.
%
%   reward probability designations:
%       Gray = 0
%       Red, orange, yellow > 0 & < 1
%       Color controls = 1
%
%   So 5 categories corresponding to balloon colors.


%% color map
% reward probability colormap
cMap(1,:) = [0.5 0.5 0.5];  % gray
cMap(2,:) = [1 0 0];        % red
cMap(4,:) = [1 0.5 0];      % orange
cMap(3,:) = [1 0.9 0];      % yellow
cMap(5,:) = [1 0 1];        % controls


%% setting up regressors for linear models.
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
rewardProbCats(balloonIDs==14) = 1;
rewardProbCats(balloonIDs==1) = 2;
rewardProbCats(balloonIDs==2) = 3;
rewardProbCats(balloonIDs==3) = 4;
rewardProbCats(balloonIDs==11 | balloonIDs==12 | balloonIDs==13) = 5;


% looping over Channels
for ch = 1:nChans
    % looping over number of units in the AP data
    nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
    for un = 1:nUnits
        fprintf('\nprocessing and plotting for unit %d of %d',un,nUnits)
        
        % getting unit times for the current channel and unit.
        unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
        if ~isempty(unitTimes)
            
            
            %% cue aligned spikes.
            % loooping over trials
            for tt = 1:nTrials
                % putting the data in a structure
                spikes.channel(ch).unit(un).trial(tt).cueTimes = unitTimes(unitTimes>balloonTimes(tt)-pre & unitTimes<balloonTimes(tt)+post) - repmat(balloonTimes(tt)-pre,length(unitTimes(unitTimes>balloonTimes(tt)-pre & unitTimes<balloonTimes(tt)+post)),1);
                
                Rst(tt,:) = psthBins(round(spikes.channel(ch).unit(un).trial(tt).cueTimes*Fspikes), binWidth, Fspikes, 1, (pre+post)*Fspikes);
            end % looping over trials
            % timing
            tSecCue = linspace(-pre, post, size(Rst,2));
            
            % setting up variables for statistics
            statTimeWin = [0.25 1.25]; % in seconds.
            
            
            %% anova for cue aligned probability.
            cueVar = squeeze(mean(Rst(:,tSecCue>pre+statTimeWin(1) & tSecCue<pre+statTimeWin(2)),2));
            [BARTstats(ch,un).pCue,BARTstats(ch,un).tblCue,BARTstats(ch,un).statsCue] = anova1(cueVar,rewardProbCats);
            
            
            %% plotting
            if ishandle(un); close(un); end
            figure(un)
            
            subplot(3,2,4)
            trials = (1:nTrials)';
            hold on
            %             scatter(trials(balloonIDs==14),rewardProbability(rewardProbCats==1),10,cMap(1,:),'filled');
            %             scatter(trials(balloonIDs==1),rewardProbability(rewardProbCats==4),10,cMap(2,:),'filled');
            %             scatter(trials(balloonIDs==2),rewardProbability(rewardProbCats==3),10,cMap(3,:),'filled');
            %             scatter(trials(balloonIDs==3),rewardProbability(rewardProbCats==2),10,cMap(4,:),'filled');
            %             scatter(trials(balloonIDs==11 | balloonIDs==12 | balloonIDs==13),rewardProbability(rewardProbCats==5),10,cMap(5,:),'filled');
            hold off
            axis tight
            
            xlabel('trials')
            ylabel('reward probability')
            
            
            %% cue aligned spikes.
            subplot(3,2,1)
            hold on
            for rr = 1:5
                patch([tSecCue fliplr(tSecCue)],[((mean(Rst(rewardProbCats==rr,:)))+((std(Rst(rewardProbCats==rr,:)))./sum(rewardProbCats==rr)))...
                    fliplr((mean(Rst(rewardProbCats==rr,:)))-((std(Rst(rewardProbCats==rr,:)))./sum(rewardProbCats==rr)))]...
                    ,cMap(rr,:),'facealpha',0.3,'edgecolor','none')
                plot(tSecCue,(mean(Rst(rewardProbCats==rr,:)))','color',cMap(rr,:))
            end
            hold off
            
            % deets
            hold off
            axis tight square
            xlim([-pre+1 post-1])
            xlabel('time relative to balloon onset (s)')
            ylabel('firing rates')
            title('cue')
            
            
            %% outcome aligned spikes.
            % loooping over trials
            for tt2 = 1:nTrials
                % putting the data in a structure
                spikes.channel(ch).unit(un).trial(tt2).outTimes = unitTimes(unitTimes>outcomeTimes(tt2)-pre & unitTimes<outcomeTimes(tt2)+post) - repmat(outcomeTimes(tt2)-pre,length(unitTimes(unitTimes>outcomeTimes(tt2)-pre & unitTimes<outcomeTimes(tt2)+post)),1);
                
                Rst(tt2,:) = psthBins(round(spikes.channel(ch).unit(un).trial(tt2).outTimes*Fspikes), binWidth, Fspikes, 1, (pre+post)*Fspikes);
            end % looping over trials
            % timing
            tSecOut = linspace(-pre, post, size(Rst,2));
            
            
            %% anova for reward
            rewardVar = squeeze(mean(Rst(outcomeType==1,tSecCue>pre+statTimeWin(1) & tSecCue<statTimeWin(2)),2));
            [BARTstats(ch,un).pReward,BARTstats(ch,un).tblReward,BARTstats(ch,un).statsReward] = anova1(rewardVar,rewardProbCats(outcomeType==1));
            
            
            %% plotting reward aligned
            figure(un)
            subplot(3,2,3)
            hold on
            for rr = 1:5
                patch([tSecOut fliplr(tSecOut)],[((mean(Rst(outcomeType==1 | rewardProbCats==rr,:)))+((std(Rst(outcomeType==1 | rewardProbCats==rr,:)))./sum(outcomeType==1 | rewardProbCats==rr)))...
                    fliplr((mean(Rst(outcomeType==1 | rewardProbCats==rr,:)))-((std(Rst(outcomeType==1 | rewardProbCats==rr,:)))./sum(outcomeType==1 | rewardProbCats==rr)))]...
                    ,cMap(rr,:),'facealpha',0.3,'edgecolor','none')
                plot(tSecOut,(mean(Rst(outcomeType==1 | rewardProbCats==rr,:)))','color',cMap(rr,:))
            end
            hold off
            
            % deets
            hold off
            axis tight square
            xlim([-pre+1 post-1])
            xlabel('time relative to positive outcome (s)')
            ylabel('firing rates')
            title('positive outcomes')
            
            
            %% anova for NO reward
            noRewardVar = squeeze(mean(Rst(outcomeType==2,tSecCue>pre+statTimeWin(1) & tSecCue<pre+statTimeWin(2)),2));
            [BARTstats(ch,un).pNoReward,BARTstats(ch,un).tblNoReward,BARTstats(ch,un).statsNoReward] = anova1(noRewardVar,rewardProbCats(outcomeType==2));
            
            
            %% plotting reward aligned
            figure(un)
            subplot(3,2,5)
            hold on
            for rr = 1:5
                patch([tSecCue fliplr(tSecCue)],[((mean(Rst(outcomeType==2 | rewardProbCats==rr,:)))+((std(Rst(outcomeType==2 | rewardProbCats==rr,:)))./sum(outcomeType==2 | rewardProbCats==rr)))...
                    fliplr((mean(Rst(outcomeType==2 | rewardProbCats==rr,:)))-((std(Rst(outcomeType==2 | rewardProbCats==rr,:)))./sum(outcomeType==2 | rewardProbCats==rr)))]...
                    ,cMap(rr,:),'facealpha',0.3,'edgecolor','none')
                plot(tSecCue,(mean(Rst(outcomeType==2 | rewardProbCats==rr,:)))','color',cMap(rr,:))
            end
            hold off
            
            % deets
            hold off
            axis tight square
            xlim([-pre+1 post-1])
            xlabel('time relative to negative outcome (s)')
            ylabel('firing rates')
            title('negative outcomes')
            
            
            %% plotting a stats summary
            subplot(3,2,6)
            hold on
            text(0,0.4,sprintf('patient: %s',ptID))
            if ch<=8
                text(0,0.2,sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{1}))
            elseif ch>8 & ch<=16
                text(0,0.2,sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{2}))
            elseif ch>16 & ch<=24
                text(0,0.2,sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{3}))
            elseif ch>24 & ch<=32
                text(0,0.2,sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,'4th BF???'))
            else
                error('trying to loop over more channels than recorded electrodes. something went wrong...')
            end
            % cue anova
            text(0,0,sprintf('Cue ANOVA: p = %.2f, t(%d) = %.2f',...
                BARTstats(ch,un).pCue,BARTstats(ch,un).tblCue{2,3},BARTstats(ch,un).tblCue{2,5}))
            % reward anova
            text(0,-0.2,sprintf('Reward ANOVA: p = %.2f, t(%d) = %.2f',...
                BARTstats(ch,un).pReward,BARTstats(ch,un).tblReward{2,3},BARTstats(ch,un).tblReward{2,5}))
            % NO reward anova
            text(0,-0.4,sprintf('No Reward ANOVA: p = %.2f, t(%d) = %.2f',...
                BARTstats(ch,un).pNoReward,BARTstats(ch,un).tblNoReward{2,3},BARTstats(ch,un).tblNoReward{2,5}))
            axis off
            hold off
            
            
            halfMaximize(un,'left')
            
            suptitle('choice prob coding -------------- reward prob coding -------------- no reward prob coding')
            
            % saving figures
            saveDir = sprintf('~/Dropbox/BART_firing/');
            saveas(un,fullfile(saveDir,sprintf('pt%s_ch%d_un%d_firingRate_reversalDiversity.pdf',ptID,ch,un)))
            close(un)
        end % if no spikes recorded for this unit.
    end % looping over units
end % looping over channels













