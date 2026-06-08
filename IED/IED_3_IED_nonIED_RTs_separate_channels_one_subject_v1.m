% finding the IEDs RTs that are significantly different from Non-IED RTs
% this is without outliers function at the very first step
% author: Nill


clear;
clc;
close all;


%% loading behavioral data

ptID = '202314';
matFile = ['\\155.100.91.44\d\Data\Nill\BART\bhvStruct\' ptID '.bhvStruct.mat'];
load(matFile)


% retrieving reaction times from bhvStruct
ReactionTimes = bhvStruct.allRTs;

%remove reaction times that are more than 10.
ReactTimeThreshold = 10;
OutlierIndices = ReactionTimes >= ReactTimeThreshold;

ReactionTimesFiltered = ReactionTimes(~OutlierIndices);


% 
% figure;
% plot(ReactionTimesFiltered);
% xlabel('Trials', FontSize= 18);
% ylabel('RTs (sec)', FontSize= 18);


%% loading neural and event data

% nevList = dir(['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);
nevList = dir(['\\155.100.91.44\d\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);

if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesBART(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));

% data parameters
selectedChans = find(isECoG);
nChans = length(selectedChans);
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;

% resampling LFP at Fnew sampling frequency


notchFilter = true;
Fnew = 500;
for ch = 1:nChans
    if notchFilter
        [b1,a1] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        tmp(ch,:) = filtfilt(b1,a1,resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs));
        
        [b2,a2] = iirnotch(120/(Fnew/2),(120/(Fnew/2))/50);
        data2K(ch,:) = filtfilt(b2,a2,tmp(ch,:));
    else
        data2K(ch,:) = resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs);
    end
end
clear NSX NEV
Fs = Fnew;

% common average re-referencing.
%   data2K = remove1stPC(data2K);

% timing parameters.
pre = 2;
post = 4;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% TF parameters
params.fpass = [1 250]; % pick a value that's 50 higher for wavelets
params.Fs = Fs;
params.dBconversion = false;
params.normalized = false;
params.theoreticalNorm = false;
params.baseline = true;

% picking pre-bird baseline period
if params.baseline
    baselineType = 'preCue';
    switch baselineType
        case {'preCue'}
            bP = [-2.5 -1.5];
        case {'preTask'}
            % normalizing based on the mean spectrum from a pre-task baseline epoch.
            secsPreTask = 50;
            if (ppData.Event.trigTimes(1)/3e4)<secsPreTask
                fprintf('only %d seconds before first trigger.',ppData.Event.trigTimes(1)./3e4);
                % TODO:: then do spectrum
            else
                % TODO:: do spectrum.
            end
    end
end

% task parameters in chronological order..
% There aren't any trigs that == 4
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
if length(balloonTimes)>length(inflateTimes)
    balloonTimes(end) = [];
end


% balloon color colormap - [yellow, orange, read, gray]
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% task parameters in chronological order..
respTimes = trigTimes(trigs==23);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
nTrials = length(outcomeType);


if length(balloonType)>nTrials; balloonType=balloonType(1:end-1); end


% TODO:: don't analyz 'NaC' trodes...
 
% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
% plotting bank/pop responses
for ch2 = 1:nChans
    % epoch the spectral data for each channel.
    for tt = 1:nTrials
        updateUser('finished spectral calculations',tt,50,nTrials);
        % epoch the data here [channels X samples X trials]
        LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
    end
end
clear ch2 tt


%% Creating new lfp based on the trials that are removed because their RTs where outliers 

LFPmatNew = LFPmat(:, :, ~OutlierIndices);

%% creating IED timepoints data:


start_channel = 1;
end_channel = size(LFPmatNew,1)-1;
nTrials = length(ReactionTimesFiltered);


IED_timepoints = nan(size(LFPmatNew));

for chz = start_channel:end_channel
    for trial = 1:nTrials
        mySig = squeeze(LFPmatNew(chz,:,trial));
        IEDStruct = detectIEDs_single_array_v3(mySig,Fnew);
        IEDsInices = IEDStruct.foundPeaks.locs;
        IED_timepoints(chz,IEDsInices,trial) = 1;
    end
end


%% Stimulus time = balloon onset time = 1001
preStim = 1:975;
periStim = 975:1025;
BalloonOnsetTime = 1001;

% the postStim commented bellow is not correct and I think the range should
% be different for each trial based on the reaction times. so the post stim
% should be a vector of ranges that has the same size as the vector of
% trials. It has to be a matrix of nTrials*2. the first element in the
% first column is the start of the range and the second elements is the end
% of the range based on the reaction time of that trial. 

% postStim = 1026:1300; 
RangeSize = 2; % explanation above
DataPointsAfterBalloonOnset = 2000;
SecondsAfterBalloonOnset = 4;
StartTimePoint = 1026;
StartofRangeColumn = 1;
EndofRangeColumn = 2;


% filling up postStim matrix
postStim = nan(nTrials, RangeSize);

postStim(:,StartofRangeColumn) = StartTimePoint;
for trial=1:nTrials
    BalloonInflationTime = StartTimePoint + round((DataPointsAfterBalloonOnset*ReactionTimesFiltered(trial))/SecondsAfterBalloonOnset);
        if BalloonInflationTime> 3001
            BalloonInflationTime = 3001;
        end
        postStim(trial, EndofRangeColumn) = BalloonInflationTime;
end



% now I'm going to find 4 different RTs. It means I am going to find
% 1: RTs related to trials without IEDs happening.
% 2: RTs related to trials with IEDs that happenned preStim onset
% 3: RTs related to trials with IEDs that happenned periStim onset
% 4: RTs related to trials with IEDs that happenned postStim onset and
% before reaction start
% of course these RTs have overlaps in terms of IEDs happening time. I mean
% there are some trials that have IEDs during the whole time including pre
% peri and post stim onset. 

% finding indices of columns where the whole column is zero it means that is a
% non IED trials


IEDtrials_preStimRT = nan(nChans-1, nTrials);
IEDtrials_periStimRT = nan(nChans-1, nTrials);
IEDtrials_postStimRT = nan(nChans-1, nTrials);




for chan = 1:nChans-1
    for trial = 1:nTrials

        % pre
        slice = IED_timepoints(chan, preStim, trial);
        if any(slice(:) == 1)
            IEDtrials_preStimRT(chan, trial) = ReactionTimesFiltered(trial);
        end


        %peri
        slice = IED_timepoints(chan, periStim, trial);
        if any(slice(:) == 1)
             IEDtrials_periStimRT(chan, trial) = ReactionTimesFiltered(trial);
        end


        %post
        
        slice = IED_timepoints(chan, postStim(trial, StartofRangeColumn):postStim(trial, EndofRangeColumn), trial);
        if any(slice(:) == 1)
             IEDtrials_postStimRT(chan, trial) = ReactionTimesFiltered(trial);
        end


    end
end

%% finding NonIED trials
% ok, non IED trials are the trials that the IED happened there after the inflation has started!
% so for finding the NonIEDTrialsRTs, I am going to find which channels and
% which trials do not have any ieds detected from the balloon onset until the reaction 

NonIEDTrialsRTs = nan(nChans-1, nTrials);


for chan = 1:nChans-1
    for trial = 1:nTrials

        if all(isnan(squeeze(IED_timepoints(chan,postStim(trial, StartofRangeColumn):postStim(trial, EndofRangeColumn),trial))))
            NonIEDTrialsRTs(chan,trial) = ReactionTimesFiltered(trial);
        end
    end
end

%% Analysis:

pVal_preStimRT = nan(start_channel,end_channel);
pVal_periStimRT = nan(start_channel,end_channel);
pVal_postStimRT = nan(start_channel,end_channel);

for chz = start_channel:end_channel

    NonIEDTrials_temp_vec = NonIEDTrialsRTs(chz,:);
    NonIEDTrials_temp_vec = NonIEDTrials_temp_vec(~isnan(NonIEDTrials_temp_vec));

    IEDtrials_preStimRT_temp_vec = IEDtrials_preStimRT(chz,:);
    IEDtrials_preStimRT_temp_vec = IEDtrials_preStimRT_temp_vec(~isnan(IEDtrials_preStimRT_temp_vec));

    IEDtrials_periStimRT_temp_vec = IEDtrials_periStimRT(chz,:);
    IEDtrials_periStimRT_temp_vec = IEDtrials_periStimRT_temp_vec(~isnan(IEDtrials_periStimRT_temp_vec));


    IEDtrials_postStimRT_temp_vec = IEDtrials_postStimRT(chz,:);
    IEDtrials_postStimRT_temp_vec = IEDtrials_postStimRT_temp_vec(~isnan(IEDtrials_postStimRT_temp_vec));


    if (size(IEDtrials_preStimRT_temp_vec)>0)
        pVal_preStimRT(chz) = ranksum(IEDtrials_preStimRT_temp_vec,NonIEDTrials_temp_vec);
    else
        pVal_preStimRT(chz) = NaN;
    end

    if (size(IEDtrials_periStimRT_temp_vec)>0)
        pVal_periStimRT(chz) = ranksum(IEDtrials_periStimRT_temp_vec,NonIEDTrials_temp_vec);
    else
        pVal_periStimRT(chz) = NaN;
    end

    if (size(IEDtrials_postStimRT_temp_vec)>0)
        pVal_postStimRT(chz) = ranksum(IEDtrials_postStimRT_temp_vec,NonIEDTrials_temp_vec);
    else
        pVal_postStimRT(chz) = NaN;
    end 
end

%% Plotting pvals less than 0.05

pVal_preStimRT_filtered = pVal_preStimRT;
pVal_preStimRT_filtered(pVal_preStimRT >= 0.05) = NaN;

pVal_periStimRT_filtered = pVal_periStimRT;
pVal_periStimRT_filtered(pVal_periStimRT >= 0.05) = NaN;

pVal_postStimRT_filtered = pVal_postStimRT;
pVal_postStimRT_filtered(pVal_postStimRT >= 0.05) = NaN;



figure;

% Create subplot for Pre Stim RT
subplot(1,3,1);
scatter(1:length(pVal_preStimRT_filtered), pVal_preStimRT_filtered, 50, 'o', 'filled', 'DisplayName', 'Pre Stim RT');
xlabel('Channel number', 'FontSize', 18);
xlim([0 nChans-1]);
ylabel('P-values', 'FontSize', 18);
title('Pre Stim RT (1:balloon onset)', 'FontSize', 20);
legend('show');
hold on;

for i = 1:length(pVal_preStimRT)
    if isnan(pVal_preStimRT_filtered(i))
        continue;
    end
    text(i, pVal_preStimRT_filtered(i), ['channel = ' anatomicalLocs(selectedChans(i))], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
hold off;

% Create subplot for Peri Stim RT
subplot(1,3,2);
scatter(1:length(pVal_periStimRT_filtered), pVal_periStimRT_filtered, 50, 's', 'filled', 'DisplayName', 'Peri Stim RT');
xlabel('Channel number', 'FontSize', 18);
xlim([0 nChans-1]);
ylabel('P-values', 'FontSize', 18);
title('Peri Stim RT (balloon onset-50ms:balloon onset+50ms)', 'FontSize', 20);
legend('show');
hold on;
for i = 1:length(pVal_periStimRT)
    if isnan(pVal_periStimRT_filtered(i))
        continue;
    end
    text(i, pVal_periStimRT_filtered(i), ['channel = ' anatomicalLocs(selectedChans(i))], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
hold off;

% Create subplot for Post Stim RT
subplot(1,3,3);
scatter(1:length(pVal_postStimRT_filtered), pVal_postStimRT_filtered, 50, 'd', 'filled', 'DisplayName', 'Post Stim RT');
xlabel('Channel number', 'FontSize', 18);
xlim([0 nChans-1]);
ylabel('P-values', 'FontSize', 18);
title('Post Stim RT (balloon onset: inflation start)', 'FontSize', 20);
legend('show');
hold on;
for i = 1:length(pVal_postStimRT)
    if isnan(pVal_postStimRT_filtered(i))
        continue;
    end
    text(i, pVal_postStimRT_filtered(i), ['channel = ' anatomicalLocs(selectedChans(i))], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
hold off;



%% debugging
