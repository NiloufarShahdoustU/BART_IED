% finding the IEDs RTs that are significantly different from Non-IED RTs
% using outliers function to find IED trials
% author: Nill



% the first two versions of this code is considering different channels.
% now, I want to consider if an IED has occured in any channels and would
% suppose that as an IED. So I won't have anything for each and every
% channel.

clear;
clc;
close all;


%% loading behavioral data

ptID = '202003';
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


% finding the times of IED trials for every channels in every trials:
start_channel = 1;
end_channel = size(LFPmatNew,1)-1;

nTrials = length(ReactionTimesFiltered);

IEDtrials = nan(nChans-1, nTrials); 

for chz = start_channel:end_channel
    IEDtrials(chz, :) = outliers(range(squeeze(LFPmatNew(chz,:,:))));
end

%% in this part, we're gonna find IED trials, but not based on the individual channels
% but based on whether an IED has occured in different channels. Based on
% finding this IED trials, we can find non IED trials.

IEDtrials_whole = any(IEDtrials, 1);


%% creating IED data:

LFP_IED_trials = nan(size(LFPmatNew));

for chz = start_channel:end_channel
    for trial = 1:nTrials
        % check if the trial is marked as 1 in IEDtrials for the current channel
        if IEDtrials(chz, trial) == 1
            % copy the data from LFPmatNew to LFP_IED_trials for this channel and trial
            LFP_IED_trials(chz, :, trial) = LFPmatNew(chz, :, trial);
        end
    end
end


% obviously some trials in LFP_IED_trials are NaN


%% finding the outliers timepoints.



IED_timepoints = nan(size(LFP_IED_trials));
for chz = start_channel:end_channel
    for trial = 1:nTrials
        if IEDtrials(chz, trial) == 1
            mySig = squeeze(LFP_IED_trials(chz,:,trial));
            IEDStruct = detectIEDs_single_array_v3(mySig,Fnew);
            IEDsInices = IEDStruct.foundPeaks.locs;
            IED_timepoints(chz,IEDsInices,trial) = 1;
        end
    end
end




%% Stimulus time = balloon onset time = 1001
preStim = 1:975;
periStim = 980:1020;
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

NonIEDTrialsRTs = nan(1, length(IEDtrials_whole));
IEDtrials_preStimRT = nan(1, length(IEDtrials_whole));
IEDtrials_periStimRT = nan(1, length(IEDtrials_whole));
IEDtrials_postStimRT = nan(1, length(IEDtrials_whole));



for col = 1:nTrials

    if IEDtrials_whole(col) == 0 
        NonIEDTrialsRTs(col) = ReactionTimes(col);
    end

    % pre
    slice = squeeze(IED_timepoints(:, preStim, col));
    slice(isnan(slice)) = 0; % Replace NaN with 0
    if any(slice(:) == 1)
        IEDtrials_preStimRT(col) = ReactionTimes(col);
    end


    %peri
    slice = squeeze(IED_timepoints(:, periStim, col));
    slice(isnan(slice)) = 0; % Replace NaN with 0
    if any(slice(:) == 1)
         IEDtrials_periStimRT(col) = ReactionTimes(col);
    end


    %post
    slice = squeeze(IED_timepoints(:, round(postStim(trial, StartofRangeColumn)):round(postStim(trial, EndofRangeColumn)), col));
    slice(isnan(slice)) = 0; % Replace NaN with 0
    if any(slice(:) == 1)
         IEDtrials_postStimRT(col) = ReactionTimes(col);
    end


end






%% Analysis:



NonIEDTrialsRTs = NonIEDTrialsRTs(~isnan(NonIEDTrialsRTs));

IEDtrials_preStimRT = IEDtrials_preStimRT(~isnan(IEDtrials_preStimRT));

IEDtrials_periStimRT = IEDtrials_periStimRT(~isnan(IEDtrials_periStimRT));

IEDtrials_postStimRT = IEDtrials_postStimRT(~isnan(IEDtrials_postStimRT));

%% P values



pVal_preStimRT  = ranksum(IEDtrials_preStimRT,NonIEDTrialsRTs);
pVal_periStimRT= ranksum(IEDtrials_periStimRT,NonIEDTrialsRTs);
pVal_postStimRT= ranksum(IEDtrials_postStimRT,NonIEDTrialsRTs);


disp(['p-value for pre-stimulus RT: ', num2str(pVal_preStimRT)]);
disp(['p-value for peri-stimulus RT: ', num2str(pVal_periStimRT)]);
disp(['p-value for post-stimulus RT: ', num2str(pVal_postStimRT)]);
