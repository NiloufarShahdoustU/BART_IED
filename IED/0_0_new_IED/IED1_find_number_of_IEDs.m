% I'm gonna find IEDs on data2k first and then save IEDs!! 

clear;
clc;
close all;


%% loading neural and event data
 
% getting the numbers of different patients

inputFolderName = 'D:\Nill\data\BART\bhvStruct_Nill_made';
fileList = dir(fullfile(inputFolderName, '*.bhvStruct.mat'));
outputFolderName = 'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\'; 

%%
 
for pt = 1:length(fileList)
% for pt = 1:1
   
    data = load(fullfile(inputFolderName, fileList(pt).name));
   
    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1}; 
    disp(' ');
    disp(['Processing patient ID: ' patientID]);

    ptID = patientID;

    nevList = dir(['D:\Nill\data\BART_preprocessed\' ptID '\Data\*.nev']);
    
    if length(nevList)>1
        error('many nev files available for this patient. Please specify...')
    elseif length(nevList)<1
        error('no nev files found...')
    else
        nevFile = fullfile(nevList.folder,nevList.name);
    end
    [trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesBART_2(ptID);
    
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


    % ReactionTimes = data.bhvStruct.allRTs;
    % here I take popped trials first and then change 0s and 1s.
    poppedTrials = data.bhvStruct.poppedTrials;
    BankedTrials = ~poppedTrials;
    % ReactTimeThreshold = 10;
    % OutlierIndices = ReactionTimes >= ReactTimeThreshold;
    % % ReactionTimesFiltered = ReactionTimes(~OutlierIndices);
    % ReactionTimesFiltered = ReactionTimes; % I want to have all the RTs
    % poppedTrials = poppedTrials(~OutlierIndices);
    % BankedTrials = BankedTrials(~OutlierIndices);


    
    
    notchFilter = true;
    for ch = 1:nChans
        if notchFilter
            [b1,a1] = iirnotch(60/(Fs/2),(60/(Fs/2))/50);
            tmp(ch,:) = filtfilt(b1,a1,resample(double(NSX.Data(selectedChans(ch),:)),Fs,Fs));
            
            [b2,a2] = iirnotch(120/(Fs/2),(120/(Fs/2))/50);
            data2K(ch,:) = filtfilt(b2,a2,tmp(ch,:));
        else
            data2K(ch,:) = resample(double(NSX.Data(selectedChans(ch),:)),Fs,Fs);
        end
        clear tmp 
         
    end
    clear NSX NEV a1 a2 b1 b2

    


    
    

    % task parameters in chronological order..
    % There aren't any trigs that == 4
    balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
    inflateTimes = trigTimes(trigs==23);
    balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
    if length(balloonTimes)>length(inflateTimes)
        balloonTimes(end) = [];
    end
    

    % task parameters in chronological order..
    respTimes = trigTimes(trigs==23);
    outcomeTimes = trigTimes(trigs==25 | trigs==26);
    outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
    nTrials = length(outcomeType);
    
    
    if length(balloonType)>nTrials; balloonType=balloonType(1:end-1); end
    


    % finding IEDs during the course of the signal:

    IED_timepoints = zeros(size(data2K));
    for chz = 1:nChans
        mySig = squeeze(data2K(chz,:));
        IEDStruct = detectIEDs_single_array_v8(mySig,Fs);
        IEDsInices = IEDStruct.foundPeaks.locs;
        IED_timepoints(chz,IEDsInices) = 1;
        
    end



    % timing parameters.
    % this time of interest should be different for each trial and also
    % should be 2 times, one RT and one IT. 



    timeOfInterestRT = data.bhvStruct.allRTs; % one RT duration per trial
    timeOfInterestIT = data.bhvStruct.allITs; % one IT duration per trial
    
    IED_timepointsRT = cell(nTrials,1);
    IED_timepointsIT = cell(nTrials,1);
    
    for tr = 1:nTrials
    
        nSamplesRT = floor(Fs * timeOfInterestRT(tr));
        nSamplesIT = floor(Fs * timeOfInterestIT(tr));
    
        IED_timepointsRT{tr} = nan(nChans, nSamplesRT);
        IED_timepointsIT{tr} = nan(nChans, nSamplesIT);
    
    end


    

    for ch2 = 1:nChans
        for tt = 1:nTrials
            threshold = 7; % in dB
            ampThreshold = 5000;

            LFPmat_temp_vec = data2K(ch2,floor(Fs*balloonTimes(tt)): floor(Fs*balloonTimes(tt))+ floor(Fs*timeOfInterestRT(tt))-1);
            IED_timepoints_temp_vec = IED_timepoints(ch2,floor(Fs*balloonTimes(tt)): floor(Fs*balloonTimes(tt))+ floor(Fs*timeOfInterestRT(tt))-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    % LFPmatRT(ch2,:,tt) = LFPmat_temp_vec;   
                     IED_timepointsRT{tt}(ch2,:) = IED_timepoints_temp_vec;
                end
            end




            LFPmat_temp_vec = data2K(ch2,floor(Fs*respTimes(tt)): floor(Fs*respTimes(tt))+ floor(Fs*timeOfInterestIT(tt))-1);
            IED_timepoints_temp_vec = IED_timepoints(ch2,floor(Fs*respTimes(tt)): floor(Fs*respTimes(tt))+ floor(Fs*timeOfInterestIT(tt))-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy  
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    % LFPmatIT(ch2,:,tt) = LFPmat_temp_vec;
                     IED_timepointsIT{tt}(ch2,:) = IED_timepoints_temp_vec;
                end
            end

        end
    end




    % tracking IED occurrence
    % Each row:
    % column 1 = trial
    % column 2 = channel index chz
    % column 3 = time index inside IED_timepointsRT{trial}(chz, :)
    
    IED_occurance_RT = [];
    IED_occurance_IT = [];

    for chz = 1:nChans
        for trial = 1:nTrials
    
            % RT IED occurrence
            tempRT = IED_timepointsRT{trial}(chz, :);
    
            % index of IEDs inside this RT window
            IED_idx_RT = find(tempRT == 1);
    
            if ~isempty(IED_idx_RT)
                tempRowsRT = [ ...
                    trial * ones(length(IED_idx_RT), 1), ...
                    chz   * ones(length(IED_idx_RT), 1), ...
                    IED_idx_RT(:) ...
                ];
    
                IED_occurance_RT = [IED_occurance_RT; tempRowsRT];
            end
    
    
            % IT IED occurrence
            tempIT = IED_timepointsIT{trial}(chz, :);
    
            % index of IEDs inside this IT window
            IED_idx_IT = find(tempIT == 1);
    
            if ~isempty(IED_idx_IT)
                tempRowsIT = [ ...
                    trial * ones(length(IED_idx_IT), 1), ...
                    chz   * ones(length(IED_idx_IT), 1), ...
                    IED_idx_IT(:) ...
                ];
    
                IED_occurance_IT = [IED_occurance_IT; tempRowsIT];
            end
    
        end
    end
   

    if ~isempty(IED_occurance_RT)
    IED_occurance_RT = sortrows(IED_occurance_RT, 1);
    end
    
    if ~isempty(IED_occurance_IT)
        IED_occurance_IT = sortrows(IED_occurance_IT, 1);
    end
    
    LFPIED.IED_occurance_RT = IED_occurance_RT;
    LFPIED.IED_occurance_IT = IED_occurance_IT;
    
    LFPIED.IED_occurance_RT = IED_occurance_RT;
    LFPIED.IED_occurance_IT = IED_occurance_IT;
    
    LFPIED.IED_occurance_RT_columns = {'trial', 'channel', 'time_index_within_RT'};
    LFPIED.IED_occurance_IT_columns = {'trial', 'channel', 'time_index_within_IT'};


    LFPIED.outcomeType = outcomeType;
    LFPIED.anatomicalLocs = anatomicalLocs;
    LFPIED.balloonType = balloonType;
    LFPIED.balloonTimes = balloonTimes;
    LFPIED.nTrials = nTrials;
    LFPIED.selectedChans = selectedChans;
    LFPIED.trodeLabels = trodeLabels;


    LFPIED.RTs = data.bhvStruct.allRTs; 
    LFPIED.ITs = data.bhvStruct.allITs;
    LFPIED.isControl = data.bhvStruct.isCtrl;
    LFPIED.poppedTrials = poppedTrials;
    LFPIED.BankedTrials = BankedTrials;



    save([outputFolderName ptID '.LFPIED.mat'],'LFPIED');

    clear ch2 tt LFPmat LFPmat_temp_vec ch data2K balloonTimes balloonType anatomicalLocs nTrials trodeLabels 
    clear fileNameParts patientID LFPmat mySig IEDdata ReactionTimes poppedTrials BankedTrials OutlierIndices
    clear poppedTrials BankedTrials ReactionTimesFiltered
    clear IED_timepointsIT IED_timepointsRT timeOfInterestRT timeOfInterestIT
end

%%  debug


trialNums = LFPIED.IED_occurance_RT(:,1);

[uniqueTrials, ~, idx] = unique(trialNums);

numIEDsPerTrial = accumarray(idx, 1);

result = [uniqueTrials numIEDsPerTrial];
