% I'm gonna create epoched LFPmat data for each patient. Also, I want to
% create a very good data structure. I am going to create the data below
% that the dimension is channels X 500(ms) sample X trials. 
% pay attention, bad channels are removed here!!

% LFPmatPreOnset1
% LFPmatPreOnset2
% LFPmatPostOnset1
% LFPmatPreResponse
% LFPmatPostInflate
% LFPmatPreOutcome

% next step:
% you need to run IED_11 for chunks to find the IED data
% author: Nill


clear;
clc;
close all;

%% loading neural and event data
 
% getting the numbers of different patients

inputFolderName = '\\155.100.91.44\d\Data\Nill\BART\bhvStruct_Nill_made';
fileList = dir(fullfile(inputFolderName, '*.bhvStruct.mat'));
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_6_chunks_500ms\LFPmat_6_chunks_500ms\'; 

%%

for pt = 1:length(fileList)
% for pt = 1:1
   
    data = load(fullfile(inputFolderName, fileList(pt).name));
   
    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1}; 
    disp(' ');
    disp(['Processing patient ID: ' patientID]);

    ptID = patientID;

    nevList = dir(['\\155.100.91.44\d\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);
    
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
    
    % resampling LFP at Fnew sampling frequency
    
    
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

    

    
    % timing parameters.
    % we are interested in 500ms of time and the sampling rate is 1000
    timeOfInterest = 0.5;

    
    

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
    
    
    % TODO:: don't analyz 'NaC' trodes...
     
    % epoching data

    LFPmatPreOnset1 = nan(nChans,Fs*timeOfInterest,nTrials); %-1000 to -500 to stim onsetr
    LFPmatPreOnset2 = nan(nChans,Fs*timeOfInterest,nTrials); %-500 to 0 to stim onset
    LFPmatPostOnset = nan(nChans,Fs*timeOfInterest,nTrials); % 0 to +500 to stim onset
    LFPmatPreResponse = nan(nChans,Fs*timeOfInterest,nTrials); % -500 to 0 to response time
    LFPmatPostResponse = nan(nChans,Fs*timeOfInterest,nTrials); % 0 to +500 after staring the inflation
    LFPmatPreOutcome = nan(nChans,Fs*timeOfInterest,nTrials); % -500 to 0 to outcome

    % plotting bank/pop responses
    for ch2 = 1:nChans
        % epoch the spectral data for each channel.
        for tt = 1:nTrials
            updateUser('finished spectral calculations',tt,50,nTrials);
            % epoch the data here [channels X samples X trials]
                        % now I need to check if this specific trial for each specific
            % channel is not a good signal. Then I won't even save it. and
            % those part of the lfpmat matrix would be null

            threshold = 7; % in dB
            ampThreshold = 5000;
   
            LFPmat_temp_vec = data2K(ch2,floor(Fs*balloonTimes(tt))-2*Fs*timeOfInterest:floor(Fs*balloonTimes(tt))-Fs*timeOfInterest-1);
            % LFPmat_mean = mean(abs(data2K(:,floor(Fs*balloonTimes(tt))-2*Fs*timeOfInterest:floor(Fs*balloonTimes(tt))-Fs*timeOfInterest-1)),1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy % skipping noisy channels
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPreOnset1(ch2,:,tt) = LFPmat_temp_vec;
                end
            end

            LFPmat_temp_vec = data2K(ch2,floor(Fs*balloonTimes(tt))-Fs*timeOfInterest:floor(Fs*balloonTimes(tt))-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPreOnset2(ch2,:,tt) = LFPmat_temp_vec;
                end
            end

            
            LFPmat_temp_vec = data2K(ch2,floor(Fs*balloonTimes(tt)): floor(Fs*balloonTimes(tt))+Fs*timeOfInterest-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPostOnset(ch2,:,tt) = LFPmat_temp_vec;    
                end
            end


            LFPmat_temp_vec = data2K(ch2,floor(Fs*respTimes(tt))-Fs*timeOfInterest:floor(Fs*respTimes(tt))-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPreResponse(ch2,:,tt) = LFPmat_temp_vec;
                end
            end


            LFPmat_temp_vec = data2K(ch2,floor(Fs*respTimes(tt)): floor(Fs*respTimes(tt))+Fs*timeOfInterest-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy  
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPostResponse(ch2,:,tt) = LFPmat_temp_vec;
                end
            end


            LFPmat_temp_vec = data2K(ch2,floor(Fs*outcomeTimes(tt))-Fs*timeOfInterest:floor(Fs*outcomeTimes(tt))-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPreOutcome(ch2,:,tt) = LFPmat_temp_vec;
                end
            end
        end
    end

    LFPmat.LFPmatPreOnset1 = LFPmatPreOnset1;
    LFPmat.LFPmatPreOnset2 = LFPmatPreOnset2;
    LFPmat.LFPmatPostOnset = LFPmatPostOnset;
    LFPmat.LFPmatPreResponse = LFPmatPreResponse;
    LFPmat.LFPmatPostResponse = LFPmatPostResponse;
    LFPmat.LFPmatPreOutcome = LFPmatPreOutcome;

    LFPmat.outcomeType = outcomeType;
    LFPmat.anatomicalLocs = anatomicalLocs;
    LFPmat.balloonType = balloonType;
    LFPmat.balloonTimes = balloonTimes;
    LFPmat.nTrials = nTrials;
    LFPmat.selectedChans = selectedChans;
    LFPmat.trodeLabels = trodeLabels;


 

    save([outputFolderName ptID '.LFPmat.mat'],'LFPmat')

    clear ch2 tt LFPmat LFPmat_temp_vec ch data2K balloonTimes balloonType anatomicalLocs nTrials trodeLabels 

end

%%  debug


% a_LFPmatPreOnset1  = squeeze(LFPmatPreOnset1(66,:,84));
% a_LFPmatPreOnset2  = squeeze(LFPmatPreOnset2(66,:,84));
% a_LFPmatPostOnset = squeeze(LFPmatPostOnset(66,:,84));
% a_LFPmatPreResponse = squeeze(LFPmatPreResponse(66,:,84));
% a_LFPmatPostResponse = squeeze(LFPmatPostResponse(66,:,84));
% a_LFPmatPreOutcome = squeeze(LFPmatPreOutcome(66,:,84));
% % % 
% 

% a_LFPmatPreOnset1  = squeeze(LFPmatPreOnset1(14,:,63));
% a_LFPmatPreOnset2  = squeeze(LFPmatPreOnset2(14,:,63));
% a_LFPmatPostOnset = squeeze(LFPmatPostOnset(14,:,63));
% a_LFPmatPreResponse = squeeze(LFPmatPreResponse(14,:,63));
% a_LFPmatPostResponse = squeeze(LFPmatPostResponse(14,:,63));
% a_LFPmatPreOutcome = squeeze(LFPmatPreOutcome(14,:,63));
% % % 



