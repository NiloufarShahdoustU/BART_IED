% I'm gonna find IEDs on data2k first and then save IEDs!! I've saved
% epochs before!


%  PreOnset1
%  PreOnset2
%  PostOnset
%  PreResponse
%  PostResponse
%  PreOutcome

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
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks\'; 

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
    


    % finding IEDs during the course of the signal:

    IED_timepoints = zeros(size(data2K));
    for chz = 1:nChans
         % disp(['Processing chan: (' int2str(chz) '/' int2str(nChans) ')']);
        mySig = squeeze(data2K(chz,:));
        IEDStruct = detectIEDs_single_array_v5(mySig,Fs);
        IEDsInices = IEDStruct.foundPeaks.locs;
        IED_timepoints(chz,IEDsInices) = 1;
        
    end

     
    % epoching data

    IED_timepointsPreOnset1 = nan(nChans,Fs*timeOfInterest,nTrials); %-1000 to -500 to stim onsetr
    IED_timepointsPreOnset2 = nan(nChans,Fs*timeOfInterest,nTrials); %-500 to 0 to stim onset
    IED_timepointsPostOnset = nan(nChans,Fs*timeOfInterest,nTrials); % 0 to +500 to stim onset
    IED_timepointsPreResponse = nan(nChans,Fs*timeOfInterest,nTrials); % -500 to 0 to response time
    IED_timepointsPostResponse = nan(nChans,Fs*timeOfInterest,nTrials); % 0 to +500 after staring the inflation
    IED_timepointsPreOutcome = nan(nChans,Fs*timeOfInterest,nTrials); % -500 to 0 to outcome



    LFPmatPreOnset1 = nan(nChans,Fs*timeOfInterest,nTrials); %-1000 to -500 to stim onsetr
    LFPmatPreOnset2 = nan(nChans,Fs*timeOfInterest,nTrials); %-500 to 0 to stim onset
    LFPmatPostOnset = nan(nChans,Fs*timeOfInterest,nTrials); % 0 to +500 to stim onset
    LFPmatPreResponse = nan(nChans,Fs*timeOfInterest,nTrials); % -500 to 0 to response time
    LFPmatPostResponse = nan(nChans,Fs*timeOfInterest,nTrials); % 0 to +500 after staring the inflation
    LFPmatPreOutcome = nan(nChans,Fs*timeOfInterest,nTrials); % -500 to 0 to outcome




    for ch2 = 1:nChans
        % epoch the spectral data for each channel.
        for tt = 1:nTrials
            % updateUser('finished spectral calculations',tt,50,nTrials);

            % epoch the data here [channels X samples X trials]
                        % now I need to check if this specific trial for each specific
            % channel is not a good signal. Then I won't even save it. and
            % those part of the lfpmat matrix would be null

            threshold = 7; % in dB
            ampThreshold = 5000;

            LFPmat_temp_vec = data2K(ch2,floor(Fs*balloonTimes(tt))-2*Fs*timeOfInterest:floor(Fs*balloonTimes(tt))-Fs*timeOfInterest-1);
            IED_timepoints_temp_vec = IED_timepoints(ch2,floor(Fs*balloonTimes(tt))-2*Fs*timeOfInterest:floor(Fs*balloonTimes(tt))-Fs*timeOfInterest-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy % skipping noisy channels
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPreOnset1(ch2,:,tt) = LFPmat_temp_vec;
                    IED_timepointsPreOnset1(ch2,:,tt) = IED_timepoints_temp_vec;

                end
            end

            LFPmat_temp_vec = data2K(ch2,floor(Fs*balloonTimes(tt))-Fs*timeOfInterest:floor(Fs*balloonTimes(tt))-1);
            IED_timepoints_temp_vec = IED_timepoints(ch2,floor(Fs*balloonTimes(tt))-Fs*timeOfInterest:floor(Fs*balloonTimes(tt))-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPreOnset2(ch2,:,tt) = LFPmat_temp_vec;
                    IED_timepointsPreOnset2(ch2,:,tt) = IED_timepoints_temp_vec;
                end
            end


            LFPmat_temp_vec = data2K(ch2,floor(Fs*balloonTimes(tt)): floor(Fs*balloonTimes(tt))+Fs*timeOfInterest-1);
            IED_timepoints_temp_vec = IED_timepoints(ch2,floor(Fs*balloonTimes(tt)): floor(Fs*balloonTimes(tt))+Fs*timeOfInterest-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPostOnset(ch2,:,tt) = LFPmat_temp_vec;   
                    IED_timepointsPostOnset(ch2,:,tt) = IED_timepoints_temp_vec;
                end
            end


            LFPmat_temp_vec = data2K(ch2,floor(Fs*respTimes(tt))-Fs*timeOfInterest:floor(Fs*respTimes(tt))-1);
            IED_timepoints_temp_vec = IED_timepoints(ch2,floor(Fs*respTimes(tt))-Fs*timeOfInterest:floor(Fs*respTimes(tt))-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPreResponse(ch2,:,tt) = LFPmat_temp_vec;
                    IED_timepointsPreResponse(ch2,:,tt) = IED_timepoints_temp_vec;
                end
            end


            LFPmat_temp_vec = data2K(ch2,floor(Fs*respTimes(tt)): floor(Fs*respTimes(tt))+Fs*timeOfInterest-1);
            IED_timepoints_temp_vec = IED_timepoints(ch2,floor(Fs*respTimes(tt)): floor(Fs*respTimes(tt))+Fs*timeOfInterest-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy  
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPostResponse(ch2,:,tt) = LFPmat_temp_vec;
                    IED_timepointsPostResponse(ch2,:,tt) = IED_timepoints_temp_vec;
                end
            end


            LFPmat_temp_vec = data2K(ch2,floor(Fs*outcomeTimes(tt))-Fs*timeOfInterest:floor(Fs*outcomeTimes(tt))-1);
            IED_timepoints_temp_vec = IED_timepoints(ch2,floor(Fs*outcomeTimes(tt))-Fs*timeOfInterest:floor(Fs*outcomeTimes(tt))-1);
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
            if ~isNoisy
                isOutofRange = detectLargeAmplitude(LFPmat_temp_vec,ampThreshold);
                if ~isOutofRange % skipping big amplitude channels
                    LFPmatPreOutcome(ch2,:,tt) = LFPmat_temp_vec;
                    IED_timepointsPreOutcome(ch2,:,tt) = IED_timepoints_temp_vec;
                end
            end
        end
    end



        % now let's create IEDtrials and channels:
    for chz = 1:nChans 
        for trial = 1:nTrials

            if any(squeeze(IED_timepointsPreOnset1(chz,:,trial))==1)
                IEDtrialsPreOnset1(chz, trial) = 1;
            else 
                IEDtrialsPreOnset1(chz, trial) = 0;
            end


            if any(squeeze(IED_timepointsPreOnset2(chz,:,trial))==1)
                IEDtrialsPreOnset2(chz, trial) = 1;
            else 
                IEDtrialsPreOnset2(chz, trial) = 0;
            end


            if any(squeeze(IED_timepointsPostOnset(chz,:,trial))==1)
                IEDtrialsPostOnset(chz, trial) = 1;
            else 
                IEDtrialsPostOnset(chz, trial) = 0;
            end                        


            if any(squeeze(IED_timepointsPreResponse(chz,:,trial))==1)
                IEDtrialsPreResponse(chz, trial) = 1;
            else 
                IEDtrialsPreResponse(chz, trial) = 0;
            end


            if any(squeeze(IED_timepointsPostResponse(chz,:,trial))==1)
                IEDtrialsPostResponse(chz, trial) = 1;
            else 
                IEDtrialsPostResponse(chz, trial) = 0;
            end

            if any(squeeze(IED_timepointsPreOutcome(chz,:,trial))==1)
                IEDtrialsPreOutcome(chz, trial) = 1;
            else 
                IEDtrialsPreOutcome(chz, trial) = 0;
            end


        end
    end   


    

    LFPIED.LFPmatPreOnset1 = LFPmatPreOnset1;
    LFPIED.LFPmatPreOnset2 = LFPmatPreOnset2;
    LFPIED.LFPmatPostOnset = LFPmatPostOnset;
    LFPIED.LFPmatPreResponse = LFPmatPreResponse;
    LFPIED.LFPmatPostResponse = LFPmatPostResponse;
    LFPIED.LFPmatPreOutcome = LFPmatPreOutcome;

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

    if (ptID == "202105")
        % this channel is a pain in the butt!!!!
        badChan = 35;
        IEDtrialsPreOnset1(badChan,:)= 0;
        IED_timepointsPreOnset1(badChan,:,:) = nan;

        IEDtrialsPreOnset2(badChan,:)= 0;
        IED_timepointsPreOnset2(badChan,:,:) = nan;

        IEDtrialsPostOnset(badChan,:)= 0;
        IED_timepointsPostOnset(badChan,:,:) = nan;

        IEDtrialsPreResponse(badChan,:)= 0;
        IED_timepointsPreResponse(badChan,:,:) = nan;

        IEDtrialsPostResponse(badChan,:)= 0;
        IED_timepointsPostResponse(badChan,:,:) = nan;

        IEDtrialsPreOutcome(badChan,:)= 0;
        IED_timepointsPreOutcome(badChan,:,:) = nan;

    end

    LFPIED.IED_timepointsPreOnset1 = IED_timepointsPreOnset1;
    LFPIED.IEDtrialsPreOnset1 = IEDtrialsPreOnset1;

    LFPIED.IED_timepointsPreOnset2 = IED_timepointsPreOnset2;
    LFPIED.IEDtrialsPreOnset2 = IEDtrialsPreOnset2;

    LFPIED.IED_timepointsPostOnset = IED_timepointsPostOnset;
    LFPIED.IEDtrialsPostOnset = IEDtrialsPostOnset;

    LFPIED.IED_timepointsPreResponse = IED_timepointsPreResponse;
    LFPIED.IEDtrialsPreResponse = IEDtrialsPreResponse;

    LFPIED.IED_timepointsPostResponse = IED_timepointsPostResponse;
    LFPIED.IEDtrialsPostResponse = IEDtrialsPostResponse;

    LFPIED.IED_timepointsPreOutcome = IED_timepointsPreOutcome;
    LFPIED.IEDtrialsPreOutcome = IEDtrialsPreOutcome;




    save([outputFolderName ptID '.LFPIED.mat'],'LFPIED');

    clear ch2 tt LFPmat LFPmat_temp_vec ch data2K balloonTimes balloonType anatomicalLocs nTrials trodeLabels 
    clear fileNameParts patientID LFPmat mySig IEDdata ReactionTimes poppedTrials BankedTrials OutlierIndices
    clear poppedTrials BankedTrials ReactionTimesFiltered
    clear LFPmatPreOnset1 IEDtrialsPreOnset1 LFP_IED_trialsPreOnset1 IED_timepointsPreOnset1
    clear LFPmatPreOnset2 IEDtrialsPreOnset2 LFP_IED_trialsPreOnset2 IED_timepointsPreOnset2
    clear LFPmatPostOnset IEDtrialsPostOnset LFP_IED_trialsPostOnset IED_timepointsPostOnset
    clear LFPmatPreResponse IEDtrialsPreResponse LFP_IED_trialsPreResponse IED_timepointsPreResponse    
    clear LFPmatPostResponse IEDtrialsPostResponse LFP_IED_trialsPostResponse IED_timepointsPostResponse  
    clear LFPmatPreOutcome IEDtrialsPreOutcome LFP_IED_trialsPreOutcome IED_timepointsPreOutcome  
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

% aaa = nansum(IED_timepoints(:));
% close all;
% aaa = squeeze(LFPIED.IED_timepointsPostOnset(2,:,:));
% imagesc(aaa);
