% Across patient and in this code I will only save different data from
% different patient. We will call them IED data and this is related to the
% accuracies
% Author: Nill

clear;
clc;
close all;
warning('off','all');

%% loading neural and event data
 
% getting the numbers of different patients 

inputFolderName_bhvStruct = '\\155.100.91.44\d\Data\Nill\BART\bhvStruct_Nill_made';
inputFolderName_LFPmat = '\\155.100.91.44\d\Data\Nill\BART\LFPmat';
fileList = dir(fullfile(inputFolderName_bhvStruct, '*.bhvStruct.mat'));
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\IEDdata_ITs\'; 
Fnew = 500;
PatientsNum = length(fileList);

%%

for pt = 1:PatientsNum
    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1}; 
    disp(' ');
    disp(['Processing patient ID: ' patientID ' (' int2str(pt) '/' int2str(PatientsNum) ')']);
    ptID = patientID;

    matFile_bhvStruct = [inputFolderName_bhvStruct '\' ptID '.bhvStruct.mat'];
    load(matFile_bhvStruct);

    matFile_bhvStruct = [inputFolderName_LFPmat '\' ptID '.LFPmat.mat'];
    load(matFile_bhvStruct);

    ReactionTimes = bhvStruct.allRTs;
    InflationTimes = bhvStruct.allITs;
    ReactTimeThreshold = 10;
    OutlierIndices = ReactionTimes >= ReactTimeThreshold;
    ReactionTimesFiltered = ReactionTimes(~OutlierIndices);
    InflationTimes = InflationTimes(~OutlierIndices);



    LFPmat = LFPmatStruct.LFPmat;
    LFPmatNew = LFPmat(:, :, ~OutlierIndices);


    start_channel = 1;
    end_channel = size(LFPmatNew,1)-1;
    nTrials = length(ReactionTimesFiltered);
    nChans = length(LFPmatStruct.selectedChans);
    IEDtrials = nan(nChans-1, nTrials); 
    for chz = start_channel:end_channel
        IEDtrials(chz, :) = outliers(range(squeeze(LFPmatNew(chz,:,:))));
    end



    LFP_IED_trials = nan(size(LFPmatNew));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrials(chz, trial) == 1
                LFP_IED_trials(chz, :, trial) = LFPmatNew(chz, :, trial);
            end
        end
    end



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


    preStim = 1001-250:1001; 
    periStim = 1001-125:1001+125; %- 250ms before and +250 ms after stim onset
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
    StartTimePoint = 1001+125;
    StartofRangeColumn = 1;
    EndofRangeColumn = 2;
    
    
    % filling up postStim matrix
    postStim = nan(nTrials, RangeSize);

    postStim(:,StartofRangeColumn) = StartTimePoint;
    postStim(:,EndofRangeColumn) = 1001+250;
    % 
    % postStim(:,StartofRangeColumn) = StartTimePoint;
    % for trial=1:nTrials
    %     BalloonInflationTime = StartTimePoint + round((DataPointsAfterBalloonOnset*ReactionTimesFiltered(trial))/SecondsAfterBalloonOnset);
    %         if BalloonInflationTime> 3001
    %             BalloonInflationTime = 3001;
    %         elseif BalloonInflationTime < 1001+125+250
    %         postStim(trial, EndofRangeColumn) = BalloonInflationTime;
    %         else
    %             postStim(trial, EndofRangeColumn) = 1001+125+250;
    %         end
    % end
    % 



% these are written as RTs but it's only the name of the varialbe
%  I'm actually looking at banked trials.

    NonIEDTrialsRTs = nan(size(IEDtrials));
    IEDtrials_preStimRT = nan(size(IEDtrials));
    IEDtrials_periStimRT = nan(size(IEDtrials));
    IEDtrials_postStimRT = nan(size(IEDtrials));
    
    
    
    
    for row = 1:nChans-1
        for col = 1:nTrials
    
            % NonIED; IED trials that are 0 are NonIED trials
            if IEDtrials(row, col) == 0 
                NonIEDTrialsRTs(row, col) = InflationTimes(col);
            end
    
            % pre
            slice = IED_timepoints(row, preStim, col);
            slice(isnan(slice)) = 0; % Replace NaN with 0
            if any(slice(:) == 1)
                IEDtrials_preStimRT(row, col) = InflationTimes(col);
            end
    
    
            %peri
            slice = IED_timepoints(row, periStim, col);
            slice(isnan(slice)) = 0; % Replace NaN with 0
            if any(slice(:) == 1)
                 IEDtrials_periStimRT(row, col) = InflationTimes(col);
            end
    
    
            %post
            slice = IED_timepoints(row, round(postStim(trial, StartofRangeColumn)):round(postStim(trial, EndofRangeColumn)), col);
            slice(isnan(slice)) = 0; % Replace NaN with 0
            if any(slice(:) == 1)
                 IEDtrials_postStimRT(row, col) = InflationTimes(col);
            end
    
    
        end
    end
    

    pVal_preStimRT = nan(start_channel,end_channel);
    pVal_periStimRT = nan(start_channel,end_channel);
    pVal_postStimRT = nan(start_channel,end_channel);
    
    RTSampleSize_NonIED = nan(1, nChans-1);
    RTSampleSize_preStim = nan(1, nChans-1);
    RTSampleSize_periStim = nan(1, nChans-1);
    RTSampleSize_postStim = nan(1, nChans-1);
    
    NumberofPermutations = 10000;
    
    
    % here I want to take a look at the number of samples that I am taking the
    % ranksum test from. So, I would need to save the length of NON ied, and
    % ied(pre, peri, post) for each channel to show them on the final data. 
    
    for chz = start_channel:end_channel
    
        NonIEDTrials_temp_vec = NonIEDTrialsRTs(chz,:);
        NonIEDTrials_temp_vec = NonIEDTrials_temp_vec(~isnan(NonIEDTrials_temp_vec));
        RTSampleSize_NonIED(chz) = length(NonIEDTrials_temp_vec);
    
        IEDtrials_preStimRT_temp_vec = IEDtrials_preStimRT(chz,:);
        IEDtrials_preStimRT_temp_vec = IEDtrials_preStimRT_temp_vec(~isnan(IEDtrials_preStimRT_temp_vec));
        RTSampleSize_preStim(chz) = length(IEDtrials_preStimRT_temp_vec);
    
        IEDtrials_periStimRT_temp_vec = IEDtrials_periStimRT(chz,:);
        IEDtrials_periStimRT_temp_vec = IEDtrials_periStimRT_temp_vec(~isnan(IEDtrials_periStimRT_temp_vec));
        RTSampleSize_periStim(chz) = length(IEDtrials_periStimRT_temp_vec);
    
    
    
        IEDtrials_postStimRT_temp_vec = IEDtrials_postStimRT(chz,:);
        IEDtrials_postStimRT_temp_vec = IEDtrials_postStimRT_temp_vec(~isnan(IEDtrials_postStimRT_temp_vec));
        RTSampleSize_postStim(chz) = length(IEDtrials_postStimRT_temp_vec);
    
    
        if (size(IEDtrials_preStimRT_temp_vec)>0)
            pVal_preStimRT(chz) = permutationTest(IEDtrials_preStimRT_temp_vec,NonIEDTrials_temp_vec, NumberofPermutations);
        else
            pVal_preStimRT(chz) = NaN;
        end
    
        if (size(IEDtrials_periStimRT_temp_vec)>0)
            pVal_periStimRT(chz) = permutationTest(IEDtrials_periStimRT_temp_vec,NonIEDTrials_temp_vec, NumberofPermutations);
        else
            pVal_periStimRT(chz) = NaN;
        end
    
        if (size(IEDtrials_postStimRT_temp_vec)>0)
            pVal_postStimRT(chz) = permutationTest(IEDtrials_postStimRT_temp_vec,NonIEDTrials_temp_vec, NumberofPermutations);
        else
            pVal_postStimRT(chz) = NaN;
        end 
    end
    
    
    
    pVal_preStimRT_filtered = pVal_preStimRT;
    pVal_preStimRT_filtered(pVal_preStimRT >= 0.05) = NaN;
    
    pVal_periStimRT_filtered = pVal_periStimRT;
    pVal_periStimRT_filtered(pVal_periStimRT >= 0.05) = NaN;
    
    pVal_postStimRT_filtered = pVal_postStimRT;
    pVal_postStimRT_filtered(pVal_postStimRT >= 0.05) = NaN;
    


    IEDdata.selectedChans = LFPmatStruct.selectedChans;
    IEDdata.anatomicalLocs = LFPmatStruct.anatomicalLocs;
    IEDdata.NonIEDTrialsRTs = NonIEDTrialsRTs;
    IEDdata.IEDtrials_preStimRT = IEDtrials_preStimRT;
    IEDdata.IEDtrials_periStimRT = IEDtrials_periStimRT;
    IEDdata.IEDtrials_postStimRT = IEDtrials_postStimRT;
    IEDdata.pVal_preStimRT = pVal_preStimRT;
    IEDdata.pVal_preStimRT_filtered = pVal_preStimRT_filtered;
    IEDdata.pVal_periStimRT = pVal_periStimRT;
    IEDdata.pVal_periStimRT_filtered = pVal_periStimRT_filtered;
    IEDdata.pVal_postStimRT = pVal_postStimRT;
    IEDdata.pVal_postStimRT_filtered = pVal_postStimRT_filtered;
    IEDdata.RTSampleSize_NonIED = RTSampleSize_NonIED;
    IEDdata.RTSampleSize_preStim = RTSampleSize_preStim;
    IEDdata.RTSampleSize_periStim = RTSampleSize_periStim;
    IEDdata.RTSampleSize_postStim = RTSampleSize_postStim;




    save([outputFolderName ptID '.IEDdata.mat'],'IEDdata')

    clear NonIEDTrials_temp_vec IEDtrials_preStimRT_temp_vec IEDtrials_periStimRT_temp_vec slice pVal_preStimRT pVal_periStimRT pVal_postStimRT
    clear LFPmatNew IEDtrials LFP_IED_trials mySig IED_timepoints NonIEDTrialsRTs IEDtrials_preStimRT IEDtrials_periStimRT IEDtrials_postStimRT
    clear pVal_preStimRT_filtered pVal_periStimRT_filtered pVal_postStimRT_filtered RTSampleSize_NonIED RTSampleSize_preStim RTSampleSize_periStim RTSampleSize_postStim
end
