% Across patient and in this code I will only save different data from
% different patient. We will call them IED data. and this is related to the
% reaction times (RTs).
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
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\IEDdata\'; 
outputIEDTrials = '\\155.100.91.44\d\Data\Nill\BART\IEDTrials\';
PatientsNum = length(fileList);
Fnew = 500;

%%

for pt = 1:PatientsNum
% for pt = 1:1
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
    % here I take popped trials first and then change 0s and 1s.
    poppedTrials = bhvStruct.poppedTrials;
    BankedTrials = ~poppedTrials;
    ReactTimeThreshold = 10;
    OutlierIndices = ReactionTimes >= ReactTimeThreshold;
    % ReactionTimesFiltered = ReactionTimes(~OutlierIndices);
    ReactionTimesFiltered = ReactionTimes; % I want to have all the RTs
    poppedTrials = poppedTrials(~OutlierIndices);
    BankedTrials = BankedTrials(~OutlierIndices);


    LFPmat = LFPmatStruct.LFPmat;
    % LFPmatNew = LFPmat(:, :, ~OutlierIndices);


    start_channel = 1;
    end_channel = size(LFPmat,1)-1;
    nTrials = length(ReactionTimesFiltered);
    nChans = length(LFPmatStruct.selectedChans);
    IEDtrials = nan(nChans-1, nTrials); 
    for chz = start_channel:end_channel
        IEDtrials(chz, :) = outliers(range(squeeze(LFPmat(chz,:,:))));
    end



    LFP_IED_trials = nan(size(LFPmat));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrials(chz, trial) == 1
                LFP_IED_trials(chz, :, trial) = LFPmat(chz, :, trial);
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

    % now let's modify IEDtrials:
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if any(squeeze(IED_timepoints(chz,:,trial))==1)
                IEDtrials(chz, trial) = 1;
            else 
                IEDtrials(chz, trial) = 0;
            end
        end
    end





    IEDdata.selectedChans = LFPmatStruct.selectedChans;
    IEDdata.IED_timepoints = IED_timepoints;
    IEDdata.anatomicalLocs = LFPmatStruct.anatomicalLocs; 
    IEDdata.IEDtrials = IEDtrials;

% I also wanna change IED trials in another folder for each patient

    IEDTrialsInfo.IEDtrials = IEDtrials;
    IEDTrialsInfo.IED_timepoints = IED_timepoints;



    save([outputIEDTrials ptID '.IEDtrials.mat'],'IEDTrialsInfo');

    save([outputFolderName ptID '.IEDdata.mat'],'IEDdata');

    clear NonIEDTrials_temp_vec IEDtrials_preStimRT_temp_vec IEDtrials_periStimRT_temp_vec slice pVal_preStimRT pVal_periStimRT pVal_postStimRT
    clear LFPmat IEDtrials LFP_IED_trials mySig IED_timepoints NonIEDTrialsRTs IEDtrials_preStimRT IEDtrials_periStimRT IEDtrials_postStimRT
    clear pVal_preStimRT_filtered pVal_periStimRT_filtered pVal_postStimRT_filtered RTSampleSize_NonIED RTSampleSize_preStim RTSampleSize_periStim RTSampleSize_postStim
end

