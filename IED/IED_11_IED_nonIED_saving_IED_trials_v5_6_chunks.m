% here I am trying to do the IED detection because I need to find the IEDs
% before starting the analysis cause I need to find out if my algorithm
% works fine!!!
% In this code I am using the newest version of IED detection code:
% which is detectIEDs_single_array_v5.m
% after this, you need to run IED_13 v4, in order to find a good IED detection  algorithm!
% and do the sanity  check of IED detection


clear;
clc;
close all;
warning('off','all');

%% loading neural and event data
 
% getting the numbers of different patients

inputFolderName_bhvStruct = '\\155.100.91.44\d\Data\Nill\BART\bhvStruct_Nill_made';
inputFolderName_LFPmat = '\\155.100.91.44\d\Data\Nill\BART\LFPmat_6_chunks_500ms';
fileList = dir(fullfile(inputFolderName_bhvStruct, '*.bhvStruct.mat'));
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_6_chunks_500ms\IEDdata_bad_chans_removed\'; 
PatientsNum = length(fileList);
Fnew = 1000;

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


%% PreOnset1
    LFPmatPreOnset1 = LFPmat.LFPmatPreOnset1;

    start_channel = 1;
    end_channel = size(LFPmatPreOnset1,1)-1;
    nTrials = length(ReactionTimesFiltered);
    nChans = length(LFPmat.selectedChans);
    IEDtrialsPreOnset1 = nan(nChans-1, nTrials); 
    for chz = start_channel:end_channel
        IEDtrialsPreOnset1(chz, :) = outliers(range(squeeze(LFPmatPreOnset1(chz,:,:))));
    end


    LFP_IED_trialsPreOnset1 = nan(size(LFPmatPreOnset1));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPreOnset1(chz, trial) == 1
                LFP_IED_trialsPreOnset1(chz, :, trial) = LFPmatPreOnset1(chz, :, trial);
            end
        end
    end



    IED_timepointsPreOnset1 = nan(size(LFP_IED_trialsPreOnset1));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPreOnset1(chz, trial) == 1
                mySig = squeeze(LFP_IED_trialsPreOnset1(chz,:,trial));
                IEDStruct = detectIEDs_single_array_v5(mySig,Fnew);
                IEDsInices = IEDStruct.foundPeaks.locs;
                IED_timepointsPreOnset1(chz,IEDsInices,trial) = 1;
            end
        end
    end

    % now let's modify IEDtrials:
    for chz = start_channel:end_channel 
        for trial = 1:nTrials
            if any(squeeze(IED_timepointsPreOnset1(chz,:,trial))==1)
                IEDtrialsPreOnset1(chz, trial) = 1;
            else 
                IEDtrialsPreOnset1(chz, trial) = 0;
            end
        end
    end                             

%% PreOnset2
    LFPmatPreOnset2 = LFPmat.LFPmatPreOnset2;

    start_channel = 1;
    end_channel = size(LFPmatPreOnset2,1)-1;
    nTrials = length(ReactionTimesFiltered);
    nChans = length(LFPmat.selectedChans);
    IEDtrialsPreOnset2 = nan(nChans-1, nTrials); 
    for chz = start_channel:end_channel
        IEDtrialsPreOnset2(chz, :) = outliers(range(squeeze(LFPmatPreOnset2(chz,:,:))));
    end


    LFP_IED_trialsPreOnset2 = nan(size(LFPmatPreOnset2));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPreOnset2(chz, trial) == 1
                LFP_IED_trialsPreOnset2(chz, :, trial) = LFPmatPreOnset2(chz, :, trial);
            end
        end
    end



    IED_timepointsPreOnset2 = nan(size(LFP_IED_trialsPreOnset2));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPreOnset2(chz, trial) == 1
                mySig = squeeze(LFP_IED_trialsPreOnset2(chz,:,trial));
                IEDStruct = detectIEDs_single_array_v5(mySig,Fnew);
                IEDsInices = IEDStruct.foundPeaks.locs;
                IED_timepointsPreOnset2(chz,IEDsInices,trial) = 1;
            end
        end
    end

    % now let's modify IEDtrials:
    for chz = start_channel:end_channel 
        for trial = 1:nTrials
            if any(squeeze(IED_timepointsPreOnset2(chz,:,trial))==1)
                IEDtrialsPreOnset2(chz, trial) = 1;
            else 
                IEDtrialsPreOnset2(chz, trial) = 0;
            end
        end
    end   




%% PostOnset
    LFPmatPostOnset = LFPmat.LFPmatPostOnset;

    start_channel = 1;
    end_channel = size(LFPmatPostOnset,1)-1;
    nTrials = length(ReactionTimesFiltered);
    nChans = length(LFPmat.selectedChans);
    IEDtrialsPostOnset = nan(nChans-1, nTrials); 
    for chz = start_channel:end_channel
        IEDtrialsPostOnset(chz, :) = outliers(range(squeeze(LFPmatPostOnset(chz,:,:))));
    end


    LFP_IED_trialsPostOnset = nan(size(LFPmatPostOnset));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPostOnset(chz, trial) == 1
                LFP_IED_trialsPostOnset(chz, :, trial) = LFPmatPostOnset(chz, :, trial);
            end
        end
    end



    IED_timepointsPostOnset = nan(size(LFP_IED_trialsPostOnset));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPostOnset(chz, trial) == 1
                mySig = squeeze(LFP_IED_trialsPostOnset(chz,:,trial));
                IEDStruct = detectIEDs_single_array_v5(mySig,Fnew);
                IEDsInices = IEDStruct.foundPeaks.locs;
                IED_timepointsPostOnset(chz,IEDsInices,trial) = 1;
            end
        end
    end

    % now let's modify IEDtrials:
    for chz = start_channel:end_channel 
        for trial = 1:nTrials
            if any(squeeze(IED_timepointsPostOnset(chz,:,trial))==1)
                IEDtrialsPostOnset(chz, trial) = 1;
            else 
                IEDtrialsPostOnset(chz, trial) = 0;
            end
        end
    end 



%% PreResponse

    LFPmatPreResponse = LFPmat.LFPmatPreResponse;

    start_channel = 1;
    end_channel = size(LFPmatPreResponse,1)-1;
    nTrials = length(ReactionTimesFiltered);
    nChans = length(LFPmat.selectedChans);
    IEDtrialsPreResponse = nan(nChans-1, nTrials); 
    for chz = start_channel:end_channel
        IEDtrialsPreResponse(chz, :) = outliers(range(squeeze(LFPmatPreResponse(chz,:,:))));
    end


    LFP_IED_trialsPreResponse = nan(size(LFPmatPreResponse));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPreResponse(chz, trial) == 1
                LFP_IED_trialsPreResponse(chz, :, trial) = LFPmatPreResponse(chz, :, trial);
            end
        end
    end



    IED_timepointsPreResponse = nan(size(LFP_IED_trialsPreResponse));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPreResponse(chz, trial) == 1
                mySig = squeeze(LFP_IED_trialsPreResponse(chz,:,trial));
                IEDStruct = detectIEDs_single_array_v5(mySig,Fnew);
                IEDsInices = IEDStruct.foundPeaks.locs;
                IED_timepointsPreResponse(chz,IEDsInices,trial) = 1;
            end
        end
    end

    % now let's modify IEDtrials:
    for chz = start_channel:end_channel 
        for trial = 1:nTrials
            if any(squeeze(IED_timepointsPreResponse(chz,:,trial))==1)
                IEDtrialsPreResponse(chz, trial) = 1;
            else 
                IEDtrialsPreResponse(chz, trial) = 0;
            end
        end
    end 



%% PostResponse

    LFPmatPostResponse = LFPmat.LFPmatPostResponse;

    start_channel = 1;
    end_channel = size(LFPmatPostResponse,1)-1;
    nTrials = length(ReactionTimesFiltered);
    nChans = length(LFPmat.selectedChans);
    IEDtrialsPostResponse = nan(nChans-1, nTrials); 
    for chz = start_channel:end_channel
        IEDtrialsPostResponse(chz, :) = outliers(range(squeeze(LFPmatPostResponse(chz,:,:))));
    end


    LFP_IED_trialsPostResponse = nan(size(LFPmatPostResponse));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPostResponse(chz, trial) == 1
                LFP_IED_trialsPostResponse(chz, :, trial) = LFPmatPostResponse(chz, :, trial);
            end
        end
    end



    IED_timepointsPostResponse = nan(size(LFP_IED_trialsPostResponse));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPostResponse(chz, trial) == 1
                mySig = squeeze(LFP_IED_trialsPostResponse(chz,:,trial));
                IEDStruct = detectIEDs_single_array_v5(mySig,Fnew);
                IEDsInices = IEDStruct.foundPeaks.locs;
                IED_timepointsPostResponse(chz,IEDsInices,trial) = 1;
            end
        end
    end

    % now let's modify IEDtrials:
    for chz = start_channel:end_channel 
        for trial = 1:nTrials
            if any(squeeze(IED_timepointsPostResponse(chz,:,trial))==1)
                IEDtrialsPostResponse(chz, trial) = 1;
            else 
                IEDtrialsPostResponse(chz, trial) = 0;
            end
        end
    end 


%% PreOutcome

    LFPmatPreOutcome = LFPmat.LFPmatPreOutcome;

    start_channel = 1;
    end_channel = size(LFPmatPreOutcome,1)-1;
    nTrials = length(ReactionTimesFiltered);
    nChans = length(LFPmat.selectedChans);
    IEDtrialsPreOutcome = nan(nChans-1, nTrials); 
    for chz = start_channel:end_channel
        IEDtrialsPreOutcome(chz, :) = outliers(range(squeeze(LFPmatPreOutcome(chz,:,:))));
    end


    LFP_IED_trialsPreOutcome = nan(size(LFPmatPreOutcome));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPreOutcome(chz, trial) == 1
                LFP_IED_trialsPreOutcome(chz, :, trial) = LFPmatPreOutcome(chz, :, trial);
            end
        end
    end



    IED_timepointsPreOutcome = nan(size(LFP_IED_trialsPreOutcome));
    for chz = start_channel:end_channel
        for trial = 1:nTrials
            if IEDtrialsPreOutcome(chz, trial) == 1
                mySig = squeeze(LFP_IED_trialsPreOutcome(chz,:,trial));
                IEDStruct = detectIEDs_single_array_v5(mySig,Fnew);
                IEDsInices = IEDStruct.foundPeaks.locs;
                IED_timepointsPreOutcome(chz,IEDsInices,trial) = 1;
            end
        end
    end

    % now let's modify IEDtrials:
    for chz = start_channel:end_channel 
        for trial = 1:nTrials
            if any(squeeze(IED_timepointsPreOutcome(chz,:,trial))==1)
                IEDtrialsPreOutcome(chz, trial) = 1;
            else 
                IEDtrialsPreOutcome(chz, trial) = 0;
            end
        end
    end 



%% save


    IEDdata.selectedChans = LFPmat.selectedChans;
    IEDdata.anatomicalLocs = LFPmat.anatomicalLocs; 
    IEDdata.RTs = bhvStruct.allRTs;
    IEDdata.ITs = bhvStruct.allITs;
    IEDdata.poppedTrials = poppedTrials;
    IEDdata.BankedTrials = BankedTrials;

    IEDdata.IED_timepointsPreOnset1 = IED_timepointsPreOnset1;
    IEDdata.IEDtrialsPreOnset1 = IEDtrialsPreOnset1;

    IEDdata.IED_timepointsPreOnset2 = IED_timepointsPreOnset2;
    IEDdata.IEDtrialsPreOnset2 = IEDtrialsPreOnset2;

    IEDdata.IED_timepointsPostOnset = IED_timepointsPostOnset;
    IEDdata.IEDtrialsPostOnset = IEDtrialsPostOnset;

    IEDdata.IED_timepointsPreResponse = IED_timepointsPreResponse;
    IEDdata.IEDtrialsPreResponse = IEDtrialsPreResponse;

    IEDdata.IED_timepointsPostResponse = IED_timepointsPostResponse;
    IEDdata.IEDtrialsPostResponse = IEDtrialsPostResponse;

    IEDdata.IED_timepointsPreOutcome = IED_timepointsPreOutcome;
    IEDdata.IEDtrialsPreOutcome = IEDtrialsPreOutcome;



    save([outputFolderName ptID '.IEDdata.mat'],'IEDdata');

    clear fileNameParts patientID LFPmat mySig IEDdata ReactionTimes poppedTrials BankedTrials OutlierIndices
    clear poppedTrials BankedTrials ReactionTimesFiltered
    clear LFPmatPreOnset1 IEDtrialsPreOnset1 LFP_IED_trialsPreOnset1 IED_timepointsPreOnset1
    clear LFPmatPreOnset2 IEDtrialsPreOnset2 LFP_IED_trialsPreOnset2 IED_timepointsPreOnset2
    clear LFPmatPostOnset IEDtrialsPostOnset LFP_IED_trialsPostOnset IED_timepointsPostOnset
    clear LFPmatPreResponse IEDtrialsPreResponse LFP_IED_trialsPreResponse IED_timepointsPreResponse    
    clear LFPmatPostResponse IEDtrialsPostResponse LFP_IED_trialsPostResponse IED_timepointsPostResponse  
    clear LFPmatPreOutcome IEDtrialsPreOutcome LFP_IED_trialsPreOutcome IED_timepointsPreOutcome  
end % pt for




%% debug
% aa = sum(IEDdata.IEDtrialsPostOnset(:));
