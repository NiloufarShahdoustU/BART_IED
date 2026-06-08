% This code is almost junk but I will keep it anyway.
% Author: Nill

clear;
clc;
close all;
warning('off','all');

%% loading neural and event data
 
% getting the numbers of different patients

inputFolderName_bhvStruct = '\\155.100.91.44\d\Data\Nill\BART\bhvStruct_Nill_made';
inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\IEDdata_RTs';
fileList = dir(fullfile(inputFolderName_bhvStruct, '*.bhvStruct.mat'));
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\IEDTrials_ONLY_pre_peri_post\'; 
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

    matFile_data = [inputFolderName_IEDdata '\' ptID '.IEDdata.mat'];
    load(matFile_data);

    ReactionTimes = bhvStruct.allRTs;


    ReactionTimesFiltered = ReactionTimes; % I want to have all the RTs


    IEDdata_pre = IEDdata.IEDtrials_preStimRT;
    IEDdata_peri = IEDdata.IEDtrials_periStimRT;
    IEDdata_post = IEDdata.IEDtrials_postStimRT;



    nTrials = size(IEDdata_pre,2); % this is different from the 
    PrePeriPostTrials = zeros(3, nTrials);

    for trial = 1:nTrials

        if sum(~isnan(IEDdata_pre(:,trial))) > 0
                PrePeriPostTrials(1,trial)= 1;
        end

        if sum(~isnan(IEDdata_peri(:,trial))) > 0
            PrePeriPostTrials(2,trial)= 1;
        end

        if sum(~isnan(IEDdata_post(:,trial)))>0
            PrePeriPostTrials(3,trial)= 1;
        end
    end


    
   
    IEDdataPrePeriPost.selectedChans = IEDdata.selectedChans;
    IEDdataPrePeriPost.anatomicalLocs = IEDdata.anatomicalLocs;
    IEDdataPrePeriPost.PrePeriPostTrials = PrePeriPostTrials;
    save([outputFolderName ptID '.IEDdata.mat'],'IEDdataPrePeriPost');
 
    clear matFile_data IEDdata_pre IEDdata_peri IEDdata_post IEDdataPrePeriPost PrePeriPostTrials

  
end

