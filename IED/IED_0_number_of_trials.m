% CCN 2026
% finding out the number of trials 

clear;
clc;
close all;
warning('off','all');

inputFolderName_IEDtrials = ...
    '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';


fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));

PatientsNum = length(fileList);



ptIDs      = strings(PatientsNum,1);
nTrialsPts = nan(PatientsNum,1);


for pt = 1:PatientsNum
% for pt = 1:1

    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    ptIDs(pt) = ptID;

    disp("Processing patient: " + ptID);
    load(fullfile(inputFolderName_IEDtrials, fileList(pt).name));

    nTrialsPts(pt) = LFPIED.nTrials;

end

%% finding out mean and std:
mean = mean(nTrialsPts);
std = std(nTrialsPts);
