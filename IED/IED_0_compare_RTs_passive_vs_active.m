% CCN 2026
% AUTHOR: Nill
% Compare ACTIVE (non-control) vs PASSIVE (control) RTs per patient

clear;
clc;
close all;
warning('off','all');

inputFolderName_IEDtrials = ...
    '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_0_compare_RTs_passive_vs_active\';

fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));

PatientsNum = length(fileList);

alpha = 0.05;

% results containers
ptIDs      = strings(PatientsNum,1);
pvals      = nan(PatientsNum,1);
isSig      = false(PatientsNum,1);
nActive    = nan(PatientsNum,1);
nPassive   = nan(PatientsNum,1);
meanActive = nan(PatientsNum,1);
meanPassive= nan(PatientsNum,1);

for pt = 1:PatientsNum

    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    ptIDs(pt) = ptID;

    disp("Processing patient: " + ptID);

    load(fullfile(inputFolderName_IEDtrials, fileList(pt).name));


    RT = LFPIED.RTs;          
    isControl = LFPIED.isControl;

    % sanity
    good = ~isnan(RT) & ~isnan(isControl);
    RT = RT(good);
    isControl = isControl(good);

    RT_active  = RT(~isControl);
    RT_passive = RT(isControl);

    nActive(pt)  = numel(RT_active);
    nPassive(pt) = numel(RT_passive);


    meanActive(pt)  = mean(RT_active);
    meanPassive(pt) = mean(RT_passive);

    % nonparametric test (robust default)
    pvals(pt) = ranksum(RT_active, RT_passive);

    isSig(pt) = pvals(pt) < alpha;

end



ResultsTable = table( ...
    ptIDs, nActive, nPassive, ...
    meanActive, meanPassive, ...
    pvals, isSig, ...
    'VariableNames', ...
    {'Patient','N_active','N_passive','MeanRT_active','MeanRT_passive','pValue','Significant'} );


%% save the table:

writetable(ResultsTable, ...
    fullfile(outputFolderName, 'Active_vs_Passive_RT_results.csv'));




