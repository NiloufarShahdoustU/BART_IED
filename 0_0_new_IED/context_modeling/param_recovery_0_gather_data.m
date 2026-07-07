
clear;
clc;
close all;
%% loading gradient data

outputFolderName = 'D:\Nill\data\BART\0_0_new_IED\context_modeling\param_recovery_0_data\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

inputFolder = 'D:\Nill\data\BART_preprocessed';

allItems = dir(inputFolder);
allItems = allItems([allItems.isdir]);
allItems = allItems(~ismember({allItems.name}, {'.', '..'}));

%  only folder names that start with '20'
folderMask = startsWith({allItems.name}, '20');
targetFolders = allItems(folderMask);
nPatients = length(targetFolders);
%%
for pt = 75:79
% for pt = 2:2

    ptID = targetFolders(pt).name;
    fprintf('processing pt: %s\n', ptID);

    dataFolder = fullfile(inputFolder, ptID, 'Data');
    matFile = dir(fullfile(dataFolder, '*TDdataGradients.mat'));

    if isempty(matFile)
        matFile = dir(fullfile(dataFolder, '*bartBHV.mat'));
        matPath = fullfile(dataFolder, matFile(1).name);
        TDdata = load(matPath);
        scoreVec = [TDdata.data.score];
        TDdataParamRecovery.scoreVec = scoreVec;
        pointsVec = [TDdata.data.points];
        Reward = [pointsVec(1), diff(scoreVec)];
        TDdataParamRecovery.Reward = Reward;
        TDdataParamRecovery.result = {TDdata.data.result};
        pointsMinusReward = pointsVec - Reward;
        TDdataParamRecovery.pointsMinusReward = pointsMinusReward;
        TDdataParamRecovery.inflate_time = [TDdata.data.inflate_time];
        TDdataParamRecovery.points = [TDdata.data.points];
        TDdataParamRecovery.is_control = [TDdata.data.is_control];
        TDdataParamRecovery.trial_type = [TDdata.data.trial_type];

    else
        matPath = fullfile(dataFolder, matFile(1).name);
        TDdata = load(matPath);

        TDdataParamRecovery.a = TDdata.a;
        TDdataParamRecovery.nTrials = TDdata.nTrials;
        scoreVec = [TDdata.data.score];
        TDdataParamRecovery.scoreVec = scoreVec;
        pointsVec = [TDdata.data.points];
        Reward = [pointsVec(1), diff(scoreVec)];
        TDdataParamRecovery.Reward = Reward;
        TDdataParamRecovery.result = {TDdata.data.result};
        pointsMinusReward = pointsVec - Reward;
        TDdataParamRecovery.pointsMinusReward = pointsMinusReward;
        TDdataParamRecovery.inflate_time = [TDdata.data.inflate_time];
        TDdataParamRecovery.points = [TDdata.data.points];
        TDdataParamRecovery.is_control = [TDdata.data.is_control];
        TDdataParamRecovery.trial_type = [TDdata.data.trial_type];
        TDdataParamRecovery.RPE = TDdata.TDdataGradients.rstdRPE;
        TDdataParamRecovery.expectedReward = TDdata.TDdataGradients.rstdV;
        TDdataParamRecovery.bestAlphaPos = TDdata.TDdataGradients.neuralFit.bestAlphaPositive;
        TDdataParamRecovery.bestAlphaNeg = TDdata.TDdataGradients.neuralFit.bestAlphaNegative;
            
    end



    save(fullfile(outputFolderName, [ptID '_TDdataParamRecovery.mat']), 'TDdataParamRecovery');

end


