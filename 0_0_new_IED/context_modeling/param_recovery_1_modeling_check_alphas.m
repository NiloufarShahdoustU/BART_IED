clear;
clc;
close all;
warning('off','all')


%%

inputFolderName = 'D:\Nill\data\BART\0_0_new_IED\context_modeling\param_recovery_1_modeling\';



matFiles = dir(fullfile(inputFolderName, '*.mat'));
nPatients = length(matFiles);

allAlphaPos = zeros(1,nPatients);
allAlphaNeg = zeros(1,nPatients);

for pt = 1:nPatients
% for pt = 1:1

    fileName = matFiles(pt).name;
    fprintf('processing pt: %s\n', fileName);

    matFilePath = fullfile(inputFolderName, fileName);
    load(matFilePath);
    nTrials = length(TDdataParamRecovery.scoreVec);

    bestAlphaPos = TDdataParamRecovery.bestAlphaPos;
    bestAlphaNeg = TDdataParamRecovery.bestAlphaNeg;

    allAlphaPos(pt) = bestAlphaPos;
    allAlphaNeg(pt) = bestAlphaNeg;

end

%%

figure;
scatter(allAlphaPos, allAlphaNeg, 20, 'filled');
xlabel('\alpha_+');
ylabel('\alpha_-');
