% Context-sensitive asymmetric RL fitting for BART
% Fixed: expectedReward and RPE are never saved as NaN

clear;
clc;
close all;
warning('off','all')

outputFolderName = 'D:\Nill\data\BART\0_0_new_IED\context_modeling\param_recovery_1_modeling\';
inputFolderName  = 'D:\Nill\data\BART\0_0_new_IED\context_modeling\param_recovery_0_data';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

matFiles = dir(fullfile(inputFolderName, '*.mat'));
nPatients = length(matFiles);

ptID_all = cell(nPatients,1);
%%
for pt = 67:71

    fileName = matFiles(pt).name;
    fprintf('processing pt: %s\n', fileName);

    matFilePath = fullfile(inputFolderName, fileName);
    load(matFilePath);

    nTrials = length(TDdataParamRecovery.scoreVec);

    [~, ptID, ~] = fileparts(fileName);
    ptID = strrep(ptID, '_TDdataParamRecovery', '');
    ptID_all{pt} = ptID;

    result    = TDdataParamRecovery.result(1:nTrials);
    isControl = TDdataParamRecovery.is_control(1:nTrials);
    trialType = TDdataParamRecovery.trial_type(1:nTrials);

    trialTypes = unique(trialType);
    nTrialTypes = max(trialType);

    outcomeBinary = nan(1,nTrials);
    outcomeBinary(strcmp(result,'banked')) = 1;
    outcomeBinary(strcmp(result,'popped')) = 0;

    pointsPerTrial = diff([0 TDdataParamRecovery.scoreVec]);

    Reward = zeros(1,nTrials);
    for t = 1:nTrials
        if strcmp(result(t),'banked')
            Reward(t) = pointsPerTrial(t);
        else
            Reward(t) = 0;
        end
    end

    a = 0.01:0.01:1;
    nAlpha = length(a);

    inverseTemperatureRSTD = zeros(nAlpha, nAlpha);
    fitScoreRSTD           = zeros(nAlpha, nAlpha);
    logLikRSTD             = zeros(nAlpha, nAlpha);
    logPriorAlphaPosMat    = zeros(nAlpha, nAlpha);
    logPriorAlphaNegMat    = zeros(nAlpha, nAlpha);

    rstdV   = zeros(nAlpha, nAlpha, nTrials);
    rstdRPE = zeros(nAlpha, nAlpha, nTrials);

    for ap = 1:nAlpha
        for an = 1:nAlpha

            RewardPE = zeros(nTrialTypes,nTrials);
            expectedReward = zeros(nTrialTypes,nTrials);

            expectedReward_byTrial = zeros(1,nTrials);
            RewardPE_byTrial       = zeros(1,nTrials);

            Ts = zeros(nTrialTypes,1);

            % -------------------------------
            % Trial 1: no NaN
            % -------------------------------
            currType = trialType(1);

            expectedReward_byTrial(1) = expectedReward(currType,1);
            RewardPE(currType,1) = Reward(1) - expectedReward(currType,1);
            RewardPE_byTrial(1) = RewardPE(currType,1);

            if ~isControl(1)
                if RewardPE(currType,1) > 0
                    expectedReward(currType,1) = expectedReward(currType,1) + ...
                        a(ap) * RewardPE(currType,1);
                elseif RewardPE(currType,1) < 0
                    expectedReward(currType,1) = expectedReward(currType,1) + ...
                        a(an) * RewardPE(currType,1);
                end
                Ts(currType) = Ts(currType) + 1;
            end

            % -------------------------------
            % Trials 2:nTrials
            % -------------------------------
            for t = 2:nTrials

                currType = trialType(t);

                expectedReward(:,t) = expectedReward(:,t-1);
                RewardPE(:,t) = RewardPE(:,t-1);

                expectedReward_byTrial(t) = expectedReward(currType,t-1);

                if isControl(t)
                    RewardPE_byTrial(t) = 0;
                    continue;
                end

                RewardPE(currType,t) = Reward(t) - expectedReward(currType,t-1);
                RewardPE_byTrial(t) = RewardPE(currType,t);

                alphaPosEff = a(ap);
                alphaNegEff = a(an);

                if RewardPE(currType,t) > 0
                    expectedReward(currType,t) = expectedReward(currType,t-1) + ...
                        alphaPosEff * RewardPE(currType,t);

                elseif RewardPE(currType,t) < 0
                    expectedReward(currType,t) = expectedReward(currType,t-1) + ...
                        alphaNegEff * RewardPE(currType,t);

                else
                    expectedReward(currType,t) = expectedReward(currType,t-1);
                end

                Ts(currType) = Ts(currType) + 1;
            end

            validIdx = false(1,nTrials);
            validIdx(2:end) = ~isnan(outcomeBinary(2:end)) & ...
                              ~isControl(2:end);

            whichLink = 'logit';

            if sum(validIdx) >= 5 && numel(unique(expectedReward_byTrial(validIdx))) > 1
                try
                    [B, dev] = glmfit(expectedReward_byTrial(validIdx)', ...
                                      outcomeBinary(validIdx)', ...
                                      'binomial', 'link', whichLink);

                    inverseTemperatureRSTD(ap,an) = B(2);

                    logPriorAlphaPos = log(max(betapdf(a(ap), 2, 2), realmin));
                    logPriorAlphaNeg = log(max(betapdf(a(an), 2, 2), realmin));

                    logPriorAlphaPosMat(ap,an) = logPriorAlphaPos;
                    logPriorAlphaNegMat(ap,an) = logPriorAlphaNeg;

                    approxLogLik = -0.5 * dev;
                    logLikRSTD(ap,an) = approxLogLik;

                    fitScoreRSTD(ap,an) = approxLogLik + ...
                                          logPriorAlphaPos + ...
                                          logPriorAlphaNeg;

                catch
                    inverseTemperatureRSTD(ap,an) = NaN;
                    fitScoreRSTD(ap,an) = NaN;
                    logLikRSTD(ap,an) = NaN;
                    logPriorAlphaPosMat(ap,an) = NaN;
                    logPriorAlphaNegMat(ap,an) = NaN;
                end
            else
                inverseTemperatureRSTD(ap,an) = NaN;
                fitScoreRSTD(ap,an) = NaN;
                logLikRSTD(ap,an) = NaN;
                logPriorAlphaPosMat(ap,an) = NaN;
                logPriorAlphaNegMat(ap,an) = NaN;
            end

            % safety check: replace any accidental NaN with 0
            expectedReward_byTrial(isnan(expectedReward_byTrial)) = 0;
            RewardPE_byTrial(isnan(RewardPE_byTrial)) = 0;

            rstdV(ap,an,:)   = expectedReward_byTrial;
            rstdRPE(ap,an,:) = RewardPE_byTrial;

        end
    end

    if all(isnan(fitScoreRSTD(:)))
        bestApIdx = NaN;
        bestAnIdx = NaN;
        bestAlphaPos = NaN;
        bestAlphaNeg = NaN;
        bestFitScore = NaN;
        bestInverseTemperature = NaN;
        bestLogLik = NaN;
    else
        [~, bestIdx] = max(fitScoreRSTD(:));
        [bestApIdx, bestAnIdx] = ind2sub(size(fitScoreRSTD), bestIdx);

        bestAlphaPos = a(bestApIdx);
        bestAlphaNeg = a(bestAnIdx);
        bestFitScore = fitScoreRSTD(bestApIdx,bestAnIdx);
        bestInverseTemperature = inverseTemperatureRSTD(bestApIdx,bestAnIdx);
        bestLogLik = logLikRSTD(bestApIdx,bestAnIdx);
    end

    TDdataParamRecovery.a = a;
    TDdataParamRecovery.expectedReward = rstdV;
    TDdataParamRecovery.RPE = rstdRPE;

    TDdataParamRecovery.bestAlphaPos = bestAlphaPos;
    TDdataParamRecovery.bestAlphaNeg = bestAlphaNeg;

    TDdataParamRecovery.result = result;
    TDdataParamRecovery.Reward = Reward;
    TDdataParamRecovery.inverseTemperatureRSTD = inverseTemperatureRSTD;

    TDdataParamRecovery.fitScoreRSTD = fitScoreRSTD;
    TDdataParamRecovery.bestFitScore = bestFitScore;
    TDdataParamRecovery.bestInverseTemperature = bestInverseTemperature;
    TDdataParamRecovery.bestApIdx = bestApIdx;
    TDdataParamRecovery.bestAnIdx = bestAnIdx;
    TDdataParamRecovery.trialTypes = trialTypes;

    TDdataParamRecovery.logLikRSTD = logLikRSTD;
    TDdataParamRecovery.bestLogLik = bestLogLik;
    TDdataParamRecovery.logPriorAlphaPos = logPriorAlphaPosMat;
    TDdataParamRecovery.logPriorAlphaNeg = logPriorAlphaNegMat;

    fprintf('best alpha+ = %.2f, alpha- = %.2f\n', bestAlphaPos, bestAlphaNeg);

    save(fullfile(outputFolderName, [ptID '_TDdataParamRecovery.mat']), ...
        'TDdataParamRecovery');
end