% Context-sensitive asymmetric RL fitting for BART
% Revised to be closer to Niv et al. (2012):
% - decaying additive learning-rate term d = 0.5/(1 + Ts)
%   where Ts is the number of previous experiences for the current trial type
%
% NOTE:
% This version still uses glmfit on banked/popped outcome
% It is therefore NOT an exact reproduction of the niv's paper's softmax
% trial-by-trial choice-likelihood fit, but it fixes the main alpha issues
% within your current framework.

clear;
clc;
close all;
warning('off','all')


inputFolderName = 'D:\Nill\data\BART\0_0_new_IED\context_modeling\param_recovery_3_simulated_fields';
outputFolderName = fullfile(pwd, 'param_recovery_4_param_recovery');

if ~exist(outputFolderName, 'dir')
    [status,msg] = mkdir(outputFolderName);
    if ~status
        error('Could not create output folder: %s', msg);
    end
end

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

matFiles = dir(fullfile(inputFolderName, '*.mat'));
nPatients = length(matFiles);

ptID_all = cell(nPatients,1);
sim_alpha_plus_all = zeros(nPatients,1);
sim_alpha_minus_all = zeros(nPatients,1);
fit_alpha_plus_all = zeros(nPatients,1);
fit_alpha_minus_all = zeros(nPatients,1);

for pt = 1:nPatients
% for pt = 1:1

    fileName = matFiles(pt).name;
    fprintf('processing pt: %s\n', fileName);

    matFilePath = fullfile(inputFolderName, fileName);
    load(matFilePath);

    nTrials = length(TDdataParamRecovery.scoreVec)-1;

    [~, ptID, ~] = fileparts(fileName);
    ptID = strrep(ptID, '_TDdataParamRecovery', '');
    ptID_all{pt} = ptID;

    result = TDdataParamRecovery.resultSimulated(1:nTrials);
    prevReward = TDdataParamRecovery.rewardSimulated(1:nTrials);
    isControl = TDdataParamRecovery.is_control(1:nTrials);

    % binary outcome for glmfit: 1 = banked, 0 = popped
    outcomeBinary = nan(1,nTrials);
    outcomeBinary(strcmp(result,'banked')) = 1;
    outcomeBinary(strcmp(result,'popped')) = 0;

    trialType = TDdataParamRecovery.trial_type(1:nTrials);
    trialTypes = unique(trialType);

    % keep indexing safe if types are 1,2,3,4
    nTrialTypes = max(trialType);

    % alpha grid
    a = 0:0.01:1;
    nAlpha = length(a);

    % store results for this patient
    inverseTemperatureRSTD = nan(nAlpha, nAlpha);   % keep B(2) from glmfit for compatibility
    fitScoreRSTD           = nan(nAlpha, nAlpha);   % posterior-like score
    logLikRSTD             = nan(nAlpha, nAlpha);   % approximate data log-likelihood from glmfit deviance
    logPriorAlphaPosMat    = nan(nAlpha, nAlpha);
    logPriorAlphaNegMat    = nan(nAlpha, nAlpha);
    rstdV                  = nan(nAlpha, nAlpha, nTrials);   % expected reward by trial
    rstdRPE                = nan(nAlpha, nAlpha, nTrials);   % reward PE by trial



    Reward = zeros(1,nTrials);
    for t = 1:nTrials
        if strcmp(result(t),'banked')
            Reward(t) = prevReward(t);
        else
            Reward(t) = 0;
        end
    end




    for ap = 1:nAlpha
        for an = 1:nAlpha

            % ============================================================
            % initialize state separately for each alpha pair
            % ============================================================
            RewardPE = zeros(nTrialTypes,nTrials);
            expectedReward = zeros(nTrialTypes,nTrials);

            % expected reward before current outcome
            expectedReward_byTrial = nan(1,nTrials);
            RewardPE_byTrial = nan(1,nTrials);

            % kept only so your structure stays similar
            Ts = zeros(nTrialTypes,1); 

            % no prior expectation on first trial
            expectedReward_byTrial(1) = NaN;
            RewardPE_byTrial(1) = NaN;

            for t = 2:nTrials
                currType = trialType(t);

                % carry forward latent states by default
                expectedReward(:,t) = expectedReward(:,t-1);
                RewardPE(:,t) = RewardPE(:,t-1);

                % do not learn from control trials
                % do not put fake repeated predictor values into trial-level fit arrays
                if isControl(t)
                    expectedReward_byTrial(t) = NaN;
                    RewardPE_byTrial(t) = NaN;
                    continue;
                end

                % predictor BEFORE observing current outcome
                expectedReward_byTrial(t) = expectedReward(currType,t-1);

                RewardPE(currType,t) = Reward(t) - expectedReward(currType,t-1);
                RewardPE_byTrial(t) = RewardPE(currType,t);

                % ========================================================
                % old:
                %   d = 0.5 / (1 + Ts(currType));
                %   alphaPosEff = min(a(ap) + d, 1);
                %   alphaNegEff = min(a(an) + d, 1);
                % ========================================================
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

            % ============================================================

            % exclude control trials from the GLM fit
            % ============================================================
            validIdx = false(1,nTrials);
            validIdx(2:end) = ~isnan(expectedReward_byTrial(2:end)) & ...
                              ~isnan(outcomeBinary(2:end)) & ...
                              ~isControl(2:end);

            whichLink = 'logit';

            if sum(validIdx) >= 5 && numel(unique(expectedReward_byTrial(validIdx))) > 1
                try
                    [B, dev] = glmfit(expectedReward_byTrial(validIdx)', ...
                                      outcomeBinary(validIdx)', ...
                                      'binomial', 'link', whichLink);

                    % keep this for compatibility with your existing outputs
                    inverseTemperatureRSTD(ap,an) = B(2);

                    % ====================================================
                    % Beta(2,2) priors on alpha+ and alpha-
                    % ====================================================
                    logPriorAlphaPos = log(max(betapdf(a(ap), 2, 2), realmin));
                    logPriorAlphaNeg = log(max(betapdf(a(an), 2, 2), realmin));

                    logPriorAlphaPosMat(ap,an) = logPriorAlphaPos;
                    logPriorAlphaNegMat(ap,an) = logPriorAlphaNeg;

                    % approximate relative log-likelihood from deviance
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



    fprintf('best alpha+ = %.2f, alpha- = %.2f\n', bestAlphaPos, bestAlphaNeg)
    sim_alpha_plus = bestAlphaPos;
    sim_alpha_minus = bestAlphaNeg;

    fit_alpha_plus = TDdataParamRecovery.bestAlphaPos;
    fit_alpha_minus = TDdataParamRecovery.bestAlphaNeg;

    sim_alpha_plus_all(pt) = sim_alpha_plus;
    sim_alpha_minus_all(pt) = sim_alpha_minus;
    fit_alpha_plus_all(pt) = fit_alpha_plus;
    fit_alpha_minus_all(pt) = fit_alpha_minus;
end


T = table(ptID_all, sim_alpha_plus_all, sim_alpha_minus_all, ...
          fit_alpha_plus_all, fit_alpha_minus_all, ...
          'VariableNames', {'ptID','sim_alpha_plus','sim_alpha_minus','fit_alpha_plus','fit_alpha_minus'});

writetable(T, fullfile(outputFolderName, 'alpha_comparison.csv'));



%% debug plots for last processed patient
% close all;
%
% if ~all(isnan(fitScoreRSTD(:)))
%     [~, bestIdx] = max(fitScoreRSTD(:));
%     [bestApIdx, bestAnIdx] = ind2sub(size(fitScoreRSTD), bestIdx);
%
%     bestExpectedReward_byTrial = squeeze(rstdV(bestApIdx,bestAnIdx,:));
%     bestRewardPE_byTrial = squeeze(rstdRPE(bestApIdx,bestAnIdx,:));
%
%     figure;
%     plot(bestExpectedReward_byTrial,'r','LineWidth',2);
%     hold on;
%     plot(prevReward,'b','LineWidth',2);
%     legend('expectedReward\_byTrial (best fit)', 'prevReward');
%     title('Best-fit expected reward by trial vs previous reward');
%
%     figure;
%     plot(abs(bestRewardPE_byTrial),'LineWidth',2);
%     title('Absolute reward prediction error by trial (best fit)');
%
%     figure;
%     imagesc(a, a, fitScoreRSTD);
%     axis xy
%     xlabel('\alpha^-')
%     ylabel('\alpha^+')
%     colorbar
%     title('Posterior-like fit surface')
%
%     figure;
%     imagesc(a, a, logLikRSTD);
%     axis xy
%     xlabel('\alpha^-')
%     ylabel('\alpha^+')
%     colorbar
%     title('Approximate log-likelihood surface')
% end