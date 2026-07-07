% Test whether longer inter-pump time (IT) is associated with
% a higher or lower probability of banking across participants.
%
% Main analysis:
%   Logistic mixed-effects model with participant-specific random intercept:
%       BR ~ IT_within + meanIT_between + (1 | patientID)
%
% IT_within:
%   Trial IT minus that participant's mean IT.
%   This is the main effect of interest:
%       beta < 0  -> longer-than-usual IT predicts LOWER banking probability
%       beta > 0  -> longer-than-usual IT predicts HIGHER banking probability
%
% meanIT_between:
%   Participant mean IT, centered around the group mean.
%   This tests whether generally slower participants bank more or less often.
%
% A random-slope model is attempted first:
%       BR ~ IT_within + meanIT_between + (1 + IT_within | patientID)
% If it fails, the code automatically uses the random-intercept model.


%%%%% Do participants with a higher average IT tend to have a higher or lower bank rate than participants with a lower average IT?

clear;
clc;
close all;

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED3_IT_vs_BR_mixed_effects\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

trialLevelData = table();
perPatientResults = table();

minimumTrialsPerPatientModel = 10;



for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    ptID = string(fileNameParts{1});

    fprintf('\nProcessing patient: %s\n', ptID);

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName));

    if ~isfield(loadedData, 'LFPIED')
        fprintf('Skipped: LFPIED structure was not found.\n');
        continue;
    end

    LFPIED = loadedData.LFPIED;

    requiredFields = {'RTs', 'ITs', 'isControl', 'BankedTrials'};

    missingField = false;

    for ff = 1:length(requiredFields)
        if ~isfield(LFPIED, requiredFields{ff})
            fprintf('Skipped: missing field %s.\n', requiredFields{ff});
            missingField = true;
        end
    end

    if missingField
        continue;
    end

    RTs = LFPIED.RTs(:);
    ITs = LFPIED.ITs(:);
    isControl = LFPIED.isControl(:);
    BRs = LFPIED.BankedTrials(:);

    if isfield(LFPIED, 'nTrials') && isfinite(LFPIED.nTrials)
        nTrials = LFPIED.nTrials;
    else
        nTrials = min([length(RTs), length(ITs), ...
            length(isControl), length(BRs)]);
    end

    nTrials = min([nTrials, length(RTs), length(ITs), ...
        length(isControl), length(BRs)]);

    RTs = RTs(1:nTrials);
    ITs = ITs(1:nTrials);
    isControl = isControl(1:nTrials);
    BRs = BRs(1:nTrials);

    trialNumber = (1:nTrials)';

    keepIdx = ...
        isControl == 0 & ...
        isfinite(RTs) & ...
        RTs > 0 & ...
        RTs <= 20 & ...
        isfinite(ITs) & ...
        ITs > 0 & ...
        isfinite(BRs) & ...
        (BRs == 0 | BRs == 1);

    patientIT = ITs(keepIdx);
    patientBR = BRs(keepIdx);
    patientTrialNumber = trialNumber(keepIdx);

    nValidTrials = length(patientIT);

    if nValidTrials == 0
        fprintf('No valid trials after filtering.\n');
        continue;
    end

    newRows = table( ...
        repmat(ptID, nValidTrials, 1), ...
        patientTrialNumber, ...
        patientIT, ...
        patientBR, ...
        'VariableNames', { ...
            'patientID', ...
            'trialNumber', ...
            'IT', ...
            'BR' ...
        });

    trialLevelData = [trialLevelData; newRows];

end

if isempty(trialLevelData)
    error('No valid trial-level data were found.');
end

% Create within-participant and between-participant IT predictors

patientList = unique(trialLevelData.patientID, 'stable');

trialLevelData.meanIT_patient = NaN(height(trialLevelData), 1);
trialLevelData.IT_within = NaN(height(trialLevelData), 1);

for pp = 1:length(patientList)

    idx = trialLevelData.patientID == patientList(pp);

    patientMeanIT = mean(trialLevelData.IT(idx), 'omitnan');

    trialLevelData.meanIT_patient(idx) = patientMeanIT;
    trialLevelData.IT_within(idx) = ...
        trialLevelData.IT(idx) - patientMeanIT;

end

grandMeanIT = mean(trialLevelData.meanIT_patient, 'omitnan');

trialLevelData.meanIT_between = ...
    trialLevelData.meanIT_patient - grandMeanIT;

trialLevelData.patientID = categorical(trialLevelData.patientID);

% Per-patient logistic regressions: BR ~ IT

patientCategories = categories(trialLevelData.patientID);

for pp = 1:length(patientCategories)

    thisPatient = patientCategories{pp};

    rows = trialLevelData( ...
        trialLevelData.patientID == thisPatient, :);

    nTrials = height(rows);
    nBanked = sum(rows.BR == 1);
    nNotBanked = sum(rows.BR == 0);
    bankRate = mean(rows.BR, 'omitnan');
    meanIT = mean(rows.IT, 'omitnan');

    betaIT = NaN;
    seIT = NaN;
    pValueIT = NaN;
    ciLowIT = NaN;
    ciHighIT = NaN;
    oddsRatioIT = NaN;
    oddsRatioCILow = NaN;
    oddsRatioCIHigh = NaN;
    modelStatus = "Not enough data";

    canFitPatientModel = ...
        nTrials >= minimumTrialsPerPatientModel && ...
        nBanked > 0 && ...
        nNotBanked > 0 && ...
        length(unique(rows.IT)) > 1;

    if canFitPatientModel

        try

            patientGLM = fitglm( ...
                rows, ...
                'BR ~ IT', ...
                'Distribution', 'binomial', ...
                'Link', 'logit');

            coefficientTable = patientGLM.Coefficients;
            confidenceIntervals = coefCI(patientGLM);

            slopeIdx = find(strcmp( ...
                coefficientTable.Properties.RowNames, 'IT'));

            betaIT = coefficientTable.Estimate(slopeIdx);
            seIT = coefficientTable.SE(slopeIdx);
            pValueIT = coefficientTable.pValue(slopeIdx);
            ciLowIT = confidenceIntervals(slopeIdx, 1);
            ciHighIT = confidenceIntervals(slopeIdx, 2);

            oddsRatioIT = exp(betaIT);
            oddsRatioCILow = exp(ciLowIT);
            oddsRatioCIHigh = exp(ciHighIT);

            modelStatus = "Per-patient logistic model fitted";

        catch
            modelStatus = "Per-patient logistic model failed";
        end

    elseif nBanked == 0 || nNotBanked == 0

        modelStatus = "No within-patient BR variation";

    end

    newPatientRow = table( ...
        string(thisPatient), ...
        nTrials, ...
        nBanked, ...
        nNotBanked, ...
        bankRate, ...
        meanIT, ...
        betaIT, ...
        seIT, ...
        pValueIT, ...
        ciLowIT, ...
        ciHighIT, ...
        oddsRatioIT, ...
        oddsRatioCILow, ...
        oddsRatioCIHigh, ...
        modelStatus, ...
        'VariableNames', { ...
            'patientID', ...
            'nTrials', ...
            'nBanked', ...
            'nNotBanked', ...
            'bankRate', ...
            'meanIT_seconds', ...
            'beta_IT', ...
            'SE_IT', ...
            'pValue_IT', ...
            'CILow_IT', ...
            'CIHigh_IT', ...
            'oddsRatio_per_1s_IT', ...
            'oddsRatioCILow', ...
            'oddsRatioCIHigh', ...
            'modelStatus' ...
        });

    perPatientResults = [perPatientResults; newPatientRow];

end

% Group logistic mixed-effects model

modelFormulaRandomSlope = ...
    'BR ~ IT_within + meanIT_between + (1 + IT_within | patientID)';

modelFormulaRandomIntercept = ...
    'BR ~ IT_within + meanIT_between + (1 | patientID)';

usedFormula = "";
modelStatus = "";
groupModel = [];

% Try the more complete random-slope model first.
try

    groupModel = fitglme( ...
        trialLevelData, ...
        modelFormulaRandomSlope, ...
        'Distribution', 'Binomial', ...
        'Link', 'logit', ...
        'FitMethod', 'Laplace');

    usedFormula = modelFormulaRandomSlope;
    modelStatus = "Random-intercept and random-slope GLME";

catch

    % Fallback to a random-intercept model if the random-slope model
    % cannot converge or is not identifiable.
    groupModel = fitglme( ...
        trialLevelData, ...
        modelFormulaRandomIntercept, ...
        'Distribution', 'Binomial', ...
        'Link', 'logit', ...
        'FitMethod', 'Laplace');

    usedFormula = modelFormulaRandomIntercept;
    modelStatus = "Random-intercept GLME fallback";

end

coefficientTable = groupModel.Coefficients;
confidenceIntervals = coefCI(groupModel);

groupResults = table();

predictorNames = ["(Intercept)", "IT_within", "meanIT_between"];

for ii = 1:length(predictorNames)

    predictorName = predictorNames(ii);

    coefficientIdx = find(strcmp( ...
        string(coefficientTable.Name), predictorName));

    if isempty(coefficientIdx)
        continue;
    end

    beta = coefficientTable.Estimate(coefficientIdx);
    se = coefficientTable.SE(coefficientIdx);
    tStat = coefficientTable.tStat(coefficientIdx);
    pValue = coefficientTable.pValue(coefficientIdx);
    ciLow = confidenceIntervals(coefficientIdx, 1);
    ciHigh = confidenceIntervals(coefficientIdx, 2);

    oddsRatio = exp(beta);
    oddsRatioCILow = exp(ciLow);
    oddsRatioCIHigh = exp(ciHigh);

    if predictorName == "IT_within"
        interpretation = ...
            "Primary: effect of being 1 second slower than the participant's own mean IT";
    elseif predictorName == "meanIT_between"
        interpretation = ...
            "Between-participant: effect of a 1-second higher participant mean IT";
    else
        interpretation = "Model intercept";
    end

    newGroupRow = table( ...
        predictorName, ...
        beta, ...
        se, ...
        tStat, ...
        pValue, ...
        ciLow, ...
        ciHigh, ...
        oddsRatio, ...
        oddsRatioCILow, ...
        oddsRatioCIHigh, ...
        interpretation, ...
        string(modelStatus), ...
        string(usedFormula), ...
        height(trialLevelData), ...
        length(patientCategories), ...
        'VariableNames', { ...
            'predictor', ...
            'beta_logOdds', ...
            'SE', ...
            'tStat', ...
            'pValue', ...
            'CILow', ...
            'CIHigh', ...
            'oddsRatio', ...
            'oddsRatioCILow', ...
            'oddsRatioCIHigh', ...
            'interpretation', ...
            'modelStatus', ...
            'modelFormula', ...
            'nTrials', ...
            'nPatients' ...
        });

    groupResults = [groupResults; newGroupRow];

end

% Print the primary answer

primaryRow = groupResults(groupResults.predictor == "IT_within", :);

fprintf('\n');


fprintf('Beta = %.6f log-odds per 1 second\n', primaryRow.beta_logOdds);
fprintf('Odds ratio = %.6f per 1 second\n', primaryRow.oddsRatio);
fprintf('95%% CI for odds ratio = [%.6f, %.6f]\n', ...
    primaryRow.oddsRatioCILow, primaryRow.oddsRatioCIHigh);
fprintf('p = %.6g\n', primaryRow.pValue);

if primaryRow.pValue < 0.05

    if primaryRow.beta_logOdds > 0

        fprintf(['Conclusion: Longer-than-usual IT is significantly ' ...
            'associated with a HIGHER probability of banking.\n']);

    elseif primaryRow.beta_logOdds < 0

        fprintf(['Conclusion: Longer-than-usual IT is significantly ' ...
            'associated with a LOWER probability of banking.\n']);

    end

else

    fprintf(['Conclusion: There is no statistically significant evidence ' ...
        'that IT is associated with banking probability.\n']);

end



% Visualization


validSlopeRows = isfinite(perPatientResults.beta_IT);
slopeValues = perPatientResults.beta_IT(validSlopeRows);

fig = figure('Visible', 'off', 'Color', 'w');
set(fig, 'Position', [100 100 620 620]);

ax = axes(fig);
hold(ax, 'on');

if ~isempty(slopeValues)

    boxplot( ...
        ax, ...
        slopeValues, ...
        ones(size(slopeValues)), ...
        'Labels', {'IT vs BR'}, ...
        'Symbol', '', ...
        'Widths', 0.32, ...
        'Colors', 'k');

    set(findobj(ax, 'Tag', 'Box'), ...
        'Color', 'k', 'LineWidth', 0.5);

    set(findobj(ax, 'Tag', 'Median'), ...
        'Color', 'k', 'LineWidth', 0.75);

    set(findobj(ax, 'Tag', 'Whisker'), ...
        'Color', 'k', 'LineWidth', 0.5);

    set(findobj(ax, 'Tag', 'Upper Adjacent Value'), ...
        'Color', 'k', 'LineWidth', 0.5);

    set(findobj(ax, 'Tag', 'Lower Adjacent Value'), ...
        'Color', 'k', 'LineWidth', 0.5);

    rng(1);
    jitter = (rand(size(slopeValues)) - 0.5) * 0.28;

    scatter( ...
        ax, ...
        1 + jitter, ...
        slopeValues, ...
        34, ...
        'MarkerFaceColor', [0.25 0.60 0.25], ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.60);

end

yline(ax, 0, '--', ...
    'Color', [0.45 0.45 0.45], ...
    'LineWidth', 1.0);

ylabel(ax, ...
    'Per-patient logistic slope: BR ~ IT (log-odds per second)', ...
    'FontSize', 11, ...
    'FontWeight', 'bold');

title(ax, ...
    sprintf('Group IT effect: OR = %.3f, p = %.3g', ...
    primaryRow.oddsRatio, primaryRow.pValue), ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

set(ax, ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'Box', 'off', ...
    'TickDir', 'out');

grid(ax, 'off');
pbaspect(ax, [1 1 1]);

hold(ax, 'off');

summaryPDF = fullfile(outputFolderName, 'IT_BR_summary.pdf');
exportgraphics(fig, summaryPDF, 'ContentType', 'vector');
close(fig);

% Save 

trialOutputFile = fullfile( ...
    outputFolderName, 'trial_level_IT_BR.csv');

patientOutputFile = fullfile( ...
    outputFolderName, 'per_patient_IT_BR_results.csv');

groupOutputFile = fullfile( ...
    outputFolderName, 'group_mixed_effects_IT_BR_results.csv');

writetable(trialLevelData, trialOutputFile);
writetable(perPatientResults, patientOutputFile);
writetable(groupResults, groupOutputFile);

fprintf('\nSaved:\n');
fprintf('%s\n', trialOutputFile);
fprintf('%s\n', patientOutputFile);
fprintf('%s\n', groupOutputFile);
fprintf('%s\n', summaryPDF);
