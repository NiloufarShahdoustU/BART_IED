% Summarize IED rate per trial vs RT, IT, and BR
% Author: Nill

clear;
clc;
close all;

inputFolder = 'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';
outputFolder = 'D:\Nill\code\BART\IED\0_0_new_IED\IED3_summary_of_IED2_without_log\';

mkdir(outputFolder);

files = dir(fullfile(inputFolder, '*.LFPIED.mat'));
nPermutations = 10000;

perPatientResults = table();
trialLevelData = table();

for pt = 1:length(files)

    fileName = files(pt).name;
    parts = strsplit(fileName, '.');
    ptID = string(parts{1});

    disp("Processing patient: " + ptID)

    data = load(fullfile(inputFolder, fileName));
    LFPIED = data.LFPIED;

    nTrials = LFPIED.nTrials;
    RT = LFPIED.RTs(1:nTrials);
    IT = LFPIED.ITs(1:nTrials);
    BR = LFPIED.BankedTrials(1:nTrials);
    control = LFPIED.isControl(1:nTrials);

    RT = RT(:);
    IT = IT(:);
    BR = BR(:);
    control = control(:);

    RTtrials = round(LFPIED.IED_occurance_RT(:, 1));
    ITtrials = round(LFPIED.IED_occurance_IT(:, 1));

    RTtrials = RTtrials(isfinite(RTtrials) & RTtrials >= 1 & RTtrials <= nTrials);
    ITtrials = ITtrials(isfinite(ITtrials) & ITtrials >= 1 & ITtrials <= nTrials);

    nIED_RT = accumarray(RTtrials, 1, [nTrials 1]);
    nIED_IT = accumarray(ITtrials, 1, [nTrials 1]);

    IEDrate_RT = nIED_RT ./ RT;
    IEDrate_IT = nIED_IT ./ IT;

    xData = {IEDrate_RT, IEDrate_IT, IEDrate_IT};
    yData = {RT, IT, BR};
    comparison = ["IED_rate_vs_RT", "IED_rate_vs_IT", "IED_rate_vs_BR"];
    xName = ["IED rate during RT (IEDs/s)", "IED rate during IT (IEDs/s)", "IED rate during IT (IEDs/s)"];
    yName = ["RT", "IT", "BR"];

    for k = 1:3

        x = xData{k};
        y = yData{k};

        keep = control == 0 & isfinite(RT) & RT <= 20 & isfinite(x) & x >= 0;

        if k == 1
            keep = keep & RT > 0;
        elseif k == 2
            keep = keep & IT > 0;
        elseif k == 3
            keep = keep & IT > 0 & (BR == 0 | BR == 1);
        end

        x = x(keep);
        y = y(keep);
        T = table(x, y);

        oddsRatio = NaN;
        oddsLow = NaN;
        oddsHigh = NaN;
        rSquared = NaN;
        adjustedRSquared = NaN;

        if k == 3
            model = fitglm(T, 'y ~ x', 'Distribution', 'binomial', 'Link', 'logit');
            modelType = "Per-patient logistic model";
        else
            model = fitlm(T, 'y ~ x');
            rSquared = model.Rsquared.Ordinary;
            adjustedRSquared = model.Rsquared.Adjusted;
            modelType = "Per-patient linear model";
        end

        coefficients = model.Coefficients;
        CI = coefCI(model);
        slopeRow = find(strcmp(coefficients.Properties.RowNames, 'x'));
        interceptRow = find(strcmp(coefficients.Properties.RowNames, '(Intercept)'));

        intercept = coefficients.Estimate(interceptRow);
        slope = coefficients.Estimate(slopeRow);
        slopeSE = coefficients.SE(slopeRow);
        slopeT = coefficients.tStat(slopeRow);
        slopeP = coefficients.pValue(slopeRow);
        slopeLow = CI(slopeRow, 1);
        slopeHigh = CI(slopeRow, 2);

        if k == 3
            oddsRatio = exp(slope);
            oddsLow = exp(slopeLow);
            oddsHigh = exp(slopeHigh);
        end

        newRow = table(ptID, comparison(k), xName(k), yName(k), length(y), ...
            mean(x), median(x), mean(y), median(y), intercept, slope, slopeSE, ...
            slopeT, slopeP, slopeLow, slopeHigh, oddsRatio, oddsLow, oddsHigh, ...
            rSquared, adjustedRSquared, modelType, ...
            'VariableNames', {'patientID', 'comparison', 'xMeasure', 'yMeasure', ...
            'nPoints', 'rawXMean', 'rawXMedian', 'yMean', 'yMedian', ...
            'modelIntercept', 'modelSlope', 'modelSlopeSE', 'modelSlopeT', ...
            'modelSlopeP', 'modelSlopeCILow', 'modelSlopeCIHigh', 'oddsRatio', ...
            'oddsRatioCILow', 'oddsRatioCIHigh', 'rSquared', 'adjRSquared', 'modelType'});

        perPatientResults = [perPatientResults; newRow];

        n = length(y);
        newTrials = table(repmat(ptID, n, 1), repmat(comparison(k), n, 1), ...
            repmat(xName(k), n, 1), repmat(yName(k), n, 1), x, y, ...
            'VariableNames', {'patientID', 'comparison', 'xMeasure', 'yMeasure', ...
            'rawX', 'y'});

        trialLevelData = [trialLevelData; newTrials];
    end
end

writetable(perPatientResults, fullfile(outputFolder, 'per_patient_results.csv'));

%% Group summary of the patient slopes

groupSlopeSummary = table();

for k = 1:3

    rows = perPatientResults(perPatientResults.comparison == comparison(k), :);
    slopes = rows.modelSlope;
    slopes = slopes(isfinite(slopes));

    nPatients = height(rows);
    nPatientsWithSlope = length(slopes);
    totalTrialPoints = sum(rows.nPoints);

    medianSlope = median(slopes);
    slopeIQR = iqr(slopes);
    meanSlope = mean(slopes);
    slopeSE = std(slopes) / sqrt(length(slopes));
    tCritical = tinv(0.975, length(slopes) - 1);
    slopeLow = meanSlope - tCritical * slopeSE;
    slopeHigh = meanSlope + tCritical * slopeSE;

    rng(1)
    signs = randi([0 1], length(slopes), nPermutations) * 2 - 1;
    permutationMeans = mean(slopes .* signs, 1);
    permutationP = (sum(abs(permutationMeans) >= abs(meanSlope)) + 1) / (nPermutations + 1);

    medianOddsRatio = NaN;
    meanOddsRatio = NaN;

    if k == 3
        medianOddsRatio = exp(medianSlope);
        meanOddsRatio = exp(meanSlope);
    end

    newRow = table(comparison(k), xName(k), yName(k), nPatients, totalTrialPoints, ...
        nPatientsWithSlope, medianSlope, slopeIQR, meanSlope, slopeLow, slopeHigh, ...
        sum(slopes > 0), sum(slopes < 0), permutationP, medianOddsRatio, meanOddsRatio, ...
        'VariableNames', {'comparison', 'xMeasure', 'yMeasure', 'nPatients', ...
        'totalTrialPoints', 'nPatientsWithSlope', 'medianSlope', 'iqrSlope', ...
        'meanSlope', 'meanSlopeCILow', 'meanSlopeCIHigh', 'nPositiveSlope', ...
        'nNegativeSlope', 'permutationP_Slope', 'medianOddsRatio', 'meanOddsRatio'});

    groupSlopeSummary = [groupSlopeSummary; newRow];
end

% Bonferroni correction
p = groupSlopeSummary.permutationP_Slope;
groupSlopeSummary.permutationP_Slope_Bonferroni = min(p * length(p), 1);

% FDR correction
[sortedP, order] = sort(p);
fdr = sortedP .* length(p) ./ (1:length(p))';

for k = length(fdr)-1:-1:1
    fdr(k) = min(fdr(k), fdr(k + 1));
end

fdr = min(fdr, 1);
correctedP = zeros(3, 1);
correctedP(order) = fdr;
groupSlopeSummary.permutationP_Slope_FDR = correctedP;

writetable(groupSlopeSummary, fullfile(outputFolder, 'group_per_patient_slope_summary.csv'));

%% Mixed-effects models

groupMixedEffectsResults = table();

for k = 1:3

    rows = trialLevelData(trialLevelData.comparison == comparison(k), :);
    rows.patientID = categorical(rows.patientID);

    if k == 3
        model = fitglme(rows, 'y ~ rawX + (1 | patientID)', ...
            'Distribution', 'Binomial', 'Link', 'logit');
        modelType = "Logistic mixed-effects: BR ~ IED rate + (1 | patientID)";
    else
        model = fitlme(rows, 'y ~ rawX + (1 | patientID)');
        modelType = "Linear mixed-effects: y ~ IED rate + (1 | patientID)";
    end

    coefficients = model.Coefficients;
    CI = coefCI(model);
    slopeRow = find(strcmp(coefficients.Name, 'rawX'));

    beta = coefficients.Estimate(slopeRow);
    SE = coefficients.SE(slopeRow);
    tStat = coefficients.tStat(slopeRow);
    pValue = coefficients.pValue(slopeRow);
    CILow = CI(slopeRow, 1);
    CIHigh = CI(slopeRow, 2);

    oddsRatio = NaN;
    oddsLow = NaN;
    oddsHigh = NaN;

    if k == 3
        oddsRatio = exp(beta);
        oddsLow = exp(CILow);
        oddsHigh = exp(CIHigh);
    end

    newRow = table(comparison(k), xName(k), yName(k), ...
        length(unique(rows.patientID)), height(rows), beta, SE, tStat, pValue, ...
        CILow, CIHigh, oddsRatio, oddsLow, oddsHigh, modelType, ...
        'VariableNames', {'comparison', 'xMeasure', 'yMeasure', 'nPatients', ...
        'nTrialPoints', 'beta_IEDrate', 'SE', 'tStat', 'pValue', ...
        'CILow', 'CIHigh', 'oddsRatio', 'oddsRatioCILow', 'oddsRatioCIHigh', 'modelType'});

    groupMixedEffectsResults = [groupMixedEffectsResults; newRow];
end

% Bonferroni correction
p = groupMixedEffectsResults.pValue;
groupMixedEffectsResults.pValue_Bonferroni = min(p * length(p), 1);

% FDR correction
[sortedP, order] = sort(p);
fdr = sortedP .* length(p) ./ (1:length(p))';

for k = length(fdr)-1:-1:1
    fdr(k) = min(fdr(k), fdr(k + 1));
end

fdr = min(fdr, 1);
correctedP = zeros(3, 1);
correctedP(order) = fdr;
groupMixedEffectsResults.pValue_FDR = correctedP;

writetable(groupMixedEffectsResults, fullfile(outputFolder, 'group_mixed_effects_results.csv'));

%% Figure

figure('Color', 'w', 'Position', [100 100 800 700]);
hold on

colors = [0.204 0.459 0.702; 0.847 0.333 0.153; 0.25 0.60 0.25];
values = perPatientResults.modelSlope;
groups = zeros(size(values));

for k = 1:3
    groups(perPatientResults.comparison == comparison(k)) = k;
end

boxplot(values, groups, 'Labels', {'RT', 'IT', 'BR'}, 'Symbol', '', ...
    'Widths', 0.5, 'Colors', 'k');

rng(1)
xPlot = groups + (rand(size(groups)) - 0.5) * 0.32;

patients = unique(perPatientResults.patientID, 'stable');

for pt = 1:length(patients)
    patientRows = find(perPatientResults.patientID == patients(pt));
    plot(xPlot(patientRows), values(patientRows), '-', ...
        'Color', [0.80 0.80 0.80], 'LineWidth', 0.5)
end

for k = 1:3
    rows = groups == k;
    scatter(xPlot(rows), values(rows), 30, colors(k, :), 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')
end

yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);

yMin = min(values);
yMax = max(values);
yRange = yMax - yMin;
ylim([yMin - 0.18*yRange, yMax + 0.30*yRange]);

for k = 1:3
    pValue = groupMixedEffectsResults.pValue_FDR(k);
    stars = "";

    if pValue < 0.001
        stars = "***";
    elseif pValue < 0.01
        stars = "**";
    elseif pValue < 0.05
        stars = "*";
    end

    rows = groups == k;
    text(k, max(values(rows)) + 0.09*yRange, stars, ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');
end

ylabel(sprintf(['Model slope\nRT/IT: seconds per IED/s; ' ...
    'BR: log-odds per IED/s']), 'FontSize', 12, 'FontWeight', 'bold');
title('Per-patient model slope', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'LineWidth', 1.1, 'Box', 'off', 'TickDir', 'out');
pbaspect([1 1 1]);

exportgraphics(gcf, fullfile(outputFolder, 'summary_boxplots.pdf'), 'ContentType', 'vector');
close(gcf)
