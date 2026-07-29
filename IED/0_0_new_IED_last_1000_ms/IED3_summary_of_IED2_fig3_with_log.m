% summarize ied count per trial vs rt, it, bank, and pop
% author: nill

clear;
clc;
close all;

inputFolder = 'D:\Nill\data\BART\0_0_new_IED_last_1000_ms\IED1_find_number_of_IEDs\';
outputFolder = 'D:\Nill\code\BART\IED\0_0_new_IED_last_1000_ms\IED3_summary_of_IED2_with_log\';

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

    % fit three unique models
    % the pr effect is created by reversing the br outcome
    % the br model also controls for it duration
    xData = {nIED_RT, nIED_IT, nIED_IT};
    yData = {RT, IT, BR};
    comparisonBase = ["IED_count_vs_RT", "IED_count_vs_IT", "IED_count_vs_BR"];
    xNameBase = ["IED count during RT", ...
        "IED count during IT", ...
        "IED count during IT"];
    yNameBase = ["RT", "IT", "BR"];

    for k = 1:3

        x = xData{k};
        y = yData{k};
        duration = zeros(size(y));

        keep = control == 0 & isfinite(RT) & RT <= 20 & isfinite(x) & x >= 0;

        if k == 1
            keep = keep & RT > 0;
        elseif k == 2
            keep = keep & IT > 0;
        elseif k == 3
            keep = keep & isfinite(IT) & IT > 0 & (BR == 0 | BR == 1);
            duration = IT;
        end

        x = x(keep);
        y = y(keep);
        duration = duration(keep);
        logX = log10(x + 1);

        oddsRatio = NaN;
        oddsLow = NaN;
        oddsHigh = NaN;
        rSquared = NaN;
        adjustedRSquared = NaN;

        if k == 3
            T = table(logX, duration, y);
            model = fitglm(T, 'y ~ logX + duration', ...
                'Distribution', 'binomial', 'Link', 'logit');
            modelType = "Per-patient logistic model adjusted for IT";
        else
            T = table(logX, y);
            model = fitlm(T, 'y ~ logX');
            rSquared = model.Rsquared.Ordinary;
            adjustedRSquared = model.Rsquared.Adjusted;
            modelType = "Per-patient linear model";
        end

        coefficients = model.Coefficients;
        CI = coefCI(model);
        slopeRow = find(strcmp(coefficients.Properties.RowNames, 'logX'));
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

        newRow = table(ptID, comparisonBase(k), xNameBase(k), yNameBase(k), length(y), ...
            mean(x), median(x), mean(y), median(y), intercept, slope, slopeSE, ...
            slopeT, slopeP, slopeLow, slopeHigh, oddsRatio, oddsLow, oddsHigh, ...
            rSquared, adjustedRSquared, modelType, ...
            'VariableNames', {'patientID', 'comparison', 'xMeasure', 'yMeasure', ...
            'nPoints', 'rawXMean', 'rawXMedian', 'yMean', 'yMedian', ...
            'modelIntercept', 'modelSlope', 'modelSlopeSE', 'modelSlopeT', ...
            'modelSlopeP', 'modelSlopeCILow', 'modelSlopeCIHigh', 'oddsRatio', ...
            'oddsRatioCILow', 'oddsRatioCIHigh', 'rSquared', 'adjRSquared', 'modelType'});

        perPatientResults = [perPatientResults; newRow];

        % pr = 1 - br, so its ied effect has the opposite sign
        if k == 3
            prRow = newRow;
            prRow.comparison = "IED_count_vs_PR";
            prRow.yMeasure = "PR";
            prRow.yMean = mean(1 - y);
            prRow.yMedian = median(1 - y);
            prRow.modelIntercept = -intercept;
            prRow.modelSlope = -slope;
            prRow.modelSlopeSE = slopeSE;
            prRow.modelSlopeT = -slopeT;
            prRow.modelSlopeP = slopeP;
            prRow.modelSlopeCILow = -slopeHigh;
            prRow.modelSlopeCIHigh = -slopeLow;
            prRow.oddsRatio = exp(-slope);
            prRow.oddsRatioCILow = exp(-slopeHigh);
            prRow.oddsRatioCIHigh = exp(-slopeLow);
            prRow.modelType = "Per-patient logistic model: PR = 1 - BR, adjusted for IT";

            perPatientResults = [perPatientResults; prRow];
        end

        n = length(y);
        newTrials = table(repmat(ptID, n, 1), repmat(comparisonBase(k), n, 1), ...
            repmat(xNameBase(k), n, 1), repmat(yNameBase(k), n, 1), ...
            x, logX, duration, y, ...
            'VariableNames', {'patientID', 'comparison', 'xMeasure', 'yMeasure', ...
            'rawX', 'logX', 'duration', 'y'});

        trialLevelData = [trialLevelData; newTrials];
    end
end

writetable(perPatientResults, fullfile(outputFolder, 'per_patient_results.csv'));

%% Group summary of the patient slopes

groupSlopeSummary = table();
comparison = ["IED_count_vs_RT", "IED_count_vs_IT", ...
    "IED_count_vs_BR", "IED_count_vs_PR"];
xName = ["IED count during RT", ...
    "IED count during IT", ...
    "IED count during IT", ...
    "IED count during IT"];
yName = ["RT", "IT", "BR", "PR"];

for k = 1:4

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

    if k == 3 || k == 4
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

% br and pr are the same test with opposite outcome coding
% therefore, correct across the three unique tests only
p = groupSlopeSummary.permutationP_Slope(1:3);
correctedBonferroni = NaN(4, 1);
correctedBonferroni(1:3) = min(p * length(p), 1);
correctedBonferroni(4) = correctedBonferroni(3);
groupSlopeSummary.permutationP_Slope_Bonferroni = correctedBonferroni;

% FDR correction
[sortedP, order] = sort(p);
fdr = sortedP .* length(p) ./ (1:length(p))';

for k = length(fdr)-1:-1:1
    fdr(k) = min(fdr(k), fdr(k + 1));
end

fdr = min(fdr, 1);
correctedP = NaN(4, 1);
correctedUnique = zeros(3, 1);
correctedUnique(order) = fdr;
correctedP(1:3) = correctedUnique;
correctedP(4) = correctedP(3);
groupSlopeSummary.permutationP_Slope_FDR = correctedP;

writetable(groupSlopeSummary, fullfile(outputFolder, 'group_per_patient_slope_summary.csv'));

%% Mixed-effects models

groupMixedEffectsResults = table();

for k = 1:3

    rows = trialLevelData(trialLevelData.comparison == comparisonBase(k), :);
    rows.patientID = categorical(rows.patientID);

    if k == 3
        model = fitglme(rows, 'y ~ logX + duration + (1 | patientID)', ...
            'Distribution', 'Binomial', 'Link', 'logit');
        modelType = "Logistic mixed-effects: BR ~ logX + IT + (1 | patientID)";
    else
        model = fitlme(rows, 'y ~ logX + (1 | patientID)');
        modelType = "Linear mixed-effects: y ~ logX + (1 | patientID)";
    end

    coefficients = model.Coefficients;
    CI = coefCI(model);
    slopeRow = find(strcmp(coefficients.Name, 'logX'));

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

    newRow = table(comparisonBase(k), xNameBase(k), yNameBase(k), ...
        length(unique(rows.patientID)), height(rows), beta, SE, tStat, pValue, ...
        CILow, CIHigh, oddsRatio, oddsLow, oddsHigh, modelType, ...
        'VariableNames', {'comparison', 'xMeasure', 'yMeasure', 'nPatients', ...
        'nTrialPoints', 'beta_log10IEDcountPlus1', 'SE', 'tStat', 'pValue', ...
        'CILow', 'CIHigh', 'oddsRatio', 'oddsRatioCILow', 'oddsRatioCIHigh', 'modelType'});

    groupMixedEffectsResults = [groupMixedEffectsResults; newRow];

    % add the same br result with pop coded as the event
    if k == 3
        prRow = newRow;
        prRow.comparison = "IED_count_vs_PR";
        prRow.yMeasure = "PR";
        prRow.beta_log10IEDcountPlus1 = -beta;
        prRow.SE = SE;
        prRow.tStat = -tStat;
        prRow.pValue = pValue;
        prRow.CILow = -CIHigh;
        prRow.CIHigh = -CILow;
        prRow.oddsRatio = exp(-beta);
        prRow.oddsRatioCILow = exp(-CIHigh);
        prRow.oddsRatioCIHigh = exp(-CILow);
        prRow.modelType = "Logistic mixed-effects: PR = 1 - BR, adjusted for IT";

        groupMixedEffectsResults = [groupMixedEffectsResults; prRow];
    end
end

% br and pr are one unique statistical test
p = groupMixedEffectsResults.pValue(1:3);
correctedBonferroni = NaN(4, 1);
correctedBonferroni(1:3) = min(p * length(p), 1);
correctedBonferroni(4) = correctedBonferroni(3);
groupMixedEffectsResults.pValue_Bonferroni = correctedBonferroni;

% FDR correction
[sortedP, order] = sort(p);
fdr = sortedP .* length(p) ./ (1:length(p))';

for k = length(fdr)-1:-1:1
    fdr(k) = min(fdr(k), fdr(k + 1));
end

fdr = min(fdr, 1);
correctedP = NaN(4, 1);
correctedUnique = zeros(3, 1);
correctedUnique(order) = fdr;
correctedP(1:3) = correctedUnique;
correctedP(4) = correctedP(3);
groupMixedEffectsResults.pValue_FDR = correctedP;

writetable(groupMixedEffectsResults, fullfile(outputFolder, 'group_mixed_effects_results.csv'));

%% Figure

figure('Color', 'w', 'Position', [100 100 950 700]);
hold on

colors = [0.204 0.459 0.702; ...
          0.847 0.333 0.153; ...
          0.25 0.60 0.25; ...
          0.498 0.184 0.553];
values = perPatientResults.modelSlope;
groups = zeros(size(values));

for k = 1:4
    groups(perPatientResults.comparison == comparison(k)) = k;
end

boxplot(values, groups, 'Labels', {'RT', 'IT', 'BR', 'PR'}, 'Symbol', '', ...
    'Widths', 0.5, 'Colors', 'k');

rng(1)
xPlot = groups + (rand(size(groups)) - 0.5) * 0.32;

patients = unique(perPatientResults.patientID, 'stable');

for pt = 1:length(patients)
    patientRows = find(perPatientResults.patientID == patients(pt));
    plot(xPlot(patientRows), values(patientRows), '-', ...
        'Color', [0.80 0.80 0.80], 'LineWidth', 0.5)
end

for k = 1:4
    rows = groups == k;
    scatter(xPlot(rows), values(rows), 30, colors(k, :), 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')
end

yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);

yMin = min(values);
yMax = max(values);
yRange = yMax - yMin;
ylim([yMin - 0.18*yRange, yMax + 0.30*yRange]);

for k = 1:4
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

ylabel(sprintf(['Model slope\nRT/IT: seconds per log10(IED count + 1); ' ...
    'BR/PR: adjusted log-odds per log10(IED count + 1)']), ...
    'FontSize', 12, 'FontWeight', 'bold');
title('Per-patient IED effects', ...
    'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'LineWidth', 1.1, 'Box', 'off', 'TickDir', 'out');
pbaspect([1 1 1]);

exportgraphics(gcf, fullfile(outputFolder, 'summary_boxplots.pdf'), 'ContentType', 'vector');
close(gcf)
