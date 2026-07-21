% Summarize number of unique IED channels per trial vs RT, IT, and BR


% Author: Nill

clear;
clc;
close all;

inputFolderName_LFPIED = 'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';
outputFolderName = 'D:\Nill\code\BART\IED\0_0_new_IED\IED3_summary_of_IED2_chans_with_log\';

summaryFolder = fullfile(outputFolderName);

if ~exist(summaryFolder, 'dir')
    mkdir(summaryFolder);
end

% BR analysis uses only the number of unique IED channels during IT.
% Trials with 0 IT-period unique IED channels are included in the BR analysis.

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

perPatientResults = table();
trialLevelData = table();

nPermutations = 10000;

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    ptID = fileNameParts{1};

    disp(' ');
    disp(['Processing patient ID: ' ptID]);

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName));
    LFPIED = loadedData.LFPIED;

    nTrials = LFPIED.nTrials;

    RTs = LFPIED.RTs(:);
    ITs = LFPIED.ITs(:);
    isControl = LFPIED.isControl(:);

    if isfield(LFPIED, 'BankedTrials')
        BRs = LFPIED.BankedTrials(:);
    else
        warning(['LFPIED.BankedTrials not found for patient ' ptID '. BR analysis will be skipped for this patient.']);
        BRs = NaN(nTrials, 1);
    end

    minLen = min([nTrials, length(RTs), length(ITs), length(isControl), length(BRs)]);

    nTrials = minLen;
    RTs = RTs(1:nTrials);
    ITs = ITs(1:nTrials);
    isControl = isControl(1:nTrials);
    BRs = BRs(1:nTrials);

    nonControlTrials = isControl == 0;
    validRT_10sec = isfinite(RTs) & RTs <= 20;
    validBR = isfinite(BRs) & (BRs == 0 | BRs == 1);

    nUniqueChans_RT = countUniqueChannelsPerTrial(LFPIED, 'IED_occurance_RT', nTrials);
    nUniqueChans_IT = countUniqueChannelsPerTrial(LFPIED, 'IED_occurance_IT', nTrials);

    % For BR, use only the unique IED channel count during IT.
    nUniqueChans_BR = nUniqueChans_IT;
    xMeasureBR = 'Unique IED channels during IT';

    keep_chans_RT = nonControlTrials & validRT_10sec & isfinite(RTs) & RTs > 0 & nUniqueChans_RT >= 0;
    keep_chans_IT = nonControlTrials & validRT_10sec & isfinite(ITs) & ITs > 0 & nUniqueChans_IT >= 0;
    keep_chans_BR = nonControlTrials & validRT_10sec & validBR & nUniqueChans_BR >= 0;

    panels = {
        'Unique_channels_vs_RT', 'Unique IED channels during RT', 'RT', nUniqueChans_RT, RTs, keep_chans_RT;
        'Unique_channels_vs_IT', 'Unique IED channels during IT', 'IT', nUniqueChans_IT, ITs, keep_chans_IT;
        'Unique_channels_vs_BR', xMeasureBR, 'BR', nUniqueChans_BR, BRs, keep_chans_BR
    };

    for pp = 1:size(panels, 1)

        comparisonName = panels{pp, 1};
        xMeasure = panels{pp, 2};
        yMeasure = panels{pp, 3};
        xAll = panels{pp, 4};
        yAll = panels{pp, 5};
        keepIdx = panels{pp, 6};

        rawX = xAll(keepIdx);
        y = yAll(keepIdx);

        [rawX, logX, y] = cleanPanelData(rawX, y, yMeasure);

        statsRow = summarizeOnePanel(rawX, logX, y, yMeasure);

        newPerPatientRow = table( ...
            string(ptID), ...
            string(comparisonName), ...
            string(xMeasure), ...
            string(yMeasure), ...
            statsRow.nPoints, ...
            statsRow.rawXMean, ...
            statsRow.rawXMedian, ...
            statsRow.yMean, ...
            statsRow.yMedian, ...
            statsRow.modelIntercept, ...
            statsRow.modelSlope, ...
            statsRow.modelSlopeSE, ...
            statsRow.modelSlopeT, ...
            statsRow.modelSlopeP, ...
            statsRow.modelSlopeCILow, ...
            statsRow.modelSlopeCIHigh, ...
            statsRow.oddsRatio, ...
            statsRow.oddsRatioCILow, ...
            statsRow.oddsRatioCIHigh, ...
            statsRow.rSquared, ...
            statsRow.adjRSquared, ...
            string(statsRow.modelType), ...
            'VariableNames', { ...
                'patientID', ...
                'comparison', ...
                'xMeasure', ...
                'yMeasure', ...
                'nPoints', ...
                'rawXMean', ...
                'rawXMedian', ...
                'yMean', ...
                'yMedian', ...
                'modelIntercept', ...
                'modelSlope', ...
                'modelSlopeSE', ...
                'modelSlopeT', ...
                'modelSlopeP', ...
                'modelSlopeCILow', ...
                'modelSlopeCIHigh', ...
                'oddsRatio', ...
                'oddsRatioCILow', ...
                'oddsRatioCIHigh', ...
                'rSquared', ...
                'adjRSquared', ...
                'modelType' ...
            });

        perPatientResults = [perPatientResults; newPerPatientRow];

        if ~isempty(logX)

            n = length(logX);

            newTrialRows = table( ...
                repmat(string(ptID), n, 1), ...
                repmat(string(comparisonName), n, 1), ...
                repmat(string(xMeasure), n, 1), ...
                repmat(string(yMeasure), n, 1), ...
                rawX(:), ...
                logX(:), ...
                y(:), ...
                'VariableNames', { ...
                    'patientID', ...
                    'comparison', ...
                    'xMeasure', ...
                    'yMeasure', ...
                    'rawX', ...
                    'logX', ...
                    'y' ...
                });

            trialLevelData = [trialLevelData; newTrialRows];

        end

    end

end

perPatientOutputFile = fullfile(summaryFolder, 'per_patient_results.csv');
writetable(perPatientResults, perPatientOutputFile);

groupSlopeSummary = makeGroupSlopeSummaryFromPerPatient(perPatientResults, nPermutations);
groupSlopeSummaryOutputFile = fullfile(summaryFolder, 'group_per_patient_slope_summary.csv');
writetable(groupSlopeSummary, groupSlopeSummaryOutputFile);

groupMixedEffectsResults = runGroupMixedEffects(trialLevelData);
groupMixedEffectsOutputFile = fullfile(summaryFolder, 'group_mixed_effects_results.csv');
writetable(groupMixedEffectsResults, groupMixedEffectsOutputFile);

summaryPDF = fullfile(summaryFolder, 'summary_boxplots.pdf');
plotSummaryBoxplots(perPatientResults, summaryPDF, groupMixedEffectsResults);

disp(' ');
disp('Summary finished.');
disp(['Saved: ' perPatientOutputFile]);
disp(['Saved: ' groupSlopeSummaryOutputFile]);
disp(['Saved: ' groupMixedEffectsOutputFile]);
disp(['Saved: ' summaryPDF]);

function [rawX, logX, y] = cleanPanelData(rawX, y, yMeasure)

    rawX = rawX(:);
    y = y(:);

    if string(yMeasure) == "BR"
        validIdx = isfinite(rawX) & isfinite(y) & rawX >= 0 & (y == 0 | y == 1);
    else
        validIdx = isfinite(rawX) & isfinite(y) & rawX >= 0 & y > 0;
    end

    rawX = rawX(validIdx);
    y = y(validIdx);

    logX = log10(rawX + 1);

end

function statsRow = summarizeOnePanel(rawX, logX, y, yMeasure)

    rawX = rawX(:);
    logX = logX(:);
    y = y(:);

    statsRow.nPoints = length(y);

    statsRow.rawXMean = NaN;
    statsRow.rawXMedian = NaN;
    statsRow.yMean = NaN;
    statsRow.yMedian = NaN;

    statsRow.modelIntercept = NaN;
    statsRow.modelSlope = NaN;
    statsRow.modelSlopeSE = NaN;
    statsRow.modelSlopeT = NaN;
    statsRow.modelSlopeP = NaN;
    statsRow.modelSlopeCILow = NaN;
    statsRow.modelSlopeCIHigh = NaN;

    statsRow.oddsRatio = NaN;
    statsRow.oddsRatioCILow = NaN;
    statsRow.oddsRatioCIHigh = NaN;

    statsRow.rSquared = NaN;
    statsRow.adjRSquared = NaN;

    statsRow.modelType = "Not enough data";

    if isempty(y)
        return;
    end

    statsRow.rawXMean = mean(rawX, 'omitnan');
    statsRow.rawXMedian = median(rawX, 'omitnan');

    if string(yMeasure) == "BR"
        statsRow.yMean = mean(y, 'omitnan');
        statsRow.yMedian = median(y, 'omitnan');
    end

    if length(logX) <= 2 || length(unique(logX)) <= 1
        return;
    end

    T = table(logX, y, 'VariableNames', {'logX', 'y'});

    if string(yMeasure) == "BR"

        try

            glm = fitglm(T, 'y ~ logX', ...
                'Distribution', 'binomial', ...
                'Link', 'logit');

            coefTable = glm.Coefficients;
            ciTable = coefCI(glm);

            slopeIdx = find(strcmp(coefTable.Properties.RowNames, 'logX'));
            interceptIdx = find(strcmp(coefTable.Properties.RowNames, '(Intercept)'));

            statsRow.modelIntercept = coefTable.Estimate(interceptIdx);
            statsRow.modelSlope = coefTable.Estimate(slopeIdx);
            statsRow.modelSlopeSE = coefTable.SE(slopeIdx);
            statsRow.modelSlopeT = coefTable.tStat(slopeIdx);
            statsRow.modelSlopeP = coefTable.pValue(slopeIdx);
            statsRow.modelSlopeCILow = ciTable(slopeIdx, 1);
            statsRow.modelSlopeCIHigh = ciTable(slopeIdx, 2);

            statsRow.oddsRatio = exp(statsRow.modelSlope);
            statsRow.oddsRatioCILow = exp(statsRow.modelSlopeCILow);
            statsRow.oddsRatioCIHigh = exp(statsRow.modelSlopeCIHigh);

            statsRow.rSquared = NaN;
            statsRow.adjRSquared = NaN;

            statsRow.modelType = "Per-patient logistic GLM: BR ~ logX";

        catch

            warning('Per-patient logistic GLM failed.');
            statsRow.modelType = "Per-patient logistic GLM failed";

        end

    else

        try

            lm = fitlm(T, 'y ~ logX');

            coefTable = lm.Coefficients;
            ciTable = coefCI(lm);

            slopeIdx = find(strcmp(coefTable.Properties.RowNames, 'logX'));
            interceptIdx = find(strcmp(coefTable.Properties.RowNames, '(Intercept)'));

            statsRow.modelIntercept = coefTable.Estimate(interceptIdx);
            statsRow.modelSlope = coefTable.Estimate(slopeIdx);
            statsRow.modelSlopeSE = coefTable.SE(slopeIdx);
            statsRow.modelSlopeT = coefTable.tStat(slopeIdx);
            statsRow.modelSlopeP = coefTable.pValue(slopeIdx);
            statsRow.modelSlopeCILow = ciTable(slopeIdx, 1);
            statsRow.modelSlopeCIHigh = ciTable(slopeIdx, 2);

            statsRow.rSquared = lm.Rsquared.Ordinary;
            statsRow.adjRSquared = lm.Rsquared.Adjusted;

            statsRow.modelType = "Per-patient linear model: y ~ logX";

        catch

            warning('Per-patient linear model failed.');
            statsRow.modelType = "Per-patient linear model failed";

        end

    end

end

function groupSlopeSummary = makeGroupSlopeSummaryFromPerPatient(perPatientResults, nPermutations)

    comparisonNames = unique(perPatientResults.comparison, 'stable');
    groupSlopeSummary = table();

    for cc = 1:length(comparisonNames)

        comparisonName = comparisonNames(cc);

        rows = perPatientResults(perPatientResults.comparison == comparisonName, :);

        validSlopeRows = isfinite(rows.modelSlope);
        slopeVals = rows.modelSlope(validSlopeRows);

        nPatients = height(rows);
        nPatientsWithSlope = length(slopeVals);
        totalTrialPoints = sum(rows.nPoints, 'omitnan');

        medianSlope = NaN;
        slopeIQR = NaN;
        meanSlope = NaN;
        slopeCILow = NaN;
        slopeCIHigh = NaN;
        permutationP_Slope = NaN;
        nPositiveSlope = NaN;
        nNegativeSlope = NaN;

        medianOddsRatio = NaN;
        meanOddsRatio = NaN;

        if nPatientsWithSlope > 0

            medianSlope = median(slopeVals, 'omitnan');
            slopeIQR = iqr(slopeVals);
            meanSlope = mean(slopeVals, 'omitnan');

            nPositiveSlope = sum(slopeVals > 0);
            nNegativeSlope = sum(slopeVals < 0);

            permutationP_Slope = permutationTestMean(slopeVals, nPermutations);

            if nPatientsWithSlope > 1
                seSlope = std(slopeVals, 'omitnan') / sqrt(nPatientsWithSlope);
                tCrit = tinv(0.975, nPatientsWithSlope - 1);

                slopeCILow = meanSlope - tCrit * seSlope;
                slopeCIHigh = meanSlope + tCrit * seSlope;
            end

            if rows.yMeasure(1) == "BR"
                medianOddsRatio = exp(medianSlope);
                meanOddsRatio = exp(meanSlope);
            end

        end

        newRow = table( ...
            comparisonName, ...
            rows.xMeasure(1), ...
            rows.yMeasure(1), ...
            nPatients, ...
            totalTrialPoints, ...
            nPatientsWithSlope, ...
            medianSlope, ...
            slopeIQR, ...
            meanSlope, ...
            slopeCILow, ...
            slopeCIHigh, ...
            nPositiveSlope, ...
            nNegativeSlope, ...
            permutationP_Slope, ...
            medianOddsRatio, ...
            meanOddsRatio, ...
            'VariableNames', { ...
                'comparison', ...
                'xMeasure', ...
                'yMeasure', ...
                'nPatients', ...
                'totalTrialPoints', ...
                'nPatientsWithSlope', ...
                'medianSlope', ...
                'iqrSlope', ...
                'meanSlope', ...
                'meanSlopeCILow', ...
                'meanSlopeCIHigh', ...
                'nPositiveSlope', ...
                'nNegativeSlope', ...
                'permutationP_Slope', ...
                'medianOddsRatio', ...
                'meanOddsRatio' ...
            });

        groupSlopeSummary = [groupSlopeSummary; newRow];

    end

    [pBonf, pFDR] = multipleHypothesisCorrection(groupSlopeSummary.permutationP_Slope);

    groupSlopeSummary.permutationP_Slope_Bonferroni = pBonf;
    groupSlopeSummary.permutationP_Slope_FDR = pFDR;

end

function groupMixedEffectsResults = runGroupMixedEffects(trialLevelData)

    comparisonNames = unique(trialLevelData.comparison, 'stable');
    groupMixedEffectsResults = table();

    for cc = 1:length(comparisonNames)

        comparisonName = comparisonNames(cc);

        rowsAll = trialLevelData(trialLevelData.comparison == comparisonName, :);
        rows = rowsAll(isfinite(rowsAll.logX) & isfinite(rowsAll.y), :);

        nTrialPoints = height(rows);
        nPatients = length(unique(rows.patientID));

        beta = NaN;
        se = NaN;
        tStat = NaN;
        pValue = NaN;
        ciLow = NaN;
        ciHigh = NaN;

        oddsRatio = NaN;
        oddsRatioCILow = NaN;
        oddsRatioCIHigh = NaN;

        modelType = "Not enough data";

        isBR = rowsAll.yMeasure(1) == "BR";

        if nTrialPoints >= 10 && nPatients >= 2 && length(unique(rows.logX)) > 1

            rows.patientID = categorical(rows.patientID);

            if isBR

                if length(unique(rows.y)) <= 1
                    modelType = "Not enough BR variation";
                else

                    try

                        glme = fitglme(rows, 'y ~ logX + (1 | patientID)', ...
                            'Distribution', 'Binomial', ...
                            'Link', 'logit');

                        coefTable = glme.Coefficients;
                        ciTable = coefCI(glme);

                        slopeIdx = find(strcmp(coefTable.Name, 'logX'));

                        beta = coefTable.Estimate(slopeIdx);
                        se = coefTable.SE(slopeIdx);
                        tStat = coefTable.tStat(slopeIdx);
                        pValue = coefTable.pValue(slopeIdx);
                        ciLow = ciTable(slopeIdx, 1);
                        ciHigh = ciTable(slopeIdx, 2);

                        oddsRatio = exp(beta);
                        oddsRatioCILow = exp(ciLow);
                        oddsRatioCIHigh = exp(ciHigh);

                        modelType = "Logistic mixed-effects: BR ~ logX + (1 | patientID)";

                    catch

                        warning(['fitglme failed for ' char(comparisonName) '. Using pooled fitglm instead.']);

                        glm = fitglm(rows, 'y ~ logX', ...
                            'Distribution', 'binomial', ...
                            'Link', 'logit');

                        coefTable = glm.Coefficients;
                        ciTable = coefCI(glm);

                        slopeIdx = find(strcmp(coefTable.Properties.RowNames, 'logX'));

                        beta = coefTable.Estimate(slopeIdx);
                        se = coefTable.SE(slopeIdx);
                        tStat = coefTable.tStat(slopeIdx);
                        pValue = coefTable.pValue(slopeIdx);
                        ciLow = ciTable(slopeIdx, 1);
                        ciHigh = ciTable(slopeIdx, 2);

                        oddsRatio = exp(beta);
                        oddsRatioCILow = exp(ciLow);
                        oddsRatioCIHigh = exp(ciHigh);

                        modelType = "Pooled logistic GLM fallback";

                    end

                end

            else

                try

                    lme = fitlme(rows, 'y ~ logX + (1 | patientID)');

                    coefTable = lme.Coefficients;
                    ciTable = coefCI(lme);

                    slopeIdx = find(strcmp(coefTable.Name, 'logX'));

                    beta = coefTable.Estimate(slopeIdx);
                    se = coefTable.SE(slopeIdx);
                    tStat = coefTable.tStat(slopeIdx);
                    pValue = coefTable.pValue(slopeIdx);
                    ciLow = ciTable(slopeIdx, 1);
                    ciHigh = ciTable(slopeIdx, 2);

                    modelType = "Linear mixed-effects: y ~ logX + (1 | patientID)";

                catch

                    warning(['fitlme failed for ' char(comparisonName) '. Using pooled fitlm instead.']);

                    lm = fitlm(rows, 'y ~ logX');

                    coefTable = lm.Coefficients;
                    ciTable = coefCI(lm);

                    slopeIdx = find(strcmp(coefTable.Properties.RowNames, 'logX'));

                    beta = coefTable.Estimate(slopeIdx);
                    se = coefTable.SE(slopeIdx);
                    tStat = coefTable.tStat(slopeIdx);
                    pValue = coefTable.pValue(slopeIdx);
                    ciLow = ciTable(slopeIdx, 1);
                    ciHigh = ciTable(slopeIdx, 2);

                    modelType = "Pooled linear model fallback";

                end

            end

        end

        newRow = table( ...
            comparisonName, ...
            rowsAll.xMeasure(1), ...
            rowsAll.yMeasure(1), ...
            nPatients, ...
            nTrialPoints, ...
            beta, ...
            se, ...
            tStat, ...
            pValue, ...
            ciLow, ...
            ciHigh, ...
            oddsRatio, ...
            oddsRatioCILow, ...
            oddsRatioCIHigh, ...
            modelType, ...
            'VariableNames', { ...
                'comparison', ...
                'xMeasure', ...
                'yMeasure', ...
                'nPatients', ...
                'nTrialPoints', ...
                'beta_log10UniqueChannelsPlus1', ...
                'SE', ...
                'tStat', ...
                'pValue', ...
                'CILow', ...
                'CIHigh', ...
                'oddsRatio', ...
                'oddsRatioCILow', ...
                'oddsRatioCIHigh', ...
                'modelType' ...
            });

        groupMixedEffectsResults = [groupMixedEffectsResults; newRow];

    end

    [pBonf, pFDR] = multipleHypothesisCorrection(groupMixedEffectsResults.pValue);

    groupMixedEffectsResults.pValue_Bonferroni = pBonf;
    groupMixedEffectsResults.pValue_FDR = pFDR;

end

function p = permutationTestMean(vals, nPermutations)

    vals = vals(:);
    vals = vals(isfinite(vals));

    if isempty(vals)
        p = NaN;
        return;
    end

    if nargin < 2 || isempty(nPermutations)
        nPermutations = 10000;
    end

    observedMean = mean(vals, 'omitnan');
    nVals = length(vals);

    rng(1);

    randomSigns = randi([0 1], nVals, nPermutations) * 2 - 1;
    permutedMeans = mean(bsxfun(@times, vals, randomSigns), 1, 'omitnan');

    p = (sum(abs(permutedMeans) >= abs(observedMean)) + 1) / (nPermutations + 1);

end

function [pBonf, pFDR] = multipleHypothesisCorrection(pVals)

    pVals = pVals(:);

    pBonf = NaN(size(pVals));
    pFDR = NaN(size(pVals));

    validIdx = isfinite(pVals) & pVals >= 0 & pVals <= 1;
    validP = pVals(validIdx);

    m = length(validP);

    if m == 0
        return;
    end

    pBonf(validIdx) = min(validP * m, 1);

    [sortedP, sortIdx] = sort(validP);
    ranks = (1:m)';

    bhVals = sortedP .* m ./ ranks;

    for ii = m-1:-1:1
        bhVals(ii) = min(bhVals(ii), bhVals(ii + 1));
    end

    bhVals = min(bhVals, 1);

    unsortedBH = NaN(size(validP));
    unsortedBH(sortIdx) = bhVals;

    pFDR(validIdx) = unsortedBH;

end

function plotSummaryBoxplots(perPatientResults, summaryPDF, groupMixedEffectsResults)

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 800 700]);

    comparisonOrder = unique(perPatientResults.comparison, 'stable');

    colors = [
        0.204   0.459   0.702
        0.847   0.333   0.153
        0.25    0.60    0.25
    ];

    if length(comparisonOrder) > size(colors, 1)
        colors = lines(length(comparisonOrder));
    end

    ax1 = axes(fig);

    plotBoxScatterOneMetric( ...
        ax1, ...
        perPatientResults, ...
        comparisonOrder, ...
        colors, ...
        'modelSlope', ...
        sprintf('Model slope\nRT/IT: seconds per log10(unique IED channels + 1); BR: log-odds per log10(unique IED channels + 1)'), ...
        'Per-patient model slope', ...
        groupMixedEffectsResults, ...
        'pValue_FDR');

    exportgraphics(fig, summaryPDF, 'ContentType', 'vector');

    close(fig);

end

function plotBoxScatterOneMetric(ax, perPatientResults, comparisonOrder, colors, valueColumn, yLabelText, titleText, sigTable, pColumn)

    hold(ax, 'on');

    values = perPatientResults.(valueColumn);
    comparisons = perPatientResults.comparison;
    patientIDs = perPatientResults.patientID;

    validRows = isfinite(values);
    values = values(validRows);
    comparisons = comparisons(validRows);
    patientIDs = patientIDs(validRows);

    groupID = NaN(size(values));

    for cc = 1:length(comparisonOrder)
        groupID(comparisons == comparisonOrder(cc)) = cc;
    end

    validGroupRows = isfinite(groupID);
    values = values(validGroupRows);
    comparisons = comparisons(validGroupRows);
    patientIDs = patientIDs(validGroupRows);
    groupID = groupID(validGroupRows);

    prettyLabels = makePrettyComparisonLabels(comparisonOrder);

    boxplot( ...
        ax, ...
        values, ...
        groupID, ...
        'Labels', cellstr(prettyLabels), ...
        'Symbol', '', ...
        'Widths', 0.5, ...
        'Colors', 'k');

    hBoxes = findobj(ax, 'Tag', 'Box');
    set(hBoxes, 'Color', 'k', 'LineWidth', 0.5);

    hMedian = findobj(ax, 'Tag', 'Median');
    set(hMedian, 'Color', 'k', 'LineWidth', 0.5);

    hWhisker = findobj(ax, 'Tag', 'Whisker');
    set(hWhisker, 'Color', 'k', 'LineWidth', 0.5);

    hUpperAdj = findobj(ax, 'Tag', 'Upper Adjacent Value');
    hLowerAdj = findobj(ax, 'Tag', 'Lower Adjacent Value');
    set([hUpperAdj; hLowerAdj], 'Color', 'k', 'LineWidth', 0.5);

    hOutliers = findobj(ax, 'Tag', 'Outliers');
    set(hOutliers, 'Marker', 'none');

    rng(1);

    xPlot = NaN(size(values));

    for cc = 1:length(comparisonOrder)

        idx = comparisons == comparisonOrder(cc);

        if any(idx)

            jitterAmount = 0.16;
            xPlot(idx) = cc + (rand(sum(idx), 1) - 0.5) * 2 * jitterAmount;

        end

    end

    lineColor = [0.65 0.65 0.65];
    lineAlpha = 0.25;

    
    uniquePatients = unique(patientIDs, 'stable');
    
    for pp = 1:length(uniquePatients)
    
        patientIdx = patientIDs == uniquePatients(pp);
    
        for cc = 1:(length(comparisonOrder) - 1)
    
            idx1 = find(patientIdx & groupID == cc, 1, 'first');
            idx2 = find(patientIdx & groupID == cc + 1, 1, 'first');
    
            if ~isempty(idx1) && ~isempty(idx2)
    
                patch( ...
                    ax, ...
                    [xPlot(idx1) xPlot(idx2)], ...
                    [values(idx1) values(idx2)], ...
                    lineColor, ...
                    'FaceColor', 'none', ...
                    'EdgeColor', lineColor, ...
                    'EdgeAlpha', lineAlpha, ...
                    'LineWidth', 0.5);
    
            end
    
        end
    
    end

    for cc = 1:length(comparisonOrder)

        idx = comparisons == comparisonOrder(cc);

        if any(idx)

            scatter( ...
                ax, ...
                xPlot(idx), ...
                values(idx), ...
                30, ...
                'MarkerFaceColor', colors(cc, :), ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 0.5);

        end

    end

    yline(ax, 0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);

    yMin = min(values, [], 'omitnan');
    yMax = max(values, [], 'omitnan');
    yRange = yMax - yMin;

    if yRange == 0 || ~isfinite(yRange)
        yRange = 1;
    end

    ylim(ax, [yMin - 0.18 * yRange, yMax + 0.30 * yRange]);

    for cc = 1:length(comparisonOrder)

        comparisonName = comparisonOrder(cc);

        pVal = NaN;

        if ~isempty(sigTable) && any(sigTable.comparison == comparisonName)

            sigRow = sigTable(sigTable.comparison == comparisonName, :);

            if ismember(pColumn, sigTable.Properties.VariableNames)
                pVal = sigRow.(pColumn)(1);
            end

        end

        sigStar = pToStars(pVal);

        if strlength(sigStar) > 0

            idx = comparisons == comparisonName;

            if any(idx)

                groupMax = max(values(idx), [], 'omitnan');
                yStar = groupMax + 0.09 * yRange;

                text( ...
                    ax, ...
                    cc, ...
                    yStar, ...
                    sigStar, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 18, ...
                    'FontWeight', 'bold', ...
                    'Color', 'k');

            end

        end

    end

    ylabel(ax, yLabelText, 'FontSize', 12, 'FontWeight', 'bold');
    title(ax, titleText, 'FontSize', 13, 'FontWeight', 'bold');

    set(ax, ...
        'FontSize', 11, ...
        'LineWidth', 1.1, ...
        'Box', 'off', ...
        'TickDir', 'out', ...
        'XTickLabelRotation', 0, ...
        'TickLabelInterpreter', 'none');

    grid(ax, 'off');

    pbaspect(ax, [1 1 1]);

    hold(ax, 'off');

end

function labels = makePrettyComparisonLabels(comparisonOrder)

    labels = strings(size(comparisonOrder));

    for ii = 1:length(comparisonOrder)

        comparisonName = comparisonOrder(ii);

        if contains(comparisonName, "RT")
            labels(ii) = "RT";
        elseif contains(comparisonName, "IT")
            labels(ii) = "IT";
        elseif contains(comparisonName, "BR")
            labels(ii) = "BR";
        else
            labels(ii) = comparisonName;
        end

    end

end

function tf = isBRComparison(comparisonName)

    tf = contains(string(comparisonName), "BR");

end

function sigStar = pToStars(pVal)

    sigStar = "";

    if ~isfinite(pVal)
        return;
    end

    if pVal < 0.001
        sigStar = "***";
    elseif pVal < 0.01
        sigStar = "**";
    elseif pVal < 0.05
        sigStar = "*";
    end

end

function nUniqueChansPerTrial = countUniqueChannelsPerTrial(LFPIED, fieldName, nTrials)

    nUniqueChansPerTrial = zeros(nTrials, 1);

    if ~isfield(LFPIED, fieldName)
        return;
    end

    IEDocc = LFPIED.(fieldName);

    if isempty(IEDocc) || size(IEDocc, 2) < 2
        return;
    end

    trialIndices = IEDocc(:, 1);
    channelIndices = IEDocc(:, 2);

    validRows = ...
        isfinite(trialIndices) & ...
        isfinite(channelIndices) & ...
        trialIndices >= 1 & ...
        trialIndices <= nTrials & ...
        channelIndices >= 1;

    trialIndices = round(trialIndices(validRows));
    channelIndices = round(channelIndices(validRows));

    if isempty(trialIndices)
        return;
    end

    uniqueTrialChannelPairs = unique([trialIndices channelIndices], 'rows');

    nUniqueChansPerTrial = accumarray( ...
        uniqueTrialChannelPairs(:, 1), ...
        1, ...
        [nTrials 1], ...
        @sum, ...
        0);

end
