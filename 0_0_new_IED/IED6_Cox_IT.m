% Cox proportional hazards model for full IT duration

%
% Author: Nill

clear;
clc;
close all;

%% Paths and analysis settings

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED6_Cox_IT\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

% If LFPIED.Fs exists, the saved value is used.
% Otherwise, the sampling frequency is silently set to 1000 Hz.
defaultSamplingFrequencyHz = 1000;

% "count" = cumulative number of IEDs observed so far.
% "any"   = 0 before the first IED and 1 afterward.
iedPredictorMode = "count";

maximumRTSeconds = 10;

% Maximum cumulative IED count shown in the model-effect visualization.
% The plotted maximum is the smaller of this value and the 95th percentile
% of trial-level IED counts, with a minimum displayed value of 1.
maximumIEDCountForVisualization = 10;

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

countingProcessData = table();
trialSummary = table();

%% Build counting-process data across all participants

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    patientID = string(fileNameParts{1});

    fprintf('\nprocessing patient: %s\n', patientID);

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName));

    if ~isfield(loadedData, 'LFPIED')
        fprintf('skipped: lfpied structure was not found.\n');
        continue;
    end

    LFPIED = loadedData.LFPIED;

    requiredFields = {'RTs', 'ITs', 'isControl', 'balloonType', 'IED_occurance_IT'};

    missingField = false;

    for ff = 1:length(requiredFields)
        if ~isfield(LFPIED, requiredFields{ff})
            fprintf('skipped: missing field %s.\n', requiredFields{ff});
            missingField = true;
        end
    end

    if missingField
        continue;
    end

    RTs = LFPIED.RTs(:);
    ITs = LFPIED.ITs(:);
    durationVector = ITs;
    isControl = LFPIED.isControl(:);
    balloonType = LFPIED.balloonType(:);

    nTrials = determineNumberOfTrials( ...
        LFPIED, RTs, ITs, isControl, balloonType);

    RTs = RTs(1:nTrials);
    ITs = ITs(1:nTrials);
    durationVector = durationVector(1:nTrials);
    isControl = isControl(1:nTrials);
    balloonType = balloonType(1:nTrials);

    balloonColorCode = mapBalloonColorCode(balloonType);

    samplingFrequencyHz = getSamplingFrequency( ...
        LFPIED, defaultSamplingFrequencyHz);

    IEDoccurrence = LFPIED.IED_occurance_IT;

    for trial = 1:nTrials

        durationSeconds = durationVector(trial);
        finalEventObserved = true;

        % include only non-control yellow, orange, and red balloon trials.
        keepTrial = ...
            isControl(trial) == 0 && ...
            isfinite(RTs(trial)) && ...
            RTs(trial) > 0 && ...
            RTs(trial) <= maximumRTSeconds && ...
            isfinite(durationSeconds) && ...
            durationSeconds > 0 && ...
            isfinite(balloonColorCode(trial)) && ...
            ismember(balloonColorCode(trial), [1 2 3]);

        if ~keepTrial
            continue;
        end

        trialRows = makeCountingProcessRows( ...
            patientID, ...
            trial, ...
            durationSeconds, ...
            IEDoccurrence, ...
            samplingFrequencyHz, ...
            balloonColorCode(trial), ...
            finalEventObserved, ...
            iedPredictorMode);

        countingProcessData = [countingProcessData; trialRows];

        if isempty(IEDoccurrence) || size(IEDoccurrence, 2) < 1
            nIEDsInTrial = 0;
        else
            nIEDsInTrial = sum( ...
                isfinite(IEDoccurrence(:, 1)) & ...
                round(IEDoccurrence(:, 1)) == trial);
        end

        newTrialSummaryRow = table( ...
            patientID, ...
            trial, ...
            durationSeconds, ...
            double(finalEventObserved), ...
            nIEDsInTrial, ...
            balloonColorCode(trial), ...
            samplingFrequencyHz, ...
            'VariableNames', { ...
                'patientID', ...
                'trialNumber', ...
                'durationSeconds', ...
                'eventObserved', ...
                'nIEDsInFullDuration', ...
                'balloonColorCode', ...
                'samplingFrequencyHz' ...
            });

        trialSummary = [trialSummary; newTrialSummaryRow];

    end

end

if isempty(countingProcessData)
    error('no valid counting-process rows were created.');
end

%% Create predictors

[countingProcessData.patientStratum, patientLevels] = ...
    findgroups(countingProcessData.patientID);

% The trial filter and color mapping above allow only the three
% non-control balloon colors: yellow, orange, and red.
if any(~ismember(countingProcessData.balloonColorCode, [1 2 3]))
    error('unexpected balloon color found. only yellow, orange, and red are allowed.');
end

% use yellow as the fixed reference balloon color.
colorLevels = [1; 2; 3];
referenceColor = 1;
referenceColorName = "yellow";

observedColorLevels = unique( ...
    countingProcessData.balloonColorCode, 'sorted');

if ~all(ismember(colorLevels, observedColorLevels))
    error('yellow, orange, and red trials must all be present in the final data.');
end

X = countingProcessData.IED_timevarying;
predictorNames = getIEDPredictorName(iedPredictorMode);

% yellow is the reference category.
X = [X, double(countingProcessData.balloonColorCode == 2)];
predictorNames(end + 1, 1) = "orange_vs_yellow";

X = [X, double(countingProcessData.balloonColorCode == 3)];
predictorNames(end + 1, 1) = "red_vs_yellow";

T = [countingProcessData.tStart, countingProcessData.tStop];
censoring = logical(countingProcessData.censored);
strata = countingProcessData.patientStratum;

validRows = ...
    all(isfinite(X), 2) & ...
    all(isfinite(T), 2) & ...
    T(:, 1) >= 0 & ...
    T(:, 2) > T(:, 1) & ...
    isfinite(strata);

X = X(validRows, :);
T = T(validRows, :);
censoring = censoring(validRows);
strata = strata(validRows);
countingProcessData = countingProcessData(validRows, :);

%% Fit participant-stratified Cox proportional hazards model

coxOptions = statset('coxphfit');
coxOptions.Display = 'final';
coxOptions.MaxIter = 1000;
coxOptions.MaxFunEvals = 5000;

[beta, logLikelihood, baselineCumulativeHazard, stats] = ...
    coxphfit( ...
        X, ...
        T, ...
        'Censoring', censoring, ...
        'Strata', strata, ...
        'Ties', 'efron', ...
        'Baseline', 0, ...
        'Options', coxOptions);

%% Participant-clustered robust inference

[clusterRobustCovariance, clusterRobustSE, ...
    clusterRobustZ, clusterRobustP] = ...
    computeClusterRobustInference( ...
        stats, ...
        beta, ...
        countingProcessData.patientStratum);

%% Create and print the results table

% Confidence intervals and primary p-values use participant-clustered
% robust standard errors. Model-based values are also saved.
ciLowBeta = beta - 1.96 .* clusterRobustSE;
ciHighBeta = beta + 1.96 .* clusterRobustSE;

hazardRatio = exp(beta);
hazardRatioCILow = exp(ciLowBeta);
hazardRatioCIHigh = exp(ciHighBeta);

interpretation = strings(length(beta), 1);

for ii = 1:length(beta)

    if ii == 1
        interpretation(ii) = ...
            "hr > 1 means a higher action hazard and therefore shorter it; " + ...
            "hr < 1 means a lower action hazard and therefore longer it.";
    else
        interpretation(ii) = ...
            "balloon color effect relative to " + referenceColorName;
    end

end

coxResults = table( ...
    predictorNames(:), ...
    beta(:), ...
    stats.se(:), ...
    stats.z(:), ...
    stats.p(:), ...
    clusterRobustSE(:), ...
    clusterRobustZ(:), ...
    clusterRobustP(:), ...
    ciLowBeta(:), ...
    ciHighBeta(:), ...
    hazardRatio(:), ...
    hazardRatioCILow(:), ...
    hazardRatioCIHigh(:), ...
    interpretation(:), ...
    'VariableNames', { ...
        'predictor', ...
        'beta_logHazard', ...
        'modelBasedSE', ...
        'modelBasedZ', ...
        'modelBasedPValue', ...
        'clusterRobustSE', ...
        'clusterRobustZ', ...
        'pValue', ...
        'betaCILow', ...
        'betaCIHigh', ...
        'hazardRatio', ...
        'hazardRatioCILow', ...
        'hazardRatioCIHigh', ...
        'interpretation' ...
    });

fprintf('\n============================================================\n');
fprintf('full-duration it cox model\n');
fprintf('============================================================\n');
fprintf('participants: %d\n', length(patientLevels));
fprintf('trials: %d\n', height(trialSummary));
fprintf('counting-process rows: %d\n', height(countingProcessData));
fprintf('observed endpoint events: %d\n', sum(~censoring));
fprintf('log likelihood: %.6f\n\n', logLikelihood);

disp(coxResults);

iedRow = coxResults(1, :);

fprintf('\nied effect:\n');
fprintf('beta = %.6f\n', iedRow.beta_logHazard);
fprintf('hazard ratio = %.6f\n', iedRow.hazardRatio);
fprintf('95%% ci = [%.6f, %.6f]\n', ...
    iedRow.hazardRatioCILow, iedRow.hazardRatioCIHigh);
fprintf('p = %.6g\n', iedRow.pValue);

if iedRow.pValue < 0.05

    if iedRow.hazardRatio > 1
        fprintf(['conclusion: more prior ieds are associated with a ' ...
            'significantly higher action hazard, corresponding to ' ...
            'shorter it.\n']);
    else
        fprintf(['conclusion: more prior ieds are associated with a ' ...
            'significantly lower action hazard, corresponding to ' ...
            'longer it.\n']);
    end

else
    fprintf(['conclusion: the association between prior ieds and it ' ...
        'action hazard is not statistically significant.\n']);
end

%% Create model-based visualizations

[visualizationFigure, IED_effect_curve] = ...
    plotCoxModelVisualization( ...
        coxResults, ...
        trialSummary, ...
        iedPredictorMode, ...
        maximumIEDCountForVisualization);

figurePDFOutputFile = fullfile( ...
    outputFolderName, 'IT_cox_model_visualization.pdf');

figurePNGOutputFile = fullfile( ...
    outputFolderName, 'IT_cox_model_visualization.png');

effectCurveOutputFile = fullfile( ...
    outputFolderName, 'IT_IED_model_effect_curve.csv');

exportgraphics( ...
    visualizationFigure, ...
    figurePDFOutputFile, ...
    'ContentType', 'vector');

exportgraphics( ...
    visualizationFigure, ...
    figurePNGOutputFile, ...
    'Resolution', 300);

writetable(IED_effect_curve, effectCurveOutputFile);

close(visualizationFigure);

%% Save numerical outputs

countingProcessOutputFile = fullfile( ...
    outputFolderName, 'IT_counting_process_data.csv');

trialSummaryOutputFile = fullfile( ...
    outputFolderName, 'IT_trial_summary.csv');

coxResultsOutputFile = fullfile( ...
    outputFolderName, 'IT_cox_results.csv');

baselineHazardOutputFile = fullfile( ...
    outputFolderName, 'IT_baseline_cumulative_hazard.csv');

modelOutputFile = fullfile( ...
    outputFolderName, 'IT_cox_model.mat');

writetable(countingProcessData, countingProcessOutputFile);
writetable(trialSummary, trialSummaryOutputFile);
writetable(coxResults, coxResultsOutputFile);
writematrix(baselineCumulativeHazard, baselineHazardOutputFile);

save( ...
    modelOutputFile, ...
    'beta', ...
    'logLikelihood', ...
    'baselineCumulativeHazard', ...
    'stats', ...
    'clusterRobustCovariance', ...
    'clusterRobustSE', ...
    'clusterRobustZ', ...
    'clusterRobustP', ...
    'predictorNames', ...
    'colorLevels', ...
    'referenceColor', ...
    'referenceColorName', ...
    'patientLevels', ...
    'iedPredictorMode', ...
    'defaultSamplingFrequencyHz', ...
    'maximumIEDCountForVisualization', ...
    'IED_effect_curve');

fprintf('\nsaved:\n');
fprintf('%s\n', countingProcessOutputFile);
fprintf('%s\n', trialSummaryOutputFile);
fprintf('%s\n', coxResultsOutputFile);
fprintf('%s\n', baselineHazardOutputFile);
fprintf('%s\n', modelOutputFile);
fprintf('%s\n', effectCurveOutputFile);
fprintf('%s\n', figurePDFOutputFile);
fprintf('%s\n', figurePNGOutputFile);

%% Local functions

function nTrials = determineNumberOfTrials( ...
    LFPIED, RTs, ITs, isControl, balloonType)

    vectorLengths = [ ...
        length(RTs), ...
        length(ITs), ...
        length(isControl), ...
        length(balloonType) ...
    ];

    if isfield(LFPIED, 'nTrials') && ...
            isscalar(LFPIED.nTrials) && ...
            isfinite(LFPIED.nTrials)

        vectorLengths(end + 1) = LFPIED.nTrials;

    end

    nTrials = floor(min(vectorLengths));

end

function samplingFrequencyHz = getSamplingFrequency( ...
    LFPIED, defaultSamplingFrequencyHz)

    if isfield(LFPIED, 'Fs') && ...
            isscalar(LFPIED.Fs) && ...
            isfinite(LFPIED.Fs) && ...
            LFPIED.Fs > 0

        samplingFrequencyHz = double(LFPIED.Fs);

    else

        samplingFrequencyHz = defaultSamplingFrequencyHz;

    end

end

function balloonColorCode = mapBalloonColorCode(balloonType)

    balloonType = double(balloonType(:));

    balloonColorCode = NaN(size(balloonType));

    validRows = ...
        isfinite(balloonType) & ...
        ismember(round(balloonType), [1 2 3 11 12 13]);

    balloonColorCode(validRows) = ...
        mod(round(balloonType(validRows)) - 1, 10) + 1;

end

function colorName = getBalloonColorName(colorCode)

    switch double(colorCode)
        case 1
            colorName = "yellow";
        case 2
            colorName = "orange";
        case 3
            colorName = "red";
        otherwise
            error('unexpected balloon color code.');
    end

end

function predictorName = getIEDPredictorName(iedPredictorMode)

    if iedPredictorMode == "count"
        predictorName = "cumulative_ied_count";
    elseif iedPredictorMode == "any"
        predictorName = "any_prior_ied";
    else
        error('iedpredictormode must be "count" or "any".');
    end

end

function [robustCovariance, robustSE, robustZ, robustP] = ...
    computeClusterRobustInference(stats, beta, clusterID)

    beta = beta(:);
    nPredictors = length(beta);

    clusterID = clusterID(:);
    uniqueClusters = unique(clusterID(isfinite(clusterID)));
    nClusters = length(uniqueClusters);

    canUseClusterScores = ...
        isfield(stats, 'scores') && ...
        size(stats.scores, 1) == length(clusterID) && ...
        size(stats.scores, 2) == nPredictors;

    if nClusters < 2 || ~canUseClusterScores

        robustCovariance = stats.covb;
        robustSE = stats.se(:);
        robustZ = stats.z(:);
        robustP = stats.p(:);
        return;

    end

    clusterScoreSums = zeros(nClusters, nPredictors);

    for gg = 1:nClusters

        rows = clusterID == uniqueClusters(gg);

        thisClusterScores = stats.scores(rows, :);
        thisClusterScores(~isfinite(thisClusterScores)) = 0;

        clusterScoreSums(gg, :) = sum(thisClusterScores, 1);

    end

    meat = clusterScoreSums' * clusterScoreSums;
    bread = stats.covb;

    robustCovariance = bread * meat * bread;

    % Finite-cluster correction.
    robustCovariance = ...
        (nClusters / (nClusters - 1)) .* robustCovariance;

    robustSE = sqrt(max(diag(robustCovariance), 0));

    zeroSERows = robustSE <= 0 | ~isfinite(robustSE);
    robustSE(zeroSERows) = stats.se(zeroSERows);

    robustZ = beta ./ robustSE;
    robustP = 2 .* normcdf(-abs(robustZ), 0, 1);

end

function trialRows = makeCountingProcessRows( ...
    patientID, ...
    trialNumber, ...
    durationSeconds, ...
    IEDoccurrence, ...
    samplingFrequencyHz, ...
    balloonColorCode, ...
    finalEventObserved, ...
    iedPredictorMode)

    % IEDoccurrence columns:
    %   1 = trial
    %   2 = channel
    %   3 = sample index within the IT interval

    if isempty(IEDoccurrence) || size(IEDoccurrence, 2) < 3

        IEDtimes = [];

    else

        validIEDRows = ...
            isfinite(IEDoccurrence(:, 1)) & ...
            isfinite(IEDoccurrence(:, 3)) & ...
            round(IEDoccurrence(:, 1)) == trialNumber & ...
            IEDoccurrence(:, 3) >= 1;

        sampleIndices = double(IEDoccurrence(validIEDRows, 3));

        % Sample index 1 is represented as 1/Fs seconds after interval
        % onset, which avoids a zero-length first interval.
        IEDtimes = sampleIndices ./ samplingFrequencyHz;

        IEDtimes = IEDtimes( ...
            isfinite(IEDtimes) & ...
            IEDtimes > 0 & ...
            IEDtimes < durationSeconds);

    end

    if isempty(IEDtimes)

        uniqueIEDtimes = [];
        numberAtEachTime = [];

    else

        [uniqueIEDtimes, ~, groupIndex] = unique(IEDtimes, 'sorted');

        numberAtEachTime = accumarray( ...
            groupIndex, ...
            1, ...
            [length(uniqueIEDtimes), 1], ...
            @sum, ...
            0);

    end

    breakTimes = [0; uniqueIEDtimes(:); durationSeconds];

    nIntervals = length(breakTimes) - 1;

    tStart = breakTimes(1:end-1);
    tStop = breakTimes(2:end);

    IEDcountSoFar = zeros(nIntervals, 1);
    IEDtimevarying = zeros(nIntervals, 1);

    cumulativeCount = 0;

    for kk = 1:nIntervals

        IEDcountSoFar(kk) = cumulativeCount;

        if iedPredictorMode == "count"
            IEDtimevarying(kk) = cumulativeCount;
        elseif iedPredictorMode == "any"
            IEDtimevarying(kk) = double(cumulativeCount > 0);
        else
            error('iedpredictormode must be "count" or "any".');
        end

        if kk <= length(numberAtEachTime)
            cumulativeCount = cumulativeCount + numberAtEachTime(kk);
        end

    end

    censored = true(nIntervals, 1);
    eventAtStop = false(nIntervals, 1);

    if finalEventObserved
        censored(end) = false;
        eventAtStop(end) = true;
    end

    trialRows = table( ...
        repmat(patientID, nIntervals, 1), ...
        repmat(trialNumber, nIntervals, 1), ...
        tStart, ...
        tStop, ...
        censored, ...
        eventAtStop, ...
        IEDtimevarying, ...
        IEDcountSoFar, ...
        repmat(balloonColorCode, nIntervals, 1), ...
        repmat(samplingFrequencyHz, nIntervals, 1), ...
        'VariableNames', { ...
            'patientID', ...
            'trialNumber', ...
            'tStart', ...
            'tStop', ...
            'censored', ...
            'eventAtStop', ...
            'IED_timevarying', ...
            'IED_count_so_far', ...
            'balloonColorCode', ...
            'samplingFrequencyHz' ...
        });

end

function [fig, effectCurveTable] = plotCoxModelVisualization( ...
    coxResults, ...
    trialSummary, ...
    iedPredictorMode, ...
    maximumIEDCountForVisualization)

    colorIT = [0.902 0.624 0.000];   % orange
    colorITLight = 0.75 .* [1 1 1] + 0.25 .* colorIT;

    fig = figure('Visible', 'off');
    set(fig, 'Position', [100 100 1250 520]);
    set(fig, 'Color', 'w');

    tiledlayout(1, 2, ...
        'Padding', 'compact', ...
        'TileSpacing', 'compact');

    %% Panel 1: hazard-ratio forest plot

    nexttile;
    hold on;

    nPredictors = height(coxResults);
    yPositions = (1:nPredictors)';

    hazardRatio = coxResults.hazardRatio;
    lowerError = hazardRatio - coxResults.hazardRatioCILow;
    upperError = coxResults.hazardRatioCIHigh - hazardRatio;

    errorbar( ...
        hazardRatio, ...
        yPositions, ...
        lowerError, ...
        upperError, ...
        'horizontal', ...
        'o', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colorIT, ...
        'MarkerEdgeColor', colorIT, ...
        'Color', colorIT, ...
        'CapSize', 7);

    xline(1, '--', 'Color', colorIT, 'LineWidth', 1);

    predictorLabels = makePrettyPredictorLabels( ...
        coxResults.predictor, ...
        iedPredictorMode);

    set(gca, ...
        'YTick', yPositions, ...
        'YTickLabel', predictorLabels, ...
        'YDir', 'reverse', ...
        'XScale', 'log', ...
        'FontSize', 10, ...
        'FontName', 'Arial', ...
        'XColor', 'k', ...
        'YColor', 'k');

    xlabel('hazard ratio with robust 95% ci', 'Color', 'k');
    title('cox model estimates', 'Color', 'k');

    xLimits = calculatePositiveLogLimits( ...
        coxResults.hazardRatioCILow, ...
        coxResults.hazardRatioCIHigh);

    xlim(xLimits);
    ylim([0.5, nPredictors + 0.5]);

    grid off;
    box off;
    hold off;

    %% Panel 2: primary IED model effect

    nexttile;
    hold on;

    betaIED = coxResults.beta_logHazard(1);
    seIED = coxResults.clusterRobustSE(1);

    betaLow = betaIED - 1.96 * seIED;
    betaHigh = betaIED + 1.96 * seIED;

    if iedPredictorMode == "count"

        trialIEDCounts = trialSummary.nIEDsInFullDuration;
        trialIEDCounts = trialIEDCounts(isfinite(trialIEDCounts));

        if isempty(trialIEDCounts)
            displayedMaximumCount = 1;
        else
            percentile95 = prctile(trialIEDCounts, 95);
            displayedMaximumCount = max(1, ceil(percentile95));
            displayedMaximumCount = min( ...
                displayedMaximumCount, ...
                maximumIEDCountForVisualization);
        end

        predictorValue = (0:displayedMaximumCount)';

        relativeHazard = exp(betaIED .* predictorValue);
        relativeHazardCILow = exp(betaLow .* predictorValue);
        relativeHazardCIHigh = exp(betaHigh .* predictorValue);

        fill( ...
            [predictorValue; flipud(predictorValue)], ...
            [relativeHazardCILow; flipud(relativeHazardCIHigh)], ...
            colorITLight, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.75);

        plot( ...
            predictorValue, ...
            relativeHazard, ...
            '-o', ...
            'LineWidth', 1.8, ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', colorIT, ...
            'MarkerEdgeColor', colorIT, ...
            'Color', colorIT);

        xlabel('cumulative prior ied count', 'Color', 'k');
        title('model-implied ied effect', 'Color', 'k');

    elseif iedPredictorMode == "any"

        predictorValue = [0; 1];

        relativeHazard = [1; exp(betaIED)];
        relativeHazardCILow = [1; exp(betaLow)];
        relativeHazardCIHigh = [1; exp(betaHigh)];

        plot( ...
            predictorValue, ...
            relativeHazard, ...
            '-', ...
            'LineWidth', 1.5, ...
            'Color', colorIT);

        errorbar( ...
            predictorValue, ...
            relativeHazard, ...
            relativeHazard - relativeHazardCILow, ...
            relativeHazardCIHigh - relativeHazard, ...
            'o', ...
            'LineStyle', 'none', ...
            'LineWidth', 1.2, ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', colorIT, ...
            'MarkerEdgeColor', colorIT, ...
            'Color', colorIT, ...
            'CapSize', 7);

        xlim([-0.25 1.25]);
        xticks([0 1]);
        xticklabels({'no prior ied', 'at least one prior ied'});
        xlabel('time-varying ied state', 'Color', 'k');
        title('model-implied ied effect', 'Color', 'k');

    else
        error('iedpredictormode must be "count" or "any".');
    end

    yline(1, '--', 'Color', colorIT, 'LineWidth', 1);

    set(gca, ...
        'YScale', 'log', ...
        'FontSize', 10, ...
        'FontName', 'Arial', ...
        'XColor', 'k', ...
        'YColor', 'k');

    ylabel('relative action hazard', 'Color', 'k');
    subtitle('reference = zero prior ieds');

    yLimits = calculatePositiveLogLimits( ...
        relativeHazardCILow, ...
        relativeHazardCIHigh);

    ylim(yLimits);

    grid off;
    box off;
    hold off;

    effectCurveTable = table( ...
        predictorValue, ...
        relativeHazard, ...
        relativeHazardCILow, ...
        relativeHazardCIHigh, ...
        'VariableNames', { ...
            'IED_predictor_value', ...
            'relativeHazard', ...
            'relativeHazardCILow', ...
            'relativeHazardCIHigh' ...
        });

end

function predictorLabels = makePrettyPredictorLabels( ...
    predictorNames, ...
    iedPredictorMode)

    predictorNames = string(predictorNames(:));
    predictorLabels = strings(size(predictorNames));

    for ii = 1:length(predictorNames)

        thisName = predictorNames(ii);

        if thisName == "cumulative_ied_count"

            predictorLabels(ii) = "cumulative prior ied count";

        elseif thisName == "any_prior_ied"

            predictorLabels(ii) = "any prior ied";

        elseif thisName == "orange_vs_yellow"

            predictorLabels(ii) = "orange vs yellow";

        elseif thisName == "red_vs_yellow"

            predictorLabels(ii) = "red vs yellow";

        elseif contains(thisName, "_vs_")

            predictorLabels(ii) = lower(replace(thisName, "_", " "));

        else

            predictorLabels(ii) = lower(replace(thisName, "_", " "));

        end

    end

    if iedPredictorMode == "count" && ~isempty(predictorLabels)
        predictorLabels(1) = "cumulative prior ied count";
    elseif iedPredictorMode == "any" && ~isempty(predictorLabels)
        predictorLabels(1) = "any prior ied";
    end

    predictorLabels = lower(predictorLabels);

end

function axisLimits = calculatePositiveLogLimits(lowerValues, upperValues)

    lowerValues = lowerValues(:);
    upperValues = upperValues(:);

    positiveLower = lowerValues(isfinite(lowerValues) & lowerValues > 0);
    positiveUpper = upperValues(isfinite(upperValues) & upperValues > 0);

    if isempty(positiveLower) || isempty(positiveUpper)
        axisLimits = [0.5 2];
        return;
    end

    minimumValue = min([positiveLower; positiveUpper; 1]);
    maximumValue = max([positiveLower; positiveUpper; 1]);

    logMinimum = log10(minimumValue);
    logMaximum = log10(maximumValue);

    if logMaximum == logMinimum
        logMinimum = logMinimum - 0.2;
        logMaximum = logMaximum + 0.2;
    else
        padding = 0.10 * (logMaximum - logMinimum);
        logMinimum = logMinimum - padding;
        logMaximum = logMaximum + padding;
    end

    axisLimits = [10^logMinimum, 10^logMaximum];

end
