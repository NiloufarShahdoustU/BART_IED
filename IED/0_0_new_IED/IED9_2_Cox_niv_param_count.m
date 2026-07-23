% Cox models for active IED count x Niv learning asymmetry
% Permutation test shuffles whole trials, not interval rows
% Author: Nill

clear;
clc;
close all;

inputFolder = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

nivFile = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\context_modeling\param_recovery_4_param_recovery\alpha_comparison.csv';

outputFolder = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED9_2_Cox_niv_param_count\';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

maximumRT = 20;
defaultSamplingFrequency = 1000;
numberOfPermutations = 500;
numberOfNivCurvePoints = 200;
randomSeed = 42;
iedMode = 'count';

colorIT = [0.847 0.333 0.153];
colorRT = [0.204 0.459 0.702];
colorBR = [0.250 0.600 0.250];

configs(1).code = 'IT';
configs(1).durationField = 'ITs';
configs(1).iedField = 'IED_occurance_IT';
configs(1).bankedOnly = false;
configs(1).windowMilliseconds = 1000;
configs(1).color = colorIT;

configs(2).code = 'RT';
configs(2).durationField = 'RTs';
configs(2).iedField = 'IED_occurance_RT';
configs(2).bankedOnly = false;
configs(2).windowMilliseconds = 500;
configs(2).color = colorRT;

configs(3).code = 'BR';
configs(3).durationField = 'ITs';
configs(3).iedField = 'IED_occurance_IT';
configs(3).bankedOnly = true;
configs(3).windowMilliseconds = 1000;
configs(3).color = colorBR;

textOutputFile = fullfile(outputFolder, ...
    'IT_RT_BR_whole_IED_x_Niv_console_output.txt');

diary off;
if exist(textOutputFile, 'file')
    delete(textOutputFile);
end
diary(textOutputFile);

%% Read Niv parameters

[nivTable, nivMean, nivSD] = readNivTable(nivFile);
writetable(nivTable, fullfile(outputFolder, 'niv_parameters_used.csv'));

fprintf('Loaded Niv parameters for %d participants\n', height(nivTable));
fprintf('Mean Niv = %.6f, SD = %.6f\n', nivMean, nivSD);

%% Plot participant Niv values

nivFigure = plotNivParameters(nivTable);
nivFigureOutputFile = fullfile(outputFolder, ...
    'Niv_alpha_plus_alpha_minus_participant_visualization.pdf');
exportgraphics(nivFigure, nivFigureOutputFile, 'ContentType', 'vector');
close(nivFigure);

%% Run IT, RT, and BR

results = cell(3, 1);

for a = 1:3
    fprintf('\n\nRUNNING %s IED x NIV ANALYSIS\n', configs(a).code);

    results{a} = runNivCox(inputFolder, configs(a), nivTable, ...
        nivMean, nivSD, maximumRT, defaultSamplingFrequency, ...
        numberOfNivCurvePoints, numberOfPermutations, ...
        randomSeed + a - 1, iedMode);

    saveOutcomeOutputs(results{a}, outputFolder);
end

IT = results{1};
RT = results{2};
BR = results{3};

%% Make one 3 x 2 figure

combinedFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [60 40 1550 1500]);

layout = tiledlayout(combinedFigure, 3, 2, ...
    'Padding', 'compact', 'TileSpacing', 'compact');

for a = 1:3
    plotForest(nexttile(layout), results{a}, configs(a));
    plotNivCurve(nexttile(layout), results{a}, configs(a), iedMode);
end

title(layout, 'Whole-IED effects moderated by Niv learning asymmetry', ...
    'FontWeight', 'bold', 'FontSize', 16);

combinedFigureOutputFile = fullfile(outputFolder, ...
    'IT_RT_BR_whole_IED_x_Niv_3x2.pdf');
exportgraphics(combinedFigure, combinedFigureOutputFile, ...
    'ContentType', 'vector');
close(combinedFigure);

combinedModelOutputFile = fullfile(outputFolder, ...
    'IT_RT_BR_whole_IED_x_Niv_all_results.mat');

save(combinedModelOutputFile, 'IT', 'RT', 'BR', 'nivTable', ...
    'nivMean', 'nivSD', 'configs', 'iedMode', ...
    'numberOfPermutations', 'randomSeed', '-v7.3');

fprintf('\ndoneeeeeee\n');
fprintf('Combined figure: %s\n', combinedFigureOutputFile);
fprintf('Combined results: %s\n', combinedModelOutputFile);
fprintf('Niv figure: %s\n', nivFigureOutputFile);
fprintf('Console output: %s\n', textOutputFile);
diary off;


%% Run one outcome

function result = runNivCox(inputFolder, config, nivTable, ...
    nivMean, nivSD, maximumRT, defaultSamplingFrequency, ...
    numberOfCurvePoints, numberOfPermutations, randomSeed, iedMode)

    code = config.code;
    windowMilliseconds = config.windowMilliseconds;
    windowSeconds = windowMilliseconds / 1000;
    files = dir(fullfile(inputFolder, '*.LFPIED.mat'));

    countingProcessData = table();
    trialSummary = table();
    permutationTrialData = table();

    for pt = 1:length(files)
        fileName = files(pt).name;
        parts = strsplit(fileName, '.');
        patientID = string(parts{1});

        nivRow = find(nivTable.normalizedPatientID == ...
            normalizeID(patientID), 1);

        if isempty(nivRow)
            fprintf('%s: no Niv values, skipped\n', patientID);
            continue;
        end

        alphaPlus = nivTable.alphaPlus(nivRow);
        alphaMinus = nivTable.alphaMinus(nivRow);
        nivRaw = nivTable.nivAsymmetryRaw(nivRow);
        nivZ = nivTable.nivAsymmetryZ(nivRow);

        fprintf('%s: processing patient %s, Niv = %.4f\n', ...
            code, patientID, nivRaw);

        data = load(fullfile(inputFolder, fileName));
        LFPIED = data.LFPIED;

        RTs = double(LFPIED.RTs(:));
        durations = LFPIED.(config.durationField);
        durations = double(durations(:));
        isControl = double(LFPIED.isControl(:));

        balloonColor = round(double(LFPIED.balloonType(:)));
        balloonColor = mod(balloonColor - 1, 10) + 1;

        IEDoccurrence = LFPIED.(config.iedField);
        nTrials = length(durations);

        if config.bankedOnly
            bankedTrials = double(LFPIED.BankedTrials(:));
        else
            bankedTrials = NaN(nTrials, 1);
        end

        samplingFrequency = defaultSamplingFrequency;
        if isfield(LFPIED, 'Fs')
            samplingFrequency = double(LFPIED.Fs);
        end

        for trial = 1:nTrials
            duration = durations(trial);

            keepTrial = isControl(trial) == 0 && ...
                RTs(trial) > 0 && RTs(trial) <= maximumRT && ...
                duration > 0 && ismember(balloonColor(trial), [1 2 3]);

            if config.bankedOnly
                keepTrial = keepTrial && ismember(bankedTrials(trial), [0 1]);
                eventObserved = bankedTrials(trial) == 1;
            else
                eventObserved = true;
            end

            if ~keepTrial
                continue;
            end

            rows = round(IEDoccurrence(:, 1)) == trial;
            IEDtimes = sort(double(IEDoccurrence(rows, 3)) / samplingFrequency);
            IEDtimes = IEDtimes(:);
            IEDtimes = IEDtimes(IEDtimes > 0 & IEDtimes < duration);

            [tStart, tStop, postIED, activeIEDCount, timeSinceIED] = ...
                makeIntervals(duration, IEDtimes, windowSeconds);

            censored = true(length(tStart), 1);
            eventAtStop = false(length(tStart), 1);
            if eventObserved
                censored(end) = false;
                eventAtStop(end) = true;
            end

            n = length(tStart);

            trialRows = table( ...
                repmat(patientID, n, 1), repmat(trial, n, 1), ...
                tStart, tStop, censored, eventAtStop, postIED, ...
                activeIEDCount, timeSinceIED, ...
                repmat(balloonColor(trial), n, 1), ...
                repmat(samplingFrequency, n, 1), ...
                repmat(alphaPlus, n, 1), repmat(alphaMinus, n, 1), ...
                repmat(nivRaw, n, 1), repmat(nivZ, n, 1), ...
                'VariableNames', {'patientID', 'trialNumber', ...
                'tStart', 'tStop', 'censored', 'eventAtStop', ...
                'postIED', 'activeIEDWindowCount', ...
                'timeSinceMostRecentIED_seconds', 'balloonColorCode', ...
                'samplingFrequencyHz', 'alphaPlus', 'alphaMinus', ...
                'nivAsymmetryRaw', 'nivAsymmetryZ'});

            countingProcessData = [countingProcessData; trialRows];

            nIEDs = length(IEDtimes);
            exposure = sum((tStop - tStart) .* postIED);
            eventInWindow = double(eventObserved && postIED(end) == 1);

            newSummary = table(patientID, trial, duration, ...
                double(eventObserved), eventInWindow, nIEDs, exposure, ...
                balloonColor(trial), samplingFrequency, alphaPlus, ...
                alphaMinus, nivRaw, nivZ, ...
                'VariableNames', {'patientID', 'trialNumber', ...
                'durationSeconds', 'eventObserved', ...
                'eventInPostIEDWindow', 'nIEDsInFullDuration', ...
                'postIEDExposureSeconds', 'balloonColorCode', ...
                'samplingFrequencyHz', 'alphaPlus', 'alphaMinus', ...
                'nivAsymmetryRaw', 'nivAsymmetryZ'});

            trialSummary = [trialSummary; newSummary];

            newPermutationTrial = table(patientID, trial, duration, ...
                double(eventObserved), balloonColor(trial), nivZ, ...
                {IEDtimes}, 'VariableNames', {'patientID', ...
                'trialNumber', 'durationSeconds', 'eventObserved', ...
                'balloonColorCode', 'nivAsymmetryZ', 'IEDtimes_seconds'});

            permutationTrialData = [permutationTrialData; newPermutationTrial];
        end
    end

    %% Fit the real model

    [countingProcessData.patientStratum, patientLevels] = ...
        findgroups(countingProcessData.patientID);

    permutationTrialData.patientStratum = ...
        findgroups(permutationTrialData.patientID);

    [X, predictorNames, displayLabels] = ...
        makePredictors(countingProcessData, iedMode);

    T = [countingProcessData.tStart, countingProcessData.tStop];
    censoring = countingProcessData.censored;
    strata = countingProcessData.patientStratum;

    options = statset('coxphfit');
    options.Display = 'final';
    options.MaxIter = 1000;
    options.MaxFunEvals = 5000;

    [beta, logLikelihood, baselineCumulativeHazard, stats] = coxphfit( ...
        X, T, 'Censoring', censoring, 'Strata', strata, ...
        'Ties', 'efron', 'Baseline', 0, 'Options', options);

    %% Participant-clustered confidence intervals

    patientNumbers = unique(strata);
    clusterScores = zeros(length(patientNumbers), length(beta));

    for p = 1:length(patientNumbers)
        clusterScores(p, :) = ...
            sum(stats.scores(strata == patientNumbers(p), :), 1);
    end

    robustCovariance = stats.covb * (clusterScores' * clusterScores) * ...
        stats.covb * length(patientNumbers) / (length(patientNumbers) - 1);
    robustSE = sqrt(max(diag(robustCovariance), 0));

    %% Permutation test

    permutationBeta = NaN(numberOfPermutations, 2);
    permutationSucceeded = false(numberOfPermutations, 1);

    groups = findgroups(permutationTrialData.patientID, ...
        permutationTrialData.balloonColorCode);

    rng(randomSeed);
    permutationOptions = options;
    permutationOptions.Display = 'off';

    fprintf('%s: running %d whole-trial permutations...\n', ...
        code, numberOfPermutations);

    for permutation = 1:numberOfPermutations
        shuffledTrials = permutationTrialData;

        for group = 1:max(groups)
            groupRows = find(groups == group);
            sourceRows = groupRows(randperm(length(groupRows)));

            % Duration and event status always move together.
            shuffledTrials.durationSeconds(groupRows) = ...
                permutationTrialData.durationSeconds(sourceRows);
            shuffledTrials.eventObserved(groupRows) = ...
                permutationTrialData.eventObserved(sourceRows);
        end

        [permX, permT, permCensoring, permStrata] = ...
            buildPermutationData(shuffledTrials, windowSeconds, iedMode);

        try
            b = coxphfit(permX, permT, 'Censoring', permCensoring, ...
                'Strata', permStrata, 'Ties', 'efron', ...
                'Baseline', 0, 'Options', permutationOptions);

            permutationBeta(permutation, :) = b(1:2);
            permutationSucceeded(permutation) = true;
        catch
            fprintf('%s permutation %d failed\n', code, permutation);
        end

        if mod(permutation, 100) == 0 || permutation == numberOfPermutations
            fprintf('%s permutation progress: %d/%d\n', ...
                code, permutation, numberOfPermutations);
        end
    end

    validPermutationBeta = permutationBeta(permutationSucceeded, :);
    permutationP = NaN(length(beta), 1);

    for term = 1:2
        permutationP(term) = ...
            (1 + sum(abs(validPermutationBeta(:, term)) >= ...
            abs(beta(term)))) / (size(validPermutationBeta, 1) + 1);
    end

    permutationResults = table((1:numberOfPermutations)', ...
        permutationBeta(:, 1), permutationBeta(:, 2), ...
        permutationSucceeded, ...
        'VariableNames', {'permutationIndex', 'permutedIEDBeta', ...
        'permutedInteractionBeta', 'fitSucceeded'});

    %% Results

    betaLow = beta - 1.96 * robustSE;
    betaHigh = beta + 1.96 * robustSE;
    hazardRatio = exp(beta);

    coxResults = table(predictorNames, displayLabels, beta(:), ...
        stats.se(:), stats.z(:), stats.p(:), robustSE(:), ...
        permutationP, betaLow(:), betaHigh(:), hazardRatio(:), ...
        exp(betaLow(:)), exp(betaHigh(:)), ...
        'VariableNames', {'predictor', 'displayLabel', ...
        'beta_logHazard', 'modelBasedSE', 'modelBasedZ', ...
        'modelBasedPValue', 'clusterRobustSE', 'permutationPValue', ...
        'betaCILow', 'betaCIHigh', 'hazardRatio', ...
        'hazardRatioCILow', 'hazardRatioCIHigh'});

    patientNiv = unique(countingProcessData(:, ...
        {'patientID', 'nivAsymmetryRaw', 'nivAsymmetryZ'}), ...
        'rows', 'stable');

    nivCurveTable = makeNivCurve(beta, robustCovariance, ...
        patientNiv.nivAsymmetryRaw, nivMean, nivSD, numberOfCurvePoints);

    rowDuration = countingProcessData.tStop - countingProcessData.tStart;
    nEventsOutside = sum(countingProcessData.eventAtStop & ...
        countingProcessData.postIED == 0);
    nEventsInside = sum(countingProcessData.eventAtStop & ...
        countingProcessData.postIED == 1);
    timeOutside = sum(rowDuration(countingProcessData.postIED == 0));
    timeInside = sum(rowDuration(countingProcessData.postIED == 1));

    fprintf('\n%s IED x Niv Cox model\n', code);
    fprintf('participants: %d\n', length(patientLevels));
    fprintf('trials: %d\n', height(trialSummary));
    fprintf('events outside/inside window: %d / %d\n', ...
        nEventsOutside, nEventsInside);
    fprintf('time outside/inside window: %.3f / %.3f s\n', ...
        timeOutside, timeInside);
    fprintf('log likelihood: %.6f\n', logLikelihood);
    disp(coxResults);
    fprintf('IED x Niv permutation p = %.6g\n', permutationP(2));

    result.analysisType = code;
    result.postIEDWindowMilliseconds = windowMilliseconds;
    result.countingProcessData = countingProcessData;
    result.trialSummary = trialSummary;
    result.coxResults = coxResults;
    result.nivCurveTable = nivCurveTable;
    result.observedPatientNiv = patientNiv;
    result.beta = beta;
    result.logLikelihood = logLikelihood;
    result.baselineCumulativeHazard = baselineCumulativeHazard;
    result.stats = stats;
    result.clusterRobustCovariance = robustCovariance;
    result.clusterRobustSE = robustSE;
    result.predictorNames = predictorNames;
    result.patientLevels = patientLevels;
    result.permutationResults = permutationResults;
    result.permutationPValue = permutationP;
    result.numberOfPermutations = numberOfPermutations;
    result.numberOfValidPermutations = size(validPermutationBeta, 1);
end


%% Make intervals for one trial

function [tStart, tStop, postIED, activeIEDCount, timeSinceIED] = ...
    makeIntervals(duration, IEDtimes, windowSeconds)

    IEDtimes = IEDtimes(:);
    windowEnds = min(IEDtimes + windowSeconds, duration);
    breakTimes = unique([0; IEDtimes; windowEnds; duration]);

    tStart = breakTimes(1:end-1);
    tStop = breakTimes(2:end);
    midpoint = (tStart + tStop) / 2;

    activeIEDCount = zeros(length(midpoint), 1);

    for ied = 1:length(IEDtimes)
        inside = midpoint >= IEDtimes(ied) & ...
            midpoint <= IEDtimes(ied) + windowSeconds;
        activeIEDCount(inside) = activeIEDCount(inside) + 1;
    end

    postIED = double(activeIEDCount > 0);
    timeSinceIED = NaN(length(midpoint), 1);

    for k = 1:length(midpoint)
        previousIEDs = IEDtimes(IEDtimes < midpoint(k));
        if ~isempty(previousIEDs)
            timeSinceIED(k) = midpoint(k) - previousIEDs(end);
        end
    end
end


%% Predictor matrix

function [X, predictorNames, displayLabels] = ...
    makePredictors(data, iedMode)

    if strcmp(iedMode, 'count')
        ied = double(data.activeIEDWindowCount);
        iedName = "active_ied_count";
        interactionName = "active_ied_count_x_niv_z";
        iedLabel = "Active IED count";
        interactionLabel = "IED count x Niv (1 SD)";
    else
        ied = double(data.postIED);
        iedName = "any_ied_occurrence";
        interactionName = "any_ied_occurrence_x_niv_z";
        iedLabel = "Any IED occurrence";
        interactionLabel = "IED occurrence x Niv (1 SD)";
    end

    X = [ied, ied .* data.nivAsymmetryZ, ...
        double(data.balloonColorCode == 2), ...
        double(data.balloonColorCode == 3)];

    predictorNames = [iedName; interactionName; ...
        "orange_vs_yellow"; "red_vs_yellow"];

    displayLabels = [iedLabel; interactionLabel; ...
        "Orange vs yellow"; "Red vs yellow"];
end


%% Rebuild Cox data after a whole-trial shuffle

function [X, T, censoring, strata] = ...
    buildPermutationData(trialData, windowSeconds, iedMode)

    nTrials = height(trialData);
    Xcells = cell(nTrials, 1);
    Tcells = cell(nTrials, 1);
    censoringCells = cell(nTrials, 1);
    strataCells = cell(nTrials, 1);

    for trial = 1:nTrials
        duration = trialData.durationSeconds(trial);
        IEDtimes = trialData.IEDtimes_seconds{trial};
        IEDtimes = IEDtimes(:);
        IEDtimes = IEDtimes(IEDtimes < duration);

        [tStart, tStop, postIED, activeIEDCount] = ...
            makeIntervals(duration, IEDtimes, windowSeconds);

        n = length(tStart);
        color = trialData.balloonColorCode(trial);
        nivZ = trialData.nivAsymmetryZ(trial);

        if strcmp(iedMode, 'count')
            ied = activeIEDCount(:);
        else
            ied = postIED(:);
        end

        Xcells{trial} = [ied, ied .* nivZ, ...
            repmat(double(color == 2), n, 1), ...
            repmat(double(color == 3), n, 1)];

        Tcells{trial} = [tStart(:), tStop(:)];

        censored = true(n, 1);
        if trialData.eventObserved(trial) == 1
            censored(end) = false;
        end

        censoringCells{trial} = censored;
        strataCells{trial} = repmat(trialData.patientStratum(trial), n, 1);
    end

    X = vertcat(Xcells{:});
    T = vertcat(Tcells{:});
    censoring = vertcat(censoringCells{:});
    strata = vertcat(strataCells{:});
end


%% Niv moderation curve

function curveTable = makeNivCurve(beta, covariance, observedNiv, ...
    nivMean, nivSD, numberOfPoints)

    nivRaw = linspace(min(observedNiv), max(observedNiv), numberOfPoints)';
    nivZ = (nivRaw - nivMean) / nivSD;

    logHR = beta(1) + beta(2) .* nivZ;
    variance = covariance(1, 1) + nivZ .^ 2 * covariance(2, 2) + ...
        2 * nivZ * covariance(1, 2);
    se = sqrt(max(variance, 0));

    curveTable = table(nivRaw, nivZ, logHR, se, ...
        exp(logHR), exp(logHR - 1.96 * se), exp(logHR + 1.96 * se), ...
        'VariableNames', {'nivAsymmetryRaw', 'nivAsymmetryZ', ...
        'logHazardRatioPerAdditionalIED', ...
        'standardErrorLogHazardRatio', ...
        'hazardRatioPerAdditionalIED', ...
        'hazardRatioCILow', 'hazardRatioCIHigh'});
end


%% Read the Niv CSV

function [nivTable, nivMean, nivSD] = readNivTable(nivFile)

    raw = readtable(nivFile, 'TextType', 'string', ...
        'VariableNamingRule', 'preserve');

    patientID = strip(string(raw.ptID));
    alphaPlus = str2double(string(raw.fit_alpha_plus));
    alphaMinus = str2double(string(raw.fit_alpha_minus));

    denominator = alphaMinus + alphaPlus;
    nivRaw = (alphaMinus - alphaPlus) ./ denominator;

    keep = strlength(patientID) > 0 & isfinite(alphaPlus) & ...
        isfinite(alphaMinus) & denominator > 0 & isfinite(nivRaw);

    patientID = patientID(keep);
    alphaPlus = alphaPlus(keep);
    alphaMinus = alphaMinus(keep);
    nivRaw = nivRaw(keep);
    normalizedPatientID = normalizeID(patientID);

    [~, firstRow] = unique(normalizedPatientID, 'stable');
    patientID = patientID(firstRow);
    normalizedPatientID = normalizedPatientID(firstRow);
    alphaPlus = alphaPlus(firstRow);
    alphaMinus = alphaMinus(firstRow);
    nivRaw = nivRaw(firstRow);

    nivMean = mean(nivRaw);
    nivSD = std(nivRaw);
    nivZ = (nivRaw - nivMean) / nivSD;

    nivTable = table(patientID, normalizedPatientID, alphaPlus, ...
        alphaMinus, nivRaw, nivZ, ...
        'VariableNames', {'patientID', 'normalizedPatientID', ...
        'alphaPlus', 'alphaMinus', 'nivAsymmetryRaw', 'nivAsymmetryZ'});
end


function normalized = normalizeID(patientID)
    normalized = lower(strip(string(patientID)));
    normalized = regexprep(normalized, '[^a-z0-9]', '');
end


%% Save one outcome

function saveOutcomeOutputs(result, outputFolder)

    code = char(result.analysisType);

    writetable(result.countingProcessData, fullfile(outputFolder, ...
        sprintf('%s_whole_IED_x_Niv_counting_process_data.csv', code)));
    writetable(result.trialSummary, fullfile(outputFolder, ...
        sprintf('%s_whole_IED_x_Niv_trial_summary.csv', code)));
    writetable(result.coxResults, fullfile(outputFolder, ...
        sprintf('%s_whole_IED_x_Niv_cox_results.csv', code)));
    writetable(result.nivCurveTable, fullfile(outputFolder, ...
        sprintf('%s_whole_IED_x_Niv_moderation_curve.csv', code)));
    writetable(result.observedPatientNiv, fullfile(outputFolder, ...
        sprintf('%s_whole_IED_x_Niv_patient_values.csv', code)));
    writetable(result.permutationResults, fullfile(outputFolder, ...
        sprintf('%s_whole_IED_x_Niv_permutation_results.csv', code)));
    writematrix(result.baselineCumulativeHazard, fullfile(outputFolder, ...
        sprintf('%s_whole_IED_x_Niv_baseline_cumulative_hazard.csv', code)));

    outcomeResult = result;
    save(fullfile(outputFolder, ...
        sprintf('%s_whole_IED_x_Niv_model.mat', code)), 'outcomeResult');
end


%% Cox forest plot

function plotForest(ax, result, config)

    data = result.coxResults;
    y = (1:height(data))';

    hold(ax, 'on');
    xline(ax, 1, '--', 'Color', [0.35 0.35 0.35]);

    for row = 1:height(data)
        plot(ax, [data.hazardRatioCILow(row), ...
            data.hazardRatioCIHigh(row)], [row row], '-', ...
            'Color', config.color, 'LineWidth', 1.6);
    end

    scatter(ax, data.hazardRatio, y, 64, config.color, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');

    set(ax, 'XScale', 'log', 'YDir', 'reverse', ...
        'YTick', y, 'YTickLabel', data.displayLabel, ...
        'TickLabelInterpreter', 'none', 'TickDir', 'out', ...
        'FontName', 'Arial', 'FontSize', 10);

    xlabel(ax, 'hazard ratio with robust 95% CI');
    title(ax, sprintf('%s: %d-ms post-IED Cox estimates', ...
        config.code, result.postIEDWindowMilliseconds), ...
        'FontWeight', 'bold');

    limits = [data.hazardRatioCILow; data.hazardRatioCIHigh; 1];
    xlim(ax, [0.8 * min(limits), 1.2 * max(limits)]);
    ylim(ax, [0.5 height(data) + 0.5]);

    pValues = data.modelBasedPValue;
    usePermutation = ~isnan(data.permutationPValue);
    pValues(usePermutation) = data.permutationPValue(usePermutation);

    for row = 1:height(data)
        if pValues(row) < 0.05
            text(ax, 1.06 * data.hazardRatioCIHigh(row), row, '*', ...
                'FontSize', 15, 'FontWeight', 'bold');
        end
    end

    box(ax, 'off');
    hold(ax, 'off');
end


%% Niv curve plot

function plotNivCurve(ax, result, config, iedMode)

    curve = result.nivCurveTable;
    interaction = result.coxResults(2, :);
    lightColor = 0.78 * [1 1 1] + 0.22 * config.color;

    hold(ax, 'on');
    fill(ax, [curve.nivAsymmetryRaw; flipud(curve.nivAsymmetryRaw)], ...
        [curve.hazardRatioCILow; flipud(curve.hazardRatioCIHigh)], ...
        lightColor, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    yline(ax, 1, '--', 'Color', [0.45 0.45 0.45]);
    xline(ax, 0, ':', 'Color', [0.55 0.55 0.55]);
    plot(ax, curve.nivAsymmetryRaw, ...
        curve.hazardRatioPerAdditionalIED, '-', ...
        'Color', config.color, 'LineWidth', 2.3);

    participantNiv = result.observedPatientNiv.nivAsymmetryRaw;
    rugY = repmat(min(curve.hazardRatioCILow), length(participantNiv), 1);
    scatter(ax, participantNiv, rugY, 14, config.color, 'filled', ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');

    set(ax, 'YScale', 'log', 'TickDir', 'out', ...
        'FontName', 'Arial', 'FontSize', 10);
    xlabel(ax, 'Niv asymmetry: (alpha- - alpha+) / (alpha- + alpha+)');

    if strcmp(iedMode, 'count')
        ylabel(ax, 'HR per additional active IED');
    else
        ylabel(ax, 'HR for any IED occurrence vs none');
    end

    title(ax, sprintf('%s: IED effect across Niv asymmetry', config.code), ...
        'FontWeight', 'bold');

    resultText = sprintf([ ...
        'IED x Niv beta = %.3f\n' ...
        'HR multiplier = %.3f\n' ...
        '95%% CI [%.3f, %.3f]\n' ...
        'permutation p = %.3g'], ...
        interaction.beta_logHazard, interaction.hazardRatio, ...
        interaction.hazardRatioCILow, interaction.hazardRatioCIHigh, ...
        interaction.permutationPValue);

    text(ax, 0.97, 0.96, resultText, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontName', 'Arial', 'FontSize', 9, 'BackgroundColor', 'w', ...
        'EdgeColor', [0.82 0.82 0.82], 'Margin', 5);

    box(ax, 'off');
    hold(ax, 'off');
end


%% Participant Niv visualization

function fig = plotNivParameters(nivTable)

    alphaPlus = nivTable.alphaPlus;
    alphaMinus = nivTable.alphaMinus;
    nivRaw = nivTable.nivAsymmetryRaw;
    patientID = nivTable.patientID;

    maximumAbsoluteNiv = max(abs(nivRaw));
    colorLimits = [-maximumAbsoluteNiv maximumAbsoluteNiv];

    blue = [0.12 0.35 0.75];
    white = [0.97 0.97 0.97];
    red = [0.78 0.18 0.18];
    colorMap = [linspace(blue(1), white(1), 128)', ...
        linspace(blue(2), white(2), 128)', ...
        linspace(blue(3), white(3), 128)'; ...
        linspace(white(1), red(1), 128)', ...
        linspace(white(2), red(2), 128)', ...
        linspace(white(3), red(3), 128)'];

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [40 40 1800 1250]);
    set(fig, 'Colormap', colorMap);
    layout = tiledlayout(fig, 1, 2, ...
        'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile(layout);
    hold(ax1, 'on');
    limit = 1.08 * max([alphaPlus; alphaMinus]);
    plot(ax1, [0 limit], [0 limit], '--', ...
        'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
    scatter(ax1, alphaPlus, alphaMinus, 92, nivRaw, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');
    set(ax1, 'CLim', colorLimits);
    xlim(ax1, [0 limit]);
    ylim(ax1, [0 limit]);
    axis(ax1, 'square');
    xlabel(ax1, 'Fitted positive learning rate, alpha+');
    ylabel(ax1, 'Fitted negative learning rate, alpha-');
    title(ax1, 'alpha+ versus alpha- by participant', ...
        'FontWeight', 'bold');
    colorbar(ax1);
    box(ax1, 'off');

    ax2 = nexttile(layout);
    [nivSorted, order] = sort(nivRaw);
    positions = (1:length(nivSorted))';
    scatter(ax2, nivSorted, positions, 62, nivSorted, 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');
    xline(ax2, 0, '--', 'Color', [0.25 0.25 0.25]);
    set(ax2, 'CLim', colorLimits, 'YDir', 'reverse', ...
        'YTick', positions, 'YTickLabel', patientID(order), ...
        'TickLabelInterpreter', 'none');
    xlabel(ax2, 'Niv learning-rate asymmetry');
    ylabel(ax2, sprintf('Participants, n = %d', height(nivTable)));
    title(ax2, 'Participant-level Niv asymmetry', ...
        'FontWeight', 'bold');
    box(ax2, 'off');

    title(layout, 'Fitted learning rates and Niv asymmetry', ...
        'FontWeight', 'bold', 'FontSize', 17);
end
