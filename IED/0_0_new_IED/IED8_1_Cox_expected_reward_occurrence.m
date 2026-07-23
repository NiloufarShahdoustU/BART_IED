% Cox models for IED occurrence x expected reward
% Permutation test shuffles whole trials, not interval rows
% Author: Nill

clear;
clc;
close all;


% these shouuld be changed if you're on your own computer not veronica's

inputFolder = ...
    'D:\Nill\data\IED1_find_number_of_IEDs\';

modelFolder = ...
    'D:\Nill\data\param_recovery_1_modeling\';

outputFolder = ...
    'D:\Nill\code\BART\IED\IED8_1_Cox_expected_reward_occurrence\';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

maximumRT = 20;
defaultSamplingFrequency = 1000;
numberOfPermutations = 500;
randomSeed = 42;
iedMode = 'occurrence';

colorIT = [0.847 0.333 0.153];
colorRT = [0.204 0.459 0.702];
colorBR = [0.250 0.600 0.250];

configs(1).code = 'IT';
configs(1).durationField = 'ITs';
configs(1).iedField = 'IED_occurance_IT';
configs(1).bankedOnly = false;
configs(1).windowMilliseconds = 1000;
configs(1).color = colorIT;
configs(1).hazardName = 'action';

configs(2).code = 'RT';
configs(2).durationField = 'RTs';
configs(2).iedField = 'IED_occurance_RT';
configs(2).bankedOnly = false;
configs(2).windowMilliseconds = 500;
configs(2).color = colorRT;
configs(2).hazardName = 'response';

configs(3).code = 'BR';
configs(3).durationField = 'ITs';
configs(3).iedField = 'IED_occurance_IT';
configs(3).bankedOnly = true;
configs(3).windowMilliseconds = 1000;
configs(3).color = colorBR;
configs(3).hazardName = 'banking';

textOutputFile = fullfile(outputFolder, ...
    'IT_RT_BR_expected_reward_overall_console_output.txt');

diary off;
if exist(textOutputFile, 'file')
    delete(textOutputFile);
end
diary(textOutputFile);

%% Run IT, RT, and BR

results = cell(3, 1);

for a = 1:3
    fprintf('\n\nRUNNING %s EXPECTED-REWARD ANALYSIS\n', configs(a).code);

    results{a} = runExpectedRewardCox(inputFolder, modelFolder, ...
        configs(a), maximumRT, defaultSamplingFrequency, ...
        numberOfPermutations, randomSeed + a - 1, iedMode);
end

analysisResults = [results{:}];
analysisConfigs = configs;

%% Save results

for a = 1:3
    code = configs(a).code;
    r = results{a};

    writetable(r.countingProcessData, fullfile(outputFolder, ...
        sprintf('%s_expected_reward_counting_process_data.csv', code)));

    writetable(r.trialSummary, fullfile(outputFolder, ...
        sprintf('%s_expected_reward_trial_summary.csv', code)));

    writetable(r.coxResults, fullfile(outputFolder, ...
        sprintf('%s_expected_reward_cox_results.csv', code)));

    writetable(r.modelEffectTable, fullfile(outputFolder, ...
        sprintf('%s_expected_reward_adjusted_effect.csv', code)));

    writetable(r.permutationResults, fullfile(outputFolder, ...
        sprintf('%s_expected_reward_permutation_results.csv', code)));

    writematrix(r.baselineCumulativeHazard, fullfile(outputFolder, ...
        sprintf('%s_expected_reward_baseline_cumulative_hazard.csv', code)));
end

combinedModelOutputFile = fullfile(outputFolder, ...
    'IT_RT_BR_expected_reward_overall_results.mat');

save(combinedModelOutputFile, 'analysisResults', 'analysisConfigs', ...
    'maximumRT', 'defaultSamplingFrequency', 'numberOfPermutations', ...
    'randomSeed', 'iedMode', '-v7.3');

%% Make one 3 x 2 figure

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [70 30 1500 1600]);

layout = tiledlayout(fig, 3, 2, ...
    'Padding', 'compact', 'TileSpacing', 'compact');

sgtitle(layout, 'Participant-stratified Cox models: IED x expected reward', ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

for a = 1:3
    plotForest(nexttile(layout), results{a}, configs(a));
    plotExpectedReward(nexttile(layout), results{a}, configs(a), iedMode);
end

combinedPDFOutputFile = fullfile(outputFolder, ...
    'IT_RT_BR_expected_reward_overall_cox_visualization_3x2.pdf');

exportgraphics(fig, combinedPDFOutputFile, 'ContentType', 'vector');
close(fig);

fprintf('\ndoneeeeeee\n');
fprintf('Figure: %s\n', combinedPDFOutputFile);
fprintf('Results: %s\n', combinedModelOutputFile);
fprintf('Console output: %s\n', textOutputFile);
diary off;


%% Run one outcome

function result = runExpectedRewardCox(inputFolder, modelFolder, config, ...
    maximumRT, defaultSamplingFrequency, numberOfPermutations, ...
    randomSeed, iedMode)

    code = config.code;
    windowMilliseconds = config.windowMilliseconds;
    windowSeconds = windowMilliseconds / 1000;
    files = dir(fullfile(inputFolder, '*.LFPIED.mat'));
    modelFiles = dir(fullfile(modelFolder, '*TDdataParamRecovery.mat'));

    countingProcessData = table();
    trialSummary = table();
    permutationTrialData = table();

    for pt = 1:length(files)
        fileName = files(pt).name;
        parts = strsplit(fileName, '.');
        patientID = string(parts{1});

        fprintf('%s: processing patient %s\n', code, patientID);

        modelFile = findModelFile(modelFiles, patientID);
        if strlength(modelFile) == 0
            fprintf('%s: no expected-reward file, skipped\n', patientID);
            continue;
        end

        data = load(fullfile(inputFolder, fileName));
        LFPIED = data.LFPIED;

        modelData = load(fullfile(modelFolder, modelFile));
        TD = modelData.TDdataParamRecovery;
        expectedRewards = squeeze(TD.expectedReward( ...
            round(TD.bestApIdx), round(TD.bestAnIdx), :));
        expectedRewards = double(expectedRewards(:));

        RTs = double(LFPIED.RTs(:));
        durations = LFPIED.(config.durationField);
        durations = double(durations(:));
        isControl = double(LFPIED.isControl(:));

        balloonColor = round(double(LFPIED.balloonType(:)));
        balloonColor = mod(balloonColor - 1, 10) + 1;

        IEDoccurrence = LFPIED.(config.iedField);
        nTrials = min(length(durations), length(expectedRewards));

        if config.bankedOnly
            bankedTrials = double(LFPIED.BankedTrials(:));
        else
            bankedTrials = NaN(length(durations), 1);
        end

        samplingFrequency = defaultSamplingFrequency;
        if isfield(LFPIED, 'Fs')
            samplingFrequency = double(LFPIED.Fs);
        end

        for trial = 1:nTrials
            duration = durations(trial);

            keepTrial = isControl(trial) == 0 && ...
                RTs(trial) > 0 && RTs(trial) <= maximumRT && ...
                duration > 0 && isfinite(expectedRewards(trial)) && ...
                ismember(balloonColor(trial), [1 2 3]);

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
                repmat(expectedRewards(trial), n, 1), ...
                repmat(samplingFrequency, n, 1), ...
                'VariableNames', {'patientID', 'trialNumber', ...
                'tStart', 'tStop', 'censored', 'eventAtStop', ...
                'postIED', 'activeIEDWindowCount', ...
                'timeSinceMostRecentIED_seconds', 'balloonColorCode', ...
                'expectedReward', 'samplingFrequencyHz'});

            countingProcessData = [countingProcessData; trialRows];

            nIEDs = length(IEDtimes);
            exposure = sum((tStop - tStart) .* postIED);
            eventInWindow = double(eventObserved && postIED(end) == 1);

            newSummary = table(patientID, trial, duration, ...
                expectedRewards(trial), bankedTrials(trial), ...
                double(eventObserved), eventInWindow, nIEDs, exposure, ...
                balloonColor(trial), samplingFrequency, ...
                'VariableNames', {'patientID', 'trialNumber', ...
                'durationSeconds', 'expectedReward', 'bankedTrial', ...
                'eventObserved', 'eventInPostIEDWindow', ...
                'nIEDsInFullDuration', 'postIEDExposureSeconds', ...
                'balloonColorCode', 'samplingFrequencyHz'});

            trialSummary = [trialSummary; newSummary];

            newPermutationTrial = table(patientID, trial, duration, ...
                double(eventObserved), balloonColor(trial), ...
                expectedRewards(trial), {IEDtimes}, ...
                'VariableNames', {'patientID', 'trialNumber', ...
                'durationSeconds', 'eventObserved', 'balloonColorCode', ...
                'expectedReward', 'IEDtimes_seconds'});

            permutationTrialData = [permutationTrialData; newPermutationTrial];
        end
    end

    %% Fit the real model

    [countingProcessData.patientStratum, patientLevels] = ...
        findgroups(countingProcessData.patientID);

    permutationTrialData.patientStratum = ...
        findgroups(permutationTrialData.patientID);

    [X, predictorNames] = makePredictors(countingProcessData, iedMode);
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

    interactionName = getInteractionName(iedMode);
    iedIndex = 1;
    interactionIndex = find(predictorNames == interactionName, 1);

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

            % Duration, event, and expected reward belong to one trial.
            shuffledTrials.durationSeconds(groupRows) = ...
                permutationTrialData.durationSeconds(sourceRows);
            shuffledTrials.eventObserved(groupRows) = ...
                permutationTrialData.eventObserved(sourceRows);
            shuffledTrials.expectedReward(groupRows) = ...
                permutationTrialData.expectedReward(sourceRows);
        end

        [permX, permT, permCensoring, permStrata] = ...
            buildPermutationData(shuffledTrials, windowSeconds, iedMode);

        try
            b = coxphfit(permX, permT, 'Censoring', permCensoring, ...
                'Strata', permStrata, 'Ties', 'efron', ...
                'Baseline', 0, 'Options', permutationOptions);

            permutationBeta(permutation, :) = ...
                [b(iedIndex), b(interactionIndex)];
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

    permutationP(iedIndex) = (1 + sum(abs(validPermutationBeta(:, 1)) >= ...
        abs(beta(iedIndex)))) / (size(validPermutationBeta, 1) + 1);

    permutationP(interactionIndex) = ...
        (1 + sum(abs(validPermutationBeta(:, 2)) >= ...
        abs(beta(interactionIndex)))) / ...
        (size(validPermutationBeta, 1) + 1);

    permutationResults = table((1:numberOfPermutations)', ...
        permutationBeta(:, 1), permutationBeta(:, 2), ...
        permutationSucceeded, ...
        'VariableNames', {'permutationIndex', 'permutedIEDBeta', ...
        'permutedInteractionBeta', 'fitSucceeded'});

    %% Results

    betaLow = beta - 1.96 * robustSE;
    betaHigh = beta + 1.96 * robustSE;
    hazardRatio = exp(beta);

    coxResults = table(predictorNames, beta(:), stats.se(:), ...
        stats.z(:), stats.p(:), robustSE(:), permutationP, ...
        betaLow(:), betaHigh(:), hazardRatio(:), ...
        exp(betaLow(:)), exp(betaHigh(:)), ...
        'VariableNames', {'predictor', 'beta_logHazard', ...
        'modelBasedSE', 'modelBasedZ', 'modelBasedPValue', ...
        'clusterRobustSE', 'permutationPValue', ...
        'betaCILow', 'betaCIHigh', 'hazardRatio', ...
        'hazardRatioCILow', 'hazardRatioCIHigh'});

    modelEffectTable = makeExpectedRewardCurve(beta, robustCovariance, ...
        predictorNames, trialSummary.expectedReward, iedMode);

    expectedRewardIndex = find(predictorNames == "expected_reward", 1);
    slopeOutside = beta(expectedRewardIndex);
    slopeInside = slopeOutside + beta(interactionIndex);

    rowDuration = countingProcessData.tStop - countingProcessData.tStart;
    nEventsOutside = sum(countingProcessData.eventAtStop & ...
        countingProcessData.postIED == 0);
    nEventsInside = sum(countingProcessData.eventAtStop & ...
        countingProcessData.postIED == 1);
    timeOutside = sum(rowDuration(countingProcessData.postIED == 0));
    timeInside = sum(rowDuration(countingProcessData.postIED == 1));

    fprintf('\n%s expected-reward Cox model\n', code);
    fprintf('participants: %d\n', length(patientLevels));
    fprintf('trials: %d\n', height(trialSummary));
    fprintf('events outside/inside window: %d / %d\n', ...
        nEventsOutside, nEventsInside);
    fprintf('time outside/inside window: %.3f / %.3f s\n', ...
        timeOutside, timeInside);
    fprintf('log likelihood: %.6f\n', logLikelihood);
    disp(coxResults);
    fprintf('IED x expected reward permutation p = %.6g\n', ...
        permutationP(interactionIndex));

    result.code = code;
    result.countingProcessData = countingProcessData;
    result.trialSummary = trialSummary;
    result.coxResults = coxResults;
    result.modelEffectTable = modelEffectTable;
    result.beta = beta;
    result.logLikelihood = logLikelihood;
    result.baselineCumulativeHazard = baselineCumulativeHazard;
    result.stats = stats;
    result.clusterRobustCovariance = robustCovariance;
    result.clusterRobustSE = robustSE;
    result.predictorNames = predictorNames;
    result.patientLevels = patientLevels;
    result.expectedRewardSlopeOutsideIED = slopeOutside;
    result.expectedRewardSlopeInsideIED = slopeInside;
    result.permutationResults = permutationResults;
    result.permutationPValue = permutationP;
    result.numberOfPermutations = numberOfPermutations;
    result.numberOfValidPermutations = size(validPermutationBeta, 1);
end


%% Find one participant's modeling file

function modelFile = findModelFile(modelFiles, patientID)

    names = string({modelFiles.name});
    matches = contains(lower(names), lower(patientID));

    if ~any(matches)
        modelFile = "";
    else
        matchingNames = names(matches);
        [~, shortest] = min(strlength(matchingNames));
        modelFile = matchingNames(shortest);
    end
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

function [X, predictorNames] = makePredictors(data, iedMode)

    if strcmp(iedMode, 'count')
        ied = double(data.activeIEDWindowCount);
        iedName = "active_ied_count";
        interactionName = "active_ied_count_x_expected_reward";
    else
        ied = double(data.postIED);
        iedName = "post_ied_indicator";
        interactionName = "post_ied_x_expected_reward";
    end

    expectedReward = double(data.expectedReward);

    X = [ied, expectedReward, ...
        double(data.balloonColorCode == 2), ...
        double(data.balloonColorCode == 3), ...
        ied .* expectedReward];

    predictorNames = [iedName; "expected_reward"; ...
        "orange_vs_yellow"; "red_vs_yellow"; interactionName];
end


function interactionName = getInteractionName(iedMode)
    if strcmp(iedMode, 'count')
        interactionName = "active_ied_count_x_expected_reward";
    else
        interactionName = "post_ied_x_expected_reward";
    end
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
        reward = trialData.expectedReward(trial);

        if strcmp(iedMode, 'count')
            ied = activeIEDCount(:);
        else
            ied = postIED(:);
        end

        Xcells{trial} = [ied, repmat(reward, n, 1), ...
            repmat(double(color == 2), n, 1), ...
            repmat(double(color == 3), n, 1), ied .* reward];

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


%% Expected-reward effect curve

function effectTable = makeExpectedRewardCurve(beta, covariance, ...
    predictorNames, observedReward, iedMode)

    rewardGrid = linspace(min(observedReward), max(observedReward), 120)';
    referenceReward = mean(observedReward);

    iedIndex = 1;
    rewardIndex = find(predictorNames == "expected_reward", 1);
    interactionIndex = find(predictorNames == getInteractionName(iedMode), 1);

    n = length(rewardGrid);
    logNoIED = zeros(n, 1);
    logOneIED = zeros(n, 1);
    seNoIED = zeros(n, 1);
    seOneIED = zeros(n, 1);

    for k = 1:n
        noIEDContrast = zeros(length(beta), 1);
        noIEDContrast(rewardIndex) = rewardGrid(k) - referenceReward;

        oneIEDContrast = noIEDContrast;
        oneIEDContrast(iedIndex) = 1;
        oneIEDContrast(interactionIndex) = rewardGrid(k);

        logNoIED(k) = noIEDContrast' * beta;
        logOneIED(k) = oneIEDContrast' * beta;
        seNoIED(k) = sqrt(max(noIEDContrast' * covariance * noIEDContrast, 0));
        seOneIED(k) = sqrt(max(oneIEDContrast' * covariance * oneIEDContrast, 0));
    end

    effectTable = table(rewardGrid, exp(logNoIED), ...
        exp(logNoIED - 1.96 * seNoIED), ...
        exp(logNoIED + 1.96 * seNoIED), exp(logOneIED), ...
        exp(logOneIED - 1.96 * seOneIED), ...
        exp(logOneIED + 1.96 * seOneIED), ...
        repmat(referenceReward, n, 1), ...
        'VariableNames', {'expectedReward', 'relativeHazardNoIED', ...
        'relativeHazardNoIEDCILow', 'relativeHazardNoIEDCIHigh', ...
        'relativeHazardPostIED', 'relativeHazardPostIEDCILow', ...
        'relativeHazardPostIEDCIHigh', 'referenceExpectedReward'});
end


%% Forest plot

function plotForest(ax, result, config)

    data = result.coxResults;
    y = (1:height(data))';

    hold(ax, 'on');
    errorbar(ax, data.beta_logHazard, y, ...
        data.beta_logHazard - data.betaCILow, ...
        data.betaCIHigh - data.beta_logHazard, ...
        'horizontal', 'o', 'LineWidth', 1.3, 'MarkerSize', 6, ...
        'MarkerFaceColor', config.color, 'MarkerEdgeColor', 'none', ...
        'Color', config.color);

    xline(ax, 0, '--', 'Color', [0.35 0.35 0.35]);

    labels = replace(data.predictor, "_", " ");
    set(ax, 'YTick', y, 'YTickLabel', labels, 'YDir', 'reverse', ...
        'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', ...
        'TickLabelInterpreter', 'none');

    xlabel(ax, 'log-hazard coefficient with robust 95% CI');
    title(ax, sprintf('%s: mechanistic Cox estimates', config.code), ...
        'FontWeight', 'bold');

    limits = [data.betaCILow; data.betaCIHigh; 0];
    padding = 0.15 * (max(limits) - min(limits));
    xlim(ax, [min(limits) - padding, max(limits) + padding]);
    ylim(ax, [0.5 height(data) + 0.5]);

    pValues = data.modelBasedPValue;
    usePermutation = ~isnan(data.permutationPValue);
    pValues(usePermutation) = data.permutationPValue(usePermutation);

    for row = 1:height(data)
        if pValues(row) < 0.05
            text(ax, data.betaCIHigh(row) + 0.03 * diff(xlim(ax)), ...
                row, '*', 'FontSize', 15, 'FontWeight', 'bold');
        end
    end

    box(ax, 'off');
    hold(ax, 'off');
end


%% Expected-reward plot

function plotExpectedReward(ax, result, config, iedMode)

    curve = result.modelEffectTable;
    x = curve.expectedReward;

    hold(ax, 'on');

    fill(ax, [x; flipud(x)], ...
        [curve.relativeHazardNoIEDCILow; ...
        flipud(curve.relativeHazardNoIEDCIHigh)], ...
        [0.78 0.78 0.78], 'EdgeColor', 'none', 'FaceAlpha', 0.30, ...
        'HandleVisibility', 'off');

    lightColor = 0.78 * [1 1 1] + 0.22 * config.color;
    fill(ax, [x; flipud(x)], ...
        [curve.relativeHazardPostIEDCILow; ...
        flipud(curve.relativeHazardPostIEDCIHigh)], ...
        lightColor, 'EdgeColor', 'none', 'FaceAlpha', 0.45, ...
        'HandleVisibility', 'off');

    yline(ax, 1, '--', 'Color', [0.45 0.45 0.45], ...
        'HandleVisibility', 'off');
    xline(ax, curve.referenceExpectedReward(1), ':', ...
        'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');

    if strcmp(iedMode, 'count')
        label0 = '0 active IEDs';
        label1 = '1 active IED';
    else
        label0 = 'outside post-IED window';
        label1 = 'inside post-IED window';
    end

    line0 = plot(ax, x, curve.relativeHazardNoIED, '--', ...
        'LineWidth', 2, 'Color', [0.25 0.25 0.25], ...
        'DisplayName', label0);
    line1 = plot(ax, x, curve.relativeHazardPostIED, '-', ...
        'LineWidth', 2.3, 'Color', config.color, 'DisplayName', label1);

    set(ax, 'YScale', 'log', 'FontName', 'Arial', ...
        'FontSize', 10, 'TickDir', 'out');
    xlabel(ax, 'expected reward');
    ylabel(ax, sprintf('relative %s hazard', config.hazardName));
    title(ax, sprintf('%s: adjusted expected-reward effect', config.code), ...
        'FontWeight', 'bold');

    interactionName = getInteractionName(iedMode);
    row = result.coxResults(result.coxResults.predictor == interactionName, :);

    resultText = sprintf([ ...
        'IED x reward beta = %.3f\n' ...
        '95%% CI [%.3f, %.3f]\n' ...
        'permutation p = %.3g\n' ...
        'slope at IED 0 = %.3f\n' ...
        'slope at IED 1 = %.3f'], ...
        row.beta_logHazard, row.betaCILow, row.betaCIHigh, ...
        row.permutationPValue, result.expectedRewardSlopeOutsideIED, ...
        result.expectedRewardSlopeInsideIED);

    text(ax, 0.97, 0.96, resultText, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
        'FontName', 'Arial', 'FontSize', 9, 'BackgroundColor', 'w', ...
        'EdgeColor', [0.8 0.8 0.8], 'Margin', 5);

    legend(ax, [line0 line1], 'Location', 'best', ...
        'Box', 'off', 'FontSize', 8);
    box(ax, 'off');
    hold(ax, 'off');
end
