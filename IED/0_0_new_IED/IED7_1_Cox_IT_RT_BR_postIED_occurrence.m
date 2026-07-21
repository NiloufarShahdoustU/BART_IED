% Participant-stratified Cox models for IT, RT, and BR
% Predictor: IED occurrence (0 = no IED, 1 = one or more IEDs)
% Author: Nill

clear;
clc;
close all;

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED7_1_Cox_IT_RT_BR_postIED_occurrence\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

textOutputFile = fullfile(outputFolderName, ...
    'IT_RT_BR_postIED_RT500ms_IT_BR1000ms_console_output.txt');

diary off;
if exist(textOutputFile, 'file')
    delete(textOutputFile);
end
diary(textOutputFile);

samplingFrequency = 1000;
maximumRT = 20;
numberOfPermutations = 5;
randomSeed = 42;

colorIT = [0.847 0.333 0.153];
colorRT = [0.204 0.459 0.702];
colorBR = [0.250 0.600 0.250];

%% Analysis settings

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

%% Run IT, RT, and BR

results = cell(3, 1);

for a = 1:3
    fprintf('\n\nRUNNING %s ANALYSIS\n', configs(a).code);

    results{a} = runCoxAnalysis( ...
        inputFolderName_LFPIED, configs(a), samplingFrequency, ...
        maximumRT, numberOfPermutations, randomSeed + a - 1);
end

%% Save results

for a = 1:3
    code = configs(a).code;
    window = configs(a).windowMilliseconds;
    r = results{a};

    writetable(r.countingProcessData, fullfile(outputFolderName, ...
        sprintf('%s_postIED%dms_counting_process_data.csv', code, window)));

    writetable(r.trialSummary, fullfile(outputFolderName, ...
        sprintf('%s_postIED%dms_trial_summary.csv', code, window)));

    writetable(r.coxResults, fullfile(outputFolderName, ...
        sprintf('%s_postIED%dms_cox_results.csv', code, window)));

    writematrix(r.baselineCumulativeHazard, fullfile(outputFolderName, ...
        sprintf('%s_postIED%dms_baseline_cumulative_hazard.csv', code, window)));

    writetable(r.patientEventRates, fullfile(outputFolderName, ...
        sprintf('%s_postIED%dms_patient_event_rates.csv', code, window)));

    writetable(r.modelEffectTable, fullfile(outputFolderName, ...
        sprintf('%s_postIED%dms_model_effect.csv', code, window)));

    writetable(r.permutationResults, fullfile(outputFolderName, ...
        sprintf('%s_postIED%dms_permutation_results.csv', code, window)));
end

combinedModelOutputFile = fullfile(outputFolderName, ...
    'IT_RT_BR_postIED_RT500ms_IT_BR1000ms_combined_results.mat');

analysisResults = [results{:}];
analysisConfigs = configs;
defaultSamplingFrequencyHz = samplingFrequency;
maximumRTSeconds = maximumRT;
permutationRandomSeed = randomSeed;

save(combinedModelOutputFile, 'analysisResults', 'analysisConfigs', ...
    'defaultSamplingFrequencyHz', 'maximumRTSeconds', ...
    'numberOfPermutations', 'permutationRandomSeed', '-v7.3');

%% Plot results

combinedFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [80 40 1450 1550]);

layout = tiledlayout(combinedFigure, 3, 2, ...
    'Padding', 'compact', 'TileSpacing', 'compact');

sgtitle(layout, ...
    'Participant-stratified Cox models: RT 500-ms; IT and BR 1000-ms post-IED windows', ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

for a = 1:3
    plotForestPanel(nexttile(layout), results{a}.coxResults, configs(a));
    plotEffectPanel(nexttile(layout), results{a}.coxResults, ...
        results{a}.modelEffectTable, configs(a));
end

combinedPDFOutputFile = fullfile(outputFolderName, ...
    'IT_RT_BR_postIED_RT500ms_IT_BR1000ms_cox_visualization_3x2.pdf');

exportgraphics(combinedFigure, combinedPDFOutputFile, 'ContentType', 'vector');
close(combinedFigure);

fprintf('\ndoneeeeeee\n');
fprintf('One visualization PDF saved:\n%s\n', combinedPDFOutputFile);
fprintf('Combined MAT file saved:\n%s\n', combinedModelOutputFile);
fprintf('Console output saved:\n%s\n', textOutputFile);
diary off;


%% Cox analysis

function result = runCoxAnalysis(inputFolder, config, samplingFrequency, ...
    maximumRT, numberOfPermutations, randomSeed)

    code = config.code;
    windowMilliseconds = config.windowMilliseconds;
    windowSeconds = windowMilliseconds / 1000;
    fileList = dir(fullfile(inputFolder, '*.LFPIED.mat'));

    countingProcessData = table();
    trialSummary = table();
    permutationTrialData = table();

    for pt = 1:length(fileList)
        fileName = fileList(pt).name;
        parts = strsplit(fileName, '.');
        patientID = string(parts{1});
        fprintf('%s: processing patient %s\n', code, char(patientID));

        loadedData = load(fullfile(inputFolder, fileName));
        LFPIED = loadedData.LFPIED;

        RTs = LFPIED.RTs(:);
        durations = LFPIED.(config.durationField);
        durations = durations(:);
        isControl = LFPIED.isControl(:);
        balloonColor = double(LFPIED.balloonType(:));
        IEDoccurrence = LFPIED.(config.iedField);
        nTrials = length(durations);

        if config.bankedOnly
            bankedTrials = LFPIED.BankedTrials(:);
        else
            bankedTrials = NaN(nTrials, 1);
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
            IEDtimes = IEDtimes(IEDtimes > 0 & IEDtimes < duration);

            windowEnds = min(IEDtimes + windowSeconds, duration);
            breakTimes = unique([0; IEDtimes; windowEnds; duration]);
            tStart = breakTimes(1:end-1);
            tStop = breakTimes(2:end);
            midpoint = (tStart + tStop) / 2;

            postIED = zeros(length(midpoint), 1);
            activeIEDWindowCount = zeros(length(midpoint), 1);

            for ied = 1:length(IEDtimes)
                insideWindow = midpoint >= IEDtimes(ied) & ...
                    midpoint <= IEDtimes(ied) + windowSeconds;
                postIED(insideWindow) = 1;
                activeIEDWindowCount(insideWindow) = ...
                    activeIEDWindowCount(insideWindow) + 1;
            end

            timeSinceIED = NaN(length(tStart), 1);
            for k = 1:length(tStart)
                previousIEDs = IEDtimes(IEDtimes < midpoint(k));
                if ~isempty(previousIEDs)
                    timeSinceIED(k) = midpoint(k) - previousIEDs(end);
                end
            end

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
                activeIEDWindowCount, timeSinceIED, ...
                repmat(balloonColor(trial), n, 1), ...
                repmat(samplingFrequency, n, 1), ...
                'VariableNames', {'patientID', 'trialNumber', 'tStart', ...
                'tStop', 'censored', 'eventAtStop', 'postIED', ...
                'activeIEDWindowCount', 'timeSinceMostRecentIED_seconds', ...
                'balloonColorCode', 'samplingFrequencyHz'});

            countingProcessData = [countingProcessData; trialRows];

            nIEDs = sum(round(IEDoccurrence(:, 1)) == trial);
            postIEDExposure = sum((tStop - tStart) .* postIED);
            eventInPostIEDWindow = double(eventObserved && postIED(end) == 1);

            newSummary = table(patientID, trial, duration, ...
                bankedTrials(trial), double(eventObserved), ...
                eventInPostIEDWindow, nIEDs, postIEDExposure, ...
                balloonColor(trial), samplingFrequency, ...
                'VariableNames', {'patientID', 'trialNumber', ...
                'durationSeconds', 'bankedTrial', 'eventObserved', ...
                'eventInPostIEDWindow', 'nIEDsInFullDuration', ...
                'postIEDExposureSeconds', 'balloonColorCode', ...
                'samplingFrequencyHz'});

            trialSummary = [trialSummary; newSummary];

            newPermutationTrial = table(patientID, trial, duration, ...
                double(eventObserved), balloonColor(trial), {IEDtimes}, ...
                'VariableNames', {'patientID', 'trialNumber', ...
                'durationSeconds', 'eventObserved', 'balloonColorCode', ...
                'IEDtimes_seconds'});

            permutationTrialData = [permutationTrialData; newPermutationTrial];
        end
    end

    %% Model predictors

    [countingProcessData.patientStratum, patientLevels] = ...
        findgroups(countingProcessData.patientID);

    permutationTrialData.patientStratum = ...
        findgroups(permutationTrialData.patientID);

    X = [countingProcessData.postIED, ...
        double(countingProcessData.balloonColorCode == 2), ...
        double(countingProcessData.balloonColorCode == 3)];

    predictorNames = [ ...
        "any_ied_occurrence_within_" + string(windowMilliseconds) + "ms"; ...
        "orange_vs_yellow"; ...
        "red_vs_yellow"];

    T = [countingProcessData.tStart, countingProcessData.tStop];
    censoring = countingProcessData.censored;
    strata = countingProcessData.patientStratum;

    nEventsOutsideWindow = sum(countingProcessData.eventAtStop == 1 & ...
        countingProcessData.postIED == 0);
    nEventsInsideWindow = sum(countingProcessData.eventAtStop == 1 & ...
        countingProcessData.postIED == 1);

    rowDuration = countingProcessData.tStop - countingProcessData.tStart;
    timeOutsideWindow = sum(rowDuration(countingProcessData.postIED == 0));
    timeInsideWindow = sum(rowDuration(countingProcessData.postIED == 1));

    %% Fit participant-stratified Cox model

    options = statset('coxphfit');
    options.Display = 'final';
    options.MaxIter = 1000;
    options.MaxFunEvals = 5000;

    [beta, logLikelihood, baselineCumulativeHazard, stats] = coxphfit( ...
        X, T, 'Censoring', censoring, 'Strata', strata, ...
        'Ties', 'efron', 'Baseline', 0, 'Options', options);

    %% Participant-clustered robust p values and confidence intervals

    patientNumbers = unique(strata);
    clusterScores = zeros(length(patientNumbers), length(beta));

    for p = 1:length(patientNumbers)
        clusterScores(p, :) = sum(stats.scores(strata == patientNumbers(p), :), 1);
    end

    robustCovariance = stats.covb * (clusterScores' * clusterScores) * ...
        stats.covb * length(patientNumbers) / (length(patientNumbers) - 1);
    robustSE = sqrt(diag(robustCovariance));

    %% Permutation test

    observedIEDBeta = beta(1);
    permutationBeta = NaN(numberOfPermutations, 1);
    permutationSucceeded = false(numberOfPermutations, 1);

    permutationGroups = findgroups(permutationTrialData.patientID, ...
        permutationTrialData.balloonColorCode);

    rng(randomSeed);
    permutationOptions = options;
    permutationOptions.Display = 'off';

    fprintf('%s: running %d permutations...\n', code, numberOfPermutations);

    for permutation = 1:numberOfPermutations
        shuffledTrials = permutationTrialData;

        for group = 1:max(permutationGroups)
            groupRows = find(permutationGroups == group);
            sourceRows = groupRows(randperm(length(groupRows)));
            shuffledTrials.durationSeconds(groupRows) = ...
                permutationTrialData.durationSeconds(sourceRows);
            shuffledTrials.eventObserved(groupRows) = ...
                permutationTrialData.eventObserved(sourceRows);
        end

        [permutedX, permutedT, permutedCensoring, permutedStrata] = ...
            buildShuffledCoxData(shuffledTrials, windowSeconds);

        try
            b = coxphfit(permutedX, permutedT, ...
                'Censoring', permutedCensoring, ...
                'Strata', permutedStrata, 'Ties', 'efron', ...
                'Baseline', 0, 'Options', permutationOptions);
            permutationBeta(permutation) = b(1);
            permutationSucceeded(permutation) = true;
        catch
            fprintf('%s permutation %d failed\n', code, permutation);
        end

        if mod(permutation, 100) == 0 || permutation == numberOfPermutations
            fprintf('%s permutation progress: %d/%d\n', ...
                code, permutation, numberOfPermutations);
        end
    end

    validPermutationBeta = permutationBeta(permutationSucceeded);
    permutationP = (1 + sum(abs(validPermutationBeta) >= ...
        abs(observedIEDBeta))) / (length(validPermutationBeta) + 1);

    permutationResults = table((1:numberOfPermutations)', ...
        permutationBeta, exp(permutationBeta), permutationSucceeded, ...
        'VariableNames', {'permutationIndex', 'permutedBeta_logHazard', ...
        'permutedHazardRatio', 'fitSucceeded'});

    %% Results table

    betaLow = beta - 1.96 * robustSE;
    betaHigh = beta + 1.96 * robustSE;
    hazardRatio = exp(beta);
    hazardRatioLow = exp(betaLow);
    hazardRatioHigh = exp(betaHigh);

    if strcmp(code, 'IT')
        primaryInterpretation = ...
            "HR > 1 means shorter IT; HR < 1 means longer IT.";
    elseif strcmp(code, 'RT')
        primaryInterpretation = ...
            "HR > 1 means shorter RT; HR < 1 means longer RT.";
    else
        primaryInterpretation = ...
            "HR > 1 means higher banking hazard; HR < 1 means lower banking hazard.";
    end

    interpretation = [primaryInterpretation; ...
        "Balloon-color effect relative to yellow."; ...
        "Balloon-color effect relative to yellow."];

    coxResults = table(predictorNames, beta(:), stats.se(:), ...
        stats.z(:), stats.p(:), robustSE(:), ...
        [permutationP; NaN; NaN], betaLow(:), betaHigh(:), ...
        hazardRatio(:), hazardRatioLow(:), hazardRatioHigh(:), ...
        interpretation, ...
        'VariableNames', {'predictor', 'beta_logHazard', 'modelBasedSE', ...
        'modelBasedZ', 'modelBasedPValue', 'clusterRobustSE', ...
        'permutationPValue', 'betaCILow', 'betaCIHigh', ...
        'hazardRatio', 'hazardRatioCILow', ...
        'hazardRatioCIHigh', 'interpretation'});

    fprintf('\n============================================================\n');
    fprintf('%s: %d-ms post-IED Cox model\n', code, windowMilliseconds);
    fprintf('============================================================\n');
    fprintf('participants: %d\n', length(patientLevels));
    fprintf('trials: %d\n', height(trialSummary));
    fprintf('counting-process rows: %d\n', height(countingProcessData));
    fprintf('events outside window: %d\n', nEventsOutsideWindow);
    fprintf('events inside window: %d\n', nEventsInsideWindow);
    fprintf('at-risk time outside window: %.3f s\n', timeOutsideWindow);
    fprintf('at-risk time inside window: %.3f s\n', timeInsideWindow);
    fprintf('log likelihood: %.6f\n\n', logLikelihood);
    disp(coxResults);

    fprintf('%s effect of any IED occurrence versus none:\n', code);
    fprintf('beta = %.6f\n', beta(1));
    fprintf('hazard ratio = %.6f\n', hazardRatio(1));
    fprintf('95%% CI = [%.6f, %.6f]\n', hazardRatioLow(1), hazardRatioHigh(1));
    fprintf('permutation p = %.6g (%d valid of %d requested)\n', ...
        permutationP, length(validPermutationBeta), numberOfPermutations);

    %% Participant event rates

    patientList = unique(countingProcessData.patientID);
    patientEventRates = table();

    for p = 1:length(patientList)
        rows = countingProcessData.patientID == patientList(p);
        patientData = countingProcessData(rows, :);
        duration = patientData.tStop - patientData.tStart;
        outside = patientData.postIED == 0;
        inside = patientData.postIED == 1;

        outsideTime = sum(duration(outside));
        insideTime = sum(duration(inside));
        outsideEvents = sum(patientData.eventAtStop(outside));
        insideEvents = sum(patientData.eventAtStop(inside));
        outsideRate = 100 * outsideEvents / outsideTime;
        insideRate = 100 * insideEvents / insideTime;
        crudeRateRatio = insideRate / outsideRate;

        newRate = table(patientList(p), outsideTime, insideTime, ...
            outsideEvents, insideEvents, outsideRate, insideRate, ...
            crudeRateRatio, ...
            'VariableNames', {'patientID', 'outsideWindowTime_seconds', ...
            'postIEDWindowTime_seconds', 'outsideWindowEvents', ...
            'postIEDWindowEvents', 'outsideWindowRate_per100s', ...
            'postIEDWindowRate_per100s', 'crudeRateRatio'});

        patientEventRates = [patientEventRates; newRate];
    end

    %% Adjusted IED effect

    exposureState = [0; 1];
    relativeHazard = hazardRatio(1) .^ exposureState;
    relativeHazardLow = hazardRatioLow(1) .^ exposureState;
    relativeHazardHigh = hazardRatioHigh(1) .^ exposureState;

    modelEffectTable = table(exposureState, relativeHazard, ...
        relativeHazardLow, relativeHazardHigh, ...
        repmat(windowMilliseconds, 2, 1), ...
        'VariableNames', {'activeIEDWindowCount', 'relativeHazard', ...
        'relativeHazardCILow', 'relativeHazardCIHigh', ...
        'postIEDWindowMilliseconds'});

    %% Return results

    result.code = code;
    result.countingProcessData = countingProcessData;
    result.trialSummary = trialSummary;
    result.coxResults = coxResults;
    result.beta = beta;
    result.logLikelihood = logLikelihood;
    result.baselineCumulativeHazard = baselineCumulativeHazard;
    result.stats = stats;
    result.clusterRobustCovariance = robustCovariance;
    result.clusterRobustSE = robustSE;
    result.predictorNames = predictorNames;
    result.patientLevels = patientLevels;
    result.nEventsOutsideWindow = nEventsOutsideWindow;
    result.nEventsInsideWindow = nEventsInsideWindow;
    result.timeOutsideWindowSeconds = timeOutsideWindow;
    result.timeInsideWindowSeconds = timeInsideWindow;
    result.patientEventRates = patientEventRates;
    result.observedIEDCounts = unique(countingProcessData.postIED);
    result.modelEffectTable = modelEffectTable;
    result.postIEDWindowMilliseconds = windowMilliseconds;
    result.permutationResults = permutationResults;
    result.permutationPValue = permutationP;
    result.numberOfPermutations = numberOfPermutations;
    result.numberOfValidPermutations = length(validPermutationBeta);
    result.permutationRandomSeed = randomSeed;
end


%% Rebuild Cox data after shuffling trial outcomes

function [X, T, censoring, strata] = buildShuffledCoxData(trialData, windowSeconds)

    nTrials = height(trialData);
    Xcells = cell(nTrials, 1);
    Tcells = cell(nTrials, 1);
    censoringCells = cell(nTrials, 1);
    strataCells = cell(nTrials, 1);

    for trial = 1:nTrials
        duration = trialData.durationSeconds(trial);
        IEDtimes = trialData.IEDtimes_seconds{trial};
        IEDtimes = IEDtimes(IEDtimes < duration);

        windowEnds = min(IEDtimes + windowSeconds, duration);
        breakTimes = unique([0; IEDtimes; windowEnds; duration]);
        tStart = breakTimes(1:end-1);
        tStop = breakTimes(2:end);
        midpoint = (tStart + tStop) / 2;

        postIED = zeros(length(midpoint), 1);

        for ied = 1:length(IEDtimes)
            insideWindow = midpoint >= IEDtimes(ied) & ...
                midpoint <= IEDtimes(ied) + windowSeconds;
            postIED(insideWindow) = 1;
        end

        color = trialData.balloonColorCode(trial);
        n = length(tStart);

        Xcells{trial} = [postIED, ...
            repmat(double(color == 2), n, 1), ...
            repmat(double(color == 3), n, 1)];

        Tcells{trial} = [tStart, tStop];

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


%% Forest plot

function plotForestPanel(ax, coxResults, config)

    hold(ax, 'on');

    y = (1:height(coxResults))';
    HR = coxResults.hazardRatio;
    low = coxResults.hazardRatioCILow;
    high = coxResults.hazardRatioCIHigh;

    errorbar(ax, HR, y, HR - low, high - HR, 'horizontal', 'o', ...
        'LineWidth', 1.2, 'MarkerSize', 6, ...
        'MarkerFaceColor', config.color, ...
        'MarkerEdgeColor', config.color, 'Color', config.color);

    xline(ax, 1, '--', 'Color', config.color);

    labels = [ ...
        "any IED occurrence within " + string(config.windowMilliseconds) + " ms"; ...
        "orange vs yellow"; ...
        "red vs yellow"];

    set(ax, 'YTick', y, 'YTickLabel', labels, 'YDir', 'reverse', ...
        'XScale', 'log', 'FontSize', 10, 'FontName', 'Arial', ...
        'TickDir', 'out');

    xlabel(ax, 'hazard ratio with robust 95% CI');
    title(ax, sprintf('%s: Cox model estimates', config.code), ...
        'FontWeight', 'bold');

    pValues = coxResults.modelBasedPValue;
    pValues(1) = coxResults.permutationPValue(1);
    significant = pValues < 0.05;

    xlim(ax, [0.8 * min([low; 1]), 1.25 * max([high; 1])]);
    ylim(ax, [0.5 height(coxResults) + 0.5]);

    for row = 1:height(coxResults)
        if significant(row)
            text(ax, 1.08 * high(row), row, '*', ...
                'FontSize', 16, 'FontWeight', 'bold');
        end
    end

    box(ax, 'off');
    hold(ax, 'off');
end


%% Adjusted effect plot

function plotEffectPanel(ax, coxResults, modelEffectTable, config)

    hold(ax, 'on');

    x = modelEffectTable.activeIEDWindowCount;
    HR = modelEffectTable.relativeHazard;
    low = modelEffectTable.relativeHazardCILow;
    high = modelEffectTable.relativeHazardCIHigh;
    lightColor = 0.75 * [1 1 1] + 0.25 * config.color;

    fill(ax, [x; flipud(x)], [low; flipud(high)], lightColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55);
    yline(ax, 1, '--', 'Color', [0.45 0.45 0.45]);
    plot(ax, x, HR, '-o', 'LineWidth', 2.2, 'Color', config.color, ...
        'MarkerFaceColor', config.color, 'MarkerEdgeColor', config.color);

    set(ax, 'YScale', 'log', 'FontSize', 10, ...
        'FontName', 'Arial', 'TickDir', 'out');
    xlim(ax, [-0.25 1.25]);
    xticks(ax, [0 1]);
    ylim(ax, [0.8 * min(low), 1.2 * max(high)]);

    xlabel(ax, 'IED occurrence in any channel (0 = no, 1 = yes)');
    ylabel(ax, sprintf('relative %s hazard', config.hazardName));
    title(ax, sprintf('%s: adjusted IED-occurrence effect', config.code), ...
        'FontWeight', 'bold');

    primary = coxResults(1, :);
    resultText = sprintf(['HR = %.3f\n95%% CI [%.3f, %.3f]\n' ...
        'permutation p = %.3g'], ...
        primary.hazardRatio, primary.hazardRatioCILow, ...
        primary.hazardRatioCIHigh, primary.permutationPValue);

    text(ax, 0.96, 0.95, resultText, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
        'FontName', 'Arial', 'FontSize', 9, 'BackgroundColor', 'w', ...
        'EdgeColor', [0.8 0.8 0.8], 'Margin', 5);

    box(ax, 'off');
    hold(ax, 'off');
end
