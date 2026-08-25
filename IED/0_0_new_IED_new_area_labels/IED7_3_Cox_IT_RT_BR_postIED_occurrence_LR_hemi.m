% Participant-stratified Cox models comparing left- and right-hemisphere IEDs
%
% The two hemisphere-specific post-IED indicators are entered together in
% the same Cox model. Therefore, the left estimate is adjusted for right
% post-IED windows and the right estimate is adjusted for left post-IED
% windows. A direct betaLeft - betaRight contrast formally compares them.
%
% Hemisphere is read from LFPIED.anatomicalLocs_wHemisphere. Labels are
% expected to end in _L or _R (for example, area_L and area_R).
%
% Author: Nill

clear;
clc;
close all;

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED_new_area_labels\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED7_3_Cox_IT_RT_BR_postIED_occurrence_LR_hemi\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

samplingFrequency = 1000;
maximumRT = 20;
numberOfPermutations = 1000;
randomSeed = 42;

colorIT = [0.847 0.333 0.153];
colorRT = [0.204 0.459 0.702];
colorBR = [0.250 0.600 0.250];
colorPR = [0.494 0.184 0.556];

leftColor = [0.121 0.466 0.705];
rightColor = [0.839 0.153 0.157];

configs(1).code = 'IT';
configs(1).durationField = 'ITs';
configs(1).iedField = 'IED_occurance_IT';
configs(1).bankedOnly = false;
configs(1).poppedOnly = false;
configs(1).windowMilliseconds = 1000;
configs(1).color = colorIT;
configs(1).hazardName = 'action';

configs(2).code = 'RT';
configs(2).durationField = 'RTs';
configs(2).iedField = 'IED_occurance_RT';
configs(2).bankedOnly = false;
configs(2).poppedOnly = false;
configs(2).windowMilliseconds = 500;
configs(2).color = colorRT;
configs(2).hazardName = 'response';

configs(3).code = 'BR';
configs(3).durationField = 'ITs';
configs(3).iedField = 'IED_occurance_IT';
configs(3).bankedOnly = true;
configs(3).poppedOnly = false;
configs(3).windowMilliseconds = 1000;
configs(3).color = colorBR;
configs(3).hazardName = 'banking';

configs(4).code = 'PR';
configs(4).durationField = 'ITs';
configs(4).iedField = 'IED_occurance_IT';
configs(4).bankedOnly = false;
configs(4).poppedOnly = true;
configs(4).windowMilliseconds = 1000;
configs(4).color = colorPR;
configs(4).hazardName = 'popping';

%% Run IT, RT, BR, and PR

results = cell(length(configs), 1);

for analysisIndex = 1:length(configs)
    fprintf('\n\nRUNNING %s HEMISPHERE ANALYSIS\n', ...
        configs(analysisIndex).code);

    results{analysisIndex} = runHemisphereCoxAnalysis( ...
        inputFolderName_LFPIED, configs(analysisIndex), ...
        samplingFrequency, maximumRT, numberOfPermutations, ...
        randomSeed + analysisIndex - 1);
end

%% Save compact results tables

hemisphereResultsAll = table();
contrastResultsAll = table();

for analysisIndex = 1:length(configs)
    thisHemisphereTable = results{analysisIndex}.hemisphereResults;
    thisHemisphereTable.outcome = repmat( ...
        string(configs(analysisIndex).code), ...
        height(thisHemisphereTable), 1);
    thisHemisphereTable = movevars( ...
        thisHemisphereTable, 'outcome', 'Before', 1);
    hemisphereResultsAll = [hemisphereResultsAll; thisHemisphereTable];

    thisContrastTable = results{analysisIndex}.contrastResults;
    thisContrastTable.outcome = string(configs(analysisIndex).code);
    thisContrastTable = movevars( ...
        thisContrastTable, 'outcome', 'Before', 1);
    contrastResultsAll = [contrastResultsAll; thisContrastTable];
end

writetable(hemisphereResultsAll, fullfile(outputFolderName, ...
    'IT_RT_BR_PR_left_right_hemisphere_cox_results.csv'));

writetable(contrastResultsAll, fullfile(outputFolderName, ...
    'IT_RT_BR_PR_left_minus_right_hemisphere_contrasts.csv'));

%% Left/right visualization

combinedFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [80 40 1500 1900]);

layout = tiledlayout(combinedFigure, 4, 2, ...
    'Padding', 'compact', 'TileSpacing', 'compact');

sgtitle(layout, [ ...
    'Left versus right hemisphere IED effects: ' ...
    'RT 500-ms; IT, BR, and PR 1000-ms post-IED windows'], ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

for analysisIndex = 1:length(configs)
    plotHemisphereForestPanel(nexttile(layout), ...
        results{analysisIndex}, configs(analysisIndex), ...
        leftColor, rightColor);

    plotHemisphereContrastPanel(nexttile(layout), ...
        results{analysisIndex}, configs(analysisIndex));
end

combinedPDFOutputFile = fullfile(outputFolderName, ...
    'IT_RT_BR_PR_left_right_hemisphere_comparison_4x2.pdf');

exportgraphics(combinedFigure, combinedPDFOutputFile, ...
    'ContentType', 'vector');
close(combinedFigure);

fprintf('\nDone.\n');
fprintf('Hemisphere results CSV:\n%s\n', fullfile(outputFolderName, ...
    'IT_RT_BR_PR_left_right_hemisphere_cox_results.csv'));
fprintf('Left-minus-right contrast CSV:\n%s\n', fullfile(outputFolderName, ...
    'IT_RT_BR_PR_left_minus_right_hemisphere_contrasts.csv'));
fprintf('Comparison PDF:\n%s\n', combinedPDFOutputFile);


%% Cox analysis for one outcome

function result = runHemisphereCoxAnalysis(inputFolder, config, ...
    defaultSamplingFrequency, maximumRT, numberOfPermutations, randomSeed)

    code = config.code;
    windowMilliseconds = config.windowMilliseconds;
    windowSeconds = windowMilliseconds / 1000;

    fileList = dir(fullfile(inputFolder, '*.LFPIED.mat'));

    if isempty(fileList)
        error('No .LFPIED.mat files were found in %s.', inputFolder);
    end

    countingProcessData = table();
    trialSummary = table();
    permutationTrialData = table();

    coverageLeftPatients = strings(0, 1);
    coverageRightPatients = strings(0, 1);
    IEDLeftPatients = strings(0, 1);
    IEDRightPatients = strings(0, 1);
    totalLeftIEDs = 0;
    totalRightIEDs = 0;

    for pt = 1:length(fileList)

        fileName = fileList(pt).name;
        fileParts = strsplit(fileName, '.');
        patientID = string(fileParts{1});
        fprintf('%s: processing patient %s\n', code, char(patientID));

        loadedData = load(fullfile(inputFolder, fileName));

        if ~isfield(loadedData, 'LFPIED')
            fprintf('Skipped %s: LFPIED structure was not found.\n', ...
                char(patientID));
            continue;
        end

        LFPIED = loadedData.LFPIED;

        requiredFields = { ...
            'RTs', config.durationField, config.iedField, ...
            'isControl', 'balloonType', 'selectedChans', ...
            'anatomicalLocs_wHemisphere'};

        if config.bankedOnly || config.poppedOnly
            requiredFields{end + 1} = 'BankedTrials';
        end

        missingFields = requiredFields( ...
            ~cellfun(@(fieldName) isfield(LFPIED, fieldName), ...
            requiredFields));

        if ~isempty(missingFields)
            fprintf('Skipped %s: missing %s.\n', char(patientID), ...
                strjoin(missingFields, ', '));
            continue;
        end

        selectedChans = round(double(LFPIED.selectedChans(:)));
        selectedChans = selectedChans( ...
            isfinite(selectedChans) & selectedChans >= 1);
        nSelectedChannels = length(selectedChans);

        if nSelectedChannels == 0
            fprintf('Skipped %s: no selected channels.\n', ...
                char(patientID));
            continue;
        end

        selectedLabels = getSelectedChannelLabels( ...
            LFPIED.anatomicalLocs_wHemisphere, ...
            selectedChans, nSelectedChannels);

        [leftChannelMask, rightChannelMask] = ...
            getHemisphereChannelMasks(selectedLabels);

        leftLocalChannels = find(leftChannelMask);
        rightLocalChannels = find(rightChannelMask);

        if ~isempty(leftLocalChannels)
            coverageLeftPatients(end + 1, 1) = patientID;
        end

        if ~isempty(rightLocalChannels)
            coverageRightPatients(end + 1, 1) = patientID;
        end

        if isempty(leftLocalChannels) && isempty(rightLocalChannels)
            fprintf(['Skipped %s: none of the selected labels ended ' ...
                'in _L or _R.\n'], char(patientID));
            continue;
        end

        RTs = double(LFPIED.RTs(:));
        durations = double(LFPIED.(config.durationField));
        durations = durations(:);
        isControl = double(LFPIED.isControl(:));
        balloonColor = double(LFPIED.balloonType(:));
        IEDoccurrence = double(LFPIED.(config.iedField));

        if isempty(IEDoccurrence)
            IEDoccurrence = zeros(0, 3);
        elseif size(IEDoccurrence, 2) < 3
            fprintf(['Skipped %s: %s must contain at least trial, ' ...
                'channel, and IED-time columns.\n'], ...
                char(patientID), config.iedField);
            continue;
        end

        vectorLengths = [length(RTs), length(durations), ...
            length(isControl), length(balloonColor)];

        if config.bankedOnly || config.poppedOnly
            bankedTrials = double(LFPIED.BankedTrials(:));
            vectorLengths(end + 1) = length(bankedTrials);
        else
            bankedTrials = NaN(length(durations), 1);
        end

        nTrials = min(vectorLengths);
        RTs = RTs(1:nTrials);
        durations = durations(1:nTrials);
        isControl = isControl(1:nTrials);
        balloonColor = balloonColor(1:nTrials);
        bankedTrials = bankedTrials(1:nTrials);

        patientHasLeftIED = false;
        patientHasRightIED = false;

        for trial = 1:nTrials

            duration = durations(trial);

            keepTrial = isControl(trial) == 0 && ...
                RTs(trial) > 0 && RTs(trial) <= maximumRT && ...
                isfinite(duration) && duration > 0 && ...
                ismember(balloonColor(trial), [1 2 3]);

            if config.bankedOnly || config.poppedOnly
                keepTrial = keepTrial && ...
                    ismember(bankedTrials(trial), [0 1]);

                if config.bankedOnly
                    eventObserved = bankedTrials(trial) == 1;
                else
                    eventObserved = bankedTrials(trial) == 0;
                end
            else
                eventObserved = true;
            end

            if ~keepTrial
                continue;
            end

            trialRows = isfinite(IEDoccurrence(:, 1)) & ...
                round(IEDoccurrence(:, 1)) == trial;

            leftRows = trialRows & ...
                isfinite(IEDoccurrence(:, 2)) & ...
                ismember(round(IEDoccurrence(:, 2)), leftLocalChannels);

            rightRows = trialRows & ...
                isfinite(IEDoccurrence(:, 2)) & ...
                ismember(round(IEDoccurrence(:, 2)), rightLocalChannels);

            leftIEDtimes = getValidIEDTimes( ...
                IEDoccurrence(leftRows, :), ...
                duration, defaultSamplingFrequency);

            rightIEDtimes = getValidIEDTimes( ...
                IEDoccurrence(rightRows, :), ...
                duration, defaultSamplingFrequency);

            totalLeftIEDs = totalLeftIEDs + length(leftIEDtimes);
            totalRightIEDs = totalRightIEDs + length(rightIEDtimes);
            patientHasLeftIED = patientHasLeftIED || ~isempty(leftIEDtimes);
            patientHasRightIED = patientHasRightIED || ~isempty(rightIEDtimes);

            [tStart, tStop, postIEDLeft, postIEDRight, ...
                activeLeftWindowCount, activeRightWindowCount] = ...
                buildHemisphereIntervals(duration, leftIEDtimes, ...
                rightIEDtimes, windowSeconds);

            nRows = length(tStart);
            censored = true(nRows, 1);
            eventAtStop = false(nRows, 1);

            if eventObserved
                censored(end) = false;
                eventAtStop(end) = true;
            end

            currentRows = table( ...
                repmat(patientID, nRows, 1), ...
                repmat(trial, nRows, 1), ...
                tStart, tStop, censored, eventAtStop, ...
                postIEDLeft, postIEDRight, ...
                activeLeftWindowCount, activeRightWindowCount, ...
                repmat(balloonColor(trial), nRows, 1), ...
                'VariableNames', { ...
                'patientID', 'trialNumber', 'tStart', 'tStop', ...
                'censored', 'eventAtStop', ...
                'postIEDLeft', 'postIEDRight', ...
                'activeLeftIEDWindowCount', ...
                'activeRightIEDWindowCount', ...
                'balloonColorCode'});

            countingProcessData = [countingProcessData; currentRows];

            rowDuration = tStop - tStart;
            currentSummary = table( ...
                patientID, trial, duration, double(eventObserved), ...
                length(leftIEDtimes), length(rightIEDtimes), ...
                sum(rowDuration .* postIEDLeft), ...
                sum(rowDuration .* postIEDRight), ...
                double(eventObserved && postIEDLeft(end) == 1), ...
                double(eventObserved && postIEDRight(end) == 1), ...
                balloonColor(trial), ...
                'VariableNames', { ...
                'patientID', 'trialNumber', 'durationSeconds', ...
                'eventObserved', 'nLeftIEDs', 'nRightIEDs', ...
                'leftPostIEDExposureSeconds', ...
                'rightPostIEDExposureSeconds', ...
                'eventInLeftPostIEDWindow', ...
                'eventInRightPostIEDWindow', ...
                'balloonColorCode'});

            trialSummary = [trialSummary; currentSummary];

            currentPermutationTrial = table( ...
                patientID, trial, duration, double(eventObserved), ...
                balloonColor(trial), {leftIEDtimes}, {rightIEDtimes}, ...
                'VariableNames', { ...
                'patientID', 'trialNumber', 'durationSeconds', ...
                'eventObserved', 'balloonColorCode', ...
                'leftIEDtimes_seconds', 'rightIEDtimes_seconds'});

            permutationTrialData = ...
                [permutationTrialData; currentPermutationTrial];
        end

        if patientHasLeftIED
            IEDLeftPatients(end + 1, 1) = patientID;
        end

        if patientHasRightIED
            IEDRightPatients(end + 1, 1) = patientID;
        end
    end

    if isempty(countingProcessData)
        error('%s: no valid counting-process data were created.', code);
    end

    if ~any(countingProcessData.postIEDLeft == 1)
        error('%s: no valid left-hemisphere post-IED exposure.', code);
    end

    if ~any(countingProcessData.postIEDRight == 1)
        error('%s: no valid right-hemisphere post-IED exposure.', code);
    end

    [countingProcessData.patientStratum, patientLevels] = ...
        findgroups(countingProcessData.patientID);

    permutationTrialData.patientStratum = ...
        findgroups(permutationTrialData.patientID);

    X = [ ...
        countingProcessData.postIEDLeft, ...
        countingProcessData.postIEDRight, ...
        double(countingProcessData.balloonColorCode == 2), ...
        double(countingProcessData.balloonColorCode == 3)];

    predictorNames = [ ...
        "left_hemisphere_ied_within_" + ...
        string(windowMilliseconds) + "ms"; ...
        "right_hemisphere_ied_within_" + ...
        string(windowMilliseconds) + "ms"; ...
        "orange_vs_yellow"; ...
        "red_vs_yellow"];

    T = [countingProcessData.tStart, countingProcessData.tStop];
    censoring = countingProcessData.censored;
    strata = countingProcessData.patientStratum;

    options = statset('coxphfit');
    options.Display = 'final';
    options.MaxIter = 1000;
    options.MaxFunEvals = 5000;

    [beta, logLikelihood, baselineCumulativeHazard, stats] = ...
        coxphfit(X, T, ...
        'Censoring', censoring, ...
        'Strata', strata, ...
        'Ties', 'efron', ...
        'Baseline', 0, ...
        'Options', options);

    patientNumbers = unique(strata);
    clusterScores = zeros(length(patientNumbers), length(beta));

    for patientIndex = 1:length(patientNumbers)
        clusterScores(patientIndex, :) = sum( ...
            stats.scores(strata == patientNumbers(patientIndex), :), 1);
    end

    robustCovariance = stats.covb * ...
        (clusterScores' * clusterScores) * stats.covb;

    if length(patientNumbers) > 1
        robustCovariance = robustCovariance * ...
            length(patientNumbers) / (length(patientNumbers) - 1);
    end

    robustSE = sqrt(max(diag(robustCovariance), 0));
    betaLow = beta - 1.96 * robustSE;
    betaHigh = beta + 1.96 * robustSE;
    hazardRatio = exp(beta);
    hazardRatioLow = exp(betaLow);
    hazardRatioHigh = exp(betaHigh);

    contrastVector = [1; -1; 0; 0];
    betaLeftMinusRight = contrastVector' * beta;
    contrastVariance = contrastVector' * ...
        robustCovariance * contrastVector;
    contrastSE = sqrt(max(contrastVariance, 0));
    contrastZ = betaLeftMinusRight / contrastSE;
    contrastModelP = 2 * normcdf(-abs(contrastZ));
    contrastLow = betaLeftMinusRight - 1.96 * contrastSE;
    contrastHigh = betaLeftMinusRight + 1.96 * contrastSE;

    %% Within-participant, within-color trial permutation test

    permutationLeftBeta = NaN(numberOfPermutations, 1);
    permutationRightBeta = NaN(numberOfPermutations, 1);
    permutationContrastBeta = NaN(numberOfPermutations, 1);
    permutationSucceeded = false(numberOfPermutations, 1);

    permutationGroups = findgroups( ...
        permutationTrialData.patientID, ...
        permutationTrialData.balloonColorCode);

    rng(randomSeed);
    permutationOptions = options;
    permutationOptions.Display = 'off';

    fprintf('%s: running %d permutations...\n', ...
        code, numberOfPermutations);

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
            buildShuffledHemisphereCoxData( ...
            shuffledTrials, windowSeconds);

        try
            permutationBeta = coxphfit( ...
                permutedX, permutedT, ...
                'Censoring', permutedCensoring, ...
                'Strata', permutedStrata, ...
                'Ties', 'efron', ...
                'Baseline', 0, ...
                'Options', permutationOptions);

            permutationLeftBeta(permutation) = permutationBeta(1);
            permutationRightBeta(permutation) = permutationBeta(2);
            permutationContrastBeta(permutation) = ...
                permutationBeta(1) - permutationBeta(2);
            permutationSucceeded(permutation) = true;
        catch permutationError
            fprintf('%s permutation %d failed: %s\n', ...
                code, permutation, permutationError.message);
        end

        if mod(permutation, 100) == 0 || ...
                permutation == numberOfPermutations
            fprintf('%s permutation progress: %d/%d\n', ...
                code, permutation, numberOfPermutations);
        end
    end

    validRows = permutationSucceeded;
    validLeftBeta = permutationLeftBeta(validRows);
    validRightBeta = permutationRightBeta(validRows);
    validContrastBeta = permutationContrastBeta(validRows);

    permutationPLeft = calculatePermutationP( ...
        validLeftBeta, beta(1));
    permutationPRight = calculatePermutationP( ...
        validRightBeta, beta(2));
    % This outcome-shuffle value is a global-null sensitivity measure.
    % The formal betaLeft = betaRight test is clusterRobustPValue below.
    permutationPContrast = calculatePermutationP( ...
        validContrastBeta, betaLeftMinusRight);

    if strcmp(code, 'IT')
        interpretationText = ...
            "HR > 1 means shorter IT; HR < 1 means longer IT.";
    elseif strcmp(code, 'RT')
        interpretationText = ...
            "HR > 1 means shorter RT; HR < 1 means longer RT.";
    elseif strcmp(code, 'BR')
        interpretationText = ...
            "HR > 1 means higher banking hazard; HR < 1 means lower banking hazard.";
    else
        interpretationText = ...
            "HR > 1 means higher popping hazard; HR < 1 means lower popping hazard.";
    end

    nCoverage = [ ...
        length(unique(coverageLeftPatients)); ...
        length(unique(coverageRightPatients))];
    nWithIED = [ ...
        length(unique(IEDLeftPatients)); ...
        length(unique(IEDRightPatients))];
    nIEDs = [totalLeftIEDs; totalRightIEDs];

    hemisphereResults = table( ...
        ["Left"; "Right"], ...
        predictorNames(1:2), ...
        beta(1:2), ...
        stats.se(1:2), ...
        stats.p(1:2), ...
        robustSE(1:2), ...
        [permutationPLeft; permutationPRight], ...
        betaLow(1:2), ...
        betaHigh(1:2), ...
        hazardRatio(1:2), ...
        hazardRatioLow(1:2), ...
        hazardRatioHigh(1:2), ...
        nCoverage, nWithIED, nIEDs, ...
        repmat(interpretationText, 2, 1), ...
        'VariableNames', { ...
        'hemisphere', 'predictor', 'beta_logHazard', ...
        'modelBasedSE', 'modelBasedPValue', ...
        'clusterRobustSE', 'permutationPValue', ...
        'betaCILow', 'betaCIHigh', ...
        'hazardRatio', 'hazardRatioCILow', ...
        'hazardRatioCIHigh', ...
        'nParticipantsWithCoverage', ...
        'nParticipantsWithIED', 'nIEDs', 'interpretation'});

    contrastResults = table( ...
        betaLeftMinusRight, contrastSE, contrastZ, contrastModelP, ...
        permutationPContrast, contrastLow, contrastHigh, ...
        exp(betaLeftMinusRight), exp(contrastLow), exp(contrastHigh), ...
        length(validContrastBeta), ...
        'VariableNames', { ...
        'betaLeftMinusRight', 'clusterRobustSE', ...
        'clusterRobustZ', 'clusterRobustPValue', ...
        'outcomeShufflePermutationPValue', 'betaCILow', 'betaCIHigh', ...
        'hazardRatioLeftDividedByRight', ...
        'hazardRatioRatioCILow', 'hazardRatioRatioCIHigh', ...
        'nValidPermutations'});

    permutationResults = table( ...
        (1:numberOfPermutations)', ...
        permutationLeftBeta, permutationRightBeta, ...
        permutationContrastBeta, permutationSucceeded, ...
        'VariableNames', { ...
        'permutationIndex', 'permutedLeftBeta_logHazard', ...
        'permutedRightBeta_logHazard', ...
        'permutedBetaLeftMinusRight', 'fitSucceeded'});

    fprintf('\n%s: %d-ms left/right post-IED Cox model\n', ...
        code, windowMilliseconds);
    fprintf('Participants analyzed: %d\n', length(patientLevels));
    fprintf('Trials: %d\n', height(trialSummary));
    fprintf('Left coverage / with IED / IEDs: %d / %d / %d\n', ...
        nCoverage(1), nWithIED(1), nIEDs(1));
    fprintf('Right coverage / with IED / IEDs: %d / %d / %d\n', ...
        nCoverage(2), nWithIED(2), nIEDs(2));
    fprintf('Log likelihood: %.6f\n', logLikelihood);
    disp(hemisphereResults);
    disp(contrastResults);

    result.code = code;
    result.countingProcessData = countingProcessData;
    result.trialSummary = trialSummary;
    result.hemisphereResults = hemisphereResults;
    result.contrastResults = contrastResults;
    result.beta = beta;
    result.logLikelihood = logLikelihood;
    result.baselineCumulativeHazard = baselineCumulativeHazard;
    result.stats = stats;
    result.clusterRobustCovariance = robustCovariance;
    result.clusterRobustSE = robustSE;
    result.predictorNames = predictorNames;
    result.patientLevels = patientLevels;
    result.permutationResults = permutationResults;
    result.numberOfPermutations = numberOfPermutations;
    result.numberOfValidPermutations = sum(validRows);
    result.permutationRandomSeed = randomSeed;
end


%% Read labels for the selected channels

function labels = getSelectedChannelLabels( ...
    anatomicalLocs, selectedChans, nSelectedChannels)

    labels = convertLabelsToString(anatomicalLocs);
    labels = labels(:);

    if length(labels) >= max(selectedChans)
        labels = labels(selectedChans);
    elseif length(labels) == nSelectedChannels
        % anatomicalLocs_wHemisphere was already selected-channel restricted.
    else
        error(['Cannot map anatomicalLocs_wHemisphere to selected channels. ' ...
            'Labels = %d, selected channels = %d, maximum channel = %d.'], ...
            length(labels), nSelectedChannels, max(selectedChans));
    end
end


function labels = convertLabelsToString(rawLabels)

    if isstring(rawLabels)
        labels = rawLabels;
    elseif iscell(rawLabels)
        labels = strings(numel(rawLabels), 1);

        for labelIndex = 1:numel(rawLabels)
            value = rawLabels{labelIndex};

            if isempty(value)
                labels(labelIndex) = "";
            elseif iscell(value) && numel(value) == 1
                labels(labelIndex) = string(value{1});
            else
                labels(labelIndex) = string(value);
            end
        end
    elseif ischar(rawLabels)
        labels = string(cellstr(rawLabels));
    elseif iscategorical(rawLabels)
        labels = string(rawLabels);
    else
        labels = string(rawLabels);
    end

    labels = strip(labels(:));
end


%% Identify _L and _R channel labels

function [leftMask, rightMask] = getHemisphereChannelMasks(labels)

    labels = upper(strip(string(labels)));
    labels(ismissing(labels)) = "";

    % The requested format is area_L or area_R. Hyphen/space endings are
    % also accepted so that small atlas-format differences do not discard
    % otherwise valid labels.
    leftMask = ~cellfun('isempty', regexp(cellstr(labels), ...
        '[_\-\s]L$', 'once'));
    rightMask = ~cellfun('isempty', regexp(cellstr(labels), ...
        '[_\-\s]R$', 'once'));

    leftMask = leftMask(:);
    rightMask = rightMask(:);
end


%% Convert one hemisphere's IED rows to valid within-trial times

function IEDtimes = getValidIEDTimes( ...
    hemisphereIEDRows, duration, samplingFrequency)

    if isempty(hemisphereIEDRows) || size(hemisphereIEDRows, 2) < 3
        IEDtimes = zeros(0, 1);
        return;
    end

    IEDtimes = double(hemisphereIEDRows(:, 3)) / samplingFrequency;
    IEDtimes = IEDtimes( ...
        isfinite(IEDtimes) & IEDtimes > 0 & IEDtimes < duration);
    IEDtimes = sort(IEDtimes(:));
end


%% Build start-stop intervals using both hemispheres

function [tStart, tStop, postIEDLeft, postIEDRight, ...
    activeLeftWindowCount, activeRightWindowCount] = ...
    buildHemisphereIntervals(duration, leftIEDtimes, ...
    rightIEDtimes, windowSeconds)

    leftWindowEnds = min(leftIEDtimes + windowSeconds, duration);
    rightWindowEnds = min(rightIEDtimes + windowSeconds, duration);

    breakTimes = unique([ ...
        0; leftIEDtimes; leftWindowEnds; ...
        rightIEDtimes; rightWindowEnds; duration]);

    tStart = breakTimes(1:end-1);
    tStop = breakTimes(2:end);
    midpoint = (tStart + tStop) / 2;

    postIEDLeft = zeros(length(midpoint), 1);
    postIEDRight = zeros(length(midpoint), 1);
    activeLeftWindowCount = zeros(length(midpoint), 1);
    activeRightWindowCount = zeros(length(midpoint), 1);

    for IEDindex = 1:length(leftIEDtimes)
        insideWindow = midpoint >= leftIEDtimes(IEDindex) & ...
            midpoint <= leftIEDtimes(IEDindex) + windowSeconds;
        postIEDLeft(insideWindow) = 1;
        activeLeftWindowCount(insideWindow) = ...
            activeLeftWindowCount(insideWindow) + 1;
    end

    for IEDindex = 1:length(rightIEDtimes)
        insideWindow = midpoint >= rightIEDtimes(IEDindex) & ...
            midpoint <= rightIEDtimes(IEDindex) + windowSeconds;
        postIEDRight(insideWindow) = 1;
        activeRightWindowCount(insideWindow) = ...
            activeRightWindowCount(insideWindow) + 1;
    end
end


%% Rebuild Cox data after shuffling trial outcomes/durations

function [X, T, censoring, strata] = ...
    buildShuffledHemisphereCoxData(trialData, windowSeconds)

    nTrials = height(trialData);
    Xcells = cell(nTrials, 1);
    Tcells = cell(nTrials, 1);
    censoringCells = cell(nTrials, 1);
    strataCells = cell(nTrials, 1);

    for trialIndex = 1:nTrials
        duration = trialData.durationSeconds(trialIndex);

        leftIEDtimes = trialData.leftIEDtimes_seconds{trialIndex};
        rightIEDtimes = trialData.rightIEDtimes_seconds{trialIndex};
        leftIEDtimes = leftIEDtimes(leftIEDtimes < duration);
        rightIEDtimes = rightIEDtimes(rightIEDtimes < duration);

        [tStart, tStop, postIEDLeft, postIEDRight] = ...
            buildHemisphereIntervals(duration, leftIEDtimes, ...
            rightIEDtimes, windowSeconds);

        color = trialData.balloonColorCode(trialIndex);
        nRows = length(tStart);

        Xcells{trialIndex} = [ ...
            postIEDLeft, postIEDRight, ...
            repmat(double(color == 2), nRows, 1), ...
            repmat(double(color == 3), nRows, 1)];

        Tcells{trialIndex} = [tStart, tStop];

        trialCensoring = true(nRows, 1);
        if trialData.eventObserved(trialIndex) == 1
            trialCensoring(end) = false;
        end
        censoringCells{trialIndex} = trialCensoring;

        strataCells{trialIndex} = repmat( ...
            trialData.patientStratum(trialIndex), nRows, 1);
    end

    X = vertcat(Xcells{:});
    T = vertcat(Tcells{:});
    censoring = vertcat(censoringCells{:});
    strata = vertcat(strataCells{:});
end


function permutationP = calculatePermutationP( ...
    validPermutationStatistics, observedStatistic)

    if isempty(validPermutationStatistics)
        permutationP = NaN;
        return;
    end

    permutationP = (1 + sum( ...
        abs(validPermutationStatistics) >= abs(observedStatistic))) / ...
        (length(validPermutationStatistics) + 1);
end


%% Paired left/right forest plot

function plotHemisphereForestPanel(ax, result, config, ...
    leftColor, rightColor)

    hold(ax, 'on');

    plotData = result.hemisphereResults;
    y = [1; 2];
    plotColors = [leftColor; rightColor];

    for hemisphereIndex = 1:2
        HR = plotData.hazardRatio(hemisphereIndex);
        low = plotData.hazardRatioCILow(hemisphereIndex);
        high = plotData.hazardRatioCIHigh(hemisphereIndex);

        errorbar(ax, HR, y(hemisphereIndex), ...
            HR - low, high - HR, 'horizontal', 'o', ...
            'LineWidth', 1.8, 'MarkerSize', 8, ...
            'MarkerFaceColor', plotColors(hemisphereIndex, :), ...
            'MarkerEdgeColor', plotColors(hemisphereIndex, :), ...
            'Color', plotColors(hemisphereIndex, :));
    end

    xline(ax, 1, '--', 'Color', [0.45 0.45 0.45]);

    set(ax, 'YTick', y, 'YTickLabel', {'Left', 'Right'}, ...
        'YDir', 'reverse', 'XScale', 'log', ...
        'FontSize', 10, 'FontName', 'Arial', 'TickDir', 'out');

    finiteLimits = [plotData.hazardRatioCILow; ...
        plotData.hazardRatioCIHigh; 1];
    finiteLimits = finiteLimits(isfinite(finiteLimits) & finiteLimits > 0);
    xlim(ax, [0.85 * min(finiteLimits), 1.18 * max(finiteLimits)]);
    ylim(ax, [0.5 2.5]);

    xlabel(ax, 'Hazard Ratio (95% CI)');
    title(ax, config.code, 'FontWeight', 'bold');

    annotationText = sprintf([ ...
        'Left:  HR %.3f, p_{perm} %.3g, N_{IED} %d\n' ...
        'Right: HR %.3f, p_{perm} %.3g, N_{IED} %d'], ...
        plotData.hazardRatio(1), plotData.permutationPValue(1), ...
        plotData.nIEDs(1), ...
        plotData.hazardRatio(2), plotData.permutationPValue(2), ...
        plotData.nIEDs(2));

    text(ax, 0.98, 0.04, annotationText, 'Units', 'normalized', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
        'FontName', 'Arial', 'FontSize', 8.5, ...
        'BackgroundColor', 'w', 'EdgeColor', [0.82 0.82 0.82], ...
        'Margin', 4);

    box(ax, 'off');
    hold(ax, 'off');
end


%% Direct left-versus-right contrast plot

function plotHemisphereContrastPanel(ax, result, config)

    hold(ax, 'on');

    contrast = result.contrastResults;
    ratio = contrast.hazardRatioLeftDividedByRight;
    low = contrast.hazardRatioRatioCILow;
    high = contrast.hazardRatioRatioCIHigh;

    xline(ax, 1, '--', 'Color', [0.45 0.45 0.45]);
    errorbar(ax, ratio, 1, ratio - low, high - ratio, ...
        'horizontal', 'd', 'LineWidth', 2, 'MarkerSize', 9, ...
        'MarkerFaceColor', config.color, ...
        'MarkerEdgeColor', config.color, 'Color', config.color);

    finiteLimits = [low; high; 1];
    finiteLimits = finiteLimits(isfinite(finiteLimits) & finiteLimits > 0);
    xlim(ax, [0.85 * min(finiteLimits), 1.18 * max(finiteLimits)]);
    ylim(ax, [0.5 1.5]);

    set(ax, 'YTick', 1, 'YTickLabel', {'Left / Right'}, ...
        'XScale', 'log', 'FontSize', 10, ...
        'FontName', 'Arial', 'TickDir', 'out');

    xlabel(ax, 'Left / Right Hazard Ratio (95% CI)');
    title(ax, config.code, 'FontWeight', 'bold');

    resultText = sprintf([ ...
        'HR_L / HR_R = %.3f\n95%% CI [%.3f, %.3f]\n' ...
        'robust contrast p = %.3g\noutcome-shuffle p = %.3g'], ...
        ratio, low, high, contrast.clusterRobustPValue, ...
        contrast.outcomeShufflePermutationPValue);

    text(ax, 0.98, 0.95, resultText, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
        'FontName', 'Arial', 'FontSize', 9, ...
        'BackgroundColor', 'w', 'EdgeColor', [0.82 0.82 0.82], ...
        'Margin', 5);

    text(ax, 0.02, 0.05, ...
        'Ratio > 1: higher post-IED hazard after left than right IEDs', ...
        'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
        'FontName', 'Arial', 'FontSize', 8.5, ...
        'Color', [0.35 0.35 0.35]);

    box(ax, 'off');
    hold(ax, 'off');
end
