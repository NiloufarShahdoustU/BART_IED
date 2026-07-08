% Brain-area-specific mechanistic Cox analyses with expected reward
%
% Model fitted separately for IT, RT, and BR within each anatomical area:
%
%   h_p(t) = h_0p(t) * exp[
%       beta_IED * IED(t)
%       + beta_V * V
%       + beta_Color * Color
%       + beta_IEDxV * IED(t) * V]
%
% IED(t) is 1 during an outcome-specific post-IED window and 0 otherwise.
% Primary windows: RT = 500 ms; IT = 1000 ms; BR = 1000 ms.
% V is the trial-wise expected reward read from TDdataParamRecovery.
% Participant-specific differences are handled by participant strata, so
% every participant has a separate baseline hazard h_0p(t).
%
% The primary coefficient is beta_IEDxV:
%   beta_IEDxV < 0: the expected-reward log-hazard slope is lower after an IED.
%   beta_IEDxV > 0: the expected-reward log-hazard slope is higher after an IED.
% Whether the absolute influence is weaker or stronger also depends on beta_V.
%
% Author: Nill

clear;
clc;
close all;

%% Paths

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

inputFolderName_modeling = ...
    'D:\Nill\data\BART\0_0_new_IED\context_modeling\param_recovery_1_modeling\';


outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED8_brain_area_IED_expected_reward\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

%% Shared settings

% Outcome-specific post-IED windows.
% RT is relatively short, so retain a temporally precise 500-ms window.
% IT is longer, so use a 1000-ms window to capture a longer behavioral effect.
% BR is defined during inflation time, so it uses the same 1000-ms window as IT.
settings.postIEDWindowMillisecondsRT = 500;
settings.postIEDWindowMillisecondsIT = 1000;
settings.postIEDWindowMillisecondsBR = 1000;

settings.maximumRTSeconds = 20;
settings.defaultSamplingFrequencyHz = 1000;
settings.useOnlyNonControlTrials = true;
settings.combineLeftAndRight = true;
settings.includeUnknownArea = false;

% Sparse-area safeguards
settings.minimumParticipantsWithCoverage = 2;
settings.minimumParticipantsWithIED = 1;
settings.minimumIEDsInArea = 1;
settings.minimumEndpointEventsInsideWindow = 1;
settings.minimumEndpointEventsOutsideWindow = 1;

% Panel colors
colorRT = [0.204  0.459  0.702];   % blue
colorIT = [0.847  0.333  0.153];   % orange
colorBR = [0.250  0.600  0.250];   % green
markerAlpha = 0.5;

%% Run IT, RT, and BR mechanistic analyses

[ITResults, ITModels] = runLoggedMechanisticAnalysis( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    outputFolderName, "IT", settings);

[RTResults, RTModels] = runLoggedMechanisticAnalysis( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    outputFolderName, "RT", settings);

[BRResults, BRModels] = runLoggedMechanisticAnalysis( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    outputFolderName, "BR", settings);

%% Create one 1x3 forest-plot figure for beta_IEDxV

combinedFigureOutputFile = fullfile(outputFolderName, ...
    'IT_RT_BR_mechanistic_IED_x_expected_reward_forest_plot.pdf');

plotCombinedInteractionForestPlots( ...
    ITResults, RTResults, BRResults, ...
    colorIT, colorRT, colorBR, markerAlpha, ...
    settings, combinedFigureOutputFile);

combinedResultsOutputFile = fullfile(outputFolderName, ...
    'IT_RT_BR_mechanistic_IED_x_expected_reward_all_results.mat');

save(combinedResultsOutputFile, ...
    'ITResults', 'RTResults', 'BRResults', ...
    'ITModels', 'RTModels', 'BRModels', ...
    'settings', 'colorIT', 'colorRT', 'colorBR', 'markerAlpha');

fprintf('\n============================================================\n');
fprintf('All three mechanistic analyses finished.\n');
fprintf('RT post-IED window: %d ms.\n', ...
    settings.postIEDWindowMillisecondsRT);
fprintf('IT post-IED window: %d ms.\n', ...
    settings.postIEDWindowMillisecondsIT);
fprintf('BR post-IED window: %d ms.\n', ...
    settings.postIEDWindowMillisecondsBR);
fprintf('Primary coefficient: post-IED x expected reward.\n');
fprintf('Combined forest plot saved: %s\n', combinedFigureOutputFile);
fprintf('Combined results saved: %s\n', combinedResultsOutputFile);

%% Local functions

function [areaResults, areaModels] = runLoggedMechanisticAnalysis( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    outputFolderName, analysisType, settings)

    analysisPrefix = analysisType;

    textOutputFile = fullfile(outputFolderName, ...
        char(analysisPrefix + ...
        "_mechanistic_IED_x_expected_reward_console_output.txt"));

    diary off;
    if exist(textOutputFile, 'file')
        delete(textOutputFile);
    end

    diary(textOutputFile);
    diaryCleanup = onCleanup(@() diary('off'));

    [areaResults, areaModels] = runMechanisticBrainAreaAnalysis( ...
        inputFolderName_LFPIED, inputFolderName_modeling, ...
        analysisType, settings);

    resultsOutputFile = fullfile(outputFolderName, ...
        char(analysisPrefix + ...
        "_mechanistic_IED_x_expected_reward_results.csv"));

    modelOutputFile = fullfile(outputFolderName, ...
        char(analysisPrefix + ...
        "_mechanistic_IED_x_expected_reward_models.mat"));

    writetable(areaResults, resultsOutputFile);

    save(modelOutputFile, ...
        'areaModels', 'areaResults', 'analysisType', 'settings');

    fittedRows = areaResults.status == "fitted" & ...
        isfinite(areaResults.pValue_IEDxExpectedReward);

    fprintf('\n============================================================\n');
    fprintf('%s mechanistic area-specific analysis finished.\n', ...
        analysisPrefix);
    fprintf('Fitted areas: %d\n', sum(fittedRows));
    fprintf('FDR-significant IED x expected-reward interactions: %d\n', ...
        sum(areaResults.significantFDR_IEDxExpectedReward));
    fprintf('Saved: %s\n', resultsOutputFile);
    fprintf('Saved: %s\n', modelOutputFile);
    fprintf('Saved: %s\n', textOutputFile);

    diary off;
    clear diaryCleanup;

end

function [areaResults, areaModels] = runMechanisticBrainAreaAnalysis( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    analysisType, settings)

    analysisPrefix = analysisType;

    % Choose the post-IED window separately for each behavioral outcome.
    switch analysisType
        case "RT"
            postIEDWindowMilliseconds = ...
                settings.postIEDWindowMillisecondsRT;
        case "IT"
            postIEDWindowMilliseconds = ...
                settings.postIEDWindowMillisecondsIT;
        case "BR"
            postIEDWindowMilliseconds = ...
                settings.postIEDWindowMillisecondsBR;
        otherwise
            error('Unknown analysis type: %s', analysisType);
    end

    postIEDWindowSeconds = postIEDWindowMilliseconds / 1000;
    maximumRTSeconds = settings.maximumRTSeconds;
    defaultSamplingFrequencyHz = settings.defaultSamplingFrequencyHz;
    useOnlyNonControlTrials = settings.useOnlyNonControlTrials;
    combineLeftAndRight = settings.combineLeftAndRight;
    includeUnknownArea = settings.includeUnknownArea;
    minimumParticipantsWithCoverage = ...
        settings.minimumParticipantsWithCoverage;
    minimumParticipantsWithIED = settings.minimumParticipantsWithIED;
    minimumIEDsInArea = settings.minimumIEDsInArea;
    minimumEndpointEventsInsideWindow = ...
        settings.minimumEndpointEventsInsideWindow;
    minimumEndpointEventsOutsideWindow = ...
        settings.minimumEndpointEventsOutsideWindow;

    fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));
    fileListModel = dir(fullfile( ...
        inputFolderName_modeling, '*TDdataParamRecovery.mat'));

    if isempty(fileList)
        error('No .LFPIED.mat files were found.');
    end

    if isempty(fileListModel)
        error('No *TDdataParamRecovery.mat files were found.');
    end

    if length(fileListModel) < length(fileList)
        error(['There are fewer TDdataParamRecovery files (%d) than ' ...
            'LFPIED files (%d).'], length(fileListModel), length(fileList));
    end

    participants = struct( ...
        'patientID', {}, ...
        'RTs', {}, ...
        'ITs', {}, ...
        'BankedTrials', {}, ...
        'isControl', {}, ...
        'balloonColorCode', {}, ...
        'expectedReward', {}, ...
        'validTrials', {}, ...
        'durationSeconds', {}, ...
        'finalEventObserved', {}, ...
        'IEDoccurrence', {}, ...
        'samplingFrequencyHz', {}, ...
        'selectedAreaLabels', {});

    allAreas = strings(0, 1);

    %% Load and organize participant data

    for pt = 1:length(fileList)

        fileName = fileList(pt).name;
        fileNameParts = strsplit(fileName, '.');
        patientID = string(fileNameParts{1});

        fprintf('\nLoading patient for %s: %s\n', ...
            analysisPrefix, patientID);

        loadedData = load(fullfile(inputFolderName_LFPIED, fileName));

        if ~isfield(loadedData, 'LFPIED')
            fprintf('Skipped: LFPIED structure was not found.\n');
            continue;
        end

        fileNameModel = fileListModel(pt).name;
        modelFilePath = fullfile(inputFolderName_modeling, fileNameModel);

        modelData = load(modelFilePath);

        if ~isfield(modelData, 'TDdataParamRecovery')
            fprintf('Skipped: TDdataParamRecovery structure was not found.\n');
            continue;
        end

        if ~isfield(modelData.TDdataParamRecovery, 'bestApIdx') || ...
                ~isfield(modelData.TDdataParamRecovery, 'bestAnIdx') || ...
                ~isfield(modelData.TDdataParamRecovery, 'expectedReward')
            fprintf(['Skipped: bestApIdx, bestAnIdx, or expectedReward ' ...
                'was missing from TDdataParamRecovery.\n']);
            continue;
        end

        bestApIdx = modelData.TDdataParamRecovery.bestApIdx;
        bestAnIdx = modelData.TDdataParamRecovery.bestAnIdx;

        expectedRewardsRaw = squeeze( ...
            modelData.TDdataParamRecovery.expectedReward( ...
            bestApIdx, bestAnIdx, :))';
        expectedRewardsRaw = double(expectedRewardsRaw(:));


        LFPIED = loadedData.LFPIED;

        switch analysisType
            case "IT"
                requiredFields = { ...
                    'selectedChans', 'anatomicalLocs', 'RTs', 'ITs', ...
                    'isControl', 'balloonType', 'IED_occurance_IT'};
            case "RT"
                requiredFields = { ...
                    'selectedChans', 'anatomicalLocs', 'RTs', ...
                    'isControl', 'balloonType', 'IED_occurance_RT'};
            case "BR"
                requiredFields = { ...
                    'selectedChans', 'anatomicalLocs', 'RTs', 'ITs', ...
                    'BankedTrials', 'isControl', 'balloonType', ...
                    'IED_occurance_IT'};
            otherwise
                error('Unknown analysis type: %s', analysisType);
        end

        missingField = false;

        for ff = 1:length(requiredFields)
            if ~isfield(LFPIED, requiredFields{ff})
                fprintf('Skipped: missing field %s.\n', ...
                    requiredFields{ff});
                missingField = true;
            end
        end

        if missingField
            continue;
        end

        selectedChans = round(LFPIED.selectedChans(:));
        nSelectedChannels = length(selectedChans);

        if nSelectedChannels == 0
            fprintf('Skipped: no selected channels.\n');
            continue;
        end

        selectedAreaLabels = getSelectedChannelLabels( ...
            LFPIED.anatomicalLocs, selectedChans, nSelectedChannels);

        selectedAreaLabels = cleanAreaLabels( ...
            selectedAreaLabels, combineLeftAndRight);

        excludedAreaRows = isExcludedAreaLabel(selectedAreaLabels);
        selectedAreaLabels(excludedAreaRows) = "Excluded";

        RTs = double(LFPIED.RTs(:));
        isControl = double(LFPIED.isControl(:));
        balloonType = double(LFPIED.balloonType(:));

        switch analysisType
            case "RT"
                ITs = NaN(size(RTs));
                BankedTrials = NaN(size(RTs));
                vectorLengths = [ ...
                    length(RTs), length(isControl), ...
                    length(balloonType), length(expectedRewardsRaw)];
            case "IT"
                ITs = double(LFPIED.ITs(:));
                BankedTrials = NaN(size(RTs));
                vectorLengths = [ ...
                    length(RTs), length(ITs), length(isControl), ...
                    length(balloonType), length(expectedRewardsRaw)];
            case "BR"
                ITs = double(LFPIED.ITs(:));
                BankedTrials = double(LFPIED.BankedTrials(:));
                vectorLengths = [ ...
                    length(RTs), length(ITs), length(BankedTrials), ...
                    length(isControl), length(balloonType), ...
                    length(expectedRewardsRaw)];
        end

        if isfield(LFPIED, 'nTrials') && ...
                isscalar(LFPIED.nTrials) && isfinite(LFPIED.nTrials)
            vectorLengths(end + 1) = LFPIED.nTrials;
        end

        nTrials = floor(min(vectorLengths));

        if nTrials < 1
            fprintf('Skipped: no aligned trials.\n');
            continue;
        end

        RTs = RTs(1:nTrials);
        isControl = isControl(1:nTrials);
        balloonType = balloonType(1:nTrials);
        balloonColorCode = mapBalloonColorCode(balloonType);
        ITs = ITs(1:nTrials);
        BankedTrials = BankedTrials(1:nTrials);
        expectedRewardsRaw = double(expectedRewardsRaw(1:nTrials));

        switch analysisType
            case "RT"
                durationSeconds = RTs;
                finalEventObserved = true(nTrials, 1);
                IEDoccurrence = LFPIED.IED_occurance_RT;
            case "IT"
                durationSeconds = ITs;
                finalEventObserved = true(nTrials, 1);
                IEDoccurrence = LFPIED.IED_occurance_IT;
            case "BR"
                durationSeconds = ITs;
                finalEventObserved = BankedTrials == 1;
                IEDoccurrence = LFPIED.IED_occurance_IT;
        end

        validTrials = ...
            isfinite(RTs) & ...
            RTs > 0 & ...
            RTs <= maximumRTSeconds & ...
            isfinite(durationSeconds) & ...
            durationSeconds > 0 & ...
            isfinite(balloonColorCode) & ...
            ismember(balloonColorCode, [1 2 3]) & ...
            isfinite(expectedRewardsRaw);

        if useOnlyNonControlTrials
            validTrials = validTrials & isControl == 0;
        end

        if analysisType == "BR"
            validTrials = validTrials & ...
                isfinite(BankedTrials) & ...
                ismember(BankedTrials, [0 1]);
        end

        if ~any(validTrials)
            fprintf('Skipped: no valid trials after filtering.\n');
            continue;
        end

        expectedReward = expectedRewardsRaw;

        samplingFrequencyHz = getSamplingFrequency( ...
            LFPIED, defaultSamplingFrequencyHz);

        newParticipant.patientID = patientID;
        newParticipant.RTs = RTs;
        newParticipant.ITs = ITs;
        newParticipant.BankedTrials = BankedTrials;
        newParticipant.isControl = isControl;
        newParticipant.balloonColorCode = balloonColorCode;
        newParticipant.expectedReward = expectedReward;
        newParticipant.validTrials = validTrials;
        newParticipant.durationSeconds = durationSeconds;
        newParticipant.finalEventObserved = finalEventObserved;
        newParticipant.IEDoccurrence = IEDoccurrence;
        newParticipant.samplingFrequencyHz = samplingFrequencyHz;
        newParticipant.selectedAreaLabels = selectedAreaLabels;

        participants(end + 1) = newParticipant; 

        participantAreas = unique(selectedAreaLabels, 'stable');

        if ~includeUnknownArea
            participantAreas = participantAreas( ...
                participantAreas ~= "Unknown");
        end

        participantAreas = participantAreas( ...
            participantAreas ~= "Excluded");

        allAreas = unique([allAreas; participantAreas], 'stable');

        fprintf('Expected reward file: %s\n', modelFilePath);
        fprintf('Aligned trials: %d; valid trials: %d\n', ...
            nTrials, sum(validTrials));

    end

    if isempty(participants)
        error('No valid participant data were loaded for %s.', ...
            analysisPrefix);
    end

    allAreas = allAreas(~isExcludedAreaLabel(allAreas));

    areaResults = initializeMechanisticResultsTable();
    areaModels = struct();

    if isempty(allAreas)
        error(['No anatomical areas were found after excluding NaC, ' ...
            'nucleus accumbens, white matter, and lateral ventricle.']);
    end

    fprintf('\n============================================================\n');
    fprintf('%s area-specific %d-ms mechanistic Cox analysis\n', ...
        analysisPrefix, postIEDWindowMilliseconds);
    fprintf('Primary test: post-IED x expected reward\n');
    fprintf('Participant effect: participant-stratified baseline hazard\n');
    fprintf('============================================================\n');
    fprintf('Loaded participants: %d\n', length(participants));
    fprintf('Anatomical areas analyzed: %d\n', length(allAreas));

    %% Fit one mechanistic Cox model per anatomical area

    for aa = 1:length(allAreas)

        anatomicalArea = allAreas(aa);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Area %d/%d: %s\n', ...
            aa, length(allAreas), anatomicalArea);

        countingProcessData = emptyMechanisticCountingProcessTable();

        nParticipantsWithCoverage = 0;
        nParticipantsWithIED = 0;
        nIEDsInArea = 0;
        nValidTrials = 0;

        for pt = 1:length(participants)

            P = participants(pt);
            areaLocalChannels = find( ...
                P.selectedAreaLabels == anatomicalArea);

            if isempty(areaLocalChannels)
                continue;
            end

            nParticipantsWithCoverage = ...
                nParticipantsWithCoverage + 1;

            areaIEDoccurrence = filterIEDsToArea( ...
                P.IEDoccurrence, areaLocalChannels);

            validAreaIEDRows = getValidAreaIEDRows( ...
                areaIEDoccurrence, P.validTrials);

            thisParticipantIEDCount = sum(validAreaIEDRows);
            nIEDsInArea = nIEDsInArea + thisParticipantIEDCount;

            if thisParticipantIEDCount > 0
                nParticipantsWithIED = nParticipantsWithIED + 1;
            end

            validTrialNumbers = find(P.validTrials);
            nValidTrials = nValidTrials + length(validTrialNumbers);

            for tt = 1:length(validTrialNumbers)

                trialNumber = validTrialNumbers(tt);

                trialRows = makeMechanisticPostIEDRows( ...
                    P.patientID, ...
                    trialNumber, ...
                    P.durationSeconds(trialNumber), ...
                    areaIEDoccurrence, ...
                    P.samplingFrequencyHz, ...
                    P.balloonColorCode(trialNumber), ...
                    P.expectedReward(trialNumber), ...
                    P.finalEventObserved(trialNumber), ...
                    postIEDWindowSeconds);

                countingProcessData = ...
                    [countingProcessData; trialRows]; 

            end

        end

        status = "fitted";

        if nParticipantsWithCoverage < minimumParticipantsWithCoverage
            status = "skipped: insufficient participants with coverage";
        elseif nParticipantsWithIED < minimumParticipantsWithIED
            status = "skipped: no participants with area IEDs";
        elseif nIEDsInArea < minimumIEDsInArea
            status = "skipped: insufficient area IEDs";
        elseif isempty(countingProcessData)
            status = "skipped: no counting-process rows";
        end

        nCountingRows = height(countingProcessData);
        nEndpointEvents = 0;
        nEventsInsideWindow = 0;
        nEventsOutsideWindow = 0;
        timeInsideWindowSeconds = 0;
        timeOutsideWindowSeconds = 0;
        nModelParticipants = 0;
        referenceBalloonColor = "";
        logLikelihood = NaN;

        IEDSummary = emptyCoefficientSummary();
        VSummary = emptyCoefficientSummary();
        interactionSummary = emptyCoefficientSummary();

        if status == "fitted"

            [countingProcessData.patientStratum, patientLevels] = ...
                findgroups(countingProcessData.patientID);

            [X, predictorNames, referenceBalloonColor] = ...
                buildMechanisticPredictorMatrix(countingProcessData);

            T = [countingProcessData.tStart, ...
                countingProcessData.tStop];
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

            nCountingRows = height(countingProcessData);
            nModelParticipants = ...
                length(unique(countingProcessData.patientID));
            nEndpointEvents = sum( ...
                countingProcessData.eventAtStop == 1);

            nEventsOutsideWindow = sum( ...
                countingProcessData.eventAtStop == 1 & ...
                countingProcessData.postIED == 0);

            nEventsInsideWindow = sum( ...
                countingProcessData.eventAtStop == 1 & ...
                countingProcessData.postIED == 1);

            timeOutsideWindowSeconds = sum( ...
                (countingProcessData.tStop - ...
                 countingProcessData.tStart) .* ...
                (countingProcessData.postIED == 0), ...
                'omitnan');

            timeInsideWindowSeconds = sum( ...
                (countingProcessData.tStop - ...
                 countingProcessData.tStart) .* ...
                (countingProcessData.postIED == 1), ...
                'omitnan');

            if nModelParticipants < minimumParticipantsWithCoverage
                status = "skipped: insufficient model participants";
            elseif timeInsideWindowSeconds <= 0
                status = "skipped: no area-specific post-IED time";
            elseif nEventsInsideWindow < ...
                    minimumEndpointEventsInsideWindow
                status = "skipped: no endpoint events inside post-IED windows";
            elseif nEventsOutsideWindow < ...
                    minimumEndpointEventsOutsideWindow
                status = "skipped: no endpoint events outside post-IED windows";
            elseif size(X, 1) <= size(X, 2)
                status = "skipped: insufficient rows for predictor matrix";
            end

        else
            patientLevels = strings(0, 1);
            predictorNames = strings(0, 1);
        end

        if status == "fitted"

            try
                coxOptions = statset('coxphfit');
                coxOptions.Display = 'off';
                coxOptions.MaxIter = 1000;
                coxOptions.MaxFunEvals = 5000;

                [beta, logLikelihood, ...
                    baselineCumulativeHazard, stats] = coxphfit( ...
                        X, ...
                        T, ...
                        'Censoring', censoring, ...
                        'Strata', strata, ...
                        'Ties', 'efron', ...
                        'Baseline', 0, ...
                        'Options', coxOptions);

                [clusterRobustCovariance, clusterRobustSE, ...
                    clusterRobustZ, clusterRobustP] = ...
                    computeClusterRobustInference( ...
                        stats, beta, ...
                        countingProcessData.patientStratum);

                IEDIndex = find( ...
                    predictorNames == "post_ied_indicator", 1);
                VIndex = find( ...
                    predictorNames == "expected_reward", 1);
                interactionIndex = find( ...
                    predictorNames == ...
                    "post_ied_x_expected_reward", 1);

                if isempty(IEDIndex) || isempty(VIndex) || ...
                        isempty(interactionIndex)
                    error('Required mechanistic predictors were not found.');
                end

                IEDSummary = summarizeCoefficient( ...
                    beta, clusterRobustSE, clusterRobustZ, ...
                    clusterRobustP, IEDIndex);

                VSummary = summarizeCoefficient( ...
                    beta, clusterRobustSE, clusterRobustZ, ...
                    clusterRobustP, VIndex);

                interactionSummary = summarizeCoefficient( ...
                    beta, clusterRobustSE, clusterRobustZ, ...
                    clusterRobustP, interactionIndex);

                safeAreaNameBase = ...
                    matlab.lang.makeValidName(char(anatomicalArea));
                safeAreaName = safeAreaNameBase;
                duplicateNumber = 1;

                while isfield(areaModels, safeAreaName)
                    duplicateNumber = duplicateNumber + 1;
                    safeAreaName = sprintf('%s_%d', ...
                        safeAreaNameBase, duplicateNumber);
                end

                areaModels.(safeAreaName) = struct( ...
                    'anatomicalArea', anatomicalArea, ...
                    'analysisType', analysisType, ...
                    'postIEDWindowMilliseconds', ...
                        postIEDWindowMilliseconds, ...
                    'beta', beta, ...
                    'predictorNames', predictorNames, ...
                    'primaryPredictorName', ...
                        "post_ied_x_expected_reward", ...
                    'logLikelihood', logLikelihood, ...
                    'baselineCumulativeHazard', ...
                        baselineCumulativeHazard, ...
                    'stats', stats, ...
                    'clusterRobustCovariance', ...
                        clusterRobustCovariance, ...
                    'clusterRobustSE', clusterRobustSE, ...
                    'clusterRobustZ', clusterRobustZ, ...
                    'clusterRobustP', clusterRobustP, ...
                    'patientLevels', patientLevels, ...
                    'referenceBalloonColor', ...
                        referenceBalloonColor, ...
                    'expectedRewardDefinition', ...
                        "trial-wise expected reward from TDdataParamRecovery", ...
                    'participantEffect', ...
                        "participant-stratified baseline hazard");

            catch modelError
                status = "failed: " + string(modelError.message);
            end

        end

        newRow = table( ...
            anatomicalArea, ...
            postIEDWindowMilliseconds, ...
            status, ...
            nParticipantsWithCoverage, ...
            nParticipantsWithIED, ...
            nModelParticipants, ...
            nValidTrials, ...
            nIEDsInArea, ...
            nCountingRows, ...
            nEndpointEvents, ...
            nEventsOutsideWindow, ...
            nEventsInsideWindow, ...
            timeOutsideWindowSeconds, ...
            timeInsideWindowSeconds, ...
            referenceBalloonColor, ...
            logLikelihood, ...
            IEDSummary.beta, ...
            IEDSummary.se, ...
            IEDSummary.z, ...
            IEDSummary.p, ...
            IEDSummary.hr, ...
            IEDSummary.hrLow, ...
            IEDSummary.hrHigh, ...
            VSummary.beta, ...
            VSummary.se, ...
            VSummary.z, ...
            VSummary.p, ...
            VSummary.hr, ...
            VSummary.hrLow, ...
            VSummary.hrHigh, ...
            interactionSummary.beta, ...
            interactionSummary.se, ...
            interactionSummary.z, ...
            interactionSummary.p, ...
            interactionSummary.betaLow, ...
            interactionSummary.betaHigh, ...
            interactionSummary.hr, ...
            interactionSummary.hrLow, ...
            interactionSummary.hrHigh, ...
            NaN, ...
            NaN, ...
            false, ...
            'VariableNames', ...
                areaResults.Properties.VariableNames);

        areaResults = [areaResults; newRow]; 

        fprintf('Status: %s\n', status);
        fprintf('Coverage participants: %d\n', ...
            nParticipantsWithCoverage);
        fprintf('Participants with area IEDs: %d\n', ...
            nParticipantsWithIED);
        fprintf('Area IEDs: %d\n', nIEDsInArea);

        if status == "fitted"
            fprintf('beta_IED = %.6f, p = %.6g\n', ...
                IEDSummary.beta, IEDSummary.p);
            fprintf('beta_V = %.6f, p = %.6g\n', ...
                VSummary.beta, VSummary.p);
            fprintf([ ...
                'PRIMARY beta_IEDxV = %.6f, interaction HR = %.6f, ' ...
                '95%% CI = [%.6f, %.6f], p = %.6g\n'], ...
                interactionSummary.beta, interactionSummary.hr, ...
                interactionSummary.hrLow, interactionSummary.hrHigh, ...
                interactionSummary.p);
        end

    end

    %% Hard removal of excluded labels

    excludedResultRows = ...
        isExcludedAreaLabel(areaResults.anatomicalArea);

    if any(excludedResultRows)
        areaResults(excludedResultRows, :) = [];
    end

    %% Multiple-comparison correction for the PRIMARY interaction

    fittedRows = areaResults.status == "fitted" & ...
        isfinite(areaResults.pValue_IEDxExpectedReward);

    if any(fittedRows)
        areaResults.pValueFDR_IEDxExpectedReward(fittedRows) = ...
            benjaminiHochberg( ...
            areaResults.pValue_IEDxExpectedReward(fittedRows));

        numberOfFittedAreas = sum(fittedRows);
        areaResults.pValueBonferroni_IEDxExpectedReward(fittedRows) = ...
            min(areaResults.pValue_IEDxExpectedReward(fittedRows) .* ...
            numberOfFittedAreas, 1);

        areaResults.significantFDR_IEDxExpectedReward(fittedRows) = ...
            areaResults.pValueFDR_IEDxExpectedReward(fittedRows) < 0.05;
    end

    areaResults = sortrows( ...
        areaResults, ...
        {'status', 'pValue_IEDxExpectedReward'}, ...
        {'ascend', 'ascend'});

end

function results = initializeMechanisticResultsTable()

    results = table( ...
        strings(0, 1), ...  % anatomicalArea
        zeros(0, 1), ...    % postIEDWindowMilliseconds
        strings(0, 1), ...  % status
        zeros(0, 1), ...    % nParticipantsWithCoverage
        zeros(0, 1), ...    % nParticipantsWithIED
        zeros(0, 1), ...    % nModelParticipants
        zeros(0, 1), ...    % nValidTrials
        zeros(0, 1), ...    % nIEDsInArea
        zeros(0, 1), ...    % nCountingRows
        zeros(0, 1), ...    % nEndpointEvents
        zeros(0, 1), ...    % nEventsOutsideWindow
        zeros(0, 1), ...    % nEventsInsideWindow
        zeros(0, 1), ...    % timeOutsideWindowSeconds
        zeros(0, 1), ...    % timeInsideWindowSeconds
        strings(0, 1), ...  % referenceBalloonColor
        zeros(0, 1), ...    % logLikelihood
        zeros(0, 1), ...    % beta_IED
        zeros(0, 1), ...    % clusterRobustSE_IED
        zeros(0, 1), ...    % clusterRobustZ_IED
        zeros(0, 1), ...    % pValue_IED
        zeros(0, 1), ...    % hazardRatio_IED
        zeros(0, 1), ...    % hazardRatioCILow_IED
        zeros(0, 1), ...    % hazardRatioCIHigh_IED
        zeros(0, 1), ...    % beta_ExpectedReward
        zeros(0, 1), ...    % clusterRobustSE_ExpectedReward
        zeros(0, 1), ...    % clusterRobustZ_ExpectedReward
        zeros(0, 1), ...    % pValue_ExpectedReward
        zeros(0, 1), ...    % hazardRatio_ExpectedReward
        zeros(0, 1), ...    % hazardRatioCILow_ExpectedReward
        zeros(0, 1), ...    % hazardRatioCIHigh_ExpectedReward
        zeros(0, 1), ...    % beta_IEDxExpectedReward
        zeros(0, 1), ...    % clusterRobustSE_IEDxExpectedReward
        zeros(0, 1), ...    % clusterRobustZ_IEDxExpectedReward
        zeros(0, 1), ...    % pValue_IEDxExpectedReward
        zeros(0, 1), ...    % betaCILow_IEDxExpectedReward
        zeros(0, 1), ...    % betaCIHigh_IEDxExpectedReward
        zeros(0, 1), ...    % interactionHazardRatio
        zeros(0, 1), ...    % interactionHazardRatioCILow
        zeros(0, 1), ...    % interactionHazardRatioCIHigh
        zeros(0, 1), ...    % pValueFDR_IEDxExpectedReward
        zeros(0, 1), ...    % pValueBonferroni_IEDxExpectedReward
        false(0, 1), ...    % significantFDR_IEDxExpectedReward
        'VariableNames', { ...
            'anatomicalArea', ...
            'postIEDWindowMilliseconds', ...
            'status', ...
            'nParticipantsWithCoverage', ...
            'nParticipantsWithIED', ...
            'nModelParticipants', ...
            'nValidTrials', ...
            'nIEDsInArea', ...
            'nCountingRows', ...
            'nEndpointEvents', ...
            'nEventsOutsideWindow', ...
            'nEventsInsideWindow', ...
            'timeOutsideWindowSeconds', ...
            'timeInsideWindowSeconds', ...
            'referenceBalloonColor', ...
            'logLikelihood', ...
            'beta_IED', ...
            'clusterRobustSE_IED', ...
            'clusterRobustZ_IED', ...
            'pValue_IED', ...
            'hazardRatio_IED', ...
            'hazardRatioCILow_IED', ...
            'hazardRatioCIHigh_IED', ...
            'beta_ExpectedReward', ...
            'clusterRobustSE_ExpectedReward', ...
            'clusterRobustZ_ExpectedReward', ...
            'pValue_ExpectedReward', ...
            'hazardRatio_ExpectedReward', ...
            'hazardRatioCILow_ExpectedReward', ...
            'hazardRatioCIHigh_ExpectedReward', ...
            'beta_IEDxExpectedReward', ...
            'clusterRobustSE_IEDxExpectedReward', ...
            'clusterRobustZ_IEDxExpectedReward', ...
            'pValue_IEDxExpectedReward', ...
            'betaCILow_IEDxExpectedReward', ...
            'betaCIHigh_IEDxExpectedReward', ...
            'interactionHazardRatio', ...
            'interactionHazardRatioCILow', ...
            'interactionHazardRatioCIHigh', ...
            'pValueFDR_IEDxExpectedReward', ...
            'pValueBonferroni_IEDxExpectedReward', ...
            'significantFDR_IEDxExpectedReward'});

end

function labels = getSelectedChannelLabels( ...
    anatomicalLocs, selectedChans, nSelectedChannels)

    labels = convertLabelsToString(anatomicalLocs);
    labels = labels(:);

    if length(labels) >= max(selectedChans)
        labels = labels(selectedChans);
    elseif length(labels) == nSelectedChannels
        % anatomicalLocs was already restricted to selected channels.
    else
        error(['Cannot map anatomicalLocs to selected channels. ' ...
            'Labels = %d, selected channels = %d, maximum channel = %d.'], ...
            length(labels), nSelectedChannels, max(selectedChans));
    end

end

function labels = convertLabelsToString(rawLabels)

    if isstring(rawLabels)
        labels = rawLabels;
    elseif iscell(rawLabels)
        labels = strings(numel(rawLabels), 1);
        for ii = 1:numel(rawLabels)
            value = rawLabels{ii};
            if isempty(value)
                labels(ii) = "";
            elseif iscell(value) && numel(value) == 1
                labels(ii) = string(value{1});
            else
                labels(ii) = string(value);
            end
        end
    elseif ischar(rawLabels)
        labels = string(cellstr(rawLabels));
    elseif iscategorical(rawLabels)
        labels = string(rawLabels);
    else
        labels = string(rawLabels);
    end

    labels = labels(:);

end

function labels = cleanAreaLabels(labels, combineLeftAndRight)

    labels = strip(string(labels));
    labels(ismissing(labels) | strlength(labels) == 0) = "Unknown";

    if combineLeftAndRight
        labels = regexprep(labels, ...
            '^\s*\(?\s*(Left|Right|LH|RH|L|R)\s*\)?[_\-\s]+', ...
            '', 'ignorecase');

        labels = regexprep(labels, ...
            '[_\-\s]+\(?\s*(Left|Right|LH|RH|L|R)\s*\)?\s*$', ...
            '', 'ignorecase');

        labels = regexprep(labels, ...
            '^\s*(Left|Right)\s+hemisphere[_\-\s]+', ...
            '', 'ignorecase');

        labels = regexprep(labels, ...
            '[_\-\s]+(Left|Right)\s+hemisphere\s*$', ...
            '', 'ignorecase');

        labels = regexprep(labels, '[_\-]+', ' ');
        labels = regexprep(labels, '\s+', ' ');
        labels = strip(labels);
        labels(strlength(labels) == 0) = "Unknown";
    end

end

function excludedRows = isExcludedAreaLabel(labels)

    normalizedLabels = lower(strip(string(labels)));
    normalizedLabels = regexprep(normalizedLabels, '[_\-]+', ' ');
    normalizedLabels = regexprep(normalizedLabels, '\s+', ' ');

    % Remove EVERY label containing NaC, including NaC60, L-NaC12, etc.
    isNaC = ...
        contains(normalizedLabels, "nac") | ...
        contains(normalizedLabels, "nucleus accumbens") | ...
        contains(normalizedLabels, "accumbens");

    % Remove white-matter labels, including WM, WM60, white_matter, etc.
    compactLabels = regexprep(normalizedLabels, '[^a-z0-9]', '');

    isWhiteMatter = ...
        contains(normalizedLabels, "white matter") | ...
        contains(compactLabels, "whitematter") | ...
        startsWith(compactLabels, "wm");

    % Remove lateral-ventricle labels, case-insensitively, including
    % variants such as lateral_ventricle and Lateral-Ventricle.
    isLateralVentricle = ...
        contains(normalizedLabels, "lateral ventricle") | ...
        contains(compactLabels, "lateralventricle");

    excludedRows = isNaC | isWhiteMatter | isLateralVentricle;

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

function areaIEDoccurrence = filterIEDsToArea( ...
    IEDoccurrence, areaLocalChannels)

    if isempty(IEDoccurrence) || size(IEDoccurrence, 2) < 3
        areaIEDoccurrence = zeros(0, 3);
        return;
    end

    localChannelIndex = round(IEDoccurrence(:, 2));

    keepRows = ...
        isfinite(localChannelIndex) & ...
        ismember(localChannelIndex, areaLocalChannels);

    areaIEDoccurrence = IEDoccurrence(keepRows, :);

end

function validRows = getValidAreaIEDRows( ...
    areaIEDoccurrence, validTrialMask)

    if isempty(areaIEDoccurrence) || size(areaIEDoccurrence, 2) < 3
        validRows = false(0, 1);
        return;
    end

    trialNumber = round(areaIEDoccurrence(:, 1));
    nTrials = length(validTrialMask);

    validRows = ...
        isfinite(trialNumber) & ...
        trialNumber >= 1 & ...
        trialNumber <= nTrials & ...
        isfinite(areaIEDoccurrence(:, 3)) & ...
        areaIEDoccurrence(:, 3) >= 1;

    validIndices = find(validRows);

    if ~isempty(validIndices)
        validRows(validIndices) = validTrialMask(trialNumber(validIndices));
    end

end



function countingProcessData = emptyMechanisticCountingProcessTable()

    countingProcessData = table( ...
        strings(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        false(0, 1), ...
        false(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        'VariableNames', { ...
            'patientID', ...
            'trialNumber', ...
            'tStart', ...
            'tStop', ...
            'censored', ...
            'eventAtStop', ...
            'postIED', ...
            'activeIEDWindowCount', ...
            'timeSinceMostRecentIED_seconds', ...
            'balloonColorCode', ...
            'expectedReward', ...
            'samplingFrequencyHz'});

end

function trialRows = makeMechanisticPostIEDRows( ...
    patientID, trialNumber, durationSeconds, IEDoccurrence, ...
    samplingFrequencyHz, balloonColorCode, expectedReward, ...
    finalEventObserved, postIEDWindowSeconds)

    if isempty(IEDoccurrence) || size(IEDoccurrence, 2) < 3
        IEDtimes = [];
    else
        validIEDRows = ...
            isfinite(IEDoccurrence(:, 1)) & ...
            isfinite(IEDoccurrence(:, 3)) & ...
            round(IEDoccurrence(:, 1)) == trialNumber & ...
            IEDoccurrence(:, 3) >= 1;

        sampleIndices = double(IEDoccurrence(validIEDRows, 3));
        IEDtimes = sampleIndices ./ samplingFrequencyHz;

        IEDtimes = IEDtimes( ...
            isfinite(IEDtimes) & ...
            IEDtimes > 0 & ...
            IEDtimes < durationSeconds);

        % Keep every IED occurrence, including duplicate times.
        IEDtimes = sort(IEDtimes);
    end

    if isempty(IEDtimes)
        breakTimes = [0; durationSeconds];
    else
        windowEndTimes = min( ...
            IEDtimes + postIEDWindowSeconds, durationSeconds);

        breakTimes = unique([ ...
            0; IEDtimes(:); windowEndTimes(:); durationSeconds], ...
            'sorted');
    end

    breakTimes = breakTimes( ...
        isfinite(breakTimes) & ...
        breakTimes >= 0 & ...
        breakTimes <= durationSeconds);

    breakTimes = unique(breakTimes, 'sorted');

    if length(breakTimes) < 2
        breakTimes = [0; durationSeconds];
    end

    tStart = breakTimes(1:end-1);
    tStop = breakTimes(2:end);

    validIntervals = tStop > tStart;
    tStart = tStart(validIntervals);
    tStop = tStop(validIntervals);

    intervalMidpoint = (tStart + tStop) ./ 2;

    postIED = zeros(length(tStart), 1);
    activeIEDWindowCount = zeros(length(tStart), 1);
    timeSinceMostRecentIED = NaN(length(tStart), 1);

    for kk = 1:length(tStart)

        if isempty(IEDtimes)
            continue;
        end

        activeWindows = ...
            intervalMidpoint(kk) >= IEDtimes & ...
            intervalMidpoint(kk) <= ...
                (IEDtimes + postIEDWindowSeconds);

        activeIEDWindowCount(kk) = sum(activeWindows);
        postIED(kk) = double(any(activeWindows));

        priorIEDtimes = IEDtimes(IEDtimes < intervalMidpoint(kk));

        if ~isempty(priorIEDtimes)
            timeSinceMostRecentIED(kk) = ...
                intervalMidpoint(kk) - max(priorIEDtimes);
        end

    end

    censored = true(length(tStart), 1);
    eventAtStop = false(length(tStart), 1);

    if finalEventObserved
        censored(end) = false;
        eventAtStop(end) = true;
    end

    trialRows = table( ...
        repmat(patientID, length(tStart), 1), ...
        repmat(trialNumber, length(tStart), 1), ...
        tStart, ...
        tStop, ...
        censored, ...
        eventAtStop, ...
        postIED, ...
        activeIEDWindowCount, ...
        timeSinceMostRecentIED, ...
        repmat(balloonColorCode, length(tStart), 1), ...
        repmat(expectedReward, length(tStart), 1), ...
        repmat(samplingFrequencyHz, length(tStart), 1), ...
        'VariableNames', { ...
            'patientID', ...
            'trialNumber', ...
            'tStart', ...
            'tStop', ...
            'censored', ...
            'eventAtStop', ...
            'postIED', ...
            'activeIEDWindowCount', ...
            'timeSinceMostRecentIED_seconds', ...
            'balloonColorCode', ...
            'expectedReward', ...
            'samplingFrequencyHz'});

end

function [X, predictorNames, referenceColorName] = ...
    buildMechanisticPredictorMatrix(countingProcessData)

    postIED = double(countingProcessData.postIED);
    expectedReward = double(countingProcessData.expectedReward);

    X = [postIED, expectedReward];
    predictorNames = [ ...
        "post_ied_indicator"; ...
        "expected_reward"];

    observedColors = unique( ...
        countingProcessData.balloonColorCode, 'sorted');

    observedColors = observedColors(ismember(observedColors, [1 2 3]));

    if isempty(observedColors)
        error('No valid balloon colors were found.');
    end

    if ismember(1, observedColors)
        referenceColor = 1;
    else
        referenceColor = observedColors(1);
    end

    colorNames = ["yellow", "orange", "red"];
    referenceColorName = colorNames(referenceColor);

    comparisonColors = observedColors(observedColors ~= referenceColor);

    for cc = 1:length(comparisonColors)
        colorCode = comparisonColors(cc);
        X(:, end + 1) = double( ...
            countingProcessData.balloonColorCode == colorCode); 
        predictorNames(end + 1, 1) = ...
            colorNames(colorCode) + "_vs_" + ...
            colorNames(referenceColor); 
    end

    X(:, end + 1) = postIED .* expectedReward;
    predictorNames(end + 1, 1) = ...
        "post_ied_x_expected_reward";

end

function summary = emptyCoefficientSummary()

    summary = struct( ...
        'beta', NaN, ...
        'se', NaN, ...
        'z', NaN, ...
        'p', NaN, ...
        'betaLow', NaN, ...
        'betaHigh', NaN, ...
        'hr', NaN, ...
        'hrLow', NaN, ...
        'hrHigh', NaN);

end

function summary = summarizeCoefficient( ...
    beta, robustSE, robustZ, robustP, coefficientIndex)

    summary = emptyCoefficientSummary();

    summary.beta = beta(coefficientIndex);
    summary.se = robustSE(coefficientIndex);
    summary.z = robustZ(coefficientIndex);
    summary.p = robustP(coefficientIndex);
    summary.betaLow = summary.beta - 1.96 * summary.se;
    summary.betaHigh = summary.beta + 1.96 * summary.se;
    summary.hr = exp(summary.beta);
    summary.hrLow = exp(summary.betaLow);
    summary.hrHigh = exp(summary.betaHigh);

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
    robustCovariance = ...
        (nClusters / (nClusters - 1)) .* robustCovariance;

    robustSE = sqrt(max(diag(robustCovariance), 0));

    zeroSERows = robustSE <= 0 | ~isfinite(robustSE);
    robustSE(zeroSERows) = stats.se(zeroSERows);

    robustZ = beta ./ robustSE;
    robustP = 2 .* normcdf(-abs(robustZ), 0, 1);

end

function adjustedP = benjaminiHochberg(pValues)

    pValues = pValues(:);
    adjustedP = NaN(size(pValues));

    validRows = isfinite(pValues);
    validP = pValues(validRows);

    if isempty(validP)
        return;
    end

    [sortedP, sortIndex] = sort(validP, 'ascend');
    m = length(sortedP);

    sortedAdjusted = sortedP .* m ./ (1:m)';

    for ii = m-1:-1:1
        sortedAdjusted(ii) = min( ...
            sortedAdjusted(ii), sortedAdjusted(ii + 1));
    end

    sortedAdjusted = min(sortedAdjusted, 1);

    unsortedAdjusted = NaN(m, 1);
    unsortedAdjusted(sortIndex) = sortedAdjusted;
    adjustedP(validRows) = unsortedAdjusted;

end



function plotCombinedInteractionForestPlots( ...
    ITResults, RTResults, BRResults, ...
    colorIT, colorRT, colorBR, markerAlpha, settings, outputPDF)

    maximumNumberOfAreas = max([ ...
        countPlottableInteractionAreas(ITResults), ...
        countPlottableInteractionAreas(RTResults), ...
        countPlottableInteractionAreas(BRResults), ...
        1]);

    figureHeight = max(650, 34 * maximumNumberOfAreas + 230);

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [80 80 2100 figureHeight]);

    layout = tiledlayout(fig, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    axIT = nexttile(layout, 1);
    plotInteractionForestPanel( ...
        axIT, ITResults, colorIT, markerAlpha, ...
        "Inflation time", "IT", ...
        settings.postIEDWindowMillisecondsIT);

    axRT = nexttile(layout, 2);
    plotInteractionForestPanel( ...
        axRT, RTResults, colorRT, markerAlpha, ...
        "Response time", "RT", ...
        settings.postIEDWindowMillisecondsRT);

    axBR = nexttile(layout, 3);
    plotInteractionForestPanel( ...
        axBR, BRResults, colorBR, markerAlpha, ...
        "Bank rate", "BR", ...
        settings.postIEDWindowMillisecondsBR);

    exportgraphics(fig, outputPDF, 'ContentType', 'vector');
    close(fig);

end

function numberOfAreas = countPlottableInteractionAreas(areaResults)

    plotRows = getInteractionForestPlotRows(areaResults);
    numberOfAreas = sum(plotRows);

end

function plotInteractionForestPanel( ...
    ax, areaResults, panelColor, markerAlpha, ...
    descriptiveTitle, abbreviation, postIEDWindowMilliseconds)

    plotRows = getInteractionForestPlotRows(areaResults);
    plotData = areaResults(plotRows, :);

    hold(ax, 'on');

    xline(ax, 0, '--', ...
        'LineWidth', 1.2, ...
        'Color', [0.25 0.25 0.25]);

    if isempty(plotData)
        ax.XLim = [-1 1];
        ax.YLim = [0 1];
        ax.YTick = [];
        text(ax, 0, 0.5, ...
            'No FDR-significant IED x value interactions', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 11);
    else
        plotData = sortrows( ...
            plotData, 'beta_IEDxExpectedReward', 'ascend');
        nAreas = height(plotData);

        for ii = 1:nAreas
            plot(ax, ...
                [plotData.betaCILow_IEDxExpectedReward(ii), ...
                 plotData.betaCIHigh_IEDxExpectedReward(ii)], ...
                [ii, ii], ...
                '-', ...
                'Color', panelColor, ...
                'LineWidth', 1.5);
        end

        scatter(ax, ...
            plotData.beta_IEDxExpectedReward, ...
            (1:nAreas)', ...
            65, ...
            panelColor, ...
            'filled', ...
            'MarkerFaceAlpha', markerAlpha, ...
            'MarkerEdgeColor', 'none');

        ax.YTick = 1:nAreas;
        ax.YTickLabel = plotData.anatomicalArea;
        ax.YLim = [0.4, nAreas + 0.6];
    end

    ax.TickLabelInterpreter = 'none';
    ax.TickDir = 'out';
    ax.FontSize = 10;

    xlabel(ax, ...
        'Interaction coefficient: post-IED x expected reward');

    windowLabel = sprintf('%s, %d-ms post-IED window', ...
        char(abbreviation), postIEDWindowMilliseconds);

    title(ax, {char(descriptiveTitle), windowLabel}, ...
        'FontWeight', 'bold');

    box(ax, 'off');
    grid(ax, 'off');

end

function plotRows = getInteractionForestPlotRows(areaResults)

    if isempty(areaResults)
        plotRows = false(0, 1);
        return;
    end

    plotRows = ...
        ~isExcludedAreaLabel(areaResults.anatomicalArea) & ...
        areaResults.status == "fitted" & ...
        areaResults.significantFDR_IEDxExpectedReward & ...
        isfinite(areaResults.beta_IEDxExpectedReward) & ...
        isfinite(areaResults.betaCILow_IEDxExpectedReward) & ...
        isfinite(areaResults.betaCIHigh_IEDxExpectedReward);

end

