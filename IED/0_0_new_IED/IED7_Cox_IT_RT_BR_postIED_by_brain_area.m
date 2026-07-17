% Combined area-specific participant-stratified Cox analyses

% Author: Nill

clear;
clc;
close all;

%% Paths

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED7_Cox_IT_RT_BR_postIED_by_brain_area\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

%% Shared settings

% Outcome-specific post-IED windows
% RT is relatively short, so use 500 ms.
% IT and BR unfold over inflation time, so use 1000 ms.
settings.postIEDWindowMillisecondsRT = 500;
settings.postIEDWindowMillisecondsIT = 1000;
settings.postIEDWindowMillisecondsBR = 1000;

settings.maximumRTSeconds = 20;
settings.defaultSamplingFrequencyHz = 1000;
settings.combineLeftAndRight = true;
settings.includeUnknownArea = false;

% Sparse-area safeguards
settings.minimumParticipantsWithCoverage = 2;
settings.minimumParticipantsWithIED = 1;
settings.minimumIEDsInArea = 1;
settings.minimumEndpointEventsInsideWindow = 1;
settings.minimumEndpointEventsOutsideWindow = 1;

% Participant-preserving permutation test
settings.nPermutations = 1000;
settings.randomSeed = 20260716;
settings.minimumSuccessfulPermutationFraction = 0.80;
rng(settings.randomSeed, 'twister');

% Requested panel colors
colorRT = [0.204  0.459  0.702];   % blue
colorIT = [0.847  0.333  0.153];   % orange
colorBR = [0.250  0.600  0.250];   % green
markerAlpha = 0.5;

%% Run IT, RT, and BR analyses

[ITResults, ITModels] = runLoggedBrainAreaAnalysis( ...
    inputFolderName_LFPIED, outputFolderName, "IT", settings);

[RTResults, RTModels] = runLoggedBrainAreaAnalysis( ...
    inputFolderName_LFPIED, outputFolderName, "RT", settings);

[BRResults, BRModels] = runLoggedBrainAreaAnalysis( ...
    inputFolderName_LFPIED, outputFolderName, "BR", settings);

%% Create one 1x3 forest-plot figure

combinedFigureOutputFile = fullfile(outputFolderName, ...
    'IT_RT_BR_brain_area_cox_forest_plot.pdf');

plotCombinedAreaForestPlots( ...
    ITResults, RTResults, BRResults, ...
    colorIT, colorRT, colorBR, markerAlpha, ...
    settings, combinedFigureOutputFile);

combinedResultsOutputFile = fullfile(outputFolderName, ...
    'IT_RT_BR_brain_area_cox_all_results.mat');

save(combinedResultsOutputFile, ...
    'ITResults', 'RTResults', 'BRResults', ...
    'ITModels', 'RTModels', 'BRModels', ...
    'settings', 'colorIT', 'colorRT', 'colorBR', 'markerAlpha');

fprintf('\n============================================================\n');
fprintf('All three analyses finished.\n');
fprintf('Post-IED windows: IT = %d ms, RT = %d ms, BR = %d ms.\n', ...
    settings.postIEDWindowMillisecondsIT, ...
    settings.postIEDWindowMillisecondsRT, ...
    settings.postIEDWindowMillisecondsBR);
fprintf('Combined forest plot saved: %s\n', combinedFigureOutputFile);
fprintf('Combined results saved: %s\n', combinedResultsOutputFile);

%% Local functions

function [areaResults, areaModels] = runLoggedBrainAreaAnalysis( ...
    inputFolderName_LFPIED, outputFolderName, analysisType, settings)

    analysisPrefix = analysisType;

    textOutputFile = fullfile(outputFolderName, ...
        char(analysisPrefix + "_brain_area_console_output.txt"));

    diary off;
    if exist(textOutputFile, 'file')
        delete(textOutputFile);
    end

    diary(textOutputFile);
    diaryCleanup = onCleanup(@() diary('off'));

    [areaResults, areaModels] = runBrainAreaAnalysis( ...
        inputFolderName_LFPIED, analysisType, settings);

    resultsOutputFile = fullfile(outputFolderName, ...
        char(analysisPrefix + "_brain_area_cox_results.csv"));

    modelOutputFile = fullfile(outputFolderName, ...
        char(analysisPrefix + "_brain_area_cox_models.mat"));

    writetable(areaResults, resultsOutputFile);

    save(modelOutputFile, ...
        'areaModels', 'areaResults', 'analysisType', 'settings');

    fittedRows = areaResults.status == "fitted" & ...
        isfinite(areaResults.pValue);

    fprintf('\n============================================================\n');
    fprintf('%s area-specific analysis finished.\n', analysisPrefix);
    fprintf('Fitted areas: %d\n', sum(fittedRows));
    fprintf('Permutation-significant areas (uncorrected p < 0.05): %d\n', ...
        sum(areaResults.significantUncorrected));
    fprintf('Saved: %s\n', resultsOutputFile);
    fprintf('Saved: %s\n', modelOutputFile);
    fprintf('Saved: %s\n', textOutputFile);

    diary off;
    clear diaryCleanup;

end

function [areaResults, areaModels] = runBrainAreaAnalysis( ...
    inputFolderName_LFPIED, analysisType, settings)

    analysisPrefix = analysisType;

    % Select the post-IED window for the current behavioral outcome.
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

    if isempty(fileList)
        error('No .LFPIED.mat files were found.');
    end

    participants = struct( ...
        'patientID', {}, ...
        'RTs', {}, ...
        'ITs', {}, ...
        'BankedTrials', {}, ...
        'isControl', {}, ...
        'balloonColorCode', {}, ...
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

        % Exclude all NaC/nucleus accumbens, white-matter, and lateral-ventricle channels.
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
                    length(RTs), length(isControl), length(balloonType)];
            case "IT"
                ITs = double(LFPIED.ITs(:));
                BankedTrials = NaN(size(RTs));
                vectorLengths = [ ...
                    length(RTs), length(ITs), ...
                    length(isControl), length(balloonType)];
            case "BR"
                ITs = double(LFPIED.ITs(:));
                BankedTrials = double(LFPIED.BankedTrials(:));
                vectorLengths = [ ...
                    length(RTs), length(ITs), length(BankedTrials), ...
                    length(isControl), length(balloonType)];
        end

        if isfield(LFPIED, 'nTrials') && ...
                isscalar(LFPIED.nTrials) && isfinite(LFPIED.nTrials)
            vectorLengths(end + 1) = LFPIED.nTrials;
        end

        nTrials = floor(min(vectorLengths));

        RTs = RTs(1:nTrials);
        isControl = isControl(1:nTrials);
        balloonType = balloonType(1:nTrials);
        balloonColorCode = mapBalloonColorCode(balloonType);
        ITs = ITs(1:nTrials);
        BankedTrials = BankedTrials(1:nTrials);

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
            isControl == 0 & ...
            isfinite(RTs) & ...
            RTs > 0 & ...
            RTs <= maximumRTSeconds & ...
            isfinite(durationSeconds) & ...
            durationSeconds > 0 & ...
            isfinite(balloonColorCode) & ...
            ismember(balloonColorCode, [1 2 3]);

        if analysisType == "BR"
            validTrials = validTrials & ...
                isfinite(BankedTrials) & ...
                ismember(BankedTrials, [0 1]);
        end

        samplingFrequencyHz = getSamplingFrequency( ...
            LFPIED, defaultSamplingFrequencyHz);

        newParticipant.patientID = patientID;
        newParticipant.RTs = RTs;
        newParticipant.ITs = ITs;
        newParticipant.BankedTrials = BankedTrials;
        newParticipant.isControl = isControl;
        newParticipant.balloonColorCode = balloonColorCode;
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

    end

    if isempty(participants)
        error('No valid participant data were loaded for %s.', ...
            analysisPrefix);
    end

    % Hard safety filter before any model is fitted.
    allAreas = allAreas(~isExcludedAreaLabel(allAreas));

    if isempty(allAreas)
        error(['No anatomical areas were found after excluding NaC, ' ...
            'white matter, and lateral ventricle.']);
    end

    fprintf('\n============================================================\n');
    fprintf('%s area-specific %d-ms post-IED Cox analysis\n', ...
        analysisPrefix, postIEDWindowMilliseconds);
    fprintf('============================================================\n');
    fprintf('Loaded participants: %d\n', length(participants));
    fprintf('Anatomical areas: %d\n', length(allAreas));

    %% Fit one Cox model per anatomical area

    areaResults = initializeAreaResultsTable();
    areaModels = struct();

    for aa = 1:length(allAreas)

        anatomicalArea = allAreas(aa);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Area %d/%d: %s\n', ...
            aa, length(allAreas), anatomicalArea);

        countingProcessData = emptyCountingProcessTable();

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

                trialRows = makePostIEDWindowRows( ...
                    P.patientID, ...
                    trialNumber, ...
                    P.durationSeconds(trialNumber), ...
                    areaIEDoccurrence, ...
                    P.samplingFrequencyHz, ...
                    P.balloonColorCode(trialNumber), ...
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

        betaPrimary = NaN;
        robustSEPrimary = NaN;
        robustZPrimary = NaN;
        pValueWald = NaN;
        pValuePermutation = NaN;
        nSuccessfulPermutations = 0;
        pValuePrimary = NaN; % Primary p value used for significance
        betaCILow = NaN;
        betaCIHigh = NaN;
        hazardRatio = NaN;
        hazardRatioCILow = NaN;
        hazardRatioCIHigh = NaN;
        logLikelihood = NaN;
        nModelParticipants = 0;
        referenceBalloonColor = "";
        safeAreaName = "";

        if status == "fitted"

            [countingProcessData.patientStratum, patientLevels] = ...
                findgroups(countingProcessData.patientID);

            [X, predictorNames, referenceBalloonColor] = ...
                buildPredictorMatrix(countingProcessData);

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
                status = "skipped: no area-specific exposed time";
            elseif nEventsInsideWindow < ...
                    minimumEndpointEventsInsideWindow
                status = "skipped: no endpoint events inside area-specific windows";
            elseif nEventsOutsideWindow < ...
                    minimumEndpointEventsOutsideWindow
                status = "skipped: no endpoint events outside area-specific windows";
            elseif size(X, 1) <= size(X, 2)
                status = "skipped: insufficient rows for predictor matrix";
            end

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

                betaPrimary = beta(1);
                robustSEPrimary = clusterRobustSE(1);
                robustZPrimary = clusterRobustZ(1);
                pValueWald = clusterRobustP(1);

                betaCILow = betaPrimary - 1.96 * robustSEPrimary;
                betaCIHigh = betaPrimary + 1.96 * robustSEPrimary;

                hazardRatio = exp(betaPrimary);
                hazardRatioCILow = exp(betaCILow);
                hazardRatioCIHigh = exp(betaCIHigh);

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
                    'beta', beta, ...
                    'predictorNames', predictorNames, ...
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
                        referenceBalloonColor);

            catch modelError
                status = "failed: " + string(modelError.message);
            end

        end

        % Shuffle complete IED trial assignments within each participant,
        % rebuild all time-varying rows, and refit the same Cox model.
        if status == "fitted" && isfinite(betaPrimary)
            [pValuePermutation, nSuccessfulPermutations] = ...
                runAreaIEDPermutationTest( ...
                    participants, anatomicalArea, ...
                    postIEDWindowSeconds, betaPrimary, settings);

            minimumSuccessfulPermutations = ceil( ...
                settings.nPermutations * ...
                settings.minimumSuccessfulPermutationFraction);

            if nSuccessfulPermutations >= minimumSuccessfulPermutations
                pValuePrimary = pValuePermutation;
            else
                status = "failed: too few successful permutation fits";
                pValuePermutation = NaN;
            end
        end

        if strlength(safeAreaName) > 0 && ...
                isfield(areaModels, char(safeAreaName))
            areaModels.(char(safeAreaName)).pValueWald = pValueWald;
            areaModels.(char(safeAreaName)).pValuePermutation = ...
                pValuePermutation;
            areaModels.(char(safeAreaName)).nSuccessfulPermutations = ...
                nSuccessfulPermutations;
            areaModels.(char(safeAreaName)).nRequestedPermutations = ...
                settings.nPermutations;
        end

        newRow = table( ...
            anatomicalArea, ...
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
            betaPrimary, ...
            robustSEPrimary, ...
            robustZPrimary, ...
            pValueWald, ...
            pValuePermutation, ...
            nSuccessfulPermutations, ...
            pValuePrimary, ...
            betaCILow, ...
            betaCIHigh, ...
            hazardRatio, ...
            hazardRatioCILow, ...
            hazardRatioCIHigh, ...
            false, ...
            'VariableNames', areaResults.Properties.VariableNames);

        areaResults = [areaResults; newRow]; 

        fprintf('Status: %s\n', status);
        fprintf('Coverage participants: %d\n', ...
            nParticipantsWithCoverage);
        fprintf('Participants with area IEDs: %d\n', ...
            nParticipantsWithIED);
        fprintf('Area IEDs: %d\n', nIEDsInArea);

        if status == "fitted"
            fprintf([ ...
                'log(HR) = %.6f, HR = %.6f, ' ...
                '95%% CI = [%.6f, %.6f], Wald p = %.6g, ' ...
                'permutation p = %.6g (%d/%d fits)\n'], ...
                betaPrimary, hazardRatio, hazardRatioCILow, ...
                hazardRatioCIHigh, pValueWald, pValuePermutation, ...
                nSuccessfulPermutations, settings.nPermutations);
        end

    end

    %% Final hard removal of excluded labels

    excludedResultRows = ...
        isExcludedAreaLabel(areaResults.anatomicalArea);

    if any(excludedResultRows)
        fprintf([ ...
            '\nHard-removing %d NaC/white-matter/lateral-ventricle result rows ' ...
            'before saving and plotting.\n'], ...
            sum(excludedResultRows));
        areaResults(excludedResultRows, :) = [];
    end

    %% Uncorrected permutation significance across fitted anatomical areas

    fittedRows = areaResults.status == "fitted" & ...
        isfinite(areaResults.pValue);

    if any(fittedRows)
        areaResults.significantUncorrected(fittedRows) = ...
            areaResults.pValue(fittedRows) < 0.05;
    end

    areaResults = sortrows( ...
        areaResults, {'status', 'pValue'}, {'ascend', 'ascend'});

end



function results = initializeAreaResultsTable()

    results = table( ...
        strings(0, 1), ...  % anatomicalArea
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
        zeros(0, 1), ...    % beta_logHazard
        zeros(0, 1), ...    % clusterRobustSE
        zeros(0, 1), ...    % clusterRobustZ
        zeros(0, 1), ...    % pValueWald
        zeros(0, 1), ...    % pValuePermutation
        zeros(0, 1), ...    % nSuccessfulPermutations
        zeros(0, 1), ...    % pValue (primary permutation p value)
        zeros(0, 1), ...    % betaCILow
        zeros(0, 1), ...    % betaCIHigh
        zeros(0, 1), ...    % hazardRatio
        zeros(0, 1), ...    % hazardRatioCILow
        zeros(0, 1), ...    % hazardRatioCIHigh
        false(0, 1), ...    % significantUncorrected
        'VariableNames', { ...
            'anatomicalArea', ...
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
            'beta_logHazard', ...
            'clusterRobustSE', ...
            'clusterRobustZ', ...
            'pValueWald', ...
            'pValuePermutation', ...
            'nSuccessfulPermutations', ...
            'pValue', ...
            'betaCILow', ...
            'betaCIHigh', ...
            'hazardRatio', ...
            'hazardRatioCILow', ...
            'hazardRatioCIHigh', ...
            'significantUncorrected'});

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

    % Create compact labels so separators and capitalization do not matter.
    compactLabels = regexprep(normalizedLabels, '[^a-z0-9]', '');

    % Remove white-matter labels, including WM, WM60, white_matter, etc.
    isWhiteMatter = ...
        contains(normalizedLabels, "white matter") | ...
        contains(compactLabels, "whitematter") | ...
        startsWith(compactLabels, "wm");

    % Remove lateral-ventricle labels without case sensitivity.
    % This includes forms such as:
    % Lateral Ventricle, lateral_ventricle, LATERAL-VENTRICLE, etc.
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

function countingProcessData = emptyCountingProcessTable()

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
            'samplingFrequencyHz'});

end

function trialRows = makePostIEDWindowRows( ...
    patientID, trialNumber, durationSeconds, IEDoccurrence, ...
    samplingFrequencyHz, balloonColorCode, finalEventObserved, ...
    postIEDWindowSeconds)

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
            'samplingFrequencyHz'});

end

function [X, predictorNames, referenceColorName] = ...
    buildPredictorMatrix(countingProcessData)

    X = countingProcessData.activeIEDWindowCount;
    predictorNames = "active_area_ied_count_within_selected_window";

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
            colorNames(colorCode) + "_vs_" + colorNames(referenceColor); 
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
    robustCovariance = ...
        (nClusters / (nClusters - 1)) .* robustCovariance;

    robustSE = sqrt(max(diag(robustCovariance), 0));

    zeroSERows = robustSE <= 0 | ~isfinite(robustSE);
    robustSE(zeroSERows) = stats.se(zeroSERows);

    robustZ = beta ./ robustSE;
    robustP = 2 .* normcdf(-abs(robustZ), 0, 1);

end

function [permutationPValue, nSuccessfulPermutations] = ...
    runAreaIEDPermutationTest( ...
        participants, anatomicalArea, postIEDWindowSeconds, ...
        observedBeta, settings)

    nPermutations = settings.nPermutations;
    permutationBetas = NaN(nPermutations, 1);

    fprintf('Running %d within-participant permutations for %s...\n', ...
        nPermutations, anatomicalArea);

    for permutationNumber = 1:nPermutations
        permutedData = emptyCountingProcessTable();

        for pt = 1:length(participants)
            P = participants(pt);
            areaLocalChannels = find( ...
                P.selectedAreaLabels == anatomicalArea);

            if isempty(areaLocalChannels)
                continue;
            end

            areaIEDoccurrence = filterIEDsToArea( ...
                P.IEDoccurrence, areaLocalChannels);

            % A trial's full IED pattern is moved as one unit. This retains
            % within-trial IED counts/timing and never mixes participants.
            permutedIEDoccurrence = permuteIEDTrialAssignments( ...
                areaIEDoccurrence, P.validTrials);

            validTrialNumbers = find(P.validTrials);

            for tt = 1:length(validTrialNumbers)
                trialNumber = validTrialNumbers(tt);

                trialRows = makePostIEDWindowRows( ...
                    P.patientID, ...
                    trialNumber, ...
                    P.durationSeconds(trialNumber), ...
                    permutedIEDoccurrence, ...
                    P.samplingFrequencyHz, ...
                    P.balloonColorCode(trialNumber), ...
                    P.finalEventObserved(trialNumber), ...
                    postIEDWindowSeconds);

                permutedData = [permutedData; trialRows]; 
            end
        end

        permutationBetas(permutationNumber) = ...
            fitPrimaryCoxCoefficient(permutedData, settings);

        if mod(permutationNumber, 50) == 0 || ...
                permutationNumber == nPermutations
            fprintf('  completed %d/%d permutations\n', ...
                permutationNumber, nPermutations);
        end
    end

    validPermutationBetas = permutationBetas( ...
        isfinite(permutationBetas));
    nSuccessfulPermutations = length(validPermutationBetas);

    if nSuccessfulPermutations == 0 || ~isfinite(observedBeta)
        permutationPValue = NaN;
        return;
    end

    % Phipson-Smyth plus-one correction prevents a zero empirical p value.
    permutationPValue = (1 + sum( ...
        abs(validPermutationBetas) >= abs(observedBeta))) ./ ...
        (nSuccessfulPermutations + 1);

end

function permutedIEDoccurrence = permuteIEDTrialAssignments( ...
    areaIEDoccurrence, validTrialMask)

    if isempty(areaIEDoccurrence) || size(areaIEDoccurrence, 2) < 3
        permutedIEDoccurrence = zeros(0, 3);
        return;
    end

    validTrialNumbers = find(validTrialMask(:));

    if isempty(validTrialNumbers)
        permutedIEDoccurrence = zeros(0, size(areaIEDoccurrence, 2));
        return;
    end

    sourceTrialNumbers = round(areaIEDoccurrence(:, 1));
    keepRows = ...
        isfinite(sourceTrialNumbers) & ...
        ismember(sourceTrialNumbers, validTrialNumbers) & ...
        isfinite(areaIEDoccurrence(:, 3)) & ...
        areaIEDoccurrence(:, 3) >= 1;

    permutedIEDoccurrence = areaIEDoccurrence(keepRows, :);

    if isempty(permutedIEDoccurrence)
        return;
    end

    % One-to-one reassignment across all eligible trials, including trials
    % with no IED. Thus the number of IED-free trials is also preserved.
    destinationTrialNumbers = validTrialNumbers( ...
        randperm(length(validTrialNumbers)));
    trialMap = NaN(length(validTrialMask), 1);
    trialMap(validTrialNumbers) = destinationTrialNumbers;

    sourceTrialNumbers = round(permutedIEDoccurrence(:, 1));
    permutedIEDoccurrence(:, 1) = trialMap(sourceTrialNumbers);

end

function betaPrimary = fitPrimaryCoxCoefficient( ...
    countingProcessData, settings)

    betaPrimary = NaN;

    if isempty(countingProcessData)
        return;
    end

    try
        [countingProcessData.patientStratum, ~] = ...
            findgroups(countingProcessData.patientID);

        [X, ~, ~] = buildPredictorMatrix(countingProcessData);
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

        nModelParticipants = length(unique( ...
            countingProcessData.patientID));
        nEventsInsideWindow = sum( ...
            countingProcessData.eventAtStop == 1 & ...
            countingProcessData.postIED == 1);
        nEventsOutsideWindow = sum( ...
            countingProcessData.eventAtStop == 1 & ...
            countingProcessData.postIED == 0);
        timeInsideWindowSeconds = sum( ...
            (countingProcessData.tStop - ...
             countingProcessData.tStart) .* ...
            (countingProcessData.postIED == 1), 'omitnan');

        invalidFit = ...
            nModelParticipants < ...
                settings.minimumParticipantsWithCoverage || ...
            nEventsInsideWindow < ...
                settings.minimumEndpointEventsInsideWindow || ...
            nEventsOutsideWindow < ...
                settings.minimumEndpointEventsOutsideWindow || ...
            timeInsideWindowSeconds <= 0 || ...
            size(X, 1) <= size(X, 2) || ...
            size(X, 2) < 1 || ...
            numel(unique(X(:, 1))) < 2 || ...
            rank(X) < size(X, 2);

        if invalidFit
            return;
        end

        coxOptions = statset('coxphfit');
        coxOptions.Display = 'off';
        coxOptions.MaxIter = 1000;
        coxOptions.MaxFunEvals = 5000;

        beta = coxphfit( ...
            X, T, ...
            'Censoring', censoring, ...
            'Strata', strata, ...
            'Ties', 'efron', ...
            'Baseline', 0, ...
            'Options', coxOptions);

        if ~isempty(beta) && isfinite(beta(1))
            betaPrimary = beta(1);
        end
    catch
        % Invalid/singular permutation fits are excluded from the empirical
        % null distribution and reported through nSuccessfulPermutations.
        betaPrimary = NaN;
    end

end

function plotCombinedAreaForestPlots( ...
    ITResults, RTResults, BRResults, ...
    colorIT, colorRT, colorBR, markerAlpha, ...
    settings, outputPDF)

    maximumNumberOfAreas = max([ ...
        countPlottableAreas(ITResults), ...
        countPlottableAreas(RTResults), ...
        countPlottableAreas(BRResults), ...
        1]);

    % Two-line area labels need a little more vertical space than ordinary
    % one-line tick labels. This keeps labels from overlapping in dense panels.
    figureHeight = max(700, 48 * maximumNumberOfAreas + 250);

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [80 80 2400 figureHeight]);

    layout = tiledlayout(fig, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    axIT = nexttile(layout, 1);
    plotAreaForestPanel( ...
        axIT, ITResults, colorIT, markerAlpha, ...
        "Inflation time", "IT", ...
        settings.postIEDWindowMillisecondsIT);

    axRT = nexttile(layout, 2);
    plotAreaForestPanel( ...
        axRT, RTResults, colorRT, markerAlpha, ...
        "Response time", "RT", ...
        settings.postIEDWindowMillisecondsRT);

    axBR = nexttile(layout, 3);
    plotAreaForestPanel( ...
        axBR, BRResults, colorBR, markerAlpha, ...
        "Bank rate", "BR", ...
        settings.postIEDWindowMillisecondsBR);

    exportgraphics(fig, outputPDF, 'ContentType', 'vector');
    close(fig);

end

function numberOfAreas = countPlottableAreas(areaResults)

    plotRows = getForestPlotRows(areaResults);
    numberOfAreas = sum(plotRows);

end

function plotAreaForestPanel( ...
    ax, areaResults, panelColor, markerAlpha, ...
    descriptiveTitle, abbreviation, postIEDWindowMilliseconds)

    plotRows = getForestPlotRows(areaResults);
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
            'No areas with uncorrected permutation p < 0.05', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 11);
    else
        plotData = sortrows( ...
            plotData, 'beta_logHazard', 'ascend');
        nAreas = height(plotData);

        for ii = 1:nAreas
            plot(ax, ...
                [plotData.betaCILow(ii), plotData.betaCIHigh(ii)], ...
                [ii, ii], ...
                '-', ...
                'Color', panelColor, ...
                'LineWidth', 1.5);
        end

        scatter(ax, ...
            plotData.beta_logHazard, ...
            (1:nAreas)', ...
            65, ...
            panelColor, ...
            'filled', ...
            'MarkerFaceAlpha', markerAlpha, ...
            'MarkerEdgeColor', 'none');

        ax.YTick = 1:nAreas;

        % Use two lines per area for a cleaner display. n IED is the number
        % of participants with >=1 valid IED in this area; n covered is the
        % number of participants with electrode coverage in this area.
        areaTickLabels = strings(nAreas, 1);
        for ii = 1:nAreas
            areaTickLabels(ii) = sprintf( ...
                '%s\nIED participants: %d   |   covered: %d', ...
                char(plotData.anatomicalArea(ii)), ...
                plotData.nParticipantsWithIED(ii), ...
                plotData.nParticipantsWithCoverage(ii));
        end
        ax.YTickLabel = areaTickLabels;
        ax.YLim = [0.4, nAreas + 0.6];
    end

    ax.TickLabelInterpreter = 'none';
    ax.TickDir = 'out';
    % Set the axes font directly for compatibility across MATLAB releases.
    ax.FontSize = 9;

    xlabel(ax, ...
        'Log hazard ratio per additional active area-specific IED');

    title(ax, {char(descriptiveTitle), ...
        sprintf('%s - %d-ms post-IED window', ...
            char(abbreviation), postIEDWindowMilliseconds), ...
        'IED participants = participants with >=1 valid area IED'}, ...
        'FontSize', 12, ...
        'FontWeight', 'bold');

    box(ax, 'off');
    grid(ax, 'off');

end

function plotRows = getForestPlotRows(areaResults)

    plotRows = ...
        ~isExcludedAreaLabel(areaResults.anatomicalArea) & ...
        areaResults.status == "fitted" & ...
        areaResults.significantUncorrected & ...
        isfinite(areaResults.beta_logHazard) & ...
        isfinite(areaResults.betaCILow) & ...
        isfinite(areaResults.betaCIHigh);

end
