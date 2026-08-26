% Left-versus-right Cox analyses for brain areas that were significant in
% the previous hemisphere-combined analyses.
%
% This script reads:
%   1) IT_RT_BR_brain_area_cox_all_results.mat
%   2) IT_RT_BR_mechanistic_IED_occurrence_x_expected_reward_all_results.mat
%
% It then:
%   - selects only areas whose saved significance flag is true;
%   - uses anatomicalLocs_wHemisphere (for example, SFG_L and SFG_R);
%   - restricts the primary L/R comparison to participants with bilateral
%     electrode coverage in the selected area;
%   - fits L and R simultaneously in one participant-stratified Cox model;
%   - tests beta_L - beta_R = 0 using the cluster-robust covariance;
%   - reports HR_L, HR_R, HR_L/HR_R, channel counts, and IED counts;
%   - creates separate 3-by-1 figures for the ordinary IED-occurrence model
%     and the IED-occurrence-by-expected-reward mechanistic model.
%
% For BR/PR, one stacked competing-transition model is used. BR is the bank
% transition and PR is the pop transition.
%
% Author: Nill

clear;
clc;
close all;

%% Paths

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED_new_area_labels\IED1_find_number_of_IEDs\';

inputFolderName_modeling = ...
    'D:\Nill\data\BART\0_0_new_IED_new_area_labels\context_modeling\param_recovery_1_modeling\';

previousOrdinaryResultsFile = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED7_Cox_IT_RT_BR_postIED_by_brain_area\IT_RT_BR_brain_area_cox_all_results.mat';

previousMechanisticResultsFile = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED8_Cox_expected_reward_by_brain_area\IT_RT_BR_mechanistic_IED_occurrence_x_expected_reward_all_results.mat';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED11_Cox_significant_brain_areas_LR\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

%% Settings

settings.postIEDWindowMillisecondsRT = 500;
settings.postIEDWindowMillisecondsIT = 1000;
settings.postIEDWindowMillisecondsBR = 1000;
settings.maximumRTSeconds = 20;
settings.defaultSamplingFrequencyHz = 1000;
settings.useOnlyNonControlTrials = true;

% Primary L/R comparison: use only participants with both L and R coverage
% for the same area. This prevents different participant samples from being
% mistaken for a hemispheric effect.
settings.requireBilateralCoverage = true;
settings.minimumBilateralParticipants = 10;
settings.minimumParticipantsWithIEDPerSide = 1;
settings.minimumIEDsPerSide = 1;

% Significance level for the new L-versus-R contrast.
settings.alpha = 0.05;

% Figure settings.
settings.figurePosition = [100 60 900 1150];
settings.leftColor = [0.204 0.459 0.702];
settings.rightColor = [0.847 0.333 0.153];
settings.markerSize = 6;

%% Read the two saved hemisphere-combined result files

if ~exist(previousOrdinaryResultsFile, 'file')
    error('Previous ordinary result file was not found:\n%s', ...
        previousOrdinaryResultsFile);
end

if ~exist(previousMechanisticResultsFile, 'file')
    error('Previous mechanistic result file was not found:\n%s', ...
        previousMechanisticResultsFile);
end

previousOrdinary = load(previousOrdinaryResultsFile, ...
    'ITResults', 'RTResults', 'BRResults');

previousMechanistic = load(previousMechanisticResultsFile, ...
    'ITResults', 'RTResults', 'BRResults');

ordinarySelectedAreas = selectSignificantAreas( ...
    previousOrdinary, "ordinary");

mechanisticSelectedAreas = selectSignificantAreas( ...
    previousMechanistic, "mechanistic");

printSelectedAreas(ordinarySelectedAreas, "Ordinary IED occurrence");
printSelectedAreas(mechanisticSelectedAreas, ...
    "IED occurrence x expected reward");

%% Run the ordinary L/R analyses

[ordinaryLRResults, ordinaryLRModels] = runAnalysisFamily( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    ordinarySelectedAreas, "ordinary", settings);

ordinaryCSV = fullfile(outputFolderName, ...
    'ordinary_significant_areas_left_right_results.csv');
ordinaryMAT = fullfile(outputFolderName, ...
    'ordinary_significant_areas_left_right_all_results.mat');
ordinaryPDF = fullfile(outputFolderName, ...
    'ordinary_significant_areas_left_right_3x1.pdf');

writetable(ordinaryLRResults, ordinaryCSV);
save(ordinaryMAT, 'ordinaryLRResults', 'ordinaryLRModels', ...
    'ordinarySelectedAreas', 'settings');
plotLRForestFigure(ordinaryLRResults, ...
    "IED occurrence", settings, ordinaryPDF, "confidenceIntervals");

%% Run the mechanistic L/R analyses

[mechanisticLRResults, mechanisticLRModels] = runAnalysisFamily( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    mechanisticSelectedAreas, "mechanistic", settings);

mechanisticCSV = fullfile(outputFolderName, ...
    'mechanistic_significant_areas_left_right_results.csv');
mechanisticMAT = fullfile(outputFolderName, ...
    'mechanistic_significant_areas_left_right_all_results.mat');
mechanisticPDF = fullfile(outputFolderName, ...
    'mechanistic_significant_areas_left_right_3x1.pdf');

writetable(mechanisticLRResults, mechanisticCSV);
save(mechanisticMAT, 'mechanisticLRResults', ...
    'mechanisticLRModels', 'mechanisticSelectedAreas', 'settings');
plotLRForestFigure(mechanisticLRResults, ...
    "IED occurrence x expected reward", settings, mechanisticPDF, ...
    "pointEstimates");

fprintf('\n============================================================\n');
fprintf('All left-versus-right analyses finished.\n');
fprintf('Ordinary results:    %s\n', ordinaryCSV);
fprintf('Ordinary figure:     %s\n', ordinaryPDF);
fprintf('Mechanistic results: %s\n', mechanisticCSV);
fprintf('Mechanistic figure:  %s\n', mechanisticPDF);
fprintf(['The pValueLeftVsRight column tests HR_L/HR_R = 1 ' ...
    '(equivalently beta_L - beta_R = 0).\n']);

%% Local functions

function selectedAreas = selectSignificantAreas(savedData, familyType)

    outcomes = ["IT", "RT", "BR"];
    selectedAreas = struct('IT', strings(0, 1), ...
        'RT', strings(0, 1), 'BR', strings(0, 1));

    for oo = 1:numel(outcomes)
        outcome = outcomes(oo);
        fieldName = char(outcome + "Results");

        if ~isfield(savedData, fieldName)
            error('Saved file does not contain %s.', fieldName);
        end

        T = savedData.(fieldName);

        if ~istable(T) || ~ismember('anatomicalArea', ...
                T.Properties.VariableNames)
            error('%s must be a table containing anatomicalArea.', ...
                fieldName);
        end

        if familyType == "ordinary"
            if outcome == "BR"
                candidateFlags = [ ...
                    "significantUncorrected", ...
                    "significantPopUncorrected", ...
                    "significantBankVsPopUncorrected"];
            else
                candidateFlags = "significantUncorrected";
            end
        else
            if outcome == "BR"
                candidateFlags = [ ...
                    "significantUncorrected_IEDxExpectedReward", ...
                    "significantPopUncorrected_IEDxExpectedReward", ...
                    "significantBankVsPopUncorrected_IEDxExpectedReward"];
            else
                candidateFlags = ...
                    "significantUncorrected_IEDxExpectedReward";
            end
        end

        selected = false(height(T), 1);
        foundFlag = false;

        for ff = 1:numel(candidateFlags)
            flagName = candidateFlags(ff);
            if ismember(char(flagName), T.Properties.VariableNames)
                flagValues = T.(char(flagName));
                selected = selected | logical(flagValues);
                foundFlag = true;
            end
        end

        if ~foundFlag
            error('No requested significance flag was found in %s.', ...
                fieldName);
        end

        areas = strip(string(T.anatomicalArea(selected)));
        areas = areas(strlength(areas) > 0 & ~ismissing(areas));
        selectedAreas.(char(outcome)) = unique(areas, 'stable');
    end
end

function printSelectedAreas(selectedAreas, familyLabel)

    fprintf('\n============================================================\n');
    fprintf('%s: significant saved areas\n', familyLabel);
    fprintf('============================================================\n');

    outcomes = ["IT", "RT", "BR"];
    for oo = 1:numel(outcomes)
        outcome = outcomes(oo);
        areas = selectedAreas.(char(outcome));
        fprintf('%s (%d): %s\n', outcome, numel(areas), ...
            strjoin(areas, ', '));
    end
end

function [allResults, allModels] = runAnalysisFamily( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    selectedAreas, familyType, settings)

    outcomes = ["IT", "RT", "BR"];
    allResults = initializeLRResultsTable();
    allModels = struct();

    for oo = 1:numel(outcomes)
        outcome = outcomes(oo);
        areas = selectedAreas.(char(outcome));

        if isempty(areas)
            fprintf('\nNo previously significant areas for %s (%s).\n', ...
                outcome, familyType);
            continue;
        end

        participants = loadParticipantsForOutcome( ...
            inputFolderName_LFPIED, inputFolderName_modeling, ...
            outcome, familyType, settings);

        for aa = 1:numel(areas)
            baseArea = areas(aa);
            [thisResults, thisModel] = fitOneBilateralArea( ...
                participants, baseArea, outcome, familyType, settings);

            allResults = [allResults; thisResults]; %#ok<AGROW>

            modelName = matlab.lang.makeValidName(char( ...
                familyType + "_" + outcome + "_" + baseArea));
            allModels.(modelName) = thisModel;
        end
    end

    % Correct the new family of direct L-versus-R tests. The uncorrected
    % Wald p-value is retained, and the BH-FDR q-value is added separately.
    allResults.qValueLeftVsRight = NaN(height(allResults), 1);
    allResults.significantLeftVsRightFDR = false(height(allResults), 1);
    validP = allResults.status == "fitted" & ...
        isfinite(allResults.pValueLeftVsRight);
    if any(validP)
        allResults.qValueLeftVsRight(validP) = ...
            benjaminiHochbergAdjustedP( ...
            allResults.pValueLeftVsRight(validP));
        allResults.significantLeftVsRightFDR(validP) = ...
            allResults.qValueLeftVsRight(validP) < settings.alpha;
    end
end

function participants = loadParticipantsForOutcome( ...
    inputFolderName_LFPIED, inputFolderName_modeling, ...
    outcome, familyType, settings)

    fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));
    if isempty(fileList)
        error('No .LFPIED.mat files were found.');
    end

    if familyType == "mechanistic"
        modelFileList = dir(fullfile( ...
            inputFolderName_modeling, '*TDdataParamRecovery.mat'));
        if isempty(modelFileList)
            error('No *TDdataParamRecovery.mat files were found.');
        end
    else
        modelFileList = [];
    end

    participants = struct( ...
        'patientID', {}, 'RTs', {}, 'ITs', {}, ...
        'BankedTrials', {}, 'balloonColorCode', {}, ...
        'expectedReward', {}, 'validTrials', {}, ...
        'durationSeconds', {}, 'finalEventObserved', {}, ...
        'IEDoccurrence', {}, 'samplingFrequencyHz', {}, ...
        'areaKeys', {}, 'hemisphere', {});

    fprintf('\nLoading %s participant data for %s...\n', ...
        familyType, outcome);

    for pt = 1:numel(fileList)
        fileName = fileList(pt).name;
        fileNameParts = strsplit(fileName, '.');
        patientID = string(fileNameParts{1});
        loadedData = load(fullfile(inputFolderName_LFPIED, fileName));

        if ~isfield(loadedData, 'LFPIED')
            fprintf('Skipped %s: LFPIED was missing.\n', patientID);
            continue;
        end

        LFPIED = loadedData.LFPIED;
        requiredFields = getRequiredLFPIEDFields(outcome);
        requiredFields{end + 1} = 'anatomicalLocs_wHemisphere';

        missingFields = requiredFields(~cellfun( ...
            @(f) isfield(LFPIED, f), requiredFields));

        if ~isempty(missingFields)
            fprintf('Skipped %s: missing %s.\n', patientID, ...
                strjoin(string(missingFields), ', '));
            continue;
        end

        selectedChans = round(double(LFPIED.selectedChans(:)));
        if isempty(selectedChans)
            continue;
        end

        rawLabels = getSelectedChannelLabels( ...
            LFPIED.anatomicalLocs_wHemisphere, selectedChans, ...
            numel(selectedChans));
        [areaKeys, hemisphere] = parseHemisphereLabels(rawLabels);

        RTs = double(LFPIED.RTs(:));
        isControl = double(LFPIED.isControl(:));
        balloonType = double(LFPIED.balloonType(:));

        switch outcome
            case "RT"
                ITs = NaN(size(RTs));
                BankedTrials = NaN(size(RTs));
                baseLengths = [numel(RTs), numel(isControl), ...
                    numel(balloonType)];
                IEDoccurrence = LFPIED.IED_occurance_RT;
            case "IT"
                ITs = double(LFPIED.ITs(:));
                BankedTrials = NaN(size(RTs));
                baseLengths = [numel(RTs), numel(ITs), ...
                    numel(isControl), numel(balloonType)];
                IEDoccurrence = LFPIED.IED_occurance_IT;
            case "BR"
                ITs = double(LFPIED.ITs(:));
                BankedTrials = double(LFPIED.BankedTrials(:));
                baseLengths = [numel(RTs), numel(ITs), ...
                    numel(BankedTrials), numel(isControl), ...
                    numel(balloonType)];
                IEDoccurrence = LFPIED.IED_occurance_IT;
            otherwise
                error('Unknown outcome: %s', outcome);
        end

        if familyType == "mechanistic"
            modelFileName = findParticipantModelFile( ...
                modelFileList, patientID);
            if strlength(modelFileName) == 0
                fprintf('Skipped %s: no unique modeling file.\n', ...
                    patientID);
                continue;
            end

            modelData = load(fullfile(inputFolderName_modeling, ...
                char(modelFileName)));
            if ~isfield(modelData, 'TDdataParamRecovery')
                continue;
            end
            TD = modelData.TDdataParamRecovery;
            if ~all(isfield(TD, ...
                    {'bestApIdx', 'bestAnIdx', 'expectedReward'}))
                continue;
            end
            expectedReward = squeeze(TD.expectedReward( ...
                TD.bestApIdx, TD.bestAnIdx, :));
            expectedReward = double(expectedReward(:));
            baseLengths(end + 1) = numel(expectedReward);
        else
            expectedReward = NaN(size(RTs));
        end

        if isfield(LFPIED, 'nTrials') && ...
                isscalar(LFPIED.nTrials) && isfinite(LFPIED.nTrials)
            baseLengths(end + 1) = double(LFPIED.nTrials);
        end

        nTrials = floor(min(baseLengths));
        if nTrials < 1
            continue;
        end

        RTs = RTs(1:nTrials);
        ITs = ITs(1:nTrials);
        BankedTrials = BankedTrials(1:nTrials);
        isControl = isControl(1:nTrials);
        balloonColorCode = mapBalloonColorCode( ...
            balloonType(1:nTrials));
        expectedReward = expectedReward(1:nTrials);

        if outcome == "RT"
            durationSeconds = RTs;
            finalEventObserved = true(nTrials, 1);
        elseif outcome == "IT"
            durationSeconds = ITs;
            finalEventObserved = true(nTrials, 1);
        else
            durationSeconds = ITs;
            finalEventObserved = BankedTrials == 1;
        end

        validTrials = isfinite(RTs) & RTs > 0 & ...
            RTs <= settings.maximumRTSeconds & ...
            isfinite(durationSeconds) & durationSeconds > 0 & ...
            isfinite(balloonColorCode) & ...
            ismember(balloonColorCode, [1 2 3]);

        if settings.useOnlyNonControlTrials
            validTrials = validTrials & isControl == 0;
        end
        if outcome == "BR"
            validTrials = validTrials & isfinite(BankedTrials) & ...
                ismember(BankedTrials, [0 1]);
        end
        if familyType == "mechanistic"
            validTrials = validTrials & isfinite(expectedReward);
        end
        if ~any(validTrials)
            continue;
        end

        newParticipant.patientID = patientID;
        newParticipant.RTs = RTs;
        newParticipant.ITs = ITs;
        newParticipant.BankedTrials = BankedTrials;
        newParticipant.balloonColorCode = balloonColorCode;
        newParticipant.expectedReward = expectedReward;
        newParticipant.validTrials = validTrials;
        newParticipant.durationSeconds = durationSeconds;
        newParticipant.finalEventObserved = finalEventObserved;
        newParticipant.IEDoccurrence = IEDoccurrence;
        newParticipant.samplingFrequencyHz = getSamplingFrequency( ...
            LFPIED, settings.defaultSamplingFrequencyHz);
        newParticipant.areaKeys = areaKeys;
        newParticipant.hemisphere = hemisphere;
        participants(end + 1) = newParticipant; %#ok<AGROW>
    end

    if isempty(participants)
        error('No valid participants were loaded for %s (%s).', ...
            outcome, familyType);
    end
end

function requiredFields = getRequiredLFPIEDFields(outcome)
    switch outcome
        case "RT"
            requiredFields = {'selectedChans', 'RTs', 'isControl', ...
                'balloonType', 'IED_occurance_RT'};
        case "IT"
            requiredFields = {'selectedChans', 'RTs', 'ITs', ...
                'isControl', 'balloonType', 'IED_occurance_IT'};
        case "BR"
            requiredFields = {'selectedChans', 'RTs', 'ITs', ...
                'BankedTrials', 'isControl', 'balloonType', ...
                'IED_occurance_IT'};
        otherwise
            error('Unknown outcome: %s', outcome);
    end
end

function [resultRows, savedModel] = fitOneBilateralArea( ...
    participants, baseArea, outcome, familyType, settings)

    baseAreaKey = normalizeAreaKey(baseArea);
    postIEDWindowSeconds = getPostIEDWindowSeconds(outcome, settings);

    countingData = emptyBilateralCountingTable();
    nCoverageLeftAll = 0;
    nCoverageRightAll = 0;
    nBilateralParticipants = 0;
    nChannelsLeft = 0;
    nChannelsRight = 0;
    nParticipantsWithIEDLeft = 0;
    nParticipantsWithIEDRight = 0;
    nIEDsLeft = 0;
    nIEDsRight = 0;
    nValidTrials = 0;

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s | %s | %s\n', familyType, outcome, baseArea);

    for pt = 1:numel(participants)
        P = participants(pt);
        leftChannels = find(P.areaKeys == baseAreaKey & ...
            P.hemisphere == "L");
        rightChannels = find(P.areaKeys == baseAreaKey & ...
            P.hemisphere == "R");

        hasLeft = ~isempty(leftChannels);
        hasRight = ~isempty(rightChannels);
        nCoverageLeftAll = nCoverageLeftAll + hasLeft;
        nCoverageRightAll = nCoverageRightAll + hasRight;

        if settings.requireBilateralCoverage && ~(hasLeft && hasRight)
            continue;
        elseif ~settings.requireBilateralCoverage && ~(hasLeft || hasRight)
            continue;
        end

        nBilateralParticipants = nBilateralParticipants + ...
            double(hasLeft && hasRight);
        nChannelsLeft = nChannelsLeft + numel(leftChannels);
        nChannelsRight = nChannelsRight + numel(rightChannels);

        leftIEDs = filterIEDsToChannels(P.IEDoccurrence, leftChannels);
        rightIEDs = filterIEDsToChannels(P.IEDoccurrence, rightChannels);
        validLeftIEDRows = getValidIEDRows(leftIEDs, P.validTrials);
        validRightIEDRows = getValidIEDRows(rightIEDs, P.validTrials);
        thisLeftCount = sum(validLeftIEDRows);
        thisRightCount = sum(validRightIEDRows);
        nIEDsLeft = nIEDsLeft + thisLeftCount;
        nIEDsRight = nIEDsRight + thisRightCount;
        nParticipantsWithIEDLeft = nParticipantsWithIEDLeft + ...
            double(thisLeftCount > 0);
        nParticipantsWithIEDRight = nParticipantsWithIEDRight + ...
            double(thisRightCount > 0);

        validTrialNumbers = find(P.validTrials);
        nValidTrials = nValidTrials + numel(validTrialNumbers);

        for tt = 1:numel(validTrialNumbers)
            trialNumber = validTrialNumbers(tt);
            expectedReward = P.expectedReward(trialNumber);

            if outcome == "BR"
                bankRows = makeBilateralRows( ...
                    P.patientID, trialNumber, ...
                    P.durationSeconds(trialNumber), ...
                    leftIEDs, rightIEDs, P.samplingFrequencyHz, ...
                    P.balloonColorCode(trialNumber), expectedReward, ...
                    P.BankedTrials(trialNumber) == 1, ...
                    postIEDWindowSeconds, "bank");
                popRows = makeBilateralRows( ...
                    P.patientID, trialNumber, ...
                    P.durationSeconds(trialNumber), ...
                    leftIEDs, rightIEDs, P.samplingFrequencyHz, ...
                    P.balloonColorCode(trialNumber), expectedReward, ...
                    P.BankedTrials(trialNumber) == 0, ...
                    postIEDWindowSeconds, "pop");
                trialRows = [bankRows; popRows];
            else
                trialRows = makeBilateralRows( ...
                    P.patientID, trialNumber, ...
                    P.durationSeconds(trialNumber), ...
                    leftIEDs, rightIEDs, P.samplingFrequencyHz, ...
                    P.balloonColorCode(trialNumber), expectedReward, ...
                    P.finalEventObserved(trialNumber), ...
                    postIEDWindowSeconds, "single");
            end
            countingData = [countingData; trialRows]; %#ok<AGROW>
        end
    end

    status = "fitted";
    if settings.requireBilateralCoverage && ...
            nBilateralParticipants < settings.minimumBilateralParticipants
        status = "skipped: insufficient bilateral participants";
    elseif nParticipantsWithIEDLeft < ...
            settings.minimumParticipantsWithIEDPerSide
        status = "skipped: insufficient participants with left IEDs";
    elseif nParticipantsWithIEDRight < ...
            settings.minimumParticipantsWithIEDPerSide
        status = "skipped: insufficient participants with right IEDs";
    elseif nIEDsLeft < settings.minimumIEDsPerSide
        status = "skipped: insufficient left IEDs";
    elseif nIEDsRight < settings.minimumIEDsPerSide
        status = "skipped: insufficient right IEDs";
    elseif isempty(countingData)
        status = "skipped: no counting-process rows";
    end

    summaries = struct('IT', emptySummaryPair(), ...
        'RT', emptySummaryPair(), 'BR', emptySummaryPair(), ...
        'PR', emptySummaryPair());
    savedModel = struct();

    if status == "fitted"
        try
            [countingData.patientStratum, patientLevels] = ...
                findgroups(countingData.patientID);
            if outcome == "BR"
                strata = findgroups(countingData.patientID, ...
                    countingData.transitionType);
            else
                strata = countingData.patientStratum;
            end

            [X, predictorNames, referenceColor] = ...
                buildBilateralPredictorMatrix( ...
                countingData, outcome, familyType);
            T = [countingData.tStart, countingData.tStop];
            censoring = logical(countingData.censored);

            validRows = all(isfinite(X), 2) & all(isfinite(T), 2) & ...
                T(:, 1) >= 0 & T(:, 2) > T(:, 1) & ...
                isfinite(strata);
            X = X(validRows, :);
            T = T(validRows, :);
            censoring = censoring(validRows);
            strata = strata(validRows);
            countingData = countingData(validRows, :);

            if size(X, 1) <= size(X, 2)
                error('Insufficient rows for the predictor matrix.');
            end

            coxOptions = statset('coxphfit');
            coxOptions.Display = 'off';
            coxOptions.MaxIter = 1000;
            coxOptions.MaxFunEvals = 5000;

            [beta, logLikelihood, baselineCumulativeHazard, stats] = ...
                coxphfit(X, T, 'Censoring', censoring, ...
                'Strata', strata, 'Ties', 'efron', 'Baseline', 0, ...
                'Options', coxOptions);

            [robustCovariance, robustSE, robustZ, robustP] = ...
                computeClusterRobustInference(stats, beta, ...
                countingData.patientStratum);

            if outcome == "BR"
                summaries.BR = getEndpointSummaries(beta, ...
                    robustCovariance, robustSE, robustZ, robustP, ...
                    predictorNames, "BR", familyType);
                summaries.PR = getEndpointSummaries(beta, ...
                    robustCovariance, robustSE, robustZ, robustP, ...
                    predictorNames, "PR", familyType);
            else
                summaries.(char(outcome)) = getEndpointSummaries( ...
                    beta, robustCovariance, robustSE, robustZ, ...
                    robustP, predictorNames, outcome, familyType);
            end

            savedModel = struct( ...
                'baseArea', baseArea, ...
                'outcome', outcome, ...
                'familyType', familyType, ...
                'beta', beta, ...
                'predictorNames', predictorNames, ...
                'logLikelihood', logLikelihood, ...
                'baselineCumulativeHazard', ...
                    baselineCumulativeHazard, ...
                'stats', stats, ...
                'clusterRobustCovariance', robustCovariance, ...
                'clusterRobustSE', robustSE, ...
                'clusterRobustZ', robustZ, ...
                'clusterRobustP', robustP, ...
                'patientLevels', patientLevels, ...
                'referenceBalloonColor', referenceColor, ...
                'nChannelsLeft', nChannelsLeft, ...
                'nChannelsRight', nChannelsRight, ...
                'nIEDsLeft', nIEDsLeft, ...
                'nIEDsRight', nIEDsRight);

        catch modelError
            status = "failed: " + string(modelError.message);
            savedModel = struct('baseArea', baseArea, ...
                'outcome', outcome, 'familyType', familyType, ...
                'errorMessage', string(modelError.message));
        end
    end

    if outcome == "BR"
        endpointNames = ["BR", "PR"];
    else
        endpointNames = outcome;
    end

    resultRows = initializeLRResultsTable();
    for ee = 1:numel(endpointNames)
        endpoint = endpointNames(ee);
        S = summaries.(char(endpoint));
        newRow = table( ...
            familyType, outcome, endpoint, baseArea, status, ...
            nCoverageLeftAll, nCoverageRightAll, ...
            nBilateralParticipants, nChannelsLeft, nChannelsRight, ...
            nParticipantsWithIEDLeft, nParticipantsWithIEDRight, ...
            nIEDsLeft, nIEDsRight, nValidTrials, ...
            height(countingData), ...
            S.left.beta, S.left.se, S.left.p, ...
            S.left.hr, S.left.hrLow, S.left.hrHigh, ...
            S.right.beta, S.right.se, S.right.p, ...
            S.right.hr, S.right.hrLow, S.right.hrHigh, ...
            S.difference.beta, S.difference.se, ...
            S.difference.z, S.difference.p, ...
            S.difference.hr, S.difference.hrLow, ...
            S.difference.hrHigh, ...
            S.left.p < settings.alpha, ...
            S.right.p < settings.alpha, ...
            S.difference.p < settings.alpha, ...
            'VariableNames', resultRows.Properties.VariableNames);
        resultRows = [resultRows; newRow]; %#ok<AGROW>
    end

    fprintf('Status: %s\n', status);
    fprintf(['Coverage L/R/all: %d/%d; bilateral: %d; ' ...
        'channels L/R: %d/%d; IEDs L/R: %d/%d\n'], ...
        nCoverageLeftAll, nCoverageRightAll, ...
        nBilateralParticipants, nChannelsLeft, nChannelsRight, ...
        nIEDsLeft, nIEDsRight);
    if status == "fitted"
        for ee = 1:numel(endpointNames)
            endpoint = endpointNames(ee);
            S = summaries.(char(endpoint));
            fprintf([ ...
                '%s: HR_L = %.4f [%.4f, %.4f], p = %.4g; ' ...
                'HR_R = %.4f [%.4f, %.4f], p = %.4g; ' ...
                'HR_L/HR_R = %.4f, L-vs-R p = %.4g\n'], ...
                endpoint, S.left.hr, S.left.hrLow, S.left.hrHigh, ...
                S.left.p, S.right.hr, S.right.hrLow, ...
                S.right.hrHigh, S.right.p, S.difference.hr, ...
                S.difference.p);
        end
    end
end

function results = initializeLRResultsTable()
    variableNames = { ...
            'analysisFamily', 'sourceOutcome', 'endpoint', ...
            'anatomicalArea', 'status', ...
            'nParticipantsCoverageLeftAll', ...
            'nParticipantsCoverageRightAll', ...
            'nBilateralParticipants', ...
            'nChannelsLeft', 'nChannelsRight', ...
            'nParticipantsWithIEDLeft', ...
            'nParticipantsWithIEDRight', ...
            'nIEDsLeft', 'nIEDsRight', 'nValidTrials', ...
            'nCountingRows', ...
            'betaLeft', 'clusterRobustSELeft', 'pValueLeft', ...
            'hazardRatioLeft', 'hazardRatioLeftCILow', ...
            'hazardRatioLeftCIHigh', ...
            'betaRight', 'clusterRobustSERight', 'pValueRight', ...
            'hazardRatioRight', 'hazardRatioRightCILow', ...
            'hazardRatioRightCIHigh', ...
            'betaLeftMinusRight', ...
            'clusterRobustSELeftMinusRight', ...
            'zLeftVsRight', 'pValueLeftVsRight', ...
            'hazardRatioLeftOverRight', ...
            'hazardRatioLeftOverRightCILow', ...
            'hazardRatioLeftOverRightCIHigh', ...
            'significantLeft', 'significantRight', ...
            'significantLeftVsRight'};

    variableTypes = [repmat("string", 1, 5), ...
        repmat("double", 1, 30), repmat("logical", 1, 3)];

    results = table('Size', [0 numel(variableNames)], ...
        'VariableTypes', cellstr(variableTypes), ...
        'VariableNames', variableNames);
end

function seconds = getPostIEDWindowSeconds(outcome, settings)
    switch outcome
        case "RT"
            milliseconds = settings.postIEDWindowMillisecondsRT;
        case "IT"
            milliseconds = settings.postIEDWindowMillisecondsIT;
        case "BR"
            milliseconds = settings.postIEDWindowMillisecondsBR;
        otherwise
            error('Unknown outcome: %s', outcome);
    end
    seconds = milliseconds / 1000;
end

function labels = getSelectedChannelLabels( ...
    anatomicalLocs, selectedChans, nSelectedChannels)

    labels = convertLabelsToString(anatomicalLocs);
    labels = labels(:);
    if numel(labels) >= max(selectedChans)
        labels = labels(selectedChans);
    elseif numel(labels) == nSelectedChannels
        % Labels were already restricted to selected channels.
    else
        error(['Cannot map anatomicalLocs_wHemisphere to selected ' ...
            'channels. Labels = %d, selected channels = %d, ' ...
            'maximum selected channel = %d.'], numel(labels), ...
            nSelectedChannels, max(selectedChans));
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
    labels = strip(labels(:));
end

function [areaKeys, hemisphere] = parseHemisphereLabels(labels)
    labels = strip(string(labels(:)));
    areaKeys = strings(size(labels));
    hemisphere = strings(size(labels));

    for ii = 1:numel(labels)
        label = labels(ii);
        if ismissing(label) || strlength(label) == 0
            continue;
        end

        textLabel = char(label);
        suffix = regexp(textLabel, ...
            '^(.*?)[_\-\s]+(L|R|LH|RH|Left|Right)$', ...
            'tokens', 'once', 'ignorecase');
        prefix = regexp(textLabel, ...
            '^(L|R|LH|RH|Left|Right)[_\-\s]+(.*)$', ...
            'tokens', 'once', 'ignorecase');

        if ~isempty(suffix)
            baseText = suffix{1};
            sideText = suffix{2};
        elseif ~isempty(prefix)
            sideText = prefix{1};
            baseText = prefix{2};
        else
            baseText = textLabel;
            sideText = '';
        end

        areaKeys(ii) = normalizeAreaKey(baseText);
        sideLower = lower(string(sideText));
        if ismember(sideLower, ["l", "lh", "left"])
            hemisphere(ii) = "L";
        elseif ismember(sideLower, ["r", "rh", "right"])
            hemisphere(ii) = "R";
        end
    end
end

function key = normalizeAreaKey(areaLabel)
    key = lower(strip(string(areaLabel)));
    key = regexprep(key, ...
        '^(left|right|lh|rh|l|r)[_\-\s]+', '');
    key = regexprep(key, ...
        '[_\-\s]+(left|right|lh|rh|l|r)$', '');
    key = regexprep(key, '[^a-z0-9]', '');
end

function samplingFrequencyHz = getSamplingFrequency( ...
    LFPIED, defaultSamplingFrequencyHz)
    if isfield(LFPIED, 'Fs') && isscalar(LFPIED.Fs) && ...
            isfinite(LFPIED.Fs) && LFPIED.Fs > 0
        samplingFrequencyHz = double(LFPIED.Fs);
    else
        samplingFrequencyHz = defaultSamplingFrequencyHz;
    end
end

function balloonColorCode = mapBalloonColorCode(balloonType)
    balloonType = double(balloonType(:));
    balloonColorCode = NaN(size(balloonType));
    validRows = isfinite(balloonType) & ...
        ismember(round(balloonType), [1 2 3 11 12 13]);
    balloonColorCode(validRows) = ...
        mod(round(balloonType(validRows)) - 1, 10) + 1;
end

function modelFileName = findParticipantModelFile(fileList, patientID)
    modelFileName = "";
    names = string({fileList.name});
    exactPattern = "^" + regexptranslate('escape', patientID) + ...
        ".*TDdataParamRecovery\\.mat$";
    matches = ~cellfun(@isempty, regexp(cellstr(names), ...
        char(exactPattern), 'once', 'ignorecase'));
    if sum(matches) == 1
        modelFileName = names(find(matches, 1));
        return;
    end

    containsMatches = contains(lower(names), lower(patientID));
    if sum(containsMatches) == 1
        modelFileName = names(find(containsMatches, 1));
    end
end

function selectedIEDs = filterIEDsToChannels( ...
    IEDoccurrence, localChannels)
    if isempty(IEDoccurrence) || size(IEDoccurrence, 2) < 3 || ...
            isempty(localChannels)
        selectedIEDs = zeros(0, 3);
        return;
    end
    localChannelIndex = round(double(IEDoccurrence(:, 2)));
    keepRows = isfinite(localChannelIndex) & ...
        ismember(localChannelIndex, localChannels);
    selectedIEDs = double(IEDoccurrence(keepRows, :));
end

function validRows = getValidIEDRows(IEDoccurrence, validTrialMask)
    if isempty(IEDoccurrence) || size(IEDoccurrence, 2) < 3
        validRows = false(0, 1);
        return;
    end
    trialNumber = round(IEDoccurrence(:, 1));
    nTrials = numel(validTrialMask);
    validRows = isfinite(trialNumber) & trialNumber >= 1 & ...
        trialNumber <= nTrials & isfinite(IEDoccurrence(:, 3)) & ...
        IEDoccurrence(:, 3) >= 1;
    validIndices = find(validRows);
    if ~isempty(validIndices)
        validRows(validIndices) = ...
            validTrialMask(trialNumber(validIndices));
    end
end

function T = emptyBilateralCountingTable()
    T = table( ...
        strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        false(0, 1), false(0, 1), zeros(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), ...
        'VariableNames', { ...
            'patientID', 'trialNumber', 'tStart', 'tStop', ...
            'censored', 'eventAtStop', 'postIEDLeft', ...
            'postIEDRight', 'balloonColorCode', 'expectedReward', ...
            'samplingFrequencyHz', 'transitionType'});
end

function trialRows = makeBilateralRows( ...
    patientID, trialNumber, durationSeconds, leftIEDs, rightIEDs, ...
    samplingFrequencyHz, balloonColorCode, expectedReward, ...
    finalEventObserved, postIEDWindowSeconds, transitionType)

    leftTimes = getTrialIEDTimes(leftIEDs, trialNumber, ...
        durationSeconds, samplingFrequencyHz);
    rightTimes = getTrialIEDTimes(rightIEDs, trialNumber, ...
        durationSeconds, samplingFrequencyHz);

    leftEnds = min(leftTimes + postIEDWindowSeconds, durationSeconds);
    rightEnds = min(rightTimes + postIEDWindowSeconds, durationSeconds);
    breakTimes = unique([0; leftTimes(:); leftEnds(:); ...
        rightTimes(:); rightEnds(:); durationSeconds], 'sorted');
    breakTimes = breakTimes(isfinite(breakTimes) & ...
        breakTimes >= 0 & breakTimes <= durationSeconds);
    if numel(breakTimes) < 2
        breakTimes = [0; durationSeconds];
    end

    tStart = breakTimes(1:end-1);
    tStop = breakTimes(2:end);
    keep = tStop > tStart;
    tStart = tStart(keep);
    tStop = tStop(keep);
    midpoint = (tStart + tStop) / 2;

    postIEDLeft = windowIndicator(midpoint, leftTimes, ...
        postIEDWindowSeconds);
    postIEDRight = windowIndicator(midpoint, rightTimes, ...
        postIEDWindowSeconds);

    censored = true(numel(tStart), 1);
    eventAtStop = false(numel(tStart), 1);
    if finalEventObserved
        censored(end) = false;
        eventAtStop(end) = true;
    end

    variableNames = { ...
        'patientID', 'trialNumber', 'tStart', 'tStop', ...
        'censored', 'eventAtStop', 'postIEDLeft', ...
        'postIEDRight', 'balloonColorCode', 'expectedReward', ...
        'samplingFrequencyHz', 'transitionType'};

    trialRows = table( ...
        repmat(patientID, numel(tStart), 1), ...
        repmat(trialNumber, numel(tStart), 1), ...
        tStart, tStop, censored, eventAtStop, ...
        postIEDLeft, postIEDRight, ...
        repmat(balloonColorCode, numel(tStart), 1), ...
        repmat(expectedReward, numel(tStart), 1), ...
        repmat(samplingFrequencyHz, numel(tStart), 1), ...
        repmat(string(transitionType), numel(tStart), 1), ...
        'VariableNames', variableNames);
end

function times = getTrialIEDTimes( ...
    IEDoccurrence, trialNumber, durationSeconds, samplingFrequencyHz)
    if isempty(IEDoccurrence) || size(IEDoccurrence, 2) < 3
        times = zeros(0, 1);
        return;
    end
    rows = isfinite(IEDoccurrence(:, 1)) & ...
        round(IEDoccurrence(:, 1)) == trialNumber & ...
        isfinite(IEDoccurrence(:, 3)) & IEDoccurrence(:, 3) >= 1;
    times = double(IEDoccurrence(rows, 3)) ./ samplingFrequencyHz;
    times = sort(times(isfinite(times) & times > 0 & ...
        times < durationSeconds));
end

function indicator = windowIndicator(midpoint, IEDtimes, windowSeconds)
    indicator = zeros(numel(midpoint), 1);
    for ii = 1:numel(midpoint)
        indicator(ii) = double(any(midpoint(ii) >= IEDtimes & ...
            midpoint(ii) <= IEDtimes + windowSeconds));
    end
end

function [X, predictorNames, referenceColorName] = ...
    buildBilateralPredictorMatrix(countingData, outcome, familyType)

    L = double(countingData.postIEDLeft > 0);
    R = double(countingData.postIEDRight > 0);
    V = double(countingData.expectedReward);
    isCompeting = outcome == "BR";

    if isCompeting
        bank = double(countingData.transitionType == "bank");
        pop = double(countingData.transitionType == "pop");
        X = [L .* bank, R .* bank, L .* pop, R .* pop];
        predictorNames = ["L_bank"; "R_bank"; "L_pop"; "R_pop"];
        if familyType == "mechanistic"
            X(:, end + 1) = V .* bank;
            X(:, end + 1) = V .* pop;
            predictorNames(end + 1:end + 2, 1) = ...
                ["V_bank"; "V_pop"];
        end
    else
        X = [L, R];
        predictorNames = ["L"; "R"];
        if familyType == "mechanistic"
            X(:, end + 1) = V;
            predictorNames(end + 1, 1) = "V";
        end
    end

    observedColors = unique(countingData.balloonColorCode, 'sorted');
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

    for cc = 1:numel(comparisonColors)
        colorCode = comparisonColors(cc);
        colorDummy = double(countingData.balloonColorCode == colorCode);
        contrastName = colorNames(colorCode) + "_vs_" + ...
            colorNames(referenceColor);
        if isCompeting
            X(:, end + 1) = colorDummy .* bank;
            predictorNames(end + 1, 1) = contrastName + "_bank";
            X(:, end + 1) = colorDummy .* pop;
            predictorNames(end + 1, 1) = contrastName + "_pop";
        else
            X(:, end + 1) = colorDummy;
            predictorNames(end + 1, 1) = contrastName;
        end
    end

    if familyType == "mechanistic"
        if isCompeting
            X(:, end + 1) = L .* V .* bank;
            X(:, end + 1) = R .* V .* bank;
            X(:, end + 1) = L .* V .* pop;
            X(:, end + 1) = R .* V .* pop;
            predictorNames(end + 1:end + 4, 1) = [ ...
                "LxV_bank"; "RxV_bank"; "LxV_pop"; "RxV_pop"];
        else
            X(:, end + 1) = L .* V;
            X(:, end + 1) = R .* V;
            predictorNames(end + 1:end + 2, 1) = ["LxV"; "RxV"];
        end
    end
end

function summaries = getEndpointSummaries( ...
    beta, robustCovariance, robustSE, robustZ, robustP, ...
    predictorNames, endpoint, familyType)

    if familyType == "ordinary"
        if endpoint == "BR"
            leftName = "L_bank";
            rightName = "R_bank";
        elseif endpoint == "PR"
            leftName = "L_pop";
            rightName = "R_pop";
        else
            leftName = "L";
            rightName = "R";
        end
    else
        if endpoint == "BR"
            leftName = "LxV_bank";
            rightName = "RxV_bank";
        elseif endpoint == "PR"
            leftName = "LxV_pop";
            rightName = "RxV_pop";
        else
            leftName = "LxV";
            rightName = "RxV";
        end
    end

    leftIndex = find(predictorNames == leftName, 1);
    rightIndex = find(predictorNames == rightName, 1);
    if isempty(leftIndex) || isempty(rightIndex)
        error('Could not find L/R predictors for endpoint %s.', endpoint);
    end

    summaries.left = summarizeCoefficient(beta, robustSE, ...
        robustZ, robustP, leftIndex);
    summaries.right = summarizeCoefficient(beta, robustSE, ...
        robustZ, robustP, rightIndex);
    summaries.difference = summarizeDifference(beta, ...
        robustCovariance, leftIndex, rightIndex);
end

function pair = emptySummaryPair()
    pair = struct('left', emptyCoefficientSummary(), ...
        'right', emptyCoefficientSummary(), ...
        'difference', emptyCoefficientSummary());
end

function summary = emptyCoefficientSummary()
    summary = struct('beta', NaN, 'se', NaN, 'z', NaN, ...
        'p', NaN, 'betaLow', NaN, 'betaHigh', NaN, ...
        'hr', NaN, 'hrLow', NaN, 'hrHigh', NaN);
end

function summary = summarizeCoefficient( ...
    beta, robustSE, robustZ, robustP, index)
    summary = emptyCoefficientSummary();
    summary.beta = beta(index);
    summary.se = robustSE(index);
    summary.z = robustZ(index);
    summary.p = robustP(index);
    summary.betaLow = summary.beta - 1.96 * summary.se;
    summary.betaHigh = summary.beta + 1.96 * summary.se;
    summary.hr = exp(summary.beta);
    summary.hrLow = exp(summary.betaLow);
    summary.hrHigh = exp(summary.betaHigh);
end

function summary = summarizeDifference( ...
    beta, robustCovariance, leftIndex, rightIndex)
    summary = emptyCoefficientSummary();
    contrast = zeros(numel(beta), 1);
    contrast(leftIndex) = 1;
    contrast(rightIndex) = -1;
    summary.beta = contrast' * beta(:);
    variance = contrast' * robustCovariance * contrast;
    summary.se = sqrt(max(variance, 0));
    if isfinite(summary.se) && summary.se > 0
        summary.z = summary.beta / summary.se;
        summary.p = 2 * normcdf(-abs(summary.z), 0, 1);
        summary.betaLow = summary.beta - 1.96 * summary.se;
        summary.betaHigh = summary.beta + 1.96 * summary.se;
        summary.hr = exp(summary.beta);
        summary.hrLow = exp(summary.betaLow);
        summary.hrHigh = exp(summary.betaHigh);
    end
end

function [robustCovariance, robustSE, robustZ, robustP] = ...
    computeClusterRobustInference(stats, beta, clusterID)
    beta = beta(:);
    nPredictors = numel(beta);
    clusterID = clusterID(:);
    uniqueClusters = unique(clusterID(isfinite(clusterID)));
    nClusters = numel(uniqueClusters);
    canUseScores = isfield(stats, 'scores') && ...
        size(stats.scores, 1) == numel(clusterID) && ...
        size(stats.scores, 2) == nPredictors;

    if nClusters < 2 || ~canUseScores
        robustCovariance = stats.covb;
        robustSE = stats.se(:);
        robustZ = stats.z(:);
        robustP = stats.p(:);
        return;
    end

    clusterScores = zeros(nClusters, nPredictors);
    for gg = 1:nClusters
        rows = clusterID == uniqueClusters(gg);
        thisScore = stats.scores(rows, :);
        thisScore(~isfinite(thisScore)) = 0;
        clusterScores(gg, :) = sum(thisScore, 1);
    end
    meat = clusterScores' * clusterScores;
    bread = stats.covb;
    robustCovariance = bread * meat * bread;
    robustCovariance = (nClusters / (nClusters - 1)) * ...
        robustCovariance;
    robustSE = sqrt(max(diag(robustCovariance), 0));
    badSE = robustSE <= 0 | ~isfinite(robustSE);
    robustSE(badSE) = stats.se(badSE);
    robustZ = beta ./ robustSE;
    robustP = 2 .* normcdf(-abs(robustZ), 0, 1);
end

function plotLRForestFigure( ...
    results, familyLabel, settings, outputFile, axisLimitMode)
    fig = figure('Color', 'w', 'Position', settings.figurePosition);
    tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    nexttile;
    plotLRPanel(results(results.endpoint == "IT", :), ...
        "IT", familyLabel, settings, axisLimitMode);
    nexttile;
    plotLRPanel(results(results.endpoint == "RT", :), ...
        "RT", familyLabel, settings, axisLimitMode);
    nexttile;
    plotLRPanel(results(ismember(results.endpoint, ["BR", "PR"]), :), ...
        "BR / PR", familyLabel, settings, axisLimitMode);

    exportgraphics(fig, outputFile, 'ContentType', 'vector');
    close(fig);
end

function plotLRPanel( ...
    T, panelTitle, familyLabel, settings, axisLimitMode)
    ax = gca;
    hold(ax, 'on');

    fitted = T.status == "fitted" & ...
        isfinite(T.hazardRatioLeft) & isfinite(T.hazardRatioRight) & ...
        T.hazardRatioLeft > 0 & T.hazardRatioRight > 0;
    T = T(fitted, :);

    if isempty(T)
        axis(ax, 'off');
        text(ax, 0.5, 0.5, ...
            'No selected area produced an estimable bilateral model', ...
            'HorizontalAlignment', 'center', 'FontSize', 11);
        title(ax, panelTitle, 'FontSize', 14, 'FontWeight', 'bold');
        return;
    end

    [~, order] = sort(T.hazardRatioLeftOverRight, 'descend');
    T = T(order, :);
    n = height(T);
    y = (1:n)';
    offset = 0.14;

    for ii = 1:n
        plot(ax, [T.hazardRatioLeftCILow(ii), ...
            T.hazardRatioLeftCIHigh(ii)], ...
            [y(ii) + offset, y(ii) + offset], '-', ...
            'Color', settings.leftColor, 'LineWidth', 1.4);
        plot(ax, T.hazardRatioLeft(ii), y(ii) + offset, 'o', ...
            'Color', settings.leftColor, ...
            'MarkerFaceColor', settings.leftColor, ...
            'MarkerSize', settings.markerSize);

        plot(ax, [T.hazardRatioRightCILow(ii), ...
            T.hazardRatioRightCIHigh(ii)], ...
            [y(ii) - offset, y(ii) - offset], '-', ...
            'Color', settings.rightColor, 'LineWidth', 1.4);
        plot(ax, T.hazardRatioRight(ii), y(ii) - offset, 's', ...
            'Color', settings.rightColor, ...
            'MarkerFaceColor', settings.rightColor, ...
            'MarkerSize', settings.markerSize);
    end

    xline(ax, 1, '--k', 'LineWidth', 1);
    set(ax, 'XScale', 'log', 'YDir', 'reverse', ...
        'YTick', y, 'FontSize', 10, 'Box', 'off');
    ylim(ax, [0.4, n + 0.6]);

    labels = T.anatomicalArea;
    brprRows = ismember(T.endpoint, ["BR", "PR"]);
    labels(brprRows) = labels(brprRows) + " (" + ...
        T.endpoint(brprRows) + ")";
    for ii = 1:n
        if T.significantLeftVsRightFDR(ii)
            labels(ii) = labels(ii) + "  L-R*";
        end
    end
    set(ax, 'YTickLabel', labels);

    if axisLimitMode == "pointEstimates"
        % For the mechanistic figure, derive the displayed range from the
        % HR point estimates. Very wide CI endpoints are clipped by the
        % axes instead of stretching the entire panel.
        finiteLimits = [T.hazardRatioLeft; ...
            T.hazardRatioRight; 1];
    else
        % Keep the original CI-based range for the ordinary figure.
        finiteLimits = [T.hazardRatioLeftCILow; ...
            T.hazardRatioLeftCIHigh; T.hazardRatioRightCILow; ...
            T.hazardRatioRightCIHigh; 1];
    end
    finiteLimits = finiteLimits(isfinite(finiteLimits) & finiteLimits > 0);
    if ~isempty(finiteLimits)
        logLimits = log10(finiteLimits);
        logLow = min(logLimits);
        logHigh = max(logLimits);
        logRange = logHigh - logLow;
        logPadding = max(0.08 * logRange, log10(1.08));
        if logRange == 0
            logPadding = log10(1.5);
        end
        low = 10^(logLow - logPadding);
        high = 10^(logHigh + logPadding);
        xlim(ax, [low high]);
    end

    formatHazardRatioTicks(ax);
    xlabel(ax, 'Hazard ratio (95% CI)');
    title(ax, panelTitle, 'FontSize', 14, 'FontWeight', 'bold');
    legend(ax, {'Left 95% CI', 'Left HR', ...
        'Right 95% CI', 'Right HR'}, ...
        'Location', 'best', 'Box', 'off');

    subtitleText = familyLabel + ...
        "; L-R* means BH-FDR q_{L vs R} < " + string(settings.alpha);
    text(ax, 0.5, 1.01, subtitleText, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9);
end

function adjustedP = benjaminiHochbergAdjustedP(pValues)
    pValues = double(pValues(:));
    adjustedP = NaN(size(pValues));
    valid = isfinite(pValues) & pValues >= 0 & pValues <= 1;
    p = pValues(valid);
    m = numel(p);
    if m == 0
        return;
    end

    [sortedP, order] = sort(p, 'ascend');
    sortedAdjusted = sortedP .* m ./ (1:m)';
    for ii = m-1:-1:1
        sortedAdjusted(ii) = min(sortedAdjusted(ii), ...
            sortedAdjusted(ii + 1));
    end
    sortedAdjusted = min(sortedAdjusted, 1);

    unsortedAdjusted = NaN(m, 1);
    unsortedAdjusted(order) = sortedAdjusted;
    adjustedP(valid) = unsortedAdjusted;
end

function formatHazardRatioTicks(ax)
    limits = xlim(ax);
    candidateTicks = [0.05 0.1 0.2 0.5 1 2 5 10 20];
    ticks = candidateTicks(candidateTicks >= limits(1) & ...
        candidateTicks <= limits(2));
    if numel(ticks) < 2
        ticks = logspace(log10(limits(1)), log10(limits(2)), 4);
    end
    xticks(ax, ticks);
    labels = arrayfun(@(x) sprintf('%.3g', x), ticks, ...
        'UniformOutput', false);
    xticklabels(ax, labels);
end
