% Summarize where RT- and IT-period IEDs occur anatomically, but plot
% only anatomical areas that are FDR-significant in EITHER analysis folder:
%
%   1) IED7_Cox_IT_RT_BR_postIED_by_brain_area
%      Significant variable: significantFDR
%
%   2) IED8_brain_area_IED_expected_reward
%      Significant variable: significantFDR_IEDxExpectedReward
%
% Definition used here
% --------------------
% An area is retained when it is significant in at least one of IT, RT, or
% BR in the post-IED folder OR significant in at least one of IT, RT, or BR
% in the expected-reward folder. Thus, this is the union across all
% significant areas from both folders.
%
% Creates:
%   1) all_IED_events_with_anatomy.csv
%   2) per_patient_brain_area_summary.csv
%   3) group_brain_area_summary.csv
%   4) channel_coverage_by_brain_area.csv
%   5) significant_area_membership_union_both_folders.csv
%   6) brain_area_IED_distribution_union_both_folders.pdf
%
% The PDF retains the original grouped-bar format:
%   - RT-period IED events per implanted channel
%   - IT-period IED events per implanted channel
%   - participant coverage counts above each bar
%
% Author: Nill

clear;
clc;
close all;

%% Paths and settings

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED9_brain_area_IED_summary\';

% Folder containing the ordinary post-IED Cox results.
postIEDOutputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED7_Cox_IT_RT_BR_postIED_by_brain_area\';

% Folder containing the mechanistic IED x expected-reward results.
expectedRewardOutputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED8_brain_area_IED_expected_reward\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

maximumRTSeconds = 20;
useOnlyNonControlTrials = true;

% Combine homologous left- and right-hemisphere labels into one area.
% Examples:
%   Left Hippocampus  + Right Hippocampus  -> Hippocampus
%   Hippocampus_L     + Hippocampus_R      -> Hippocampus
%   L-Amygdala        + R-Amygdala         -> Amygdala
combineLeftAndRight = true;

% Colors used in the existing analyses.
colorRT = [0.204 0.459 0.702];
colorIT = [0.847 0.333 0.153];

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

if isempty(fileList)
    error('No .LFPIED.mat files were found in the input folder.');
end

allEvents = emptyIEDAnatomyTable();
coverageTable = table();
perPatientSummary = table();

%% Process participants

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    patientID = string(fileNameParts{1});

    fprintf('\nProcessing patient: %s\n', patientID);

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName));

    if ~isfield(loadedData, 'LFPIED')
        fprintf('Skipped: LFPIED structure was not found.\n');
        continue;
    end

    LFPIED = loadedData.LFPIED;

    requiredFields = { ...
        'selectedChans', ...
        'anatomicalLocs', ...
        'RTs', ...
        'ITs', ...
        'isControl', ...
        'IED_occurance_RT', ...
        'IED_occurance_IT'};

    missingField = false;

    for ff = 1:length(requiredFields)
        if ~isfield(LFPIED, requiredFields{ff})
            fprintf('Skipped: missing field %s.\n', requiredFields{ff});
            missingField = true;
        end
    end

    if missingField
        continue;
    end

    selectedChans = round(LFPIED.selectedChans(:));
    nSelectedChannels = length(selectedChans);

    if nSelectedChannels == 0
        fprintf('Skipped: no selected ECoG channels.\n');
        continue;
    end

    selectedAreaLabels = getSelectedChannelLabels( ...
        LFPIED.anatomicalLocs, selectedChans, nSelectedChannels, ...
        'anatomicalLocs');

    if isfield(LFPIED, 'trodeLabels')
        selectedElectrodeLabels = getSelectedChannelLabels( ...
            LFPIED.trodeLabels, selectedChans, nSelectedChannels, ...
            'trodeLabels');
    else
        selectedElectrodeLabels = "Channel_" + string(selectedChans);
    end

    selectedAreaLabels = cleanAreaLabels( ...
        selectedAreaLabels, combineLeftAndRight);

    % Save channel coverage even if a channel had no IED.
    participantCoverage = table( ...
        repmat(patientID, nSelectedChannels, 1), ...
        (1:nSelectedChannels)', ...
        selectedChans, ...
        selectedElectrodeLabels, ...
        selectedAreaLabels, ...
        'VariableNames', { ...
            'patientID', ...
            'localChannelIndex', ...
            'originalChannelIndex', ...
            'electrodeLabel', ...
            'anatomicalArea'});

    coverageTable = [coverageTable; participantCoverage];

    RTs = LFPIED.RTs(:);
    ITs = LFPIED.ITs(:);
    isControl = LFPIED.isControl(:);

    if isfield(LFPIED, 'nTrials') && isscalar(LFPIED.nTrials) && ...
            isfinite(LFPIED.nTrials)
        nTrials = round(LFPIED.nTrials);
    else
        nTrials = min([length(RTs), length(ITs), length(isControl)]);
    end

    nTrials = min([nTrials, length(RTs), length(ITs), length(isControl)]);

    RTs = RTs(1:nTrials);
    ITs = ITs(1:nTrials);
    isControl = isControl(1:nTrials);

    validRTTrials = ...
        isfinite(RTs) & ...
        RTs > 0 & ...
        RTs <= maximumRTSeconds;

    validITTrials = ...
        validRTTrials & ...
        isfinite(ITs) & ...
        ITs > 0;

    if useOnlyNonControlTrials
        validRTTrials = validRTTrials & isControl == 0;
        validITTrials = validITTrials & isControl == 0;
    end

    eventsRT = makeIEDAnatomyTable( ...
        patientID, ...
        "RT", ...
        LFPIED.IED_occurance_RT, ...
        validRTTrials, ...
        selectedChans, ...
        selectedElectrodeLabels, ...
        selectedAreaLabels);

    eventsIT = makeIEDAnatomyTable( ...
        patientID, ...
        "IT", ...
        LFPIED.IED_occurance_IT, ...
        validITTrials, ...
        selectedChans, ...
        selectedElectrodeLabels, ...
        selectedAreaLabels);

    participantEvents = [eventsRT; eventsIT];
    allEvents = [allEvents; participantEvents];

    % Include zero-IED rows for every covered area in both phases.
    participantAreas = unique(selectedAreaLabels, 'stable');

    for aa = 1:length(participantAreas)

        thisArea = participantAreas(aa);
        areaChannelMask = selectedAreaLabels == thisArea;
        nCoveredChannels = sum(areaChannelMask);

        for phase = ["RT", "IT"]

            phaseRows = participantEvents.phase == phase & ...
                participantEvents.anatomicalArea == thisArea;

            nIEDs = sum(phaseRows);

            if nIEDs > 0
                nUniqueIEDChannels = length(unique( ...
                    participantEvents.localChannelIndex(phaseRows)));
            else
                nUniqueIEDChannels = 0;
            end

            allPhaseRows = participantEvents.phase == phase;
            nAllIEDsInPhase = sum(allPhaseRows);

            if nAllIEDsInPhase > 0
                percentOfPatientIEDs = 100 * nIEDs / nAllIEDsInPhase;
            else
                percentOfPatientIEDs = 0;
            end

            IEDsPerCoveredChannel = nIEDs / nCoveredChannels;

            newSummaryRow = table( ...
                patientID, ...
                phase, ...
                thisArea, ...
                nCoveredChannels, ...
                nUniqueIEDChannels, ...
                nIEDs, ...
                percentOfPatientIEDs, ...
                IEDsPerCoveredChannel, ...
                'VariableNames', { ...
                    'patientID', ...
                    'phase', ...
                    'anatomicalArea', ...
                    'nCoveredChannels', ...
                    'nUniqueIEDChannels', ...
                    'nIEDs', ...
                    'percentOfPatientIEDsInPhase', ...
                    'IEDsPerCoveredChannel'});

            perPatientSummary = [perPatientSummary; newSummaryRow];

        end
    end

    fprintf('RT IED events retained: %d\n', height(eventsRT));
    fprintf('IT IED events retained: %d\n', height(eventsIT));

end

if isempty(coverageTable)
    error('No valid channel coverage information was created.');
end

%% Make group summary

groupSummary = makeGroupSummary(perPatientSummary);

%% Save tables

allEventsOutputFile = fullfile( ...
    outputFolderName, 'all_IED_events_with_anatomy.csv');

perPatientOutputFile = fullfile( ...
    outputFolderName, 'per_patient_brain_area_summary.csv');

groupOutputFile = fullfile( ...
    outputFolderName, 'group_brain_area_summary.csv');

coverageOutputFile = fullfile( ...
    outputFolderName, 'channel_coverage_by_brain_area.csv');

writetable(allEvents, allEventsOutputFile);
writetable(perPatientSummary, perPatientOutputFile);
writetable(groupSummary, groupOutputFile);
writetable(coverageTable, coverageOutputFile);

%% Read both result folders and identify the union of significant areas

% Post-IED Cox result files.
postIED_IT_File = fullfile( ...
    postIEDOutputFolderName, 'IT_brain_area_cox_results.csv');

postIED_RT_File = fullfile( ...
    postIEDOutputFolderName, 'RT_brain_area_cox_results.csv');

postIED_BR_File = fullfile( ...
    postIEDOutputFolderName, 'BR_brain_area_cox_results.csv');

% Mechanistic IED x expected-reward result files.
expectedReward_IT_File = fullfile( ...
    expectedRewardOutputFolderName, ...
    'IT_mechanistic_IED_x_expected_reward_results.csv');

expectedReward_RT_File = fullfile( ...
    expectedRewardOutputFolderName, ...
    'RT_mechanistic_IED_x_expected_reward_results.csv');

expectedReward_BR_File = fullfile( ...
    expectedRewardOutputFolderName, ...
    'BR_mechanistic_IED_x_expected_reward_results.csv');

[significantAreas, significantAreaMembership] = ...
    getSignificantAreaUnionAcrossFolders( ...
        postIED_IT_File, postIED_RT_File, postIED_BR_File, ...
        expectedReward_IT_File, expectedReward_RT_File, ...
        expectedReward_BR_File, combineLeftAndRight);

membershipOutputFile = fullfile( ...
    outputFolderName, ...
    'significant_area_membership_union_both_folders.csv');

writetable(significantAreaMembership, membershipOutputFile);

%% Plot the union of FDR-significant areas from both folders

figureOutputFile = fullfile( ...
    outputFolderName, ...
    'brain_area_IED_distribution_union_both_folders.pdf');

plotBrainAreaSummaryForSelectedAreas( ...
    groupSummary, ...
    significantAreas, ...
    significantAreaMembership, ...
    figureOutputFile, ...
    colorRT, ...
    colorIT);

%% Display completion information

fprintf('\nBrain-area summary finished.\n');
fprintf('Union of FDR-significant areas across both folders: %d\n', ...
    length(significantAreas));
fprintf('Saved: %s\n', allEventsOutputFile);
fprintf('Saved: %s\n', perPatientOutputFile);
fprintf('Saved: %s\n', groupOutputFile);
fprintf('Saved: %s\n', coverageOutputFile);
fprintf('Saved: %s\n', membershipOutputFile);
fprintf('Saved: %s\n', figureOutputFile);

%% Local functions

function labels = getSelectedChannelLabels( ...
    allLabels, selectedChans, nSelectedChannels, fieldName)

    allLabels = convertLabelsToString(allLabels);
    allLabels = allLabels(:);

    % Most likely structure: labels correspond to all original channels.
    if length(allLabels) >= max(selectedChans)
        labels = allLabels(selectedChans);

    % Alternative structure: labels were already restricted to selectedChans.
    elseif length(allLabels) == nSelectedChannels
        labels = allLabels;

    else
        error([ ...
            'Cannot map LFPIED.%s to selected channels. ' ...
            'Number of labels = %d, number of selected channels = %d, ' ...
            'maximum selected channel index = %d.'], ...
            fieldName, length(allLabels), nSelectedChannels, ...
            max(selectedChans));
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
            elseif ischar(value) || isstring(value)
                labels(ii) = string(value);
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
        % Remove common hemisphere markers from the beginning or end of
        % each label. The delimiter requirement prevents removal of letters
        % that are part of the anatomical name itself.
        %
        % Examples handled:
        %   Left Hippocampus, Right Hippocampus
        %   L Hippocampus, R-Hippocampus
        %   Hippocampus Left, Hippocampus_R
        %   (L) Hippocampus, Hippocampus (Right)
        labels = regexprep(labels, ...
            '^\s*\(?\s*(Left|Right|LH|RH|L|R)\s*\)?[_\-\s]+', ...
            '', 'ignorecase');

        labels = regexprep(labels, ...
            '[_\-\s]+\(?\s*(Left|Right|LH|RH|L|R)\s*\)?\s*$', ...
            '', 'ignorecase');

        % Remove optional words such as "hemisphere" when they accompany
        % the side designation.
        labels = regexprep(labels, ...
            '^\s*(Left|Right)\s+hemisphere[_\-\s]+', ...
            '', 'ignorecase');

        labels = regexprep(labels, ...
            '[_\-\s]+(Left|Right)\s+hemisphere\s*$', ...
            '', 'ignorecase');

        % Clean separators and repeated whitespace left after removal.
        labels = regexprep(labels, '[_\-]+', ' ');
        labels = regexprep(labels, '\s+', ' ');
        labels = strip(labels);
        labels(strlength(labels) == 0) = "Unknown";
    end

end

function eventTable = emptyIEDAnatomyTable()

    eventTable = table( ...
        strings(0, 1), ...
        strings(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        zeros(0, 1), ...
        'VariableNames', { ...
            'patientID', ...
            'phase', ...
            'trialNumber', ...
            'localChannelIndex', ...
            'originalChannelIndex', ...
            'electrodeLabel', ...
            'anatomicalArea', ...
            'timeIndexWithinWindow'});

end

function eventTable = makeIEDAnatomyTable( ...
    patientID, phase, IEDocc, validTrialMask, selectedChans, ...
    selectedElectrodeLabels, selectedAreaLabels)

    eventTable = emptyIEDAnatomyTable();

    if isempty(IEDocc) || size(IEDocc, 2) < 3
        return;
    end

    trialNumber = round(IEDocc(:, 1));
    localChannelIndex = round(IEDocc(:, 2));
    timeIndexWithinWindow = IEDocc(:, 3);

    nTrials = length(validTrialMask);
    nSelectedChannels = length(selectedChans);

    validRows = ...
        isfinite(trialNumber) & ...
        trialNumber >= 1 & ...
        trialNumber <= nTrials & ...
        isfinite(localChannelIndex) & ...
        localChannelIndex >= 1 & ...
        localChannelIndex <= nSelectedChannels & ...
        isfinite(timeIndexWithinWindow);

    if ~any(validRows)
        return;
    end

    trialNumber = trialNumber(validRows);
    localChannelIndex = localChannelIndex(validRows);
    timeIndexWithinWindow = timeIndexWithinWindow(validRows);

    keepByTrial = validTrialMask(trialNumber);

    trialNumber = trialNumber(keepByTrial);
    localChannelIndex = localChannelIndex(keepByTrial);
    timeIndexWithinWindow = timeIndexWithinWindow(keepByTrial);

    if isempty(trialNumber)
        return;
    end

    originalChannelIndex = selectedChans(localChannelIndex);
    electrodeLabel = selectedElectrodeLabels(localChannelIndex);
    anatomicalArea = selectedAreaLabels(localChannelIndex);

    nEvents = length(trialNumber);

    eventTable = table( ...
        repmat(patientID, nEvents, 1), ...
        repmat(phase, nEvents, 1), ...
        trialNumber, ...
        localChannelIndex, ...
        originalChannelIndex, ...
        electrodeLabel, ...
        anatomicalArea, ...
        timeIndexWithinWindow, ...
        'VariableNames', { ...
            'patientID', ...
            'phase', ...
            'trialNumber', ...
            'localChannelIndex', ...
            'originalChannelIndex', ...
            'electrodeLabel', ...
            'anatomicalArea', ...
            'timeIndexWithinWindow'});

end

function groupSummary = makeGroupSummary(perPatientSummary)

    groupSummary = table();

    if isempty(perPatientSummary)
        return;
    end

    phases = unique(perPatientSummary.phase, 'stable');
    areas = unique(perPatientSummary.anatomicalArea, 'stable');

    for pp = 1:length(phases)
        for aa = 1:length(areas)

            thisPhase = phases(pp);
            thisArea = areas(aa);

            rows = ...
                perPatientSummary.phase == thisPhase & ...
                perPatientSummary.anatomicalArea == thisArea;

            if ~any(rows)
                continue;
            end

            areaRows = perPatientSummary(rows, :);

            nParticipantsWithCoverage = height(areaRows);
            nParticipantsWithIED = sum(areaRows.nIEDs > 0);
            nCoveredChannels = sum(areaRows.nCoveredChannels);
            nUniqueIEDChannels = sum(areaRows.nUniqueIEDChannels);
            nIEDs = sum(areaRows.nIEDs);

            phaseRows = perPatientSummary.phase == thisPhase;
            totalIEDsInPhase = sum(perPatientSummary.nIEDs(phaseRows));

            if totalIEDsInPhase > 0
                percentOfAllIEDsInPhase = 100 * nIEDs / totalIEDsInPhase;
            else
                percentOfAllIEDsInPhase = 0;
            end

            if nCoveredChannels > 0
                IEDsPerCoveredChannel = nIEDs / nCoveredChannels;
            else
                IEDsPerCoveredChannel = NaN;
            end

            if nUniqueIEDChannels > 0
                IEDsPerIEDChannel = nIEDs / nUniqueIEDChannels;
            else
                IEDsPerIEDChannel = 0;
            end

            meanIEDsPerPatient = mean(areaRows.nIEDs, 'omitnan');
            medianIEDsPerPatient = median(areaRows.nIEDs, 'omitnan');

            newRow = table( ...
                thisPhase, ...
                thisArea, ...
                nParticipantsWithCoverage, ...
                nParticipantsWithIED, ...
                nCoveredChannels, ...
                nUniqueIEDChannels, ...
                nIEDs, ...
                percentOfAllIEDsInPhase, ...
                IEDsPerCoveredChannel, ...
                IEDsPerIEDChannel, ...
                meanIEDsPerPatient, ...
                medianIEDsPerPatient, ...
                'VariableNames', { ...
                    'phase', ...
                    'anatomicalArea', ...
                    'nParticipantsWithCoverage', ...
                    'nParticipantsWithIED', ...
                    'nCoveredChannels', ...
                    'nUniqueIEDChannels', ...
                    'nIEDs', ...
                    'percentOfAllIEDsInPhase', ...
                    'IEDsPerCoveredChannel', ...
                    'IEDsPerIEDChannel', ...
                    'meanIEDsPerPatient', ...
                    'medianIEDsPerPatient'});

            groupSummary = [groupSummary; newRow];

        end
    end

    if ~isempty(groupSummary)
        groupSummary = sortrows( ...
            groupSummary, {'phase', 'nIEDs'}, {'ascend', 'descend'});
    end

end

function [significantAreas, membershipTable] = ...
    getSignificantAreaUnionAcrossFolders( ...
        postIED_IT_File, postIED_RT_File, postIED_BR_File, ...
        expectedReward_IT_File, expectedReward_RT_File, ...
        expectedReward_BR_File, combineLeftAndRight)

    % Within each folder, take the union across IT, RT, and BR.
    % Then take the union between the two folder-level unions.

    postIED_IT_Areas = readSignificantAreas( ...
        postIED_IT_File, "significantFDR", ...
        combineLeftAndRight, "Post-IED IT");

    postIED_RT_Areas = readSignificantAreas( ...
        postIED_RT_File, "significantFDR", ...
        combineLeftAndRight, "Post-IED RT");

    postIED_BR_Areas = readSignificantAreas( ...
        postIED_BR_File, "significantFDR", ...
        combineLeftAndRight, "Post-IED BR");

    expectedReward_IT_Areas = readSignificantAreas( ...
        expectedReward_IT_File, ...
        "significantFDR_IEDxExpectedReward", ...
        combineLeftAndRight, "Expected-reward IT");

    expectedReward_RT_Areas = readSignificantAreas( ...
        expectedReward_RT_File, ...
        "significantFDR_IEDxExpectedReward", ...
        combineLeftAndRight, "Expected-reward RT");

    expectedReward_BR_Areas = readSignificantAreas( ...
        expectedReward_BR_File, ...
        "significantFDR_IEDxExpectedReward", ...
        combineLeftAndRight, "Expected-reward BR");

    postIEDUnion = unique( ...
        [postIED_IT_Areas; postIED_RT_Areas; postIED_BR_Areas], ...
        'stable');

    expectedRewardUnion = unique( ...
        [expectedReward_IT_Areas; expectedReward_RT_Areas; ...
         expectedReward_BR_Areas], ...
        'stable');

    allCandidateAreas = unique( ...
        [postIEDUnion; expectedRewardUnion], 'stable');

    allCandidateAreas = allCandidateAreas( ...
        ~isExcludedAreaLabel(allCandidateAreas) & ...
        allCandidateAreas ~= "Unknown" & ...
        strlength(strip(allCandidateAreas)) > 0);

    nAreas = length(allCandidateAreas);

    significantPostIED_IT = ...
        ismember(allCandidateAreas, postIED_IT_Areas);
    significantPostIED_RT = ...
        ismember(allCandidateAreas, postIED_RT_Areas);
    significantPostIED_BR = ...
        ismember(allCandidateAreas, postIED_BR_Areas);

    significantExpectedReward_IT = ...
        ismember(allCandidateAreas, expectedReward_IT_Areas);
    significantExpectedReward_RT = ...
        ismember(allCandidateAreas, expectedReward_RT_Areas);
    significantExpectedReward_BR = ...
        ismember(allCandidateAreas, expectedReward_BR_Areas);

    significantInPostIEDAny = ...
        significantPostIED_IT | ...
        significantPostIED_RT | ...
        significantPostIED_BR;

    significantInExpectedRewardAny = ...
        significantExpectedReward_IT | ...
        significantExpectedReward_RT | ...
        significantExpectedReward_BR;

    significantInEitherFolder = ...
        significantInPostIEDAny | significantInExpectedRewardAny;

    postIEDSignificantIn = strings(nAreas, 1);
    expectedRewardSignificantIn = strings(nAreas, 1);

    for aa = 1:nAreas
        postIEDLabels = strings(0, 1);
        expectedRewardLabels = strings(0, 1);

        if significantPostIED_IT(aa)
            postIEDLabels(end + 1, 1) = "IT";
        end
        if significantPostIED_RT(aa)
            postIEDLabels(end + 1, 1) = "RT";
        end
        if significantPostIED_BR(aa)
            postIEDLabels(end + 1, 1) = "BR";
        end

        if significantExpectedReward_IT(aa)
            expectedRewardLabels(end + 1, 1) = "IT";
        end
        if significantExpectedReward_RT(aa)
            expectedRewardLabels(end + 1, 1) = "RT";
        end
        if significantExpectedReward_BR(aa)
            expectedRewardLabels(end + 1, 1) = "BR";
        end

        postIEDSignificantIn(aa) = strjoin(postIEDLabels, ", ");
        expectedRewardSignificantIn(aa) = ...
            strjoin(expectedRewardLabels, ", ");
    end

    membershipTable = table( ...
        allCandidateAreas, ...
        significantPostIED_IT, ...
        significantPostIED_RT, ...
        significantPostIED_BR, ...
        significantInPostIEDAny, ...
        postIEDSignificantIn, ...
        significantExpectedReward_IT, ...
        significantExpectedReward_RT, ...
        significantExpectedReward_BR, ...
        significantInExpectedRewardAny, ...
        expectedRewardSignificantIn, ...
        significantInEitherFolder, ...
        'VariableNames', { ...
            'anatomicalArea', ...
            'significantPostIED_IT', ...
            'significantPostIED_RT', ...
            'significantPostIED_BR', ...
            'significantInPostIEDAny', ...
            'postIEDSignificantIn', ...
            'significantExpectedReward_IT', ...
            'significantExpectedReward_RT', ...
            'significantExpectedReward_BR', ...
            'significantInExpectedRewardAny', ...
            'expectedRewardSignificantIn', ...
            'significantInEitherFolder'});

    significantAreas = allCandidateAreas(significantInEitherFolder);

    if isempty(significantAreas)
        error([ ...
            'No FDR-significant anatomical area was found in either folder. ' ...
            'The code checks the IT, RT, and BR result files from both ' ...
            'the post-IED and expected-reward folders.']);
    end

    fprintf('\nUnion of FDR-significant areas across both folders:\n');
    disp(membershipTable(membershipTable.significantInEitherFolder, :));

end

function significantAreas = readSignificantAreas( ...
    resultsFile, significanceVariable, ...
    combineLeftAndRight, analysisLabel)

    if ~isfile(resultsFile)
        error([ ...
            'Result file was not found:\n%s\n' ...
            'Check the two analysis-folder paths and confirm that the ' ...
            'IT, RT, and BR analyses were run first.'], resultsFile);
    end

    results = readtable(resultsFile, 'TextType', 'string');

    significanceVariable = char(significanceVariable);
    requiredVariables = {'anatomicalArea', significanceVariable};

    for vv = 1:length(requiredVariables)
        if ~ismember(requiredVariables{vv}, ...
                results.Properties.VariableNames)
            error('Missing variable %s in %s.', ...
                requiredVariables{vv}, resultsFile);
        end
    end

    results.anatomicalArea = cleanAreaLabels( ...
        results.anatomicalArea, combineLeftAndRight);

    significantRows = convertToLogical( ...
        results.(significanceVariable));

    if ismember('status', results.Properties.VariableNames)
        significantRows = significantRows & ...
            lower(strip(string(results.status))) == "fitted";
    end

    significantRows = significantRows & ...
        ~isExcludedAreaLabel(results.anatomicalArea) & ...
        results.anatomicalArea ~= "Unknown" & ...
        strlength(strip(results.anatomicalArea)) > 0;

    significantAreas = unique( ...
        results.anatomicalArea(significantRows), 'stable');

    fprintf('%s FDR-significant areas loaded: %d\n', ...
        analysisLabel, length(significantAreas));

end

function logicalValues = convertToLogical(values)

    if islogical(values)
        logicalValues = values;

    elseif isnumeric(values)
        logicalValues = isfinite(values) & values ~= 0;

    else
        normalizedValues = lower(strip(string(values)));
        logicalValues = ismember( ...
            normalizedValues, ["true", "1", "yes"]);
    end

    logicalValues = logicalValues(:);

end

function excludedRows = isExcludedAreaLabel(labels)

    normalizedLabels = lower(strip(string(labels)));
    normalizedLabels = regexprep(normalizedLabels, '[_\-]+', ' ');
    normalizedLabels = regexprep(normalizedLabels, '\s+', ' ');

    % Remove every NaC/nucleus-accumbens label, including NaC60.
    isNaC = ...
        contains(normalizedLabels, "nac") | ...
        contains(normalizedLabels, "nucleus accumbens") | ...
        contains(normalizedLabels, "accumbens");

    % Remove white-matter labels, including WM, WM60, and white_matter.
    compactLabels = regexprep(normalizedLabels, '[^a-z0-9]', '');

    isWhiteMatter = ...
        contains(normalizedLabels, "white matter") | ...
        contains(compactLabels, "whitematter") | ...
        startsWith(compactLabels, "wm");

    % Remove lateral-ventricle labels regardless of case or separators.
    isLateralVentricle = ...
        contains(normalizedLabels, "lateral ventricle") | ...
        contains(compactLabels, "lateralventricle");

    excludedRows = isNaC | isWhiteMatter | isLateralVentricle;

end

function plotBrainAreaSummaryForSelectedAreas( ...
    groupSummary, selectedAreas, membershipTable, ...
    outputPDF, colorRT, colorIT)

    if isempty(groupSummary)
        warning('Group summary is empty. Figure was not created.');
        return;
    end

    selectedAreas = string(selectedAreas(:));

    availableAreas = unique( ...
        string(groupSummary.anatomicalArea), 'stable');

    missingAreas = selectedAreas(~ismember(selectedAreas, availableAreas));

    if ~isempty(missingAreas)
        warning([ ...
            'These union-selected areas were not found in the ' ...
            'IED summary and will not be plotted: %s'], ...
            strjoin(missingAreas, ', '));
    end

    selectedAreas = selectedAreas( ...
        ismember(selectedAreas, availableAreas));

    if isempty(selectedAreas)
        error([ ...
            'None of the union-selected areas were found in the ' ...
            'IED group summary. Check that all scripts use the same ' ...
            'input data and anatomical-label cleaning rules.']);
    end

    nAreas = length(selectedAreas);
    rawCounts = zeros(nAreas, 2);
    normalizedCounts = zeros(nAreas, 2);
    participantCounts = zeros(nAreas, 2);

    for aa = 1:nAreas
        for pp = 1:2

            if pp == 1
                phase = "RT";
            else
                phase = "IT";
            end

            row = ...
                groupSummary.anatomicalArea == selectedAreas(aa) & ...
                groupSummary.phase == phase;

            if any(row)
                rowIndex = find(row, 1);

                rawCounts(aa, pp) = ...
                    groupSummary.nIEDs(rowIndex);

                normalizedCounts(aa, pp) = ...
                    groupSummary.IEDsPerCoveredChannel(rowIndex);

                % Number of participants who had at least one implanted
                % channel covering this anatomical area. These are the
                % participants contributing to the bar denominator.
                participantCounts(aa, pp) = ...
                    groupSummary.nParticipantsWithCoverage(rowIndex);
            end
        end
    end

    % Preserve the original IED7 ordering rule: descending total raw IEDs
    % across RT and IT, now applied to the union-selected areas.

    [~, sortIndex] = sort(normalizedCounts(:, 2), 'descend');

    selectedAreas = selectedAreas(sortIndex);
    rawCounts = rawCounts(sortIndex, :); 
    normalizedCounts = normalizedCounts(sortIndex, :);
    participantCounts = participantCounts(sortIndex, :);

    % Keep membership information aligned with the plotted order. It is
    % saved in the CSV, but it is intentionally not added as a second line
    % to the x-axis labels because MATLAB was treating those second lines
    % as extra tick labels in the exported PDF.
    membershipRows = ismember( ...
        membershipTable.anatomicalArea, selectedAreas);
    membershipForPlot = membershipTable(membershipRows, :); 

    % Clean x-axis: one anatomical area per tick.
    xTickLabels = cellstr(selectedAreas);

    figureWidth = max(1100, 170 * nAreas);

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100 100 figureWidth 820]);

    b = bar(normalizedCounts, 'grouped');
    b(1).FaceColor = colorRT;
    b(2).FaceColor = colorIT;
    b(1).EdgeColor = 'none';
    b(2).EdgeColor = 'none';
    b(1).FaceAlpha = 0.5;
    b(2).FaceAlpha = 0.5;

    ax = gca;
    ax.XTick = 1:nAreas;
    ax.XTickLabel = xTickLabels;
    ax.XTickLabelRotation = 35;
    ax.TickLabelInterpreter = 'none';
    ax.FontSize = 11;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    ax.Layer = 'top';

    xlabel('Anatomical area');
    ylabel('IED events per implanted channel');
    title('IED distribution in the union of significant brain areas');
    legend({'RT', 'IT'}, 'Location', 'northeast', 'Box', 'off');

    xlim([0.5, nAreas + 0.5]);

    maximumY = max(normalizedCounts(:));
    if ~isfinite(maximumY) || maximumY <= 0
        maximumY = 1;
    end

    % Leave room above the bars for participant-count labels.
    labelOffset = 0.025 * maximumY;
    ylim([0, maximumY * 1.18]);

    % Write participant counts above every RT and IT bar.
    % n = number of participants with channel coverage in that area.
    drawnow;

    for pp = 1:2
        xPositions = b(pp).XEndPoints;
        yPositions = b(pp).YEndPoints;
        countLabels = compose('n=%d', participantCounts(:, pp));

        text( ...
            xPositions(:), ...
            yPositions(:) + labelOffset, ...
            countLabels(:), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 9, ...
            'Color', [0.15 0.15 0.15], ...
            'Interpreter', 'none');
    end

    box off;
    grid off;

    exportgraphics(fig, outputPDF, 'ContentType', 'vector');
    close(fig);

end
