% Summarize where RT- and IT-period IEDs occur anatomically, but plot
% only anatomical areas that are FDR-significant in the combined IED8
% area-specific Cox analysis.
%
% IMPORTANT
% ---------
% Run IED8_Cox_IT_RT_BR_postIED_by_brain_area.m first. This script reads:
%   IT_brain_area_cox_results.csv
%   RT_brain_area_cox_results.csv
%   BR_brain_area_cox_results.csv
%
% An anatomical area is included in the figure when significantFDR is true
% in at least one of the IT, RT, or BR IED8 result tables. The plotted bars
% retain the IED7 definition:
%   - RT-period IED events per implanted channel
%   - IT-period IED events per implanted channel
%
% The BR model in IED8 uses IT-period IEDs. Therefore, an area significant
% only for BR is still shown with its RT and IT descriptive IED counts in
% this IED7-style grouped bar figure.
%
% IMPORTANT CHANNEL MAPPING
% -------------------------
% LFPIED.IED_occurance_RT(:,2) and LFPIED.IED_occurance_IT(:,2)
% contain the LOCAL ECoG channel index (chz = 1:length(selectedChans)).
% Therefore, the anatomical label is obtained safely as:
%
%   localChannelIndex    = IEDocc(:,2);
%   originalChannelIndex = LFPIED.selectedChans(localChannelIndex);
%   anatomicalArea       = LFPIED.anatomicalLocs(originalChannelIndex);
%
% Default IED7 filtering is retained:
%   - only non-control trials
%   - trials with RT > 10 s are excluded from both RT and IT
%   - RT must be finite and > 0
%   - IT must be finite and > 0 for the IT analysis
%
% Creates:
%   1) all_IED_events_with_anatomy.csv
%   2) per_patient_brain_area_summary.csv
%   3) group_brain_area_summary.csv
%   4) channel_coverage_by_brain_area.csv
%   5) IED8_significant_area_membership.csv
%   6) brain_area_IED_distribution_IED8_significant_areas.pdf
%
% The PDF is a single grouped bar panel and contains only the union of
% FDR-significant brain areas from the IT, RT, and BR IED8 analyses.
%
% Author: Nill

clear;
clc;
close all;

%% Paths and settings

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED8_brain_area_IED_summary\';

% Folder containing the CSV outputs created by the combined IED8 code.
ied8OutputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED7_Cox_IT_RT_BR_postIED_by_brain_area\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

maximumRTSeconds = 10;
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

%% Read IED8 results and identify the significant-area union

ITResultsFile = fullfile( ...
    ied8OutputFolderName, 'IT_brain_area_cox_results.csv');

RTResultsFile = fullfile( ...
    ied8OutputFolderName, 'RT_brain_area_cox_results.csv');

BRResultsFile = fullfile( ...
    ied8OutputFolderName, 'BR_brain_area_cox_results.csv');

[significantAreas, significantAreaMembership] = ...
    getIED8SignificantAreaUnion( ...
        ITResultsFile, RTResultsFile, BRResultsFile, ...
        combineLeftAndRight);

membershipOutputFile = fullfile( ...
    outputFolderName, 'IED8_significant_area_membership.csv');

writetable(significantAreaMembership, membershipOutputFile);

%% Plot only the IED8 FDR-significant anatomical areas

figureOutputFile = fullfile( ...
    outputFolderName, ...
    'brain_area_IED_distribution_IED8_significant_areas.pdf');

plotBrainAreaSummaryForSelectedAreas( ...
    groupSummary, ...
    significantAreas, ...
    significantAreaMembership, ...
    figureOutputFile, ...
    colorRT, ...
    colorIT);

%% Display completion information

fprintf('\nBrain-area summary finished.\n');
fprintf('IED8 FDR-significant areas in union: %d\n', ...
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
    getIED8SignificantAreaUnion( ...
        ITResultsFile, RTResultsFile, BRResultsFile, ...
        combineLeftAndRight)

    ITAreas = readIED8SignificantAreas( ...
        ITResultsFile, combineLeftAndRight, "IT");

    RTAreas = readIED8SignificantAreas( ...
        RTResultsFile, combineLeftAndRight, "RT");

    BRAreas = readIED8SignificantAreas( ...
        BRResultsFile, combineLeftAndRight, "BR");

    significantAreas = unique( ...
        [ITAreas; RTAreas; BRAreas], 'stable');

    significantAreas = significantAreas( ...
        ~isExcludedAreaLabel(significantAreas));

    significantAreas = significantAreas( ...
        significantAreas ~= "Unknown" & ...
        strlength(strip(significantAreas)) > 0);

    if isempty(significantAreas)
        error([ ...
            'No FDR-significant anatomical areas were found in the ' ...
            'IED8 IT, RT, or BR result files.']);
    end

    nAreas = length(significantAreas);

    membershipTable = table( ...
        significantAreas, ...
        ismember(significantAreas, ITAreas), ...
        ismember(significantAreas, RTAreas), ...
        ismember(significantAreas, BRAreas), ...
        strings(nAreas, 1), ...
        'VariableNames', { ...
            'anatomicalArea', ...
            'significantIT', ...
            'significantRT', ...
            'significantBR', ...
            'significantIn'});

    for aa = 1:nAreas
        labels = strings(0, 1);

        if membershipTable.significantIT(aa)
            labels(end + 1, 1) = "IT"; %#ok<AGROW>
        end

        if membershipTable.significantRT(aa)
            labels(end + 1, 1) = "RT"; %#ok<AGROW>
        end

        if membershipTable.significantBR(aa)
            labels(end + 1, 1) = "BR"; %#ok<AGROW>
        end

        membershipTable.significantIn(aa) = strjoin(labels, ", ");
    end

    fprintf('\nIED8 FDR-significant areas:\n');
    disp(membershipTable);

end

function significantAreas = readIED8SignificantAreas( ...
    resultsFile, combineLeftAndRight, analysisLabel)

    if ~isfile(resultsFile)
        error([ ...
            'IED8 result file was not found:\n%s\n' ...
            'Run IED8_Cox_IT_RT_BR_postIED_by_brain_area.m first, ' ...
            'or correct ied8OutputFolderName.'], resultsFile);
    end

    results = readtable(resultsFile, 'TextType', 'string');

    requiredVariables = {'anatomicalArea', 'significantFDR'};

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
        results.significantFDR);

    if ismember('status', results.Properties.VariableNames)
        significantRows = significantRows & ...
            lower(strip(string(results.status))) == "fitted";
    end

    significantRows = significantRows & ...
        ~isExcludedAreaLabel(results.anatomicalArea) & ...
        results.anatomicalArea ~= "Unknown";

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

    excludedRows = isNaC | isWhiteMatter;

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
            'These IED8-significant areas were not found in the IED7 ' ...
            'group summary and will not be plotted: %s'], ...
            strjoin(missingAreas, ', '));
    end

    selectedAreas = selectedAreas( ...
        ismember(selectedAreas, availableAreas));

    if isempty(selectedAreas)
        error([ ...
            'None of the IED8-significant areas were found in the ' ...
            'IED7 group summary. Check that both scripts use the same ' ...
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
    % across RT and IT, now applied only to IED8-significant areas.
    [~, sortIndex] = sort(sum(rawCounts, 2), 'descend');

    selectedAreas = selectedAreas(sortIndex);
    rawCounts = rawCounts(sortIndex, :); %#ok<NASGU>
    normalizedCounts = normalizedCounts(sortIndex, :);
    participantCounts = participantCounts(sortIndex, :);

    % Keep membership information aligned with the plotted order. It is
    % saved in the CSV, but it is intentionally not added as a second line
    % to the x-axis labels because MATLAB was treating those second lines
    % as extra tick labels in the exported PDF.
    membershipRows = ismember( ...
        membershipTable.anatomicalArea, selectedAreas);
    membershipForPlot = membershipTable(membershipRows, :); %#ok<NASGU>

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
    title('IED distribution in IED8 FDR-significant brain areas');
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
