% Summarize where RT- and IT-period IEDs occur anatomically.
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
% If anatomicalLocs was already saved only for selected channels, this
% script detects that situation and maps the labels directly.
%
% Default filtering matches the behavioral analyses:
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
%   5) brain_area_IED_distribution.pdf (single normalized panel)
%
% The CSV summaries include raw and normalized counts. The PDF shows only
% IEDs per implanted channel because raw counts depend on electrode coverage.
%
% Author: Nill

clear;
clc;
close all;

%% Paths and settings

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED7_brain_area_IED_summary\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

maximumRTSeconds = 10;
useOnlyNonControlTrials = true;
maximumNumberOfAreasInFigure = 20;

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

%% Plot group distribution

figureOutputFile = fullfile( ...
    outputFolderName, 'brain_area_IED_distribution.pdf');

plotBrainAreaSummary( ...
    groupSummary, ...
    figureOutputFile, ...
    maximumNumberOfAreasInFigure, ...
    colorRT, ...
    colorIT);

%% Display completion information

fprintf('\nBrain-area summary finished.\n');
fprintf('Saved: %s\n', allEventsOutputFile);
fprintf('Saved: %s\n', perPatientOutputFile);
fprintf('Saved: %s\n', groupOutputFile);
fprintf('Saved: %s\n', coverageOutputFile);
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

function plotBrainAreaSummary( ...
    groupSummary, outputPDF, maximumNumberOfAreas, colorRT, colorIT)

    if isempty(groupSummary)
        warning('Group summary is empty. Figure was not created.');
        return;
    end

    areas = unique(groupSummary.anatomicalArea, 'stable');

    totalAcrossPhases = zeros(length(areas), 1);

    for aa = 1:length(areas)
        rows = groupSummary.anatomicalArea == areas(aa);
        totalAcrossPhases(aa) = sum(groupSummary.nIEDs(rows));
    end

    [~, sortIdx] = sort(totalAcrossPhases, 'descend');
    sortIdx = sortIdx(1:min(maximumNumberOfAreas, length(sortIdx)));
    selectedAreas = areas(sortIdx);

    nAreas = length(selectedAreas);
    counts = zeros(nAreas, 2);
    normalizedCounts = zeros(nAreas, 2);

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
                counts(aa, pp) = groupSummary.nIEDs(find(row, 1));
                normalizedCounts(aa, pp) = ...
                    groupSummary.IEDsPerCoveredChannel(find(row, 1));
            end
        end
    end

    % Single coverage-normalized panel.
    % Anatomical areas are shown on the x-axis and normalized IED counts
    % are shown on the y-axis. Areas remain ordered from highest to lowest
    % total raw IED count across RT and IT.

    fig = figure('Visible', 'off');
    set(fig, 'Position', [100 100 1550 850]);

    b = bar(normalizedCounts, 'grouped');
    b(1).FaceColor = colorRT;
    b(2).FaceColor = colorIT;
    b(1).EdgeColor = 'none';
    b(2).EdgeColor = 'none';
    b(1).FaceAlpha = 0.5;
    b(2).FaceAlpha = 0.5;

    ax = gca;
    ax.XTick = 1:nAreas;
    ax.XTickLabel = selectedAreas;
    ax.XTickLabelRotation = 45;
    ax.TickLabelInterpreter = 'none';
    ax.FontSize = 10;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    ax.Layer = 'top';

    xlabel('Anatomical area');
    ylabel('IED events per implanted channel');
    title('Coverage-normalized IED distribution by brain area');
    legend({'RT', 'IT'}, 'Location', 'best');

    xlim([0.5, nAreas + 0.5]);

    maximumY = max(normalizedCounts(:));
    if ~isfinite(maximumY) || maximumY <= 0
        maximumY = 1;
    end
    ylim([0, maximumY * 1.10]);
    box off;
    grid off;

    exportgraphics(fig, outputPDF, 'ContentType', 'vector');
    close(fig);

end
