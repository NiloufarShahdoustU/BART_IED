% Save participant-level channel counts for significant anatomical-area models.
%
% Run IED7_Cox_IT_RT_BR_postIED_by_brain_area first. This script then reads
% its combined results and creates one row for every participant who has
% IED channels in an FDR-significant anatomical area. Channels without an
% IED in the analyzed outcome are not counted or saved.
%
% IMPORTANT: p-values are from the pooled, participant-stratified model for
% each anatomical area. They are not participant-specific p-values.

clear;
clc;

%% Paths

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

brainAreaResultsFolder = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED7_Cox_IT_RT_BR_postIED_by_brain_area\';

combinedResultsFile = fullfile(brainAreaResultsFolder, ...
    'IT_RT_BR_brain_area_cox_all_results.mat');

outputFile = fullfile(brainAreaResultsFolder, ...
    'significant_channels_by_participant.csv');

barPlotFile = fullfile(brainAreaResultsFolder, ...
    'significant_channels_IT_RT_BR_barplot.pdf');

%% Settings matching the brain-area analysis

combineLeftAndRight = true;
includeParticipantsWithoutAreaIEDs = false;

%% Load area-level results

if ~isfile(combinedResultsFile)
    error('Combined brain-area results were not found: %s', ...
        combinedResultsFile);
end

S = load(combinedResultsFile, 'ITResults', 'RTResults', 'BRResults');

analysisNames = ["IT", "RT", "BR"];
analysisTables = {S.ITResults, S.RTResults, S.BRResults};
iedFields = ["IED_occurance_IT", "IED_occurance_RT", ...
    "IED_occurance_IT"];

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));
if isempty(fileList)
    error('No .LFPIED.mat files were found in: %s', ...
        inputFolderName_LFPIED);
end

outputTable = initializeOutputTable();

%% Match every participant to every significant area

for aa = 1:numel(analysisNames)

    analysisName = analysisNames(aa);
    areaResults = analysisTables{aa};

    significantRows = ...
        areaResults.status == "fitted" & ...
        areaResults.significantFDR & ...
        isfinite(areaResults.pValue) & ...
        isfinite(areaResults.pValueFDR);

    significantResults = areaResults(significantRows, :);

    fprintf('\n%s: %d FDR-significant anatomical areas.\n', ...
        analysisName, height(significantResults));

    for pp = 1:numel(fileList)

        fileName = fileList(pp).name;
        fileNameParts = strsplit(fileName, '.');
        patientID = string(fileNameParts{1});

        loadedData = load(fullfile(inputFolderName_LFPIED, fileName));
        if ~isfield(loadedData, 'LFPIED')
            fprintf('%s %s: skipped; LFPIED was not found.\n', ...
                analysisName, patientID);
            continue;
        end

        LFPIED = loadedData.LFPIED;
        iedField = char(iedFields(aa));

        requiredFields = {'selectedChans', 'anatomicalLocs', iedField};
        if ~all(isfield(LFPIED, requiredFields))
            fprintf('%s %s: skipped; required fields are missing.\n', ...
                analysisName, patientID);
            continue;
        end

        selectedChans = round(double(LFPIED.selectedChans(:)));
        if isempty(selectedChans)
            continue;
        end

        areaLabels = getSelectedChannelLabels( ...
            LFPIED.anatomicalLocs, selectedChans, numel(selectedChans));
        areaLabels = cleanAreaLabels(areaLabels, combineLeftAndRight);
        areaLabels(isExcludedAreaLabel(areaLabels)) = "Excluded";

        IEDoccurrence = double(LFPIED.(iedField));

        if isempty(IEDoccurrence) || size(IEDoccurrence, 2) < 2
            localIEDChannels = zeros(0, 1);
        else
            localIEDChannels = round(IEDoccurrence(:, 2));
            localIEDChannels = unique(localIEDChannels( ...
                isfinite(localIEDChannels) & ...
                localIEDChannels >= 1 & ...
                localIEDChannels <= numel(selectedChans)));
        end

        for ss = 1:height(significantResults)

            resultRow = significantResults(ss, :);
            areaLocalChannels = find( ...
                areaLabels == string(resultRow.anatomicalArea));

            if isempty(areaLocalChannels)
                continue;
            end

            areaLocalChannelsWithIED = intersect( ...
                areaLocalChannels(:), localIEDChannels(:), 'stable');

            if ~includeParticipantsWithoutAreaIEDs && ...
                    isempty(areaLocalChannelsWithIED)
                continue;
            end

            % Keep ONLY channels that contributed IEDs to a significant
            % anatomical-area result.
            significantLocalChannels = areaLocalChannelsWithIED;
            significantOriginalChannels = ...
                selectedChans(significantLocalChannels);

            newRow = table( ...
                patientID, ...
                analysisName, ...
                string(resultRow.anatomicalArea), ...
                resultRow.pValue, ...
                resultRow.pValueFDR, ...
                resultRow.hazardRatio, ...
                numel(significantLocalChannels), ...
                join(string(significantLocalChannels(:)'), ","), ...
                join(string(significantOriginalChannels(:)'), ","), ...
                'VariableNames', outputTable.Properties.VariableNames);

            outputTable = [outputTable; newRow]; %#ok<AGROW>
        end
    end
end

%% Sort and save

if ~isempty(outputTable)
    outputTable = sortrows(outputTable, ...
        {'analysis', 'pValueFDR', 'anatomicalArea', 'patientID'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend'});
end

writetable(outputTable, outputFile);

createSignificantAreaBarPlot( ...
    outputTable, analysisNames, barPlotFile);

fprintf('\n============================================================\n');
fprintf('Finished. Saved %d participant-area rows:\n%s\n', ...
    height(outputTable), outputFile);
fprintf('Combined IT/RT/BR bar plot saved:\n%s\n', barPlotFile);

%% Local functions

function T = initializeOutputTable()
    T = table( ...
        strings(0,1), strings(0,1), strings(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), ...
        strings(0,1), strings(0,1), ...
        'VariableNames', { ...
            'patientID', 'analysis', 'anatomicalArea', ...
            'pValue', 'pValueFDR', 'hazardRatio', ...
            'nSignificantChannels', ...
            'significantLocalChannelIndices', ...
            'significantOriginalChannelNumbers'});
end

function labels = getSelectedChannelLabels( ...
        anatomicalLocs, selectedChans, nSelectedChannels)
    labels = convertLabelsToString(anatomicalLocs);
    labels = labels(:);
    if numel(labels) >= max(selectedChans)
        labels = labels(selectedChans);
    elseif numel(labels) ~= nSelectedChannels
        error(['Cannot map anatomicalLocs to selected channels. ', ...
            'Labels=%d, selected channels=%d, maximum channel=%d.'], ...
            numel(labels), nSelectedChannels, max(selectedChans));
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
    normalized = lower(strip(string(labels)));
    normalized = regexprep(normalized, '[_\-]+', ' ');
    normalized = regexprep(normalized, '\s+', ' ');
    compact = regexprep(normalized, '[^a-z0-9]', '');
    isNaC = contains(normalized, "nac") | ...
        contains(normalized, "nucleus accumbens") | ...
        contains(normalized, "accumbens");
    isWhiteMatter = contains(normalized, "white matter") | ...
        contains(compact, "whitematter") | startsWith(compact, "wm");
    isLateralVentricle = contains(normalized, "lateral ventricle") | ...
        contains(compact, "lateralventricle");
    excludedRows = isNaC | isWhiteMatter | isLateralVentricle;
end

function createSignificantAreaBarPlot(T, analysisNames, outputFile)

    panelColors = [ ...
        0.847 0.333 0.153; ... % IT: orange
        0.204 0.459 0.702; ... % RT: blue
        0.250 0.600 0.250];    % BR: green

    fig = figure( ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 1800 800]);

    layout = tiledlayout(fig, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    for aa = 1:numel(analysisNames)

        ax = nexttile(layout, aa);
        outcomeRows = T.analysis == analysisNames(aa);
        outcomeTable = T(outcomeRows, :);

        if isempty(outcomeTable)
            axis(ax, 'off');
            text(ax, 0.5, 0.5, ...
                sprintf('No FDR-significant %s areas', analysisNames(aa)), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 12, ...
                'FontWeight', 'bold');
            title(ax, analysisNames(aa), ...
                'FontSize', 14, 'FontWeight', 'bold');
            continue;
        end

        [areaGroup, areaNames] = findgroups(outcomeTable.anatomicalArea);
        totalChannels = splitapply(@sum, ...
            outcomeTable.nSignificantChannels, areaGroup);
        participants = splitapply(@numel, ...
            outcomeTable.patientID, areaGroup);
        pValueFDR = splitapply(@(x) x(1), ...
            outcomeTable.pValueFDR, areaGroup);

        summaryTable = table( ...
            areaNames, totalChannels, participants, pValueFDR, ...
            'VariableNames', { ...
                'anatomicalArea', 'totalChannels', ...
                'nParticipants', 'pValueFDR'});

        summaryTable = sortrows(summaryTable, ...
            {'totalChannels', 'anatomicalArea'}, ...
            {'ascend', 'ascend'});

        y = 1:height(summaryTable);
        barh(ax, y, summaryTable.totalChannels, ...
            0.72, ...
            'FaceColor', panelColors(aa, :), ...
            'EdgeColor', 'none');

        ax.YTick = y;
        ax.YTickLabel = summaryTable.anatomicalArea;
        ax.TickLabelInterpreter = 'none';
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
        ax.Box = 'off';


        xlabel(ax, 'significant IED channels across participants', ...
            'FontSize', 12, 'FontWeight', 'bold');
        title(ax, sprintf('%s significant areas', analysisNames(aa)), ...
            'FontSize', 14, 'FontWeight', 'bold');

        xPadding = max(0.5, 0.02 * max(summaryTable.totalChannels));
        xlim(ax, [0, max(summaryTable.totalChannels) + ...
            max(1, 7 * xPadding)]);

        for rr = 1:height(summaryTable)
            label = sprintf('%d channels; n=%d', ...
                summaryTable.totalChannels(rr), ...
                summaryTable.nParticipants(rr));
            text(ax, summaryTable.totalChannels(rr) + xPadding, rr, ...
                label, ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, ...
                'FontWeight', 'bold');
        end
    end

    title(layout, ...
        'significant anatomical areas', ...
        'FontSize', 16, ...
        'FontWeight', 'bold');

    exportgraphics(fig, outputFile, 'ContentType', 'vector');
    close(fig);
end