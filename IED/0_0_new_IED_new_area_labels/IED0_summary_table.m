%% create regional summary tables, with and without hemisphere

clear;
clc;
close all;

%% paths

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED_new_area_labels\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED0_summary_table\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

%% settings

combineLeftRight = true;

IEDFieldNames = { ...
    'IED_occurance_RT', ...
    'IED_occurance_IT'};

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

[~, sortIndex] = sort({fileList.name});
fileList = fileList(sortIndex);

%% empty tables: original regions without hemisphere

channelTable = table( ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames', { ...
    'ParticipantID', ...
    'Region', ...
    'ElectrodeName', ...
    'LocalChan', ...
    'OriginalChan'});

IEDTable = table( ...
    strings(0,1), ...
    strings(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames', { ...
    'ParticipantID', ...
    'Region', ...
    'LocalChan', ...
    'OriginalChan'});

%% empty tables: regions separated by left/right hemisphere

channelTable_wHemisphere = table( ...
    strings(0,1), ...
    strings(0,1), ...
    strings(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames', { ...
    'ParticipantID', ...
    'Region', ...
    'ElectrodeName', ...
    'LocalChan', ...
    'OriginalChan'});

IEDTable_wHemisphere = table( ...
    strings(0,1), ...
    strings(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames', { ...
    'ParticipantID', ...
    'Region', ...
    'LocalChan', ...
    'OriginalChan'});

participantTrials = table( ...
    strings(0,1), ...
    zeros(0,1), ...
    'VariableNames', { ...
    'ParticipantID', ...
    'NTrials'});

%% loop over participants

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    ptID = string(fileNameParts{1});

    disp(' ');
    disp(['Processing patient ID: ' char(ptID)]);

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName));

    if ~isfield(loadedData, 'LFPIED')
        disp('No LFPIED variable found, skipped');
        continue;
    end

    LFPIED = loadedData.LFPIED;

    if ~isfield(LFPIED, 'selectedChans') || ...
            ~isfield(LFPIED, 'anatomicalLocs') || ...
            ~isfield(LFPIED, 'anatomicalLocs_wHemisphere')

        disp(['Missing selectedChans, anatomicalLocs, or ' ...
            'anatomicalLocs_wHemisphere, skipped']);
        continue;
    end

    %% selected channels

    selectedChans = LFPIED.selectedChans;

    if islogical(selectedChans)
        selectedChans = find(selectedChans);
    end

    selectedChans = selectedChans(:);
    selectedChans = double(selectedChans);
    selectedChans = selectedChans(isfinite(selectedChans));
    selectedChans = round(selectedChans);
    selectedChans = selectedChans(selectedChans > 0);
    selectedChans = unique(selectedChans, 'stable');

    nSelectedChans = length(selectedChans);

    if nSelectedChans == 0
        disp('No selected channels, skipped');
        continue;
    end

    %% regions without hemisphere

    anatomicalLocs = convert_to_string(LFPIED.anatomicalLocs);

    if length(anatomicalLocs) == nSelectedChans

        chanRegions = anatomicalLocs;

    elseif max(selectedChans) <= length(anatomicalLocs)

        chanRegions = anatomicalLocs(selectedChans);

    else

        disp('Could not match anatomicalLocs with selectedChans, skipped');
        continue;

    end

    chanRegions = clean_regions(chanRegions, combineLeftRight);

    %% regions with hemisphere

    anatomicalLocs_wHemisphere = ...
        convert_to_string(LFPIED.anatomicalLocs_wHemisphere);

    if length(anatomicalLocs_wHemisphere) == nSelectedChans

        chanRegions_wHemisphere = anatomicalLocs_wHemisphere;

    elseif max(selectedChans) <= length(anatomicalLocs_wHemisphere)

        chanRegions_wHemisphere = ...
            anatomicalLocs_wHemisphere(selectedChans);

    else

        disp(['Could not match anatomicalLocs_wHemisphere with ' ...
            'selectedChans, skipped']);
        continue;

    end

    % This preserves/standardizes labels as area_L or area_R.
    chanRegions_wHemisphere = ...
        clean_regions_with_hemisphere(chanRegions_wHemisphere);

    %% remove bad/unknown regions from original region table

    keepChan = true(nSelectedChans, 1);

    keepChan = keepChan & chanRegions ~= "Unknown";
    keepChan = keepChan & ~contains(lower(chanRegions), "white matter");
    keepChan = keepChan & ~contains(lower(chanRegions), "ventral dc");
    keepChan = keepChan & ~contains(lower(chanRegions), "lateral ventricle");
    keepChan = keepChan & ~contains(lower(chanRegions), "nac");
    keepChan = keepChan & ~contains(lower(chanRegions), "nucleus accumbens");

    %% remove bad/unknown regions from hemisphere-specific table

    keepChan_wHemisphere = true(nSelectedChans, 1);

    keepChan_wHemisphere = keepChan_wHemisphere & ...
        chanRegions_wHemisphere ~= "Unknown";
    keepChan_wHemisphere = keepChan_wHemisphere & ...
        (endsWith(chanRegions_wHemisphere, "_L") | ...
        endsWith(chanRegions_wHemisphere, "_R"));
    keepChan_wHemisphere = keepChan_wHemisphere & ...
        ~contains(lower(chanRegions_wHemisphere), "white matter");
    keepChan_wHemisphere = keepChan_wHemisphere & ...
        ~contains(lower(chanRegions_wHemisphere), "ventral dc");
    keepChan_wHemisphere = keepChan_wHemisphere & ...
        ~contains(lower(chanRegions_wHemisphere), "lateral ventricle");
    keepChan_wHemisphere = keepChan_wHemisphere & ...
        ~contains(lower(chanRegions_wHemisphere), "nac");
    keepChan_wHemisphere = keepChan_wHemisphere & ...
        ~contains(lower(chanRegions_wHemisphere), "nucleus accumbens");

    %% electrode names

    if isfield(LFPIED, 'trodeLabels')

        trodeLabels = convert_to_string(LFPIED.trodeLabels);

        if length(trodeLabels) == nSelectedChans

            chanLabels = trodeLabels;

        elseif max(selectedChans) <= length(trodeLabels)

            chanLabels = trodeLabels(selectedChans);

        else

            chanLabels = "Chan" + string(selectedChans);

        end

        electrodeNames = get_electrode_names(chanLabels);

    else

        % if electrode labels do not exist, each channel is counted as
        % one electrode
        electrodeNames = "Chan" + string(selectedChans);

    end

    %% save channel coverage rows without hemisphere

    localChanIndex = find(keepChan);

    if ~isempty(localChanIndex)

        tempChannelTable = table( ...
            repmat(ptID, length(localChanIndex), 1), ...
            chanRegions(localChanIndex), ...
            electrodeNames(localChanIndex), ...
            localChanIndex(:), ...
            selectedChans(localChanIndex), ...
            'VariableNames', { ...
            'ParticipantID', ...
            'Region', ...
            'ElectrodeName', ...
            'LocalChan', ...
            'OriginalChan'});

        channelTable = [channelTable; tempChannelTable];

    end

    %% save channel coverage rows with hemisphere

    localChanIndex_wHemisphere = find(keepChan_wHemisphere);

    if ~isempty(localChanIndex_wHemisphere)

        tempChannelTable_wHemisphere = table( ...
            repmat(ptID, length(localChanIndex_wHemisphere), 1), ...
            chanRegions_wHemisphere(localChanIndex_wHemisphere), ...
            electrodeNames(localChanIndex_wHemisphere), ...
            localChanIndex_wHemisphere(:), ...
            selectedChans(localChanIndex_wHemisphere), ...
            'VariableNames', { ...
            'ParticipantID', ...
            'Region', ...
            'ElectrodeName', ...
            'LocalChan', ...
            'OriginalChan'});

        channelTable_wHemisphere = ...
            [channelTable_wHemisphere; tempChannelTable_wHemisphere];

    end

    %% nTrials

    nTrials = NaN;

    if isfield(LFPIED, 'nTrials')
        if isnumeric(LFPIED.nTrials) && isscalar(LFPIED.nTrials)
            nTrials = double(LFPIED.nTrials);
        end
    end

    participantTrials = [participantTrials; ...
        table(ptID, nTrials, ...
        'VariableNames', {'ParticipantID', 'NTrials'})];

    %% count IEDs

    for ff = 1:length(IEDFieldNames)

        fieldName = IEDFieldNames{ff};

        if ~isfield(LFPIED, fieldName)
            continue;
        end

        IEDoccurance = LFPIED.(fieldName);

        if istable(IEDoccurance)
            IEDoccurance = table2array(IEDoccurance);
        end

        if isempty(IEDoccurance) || size(IEDoccurance, 2) < 2
            continue;
        end

        % column 2 is local channel index inside selectedChans
        localIEDChan = round(double(IEDoccurance(:, 2)));

        validRows = ...
            isfinite(localIEDChan) & ...
            localIEDChan >= 1 & ...
            localIEDChan <= nSelectedChans;

        localIEDChan = localIEDChan(validRows);

        if isempty(localIEDChan)
            continue;
        end

        %% IED rows without hemisphere

        keepIED = keepChan(localIEDChan);
        localIEDChan_noHemisphere = localIEDChan(keepIED);

        if ~isempty(localIEDChan_noHemisphere)

            tempIEDTable = table( ...
                repmat(ptID, length(localIEDChan_noHemisphere), 1), ...
                chanRegions(localIEDChan_noHemisphere), ...
                localIEDChan_noHemisphere(:), ...
                selectedChans(localIEDChan_noHemisphere), ...
                'VariableNames', { ...
                'ParticipantID', ...
                'Region', ...
                'LocalChan', ...
                'OriginalChan'});

            IEDTable = [IEDTable; tempIEDTable];

        end

        %% IED rows with hemisphere

        keepIED_wHemisphere = keepChan_wHemisphere(localIEDChan);
        localIEDChan_wHemisphere = ...
            localIEDChan(keepIED_wHemisphere);

        if ~isempty(localIEDChan_wHemisphere)

            tempIEDTable_wHemisphere = table( ...
                repmat(ptID, length(localIEDChan_wHemisphere), 1), ...
                chanRegions_wHemisphere(localIEDChan_wHemisphere), ...
                localIEDChan_wHemisphere(:), ...
                selectedChans(localIEDChan_wHemisphere), ...
                'VariableNames', { ...
                'ParticipantID', ...
                'Region', ...
                'LocalChan', ...
                'OriginalChan'});

            IEDTable_wHemisphere = ...
                [IEDTable_wHemisphere; tempIEDTable_wHemisphere];

        end

    end

end

%% make regional summary table without hemisphere

regionalSummary = make_regional_summary(channelTable, IEDTable);

%% make regional summary table with left/right hemisphere

% The Region column contains labels such as area_L and area_R, so the two
% hemispheres are counted as separate regional rows.
regional_wHemisphere_summary = make_regional_summary( ...
    channelTable_wHemisphere, IEDTable_wHemisphere);

%% nTrials summary

validTrials = participantTrials.NTrials;
validTrials = validTrials(isfinite(validTrials) & validTrials > 0);

nParticipants = length(validTrials);
meanNTrials = mean(validTrials);
stdNTrials = std(validTrials);

nTrialsSummary = table( ...
    nParticipants, ...
    meanNTrials, ...
    stdNTrials, ...
    min(validTrials), ...
    max(validTrials), ...
    'VariableNames', { ...
    'NParticipants', ...
    'MeanNTrials', ...
    'SDNTrials', ...
    'MinNTrials', ...
    'MaxNTrials'});

%% display

disp(' ');
disp('Regional summary table without hemisphere:');
disp(regionalSummary);

disp(' ');
disp('Regional summary table with hemisphere:');
disp(regional_wHemisphere_summary);

disp(' ');
disp('nTrials summary:');
disp(nTrialsSummary);

fprintf('\nnTrials = %.2f +/- %.2f\n', meanNTrials, stdNTrials);

%% save these three csv files

writetable( ...
    regionalSummary, ...
    fullfile(outputFolderName, 'regional_summary_table.csv'));

writetable( ...
    regional_wHemisphere_summary, ...
    fullfile(outputFolderName, ...
    'regional_wHemisphere_summary_table.csv'));

writetable( ...
    nTrialsSummary, ...
    fullfile(outputFolderName, 'nTrials_summary.csv'));

%% create and save 3 x 1 bar-plot PDFs

create_hemisphere_bar_figure( ...
    regionalSummary, ...
    regional_wHemisphere_summary, ...
    'TotalIEDs', ...
    'Number of IEDs', ...
    'IED counts by region and hemisphere', ...
    fullfile(outputFolderName, 'IEDs.pdf'));

create_hemisphere_bar_figure( ...
    regionalSummary, ...
    regional_wHemisphere_summary, ...
    'TotalParticipantsWithElectrodes', ...
    'Number of participants', ...
    'Participant counts by region and hemisphere', ...
    fullfile(outputFolderName, 'pts.pdf'));

create_hemisphere_bar_figure( ...
    regionalSummary, ...
    regional_wHemisphere_summary, ...
    'TotalChannels', ...
    'Number of channels', ...
    'Channel counts by region and hemisphere', ...
    fullfile(outputFolderName, 'chans.pdf'));

disp(' ');
disp('Saved:');
disp(fullfile(outputFolderName, 'regional_summary_table.csv'));
disp(fullfile(outputFolderName, ...
    'regional_wHemisphere_summary_table.csv'));
disp(fullfile(outputFolderName, 'nTrials_summary.csv'));
disp(fullfile(outputFolderName, 'IEDs.pdf'));
disp(fullfile(outputFolderName, 'pts.pdf'));
disp(fullfile(outputFolderName, 'chans.pdf'));

%% local functions

function summaryTable = make_regional_summary(channelTable, IEDTable)

    regions = unique(channelTable.Region);
    regions = sort(regions);

    summaryTable = table( ...
        strings(length(regions),1), ...
        zeros(length(regions),1), ...
        zeros(length(regions),1), ...
        zeros(length(regions),1), ...
        zeros(length(regions),1), ...
        'VariableNames', { ...
        'Region', ...
        'TotalParticipantsWithElectrodes', ...
        'TotalElectrodes', ...
        'TotalChannels', ...
        'TotalIEDs'});

    for rr = 1:length(regions)

        thisRegion = regions(rr);

        regionChannelRows = channelTable.Region == thisRegion;
        regionChannelTable = channelTable(regionChannelRows, :);

        summaryTable.Region(rr) = thisRegion;

        summaryTable.TotalParticipantsWithElectrodes(rr) = ...
            length(unique(regionChannelTable.ParticipantID));

        participantElectrodeKey = ...
            regionChannelTable.ParticipantID + "_" + ...
            regionChannelTable.ElectrodeName;

        summaryTable.TotalElectrodes(rr) = ...
            length(unique(participantElectrodeKey));

        summaryTable.TotalChannels(rr) = ...
            height(regionChannelTable);

        summaryTable.TotalIEDs(rr) = ...
            sum(IEDTable.Region == thisRegion);

    end

    %% sort based on decreasing number of participants

    summaryTable = sortrows( ...
        summaryTable, ...
        {'TotalParticipantsWithElectrodes', ...
        'TotalElectrodes', 'TotalChannels', 'TotalIEDs'}, ...
        {'descend', 'descend', 'descend', 'descend'});

end

function labels = convert_to_string(rawLabels)

    if isstring(rawLabels)

        labels = rawLabels(:);

    elseif iscell(rawLabels)

        labels = strings(numel(rawLabels), 1);

        for ii = 1:numel(rawLabels)

            if isempty(rawLabels{ii})
                labels(ii) = "Unknown";
            else
                labels(ii) = string(rawLabels{ii});
            end

        end

    elseif iscategorical(rawLabels)

        labels = string(rawLabels(:));

    elseif ischar(rawLabels)

        labels = string(cellstr(rawLabels));

    else

        labels = string(rawLabels(:));

    end

    labels = strip(labels);

    badRows = ...
        ismissing(labels) | ...
        strlength(labels) == 0 | ...
        lower(labels) == "nan" | ...
        lower(labels) == "<missing>";

    labels(badRows) = "Unknown";

end

function labels = clean_regions(labels, combineLeftRight)

    labels = string(labels);
    labels = strip(labels);

    labels = regexprep(labels, '_', ' ');
    labels = regexprep(labels, '\s+', ' ');

    if combineLeftRight

        labels = regexprep(labels, ...
            '^(Left|Right)\s+', '', 'ignorecase');
        labels = regexprep(labels, ...
            '\s+(Left|Right)$', '', 'ignorecase');

        labels = regexprep(labels, ...
            '^[LR][\s\-_]+', '', 'ignorecase');
        labels = regexprep(labels, ...
            '[\s\-_]+[LR]$', '', 'ignorecase');

    end

    labels = strip(labels);

    labels(strlength(labels) == 0) = "Unknown";

end

function labels = clean_regions_with_hemisphere(labels)

    % Convert possible formats such as area_L, area L, L_area,
    % Left area, area_Left, and area Left into area_L or area_R.

    labels = string(labels);
    labels = strip(labels);

    labels = regexprep(labels, '_', ' ');
    labels = regexprep(labels, '-', ' ');
    labels = regexprep(labels, '\s+', ' ');

    lowerLabels = lower(labels);

    isLeft = ...
        startsWith(lowerLabels, "left ") | ...
        startsWith(lowerLabels, "l ") | ...
        endsWith(lowerLabels, " left") | ...
        endsWith(lowerLabels, " l");

    isRight = ...
        startsWith(lowerLabels, "right ") | ...
        startsWith(lowerLabels, "r ") | ...
        endsWith(lowerLabels, " right") | ...
        endsWith(lowerLabels, " r");

    baseLabels = regexprep(labels, ...
        '^(Left|Right|L|R)\s+', '', 'ignorecase');
    baseLabels = regexprep(baseLabels, ...
        '\s+(Left|Right|L|R)$', '', 'ignorecase');
    baseLabels = strip(baseLabels);

    labels = strings(size(baseLabels));
    labels(isLeft) = baseLabels(isLeft) + "_L";
    labels(isRight) = baseLabels(isRight) + "_R";

    unknownRows = ...
        (~isLeft & ~isRight) | ...
        strlength(baseLabels) == 0 | ...
        lower(baseLabels) == "unknown";

    labels(unknownRows) = "Unknown";

end

function electrodeNames = get_electrode_names(chanLabels)

    chanLabels = string(chanLabels);
    chanLabels = strip(chanLabels);

    electrodeNames = regexprep(chanLabels, ...
        '[\s\-_]*\d+[A-Za-z]?$', '');
    electrodeNames = strip(electrodeNames);

    emptyRows = strlength(electrodeNames) == 0;
    electrodeNames(emptyRows) = chanLabels(emptyRows);

end

function create_hemisphere_bar_figure( ...
        regionalSummary, regionalSummary_wHemisphere, ...
        metricField, yAxisLabel, mainTitle, outputFile)

    %% Sort regions using the TOTAL counts only

    totalValuesUnsorted = regionalSummary.(metricField);

    [totalValues, sortIndex] = sort( ...
        totalValuesUnsorted, 'descend');

    orderedRegions = regionalSummary.Region(sortIndex);
    nRegions = length(orderedRegions);

    %% Match left and right values to the TOTAL region order

    % NaN is intentional: when a left or right region does not exist,
    % MATLAB leaves that bar position empty.
    leftValues = NaN(nRegions, 1);
    rightValues = NaN(nRegions, 1);

    hemisphereMetricValues = ...
        regionalSummary_wHemisphere.(metricField);

    for rr = 1:nRegions

        leftRegion = orderedRegions(rr) + "_L";
        rightRegion = orderedRegions(rr) + "_R";

        leftRow = find( ...
            regionalSummary_wHemisphere.Region == leftRegion, 1);

        rightRow = find( ...
            regionalSummary_wHemisphere.Region == rightRegion, 1);

        if ~isempty(leftRow)
            leftValues(rr) = hemisphereMetricValues(leftRow);
        end

        if ~isempty(rightRow)
            rightValues(rr) = hemisphereMetricValues(rightRow);
        end

    end

    %% Set each subplot y-axis limit based on its own tallest bar

    subplotValues = {totalValues, leftValues, rightValues};
    yMaximum = ones(3, 1);

    for pp = 1:3

        finiteValues = subplotValues{pp}( ...
            isfinite(subplotValues{pp}));

        if ~isempty(finiteValues) && max(finiteValues) > 0
            largestValue = max(finiteValues);
            yMaximum(pp) = ...
                largestValue + max(1, 0.12 * largestValue);
        end

    end

    %% Figure

    figureHandle = figure( ...
        'Visible', 'off', ...
        'Color', 'w', ...
        'Position', [50 50 2400 1700]);

    % Total
    subplot(3, 1, 1);
    plot_count_bars( ...
        orderedRegions, totalValues, ...
        'Total', yAxisLabel, yMaximum(1), ...
        [0.35 0.35 0.35]);

    % Left hemisphere
    subplot(3, 1, 2);
    plot_count_bars( ...
        orderedRegions, leftValues, ...
        'Left hemisphere', yAxisLabel, yMaximum(2), ...
        [0.204 0.459 0.702]);

    % Right hemisphere
    subplot(3, 1, 3);
    plot_count_bars( ...
        orderedRegions, rightValues, ...
        'Right hemisphere', yAxisLabel, yMaximum(3), ...
        [0.847 0.333 0.153]);

    sgtitle(mainTitle, ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');

    exportgraphics( ...
        figureHandle, ...
        outputFile, ...
        'ContentType', 'vector');

    close(figureHandle);

end

function plot_count_bars( ...
        orderedRegions, values, panelTitle, ...
        yAxisLabel, yMaximum, barColor)

    nRegions = length(orderedRegions);
    xPositions = (1:nRegions)';

    bar( ...
        xPositions, ...
        values, ...
        0.80, ...
        'FaceColor', barColor, ...
        'EdgeColor', 'none');

    hold on;

    %% Put the count above every existing bar

    validBars = isfinite(values);
    textOffset = 0.015 * yMaximum;

    text( ...
        xPositions(validBars), ...
        values(validBars) + textOffset, ...
        compose('%.0f', values(validBars)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 7, ...
        'Color', [0.1 0.1 0.1], ...
        'Clipping', 'off');

    hold off;

    xlim([0.25, nRegions + 0.75]);
    ylim([0, yMaximum]);

    xticks(xPositions);
    xticklabels(orderedRegions);
    xtickangle(0);

    ylabel(yAxisLabel, 'FontSize', 11);
    title(panelTitle, 'FontSize', 13, 'FontWeight', 'bold');

    set(gca, ...
        'FontSize', 7, ...
        'TickLabelInterpreter', 'none', ...
        'TickDir', 'out', ...
        'TickLength', [0.003 0.003], ...
        'Box', 'off', ...
        'LineWidth', 1);

    % grid on;
    % set(gca, 'XGrid', 'off', 'YGrid', 'on');

end