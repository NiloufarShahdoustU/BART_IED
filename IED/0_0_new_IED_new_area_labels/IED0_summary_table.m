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

%% Left-versus-right proportion tests with BH-FDR correction

% Participant denominators are the numbers of unique participants with at
% least one retained channel in the corresponding hemisphere. Channel and
% IED denominators are the total retained channel and IED counts in that
% hemisphere. For every region, the two-proportion z-test compares:
%
%        left regional count / left hemisphere denominator
%    vs  right regional count / right hemisphere denominator
%
% BH-FDR correction is performed separately for participants, channels,
% and IEDs across all plotted regions.

leftChannelRows = endsWith(channelTable_wHemisphere.Region, "_L");
rightChannelRows = endsWith(channelTable_wHemisphere.Region, "_R");
leftIEDRows = endsWith(IEDTable_wHemisphere.Region, "_L");
rightIEDRows = endsWith(IEDTable_wHemisphere.Region, "_R");

participantDenominators = [ ...
    length(unique(channelTable_wHemisphere.ParticipantID(leftChannelRows))), ...
    length(unique(channelTable_wHemisphere.ParticipantID(rightChannelRows)))];

channelDenominators = [sum(leftChannelRows), sum(rightChannelRows)];
IEDDenominators = [sum(leftIEDRows), sum(rightIEDRows)];

participantStats = calculate_hemisphere_proportion_stats( ...
    regionalSummary, regional_wHemisphere_summary, ...
    'TotalParticipantsWithElectrodes', 'Participants', ...
    participantDenominators);

channelStats = calculate_hemisphere_proportion_stats( ...
    regionalSummary, regional_wHemisphere_summary, ...
    'TotalChannels', 'Channels', channelDenominators);

IEDStats = calculate_hemisphere_proportion_stats( ...
    regionalSummary, regional_wHemisphere_summary, ...
    'TotalIEDs', 'IEDs', IEDDenominators);

hemisphereProportionStats = [ ...
    participantStats; channelStats; IEDStats];

writetable( ...
    hemisphereProportionStats, ...
    fullfile(outputFolderName, ...
    'left_right_proportion_tests_BH_FDR.csv'));

%% Create one 6 x 1 figure

% Rows 1-2: participants (overall, then left/right)
% Rows 3-4: channels     (overall, then left/right)
% Rows 5-6: IEDs         (overall, then left/right)

create_six_panel_hemisphere_figure( ...
    regionalSummary, ...
    regional_wHemisphere_summary, ...
    participantStats, channelStats, IEDStats, ...
    fullfile(outputFolderName, ...
    'participants_channels_IEDs_6x1.pdf'));

disp(' ');
disp('Saved:');
disp(fullfile(outputFolderName, 'regional_summary_table.csv'));
disp(fullfile(outputFolderName, ...
    'regional_wHemisphere_summary_table.csv'));
disp(fullfile(outputFolderName, 'nTrials_summary.csv'));
disp(fullfile(outputFolderName, ...
    'left_right_proportion_tests_BH_FDR.csv'));
disp(fullfile(outputFolderName, ...
    'participants_channels_IEDs_6x1.pdf'));

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

function statsTable = calculate_hemisphere_proportion_stats( ...
        regionalSummary, regionalSummary_wHemisphere, ...
        metricField, metricName, denominators)

    %% Test only the regions displayed in the figure

    keepRegion = ...
        regionalSummary.TotalParticipantsWithElectrodes >= 10;
    regions = regionalSummary.Region(keepRegion);
    nRegions = length(regions);

    leftCount = zeros(nRegions, 1);
    rightCount = zeros(nRegions, 1);
    leftDenominator = repmat(denominators(1), nRegions, 1);
    rightDenominator = repmat(denominators(2), nRegions, 1);
    leftProportion = NaN(nRegions, 1);
    rightProportion = NaN(nRegions, 1);
    zStatistic = NaN(nRegions, 1);
    pValue = NaN(nRegions, 1);

    hemisphereValues = regionalSummary_wHemisphere.(metricField);

    for rr = 1:nRegions

        leftRow = find( ...
            regionalSummary_wHemisphere.Region == regions(rr) + "_L", 1);
        rightRow = find( ...
            regionalSummary_wHemisphere.Region == regions(rr) + "_R", 1);

        if ~isempty(leftRow)
            leftCount(rr) = hemisphereValues(leftRow);
        end

        if ~isempty(rightRow)
            rightCount(rr) = hemisphereValues(rightRow);
        end

        nLeft = leftDenominator(rr);
        nRight = rightDenominator(rr);

        if nLeft <= 0 || nRight <= 0
            continue;
        end

        leftProportion(rr) = leftCount(rr) / nLeft;
        rightProportion(rr) = rightCount(rr) / nRight;

        pooledProportion = ...
            (leftCount(rr) + rightCount(rr)) / (nLeft + nRight);

        standardError = sqrt( ...
            pooledProportion * (1 - pooledProportion) * ...
            (1 / nLeft + 1 / nRight));

        if standardError > 0
            zStatistic(rr) = ...
                (leftProportion(rr) - rightProportion(rr)) / ...
                standardError;

            % Two-sided standard-normal p-value. erfc avoids requiring
            % the Statistics and Machine Learning Toolbox for normcdf.
            pValue(rr) = erfc(abs(zStatistic(rr)) / sqrt(2));
        else
            zStatistic(rr) = 0;
            pValue(rr) = 1;
        end

    end

    qValue = benjamini_hochberg(pValue);
    significantFDR = qValue < 0.05;
    metricColumn = repmat(string(metricName), nRegions, 1);

    statsTable = table( ...
        metricColumn, regions, ...
        leftCount, rightCount, ...
        leftDenominator, rightDenominator, ...
        leftProportion, rightProportion, ...
        zStatistic, pValue, qValue, significantFDR, ...
        'VariableNames', { ...
        'Metric', 'Region', ...
        'LeftCount', 'RightCount', ...
        'LeftDenominator', 'RightDenominator', ...
        'LeftProportion', 'RightProportion', ...
        'ZStatistic', 'PValue', 'QValue_BH_FDR', ...
        'Significant_BH_FDR_0_05'});

end

function qValues = benjamini_hochberg(pValues)

    %% Benjamini-Hochberg FDR-adjusted p-values

    qValues = NaN(size(pValues));
    validIndex = find(isfinite(pValues));

    if isempty(validIndex)
        return;
    end

    validP = pValues(validIndex);
    [sortedP, sortIndex] = sort(validP, 'ascend');
    m = length(sortedP);

    adjustedSorted = sortedP .* m ./ (1:m)';

    % Enforce monotonicity from the largest p-value downward.
    for ii = m-1:-1:1
        adjustedSorted(ii) = min( ...
            adjustedSorted(ii), adjustedSorted(ii + 1));
    end

    adjustedSorted = min(adjustedSorted, 1);

    unsortedAdjusted = NaN(m, 1);
    unsortedAdjusted(sortIndex) = adjustedSorted;
    qValues(validIndex) = unsortedAdjusted;

end

function create_six_panel_hemisphere_figure( ...
        regionalSummary, regionalSummary_wHemisphere, ...
        participantStats, channelStats, IEDStats, outputFile)

    %% Keep only regions covered by at least 10 participants in total

    keepRegion = ...
        regionalSummary.TotalParticipantsWithElectrodes >= 10;
    plottingSummary = regionalSummary(keepRegion, :);

    %% One common region order, based on decreasing participant coverage

    plottingSummary = sortrows( ...
        plottingSummary, ...
        {'TotalParticipantsWithElectrodes', ...
        'TotalChannels', 'TotalIEDs'}, ...
        {'descend', 'descend', 'descend'});

    orderedRegions = plottingSummary.Region;

    metricFields = { ...
        'TotalParticipantsWithElectrodes', ...
        'TotalChannels', ...
        'TotalIEDs'};

    metricTitles = { ...
        'Participants', ...
        'Channels', ...
        'IEDs'};

    yAxisLabels = { ...
        'Number of participants', ...
        'Number of channels', ...
        'Number of IEDs'};

    statsTables = {participantStats, channelStats, IEDStats};

    totalColor = [0.35 0.35 0.35];
    leftColor = [0.204 0.459 0.702];
    rightColor = [0.847 0.333 0.153];

    % Requested widths: each left/right bar is half as thick as an
    % overall bar.
    overallBarWidth = 0.42;
    hemisphereBarWidth = overallBarWidth / 2;

    figureHandle = figure( ...
        'Visible', 'off', ...
        'Color', 'w', ...
        'Position', [50 50 700 3300]);

    tiledlayout(6, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    for mm = 1:3

        metricField = metricFields{mm};
        totalValues = get_total_values( ...
            plottingSummary, orderedRegions, metricField);

        [leftValues, rightValues] = get_hemisphere_values( ...
            regionalSummary_wHemisphere, orderedRegions, metricField);

        thisStats = order_stats_table( ...
            statsTables{mm}, orderedRegions);

        %% Overall panel

        nexttile;
        totalMaximum = calculate_y_maximum(totalValues, false);

        plot_overall_bars( ...
            orderedRegions, totalValues, ...
            [metricTitles{mm} ' - overall'], ...
            yAxisLabels{mm}, totalMaximum, ...
            totalColor, overallBarWidth);

        %% Paired left/right panel

        nexttile;
        pairedMaximum = calculate_y_maximum( ...
            [leftValues; rightValues], true);

        plot_paired_hemisphere_bars( ...
            orderedRegions, leftValues, rightValues, ...
            thisStats.QValue_BH_FDR, ...
            thisStats.Significant_BH_FDR_0_05, ...
            [metricTitles{mm} ' - left versus right'], ...
            yAxisLabels{mm}, pairedMaximum, ...
            leftColor, rightColor, hemisphereBarWidth);

    end

    sgtitle( ...
        {'Regional coverage and IED counts', ...
        'Left-right statistics: two-proportion z-tests; ' ...
        'BH-FDR within each metric'}, ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');

    exportgraphics( ...
        figureHandle, outputFile, 'ContentType', 'vector');

    close(figureHandle);

end

function values = get_total_values( ...
        plottingSummary, orderedRegions, metricField)

    values = NaN(length(orderedRegions), 1);

    for rr = 1:length(orderedRegions)
        rowIndex = find( ...
            plottingSummary.Region == orderedRegions(rr), 1);
        if ~isempty(rowIndex)
            values(rr) = plottingSummary.(metricField)(rowIndex);
        end
    end

end


function [leftValues, rightValues] = get_hemisphere_values( ...
        regionalSummary_wHemisphere, orderedRegions, metricField)

    nRegions = length(orderedRegions);
    leftValues = NaN(nRegions, 1);
    rightValues = NaN(nRegions, 1);
    hemisphereValues = regionalSummary_wHemisphere.(metricField);

    for rr = 1:nRegions

        leftRow = find( ...
            regionalSummary_wHemisphere.Region == ...
            orderedRegions(rr) + "_L", 1);
        rightRow = find( ...
            regionalSummary_wHemisphere.Region == ...
            orderedRegions(rr) + "_R", 1);

        if ~isempty(leftRow)
            leftValues(rr) = hemisphereValues(leftRow);
        end

        if ~isempty(rightRow)
            rightValues(rr) = hemisphereValues(rightRow);
        end

    end

end

function orderedStats = order_stats_table(statsTable, orderedRegions)

    orderedStats = statsTable([],:);

    for rr = 1:length(orderedRegions)
        rowIndex = find( ...
            statsTable.Region == orderedRegions(rr), 1);
        if ~isempty(rowIndex)
            orderedStats = [orderedStats; statsTable(rowIndex,:)];
        end
    end

end

function yMaximum = calculate_y_maximum(values, reserveStatisticsSpace)

    finiteValues = values(isfinite(values));

    if isempty(finiteValues) || max(finiteValues) <= 0
        yMaximum = 1;
        return;
    end

    largestValue = max(finiteValues);

    if reserveStatisticsSpace
        yMaximum = largestValue + max(1, 0.22 * largestValue);
    else
        yMaximum = largestValue + max(1, 0.12 * largestValue);
    end

end

function plot_overall_bars( ...
        orderedRegions, values, panelTitle, yAxisLabel, ...
        yMaximum, barColor, barWidth)

    nRegions = length(orderedRegions);
    xPositions = (1:nRegions)';

    bar( ...
        xPositions, values, barWidth, ...
        'FaceColor', barColor, ...
        'EdgeColor', 'none');

    hold on;
    add_count_labels(xPositions, values, yMaximum);
    hold off;

    format_bar_axes( ...
        orderedRegions, panelTitle, yAxisLabel, yMaximum);

end

function plot_paired_hemisphere_bars( ...
        orderedRegions, leftValues, rightValues, ...
        qValues, significantFDR, panelTitle, yAxisLabel, ...
        yMaximum, leftColor, rightColor, barWidth)

    nRegions = length(orderedRegions);
    xPositions = (1:nRegions)';

    % The center-to-center separation keeps the two half-width bars next
    % to each other without overlap.
    xOffset = barWidth / 2;
    leftX = xPositions - xOffset;
    rightX = xPositions + xOffset;

    leftBars = bar( ...
        leftX, leftValues, barWidth, ...
        'FaceColor', leftColor, ...
        'EdgeColor', 'none');

    hold on;

    rightBars = bar( ...
        rightX, rightValues, barWidth, ...
        'FaceColor', rightColor, ...
        'EdgeColor', 'none');

    % add_count_labels(leftX, leftValues, yMaximum);
    % add_count_labels(rightX, rightValues, yMaximum);

    %% Mark significant FDR-adjusted left-right comparisons

    for rr = 1:nRegions

        if significantFDR(rr) && isfinite(qValues(rr))

            pairValues = [leftValues(rr), rightValues(rr)];
            pairValues = pairValues(isfinite(pairValues));

            if isempty(pairValues)
                continue;
            end

            markerY = max(pairValues) + 0.09 * yMaximum;
            bracketHeight = 0.018 * yMaximum;

            plot( ...
                [leftX(rr), leftX(rr), rightX(rr), rightX(rr)], ...
                [markerY - bracketHeight, markerY, ...
                markerY, markerY - bracketHeight], ...
                'k-', 'LineWidth', 0.8);

            if qValues(rr) < 0.001
                significanceLabel = '***';
            elseif qValues(rr) < 0.01
                significanceLabel = '**';
            else
                significanceLabel = '*';
            end

            text( ...
                xPositions(rr), markerY + 0.015 * yMaximum, ...
                significanceLabel, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 8, ...
                'FontWeight', 'bold', ...
                'Color', [0.05 0.05 0.05], ...
                'Clipping', 'off');

        end

    end

    legend( ...
        [leftBars, rightBars], ...
        {'Left', 'Right'}, ...
        'Location', 'northeast', ...
        'Box', 'off', ...
        'FontSize', 8);

    hold off;

    format_bar_axes( ...
        orderedRegions, panelTitle, yAxisLabel, yMaximum);

end

function add_count_labels(xPositions, values, yMaximum)

    validBars = isfinite(values);
    textOffset = 0.012 * yMaximum;

    text( ...
        xPositions(validBars), ...
        values(validBars) + textOffset, ...
        compose('%.0f', values(validBars)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 6.5, ...
        'Color', [0.1 0.1 0.1], ...
        'Clipping', 'off');

end

function format_bar_axes( ...
        orderedRegions, panelTitle, yAxisLabel, yMaximum)

    nRegions = length(orderedRegions);
    xPositions = (1:nRegions)';

    xlim([0.45, nRegions + 0.55]);
    ylim([0, yMaximum]);

    xticks(xPositions);
    xticklabels(orderedRegions);
    xtickangle(0);

    ylabel(yAxisLabel, 'FontSize', 10);
    title(panelTitle, 'FontSize', 12, 'FontWeight', 'bold');

    set(gca, ...
        'FontSize', 7, ...
        'TickLabelInterpreter', 'none', ...
        'TickDir', 'out', ...
        'TickLength', [0.003 0.003], ...
        'Box', 'off', ...
        'LineWidth', 1);

end
