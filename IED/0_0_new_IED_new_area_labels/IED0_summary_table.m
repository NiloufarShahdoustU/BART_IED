%% create regional summary table

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

%% empty tables

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
            ~isfield(LFPIED, 'anatomicalLocs')

        disp('Missing selectedChans or anatomicalLocs, skipped');
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

    %% regions

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

    %% remove bad/unknown regions

    keepChan = true(nSelectedChans, 1);

    keepChan = keepChan & chanRegions ~= "Unknown";
    keepChan = keepChan & ~contains(lower(chanRegions), "white matter");
    keepChan = keepChan & ~contains(lower(chanRegions), "ventral dc");
    keepChan = keepChan & ~contains(lower(chanRegions), "lateral ventricle");
    keepChan = keepChan & ~contains(lower(chanRegions), "nac");
    keepChan = keepChan & ~contains(lower(chanRegions), "nucleus accumbens");

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

    %% save channel coverage rows

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

        % keep only valid regions
        keepIED = keepChan(localIEDChan);

        localIEDChan = localIEDChan(keepIED);

        if isempty(localIEDChan)
            continue;
        end

        tempIEDTable = table( ...
            repmat(ptID, length(localIEDChan), 1), ...
            chanRegions(localIEDChan), ...
            localIEDChan(:), ...
            selectedChans(localIEDChan), ...
            'VariableNames', { ...
            'ParticipantID', ...
            'Region', ...
            'LocalChan', ...
            'OriginalChan'});

        IEDTable = [IEDTable; tempIEDTable];

    end

end

%% make regional summary table

regions = unique(channelTable.Region);
regions = sort(regions);

regionalSummary = table( ...
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

    regionalSummary.Region(rr) = thisRegion;

    regionalSummary.TotalParticipantsWithElectrodes(rr) = ...
        length(unique(regionChannelTable.ParticipantID));

    participantElectrodeKey = ...
        regionChannelTable.ParticipantID + "_" + ...
        regionChannelTable.ElectrodeName;

    regionalSummary.TotalElectrodes(rr) = ...
        length(unique(participantElectrodeKey));

    regionalSummary.TotalChannels(rr) = ...
        height(regionChannelTable);

    regionalSummary.TotalIEDs(rr) = ...
        sum(IEDTable.Region == thisRegion);

end

%% sort based on decreasing number of participants

regionalSummary = sortrows( ...
    regionalSummary, ...
    {'TotalParticipantsWithElectrodes', 'TotalElectrodes', 'TotalChannels', 'TotalIEDs'}, ...
    {'descend', 'descend', 'descend', 'descend'});

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
disp('Regional summary table:');
disp(regionalSummary);

disp(' ');
disp('nTrials summary:');
disp(nTrialsSummary);

fprintf('\nnTrials = %.2f +/- %.2f\n', meanNTrials, stdNTrials);

%% save only these two csv files

writetable( ...
    regionalSummary, ...
    fullfile(outputFolderName, 'regional_summary_table.csv'));

writetable( ...
    nTrialsSummary, ...
    fullfile(outputFolderName, 'nTrials_summary.csv'));

disp(' ');
disp('Saved only:');
disp(fullfile(outputFolderName, 'regional_summary_table.csv'));
disp(fullfile(outputFolderName, 'nTrials_summary.csv'));

%% local functions

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

        labels = regexprep(labels, '^(Left|Right)\s+', '', 'ignorecase');
        labels = regexprep(labels, '\s+(Left|Right)$', '', 'ignorecase');

        labels = regexprep(labels, '^[LR][\s\-_]+', '', 'ignorecase');
        labels = regexprep(labels, '[\s\-_]+[LR]$', '', 'ignorecase');

    end

    labels = strip(labels);

    labels(strlength(labels) == 0) = "Unknown";

end

function electrodeNames = get_electrode_names(chanLabels)

    chanLabels = string(chanLabels);
    chanLabels = strip(chanLabels);

    electrodeNames = regexprep(chanLabels, '[\s\-_]*\d+[A-Za-z]?$', '');
    electrodeNames = strip(electrodeNames);

    emptyRows = strlength(electrodeNames) == 0;
    electrodeNames(emptyRows) = chanLabels(emptyRows);

end