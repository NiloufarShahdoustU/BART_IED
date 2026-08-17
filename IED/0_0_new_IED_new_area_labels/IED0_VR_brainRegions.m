%% Brain Regions for VR Project %% only for those that did the full sessions

% Author: RC - 04/02/2026
% new edits: Niloufar Shahdoust 08/2026


%%

clc;
clear;
close all;

% instead of stuff above, I'm gonna read pts in the folder below so that we
% only read the data from the patients in our IED algorithm:
inputFolder = 'D:\Nill\data\BART\bhvStruct_Nill_made';
fileList = dir(fullfile(inputFolder, '*.bhvStruct.mat'));

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED0_VR_brainRegions\';
if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end
electrodePath = 'D:\Nill\data\BART_preprocessed\';


% (2) Build channel-indexed electrode/atlas matrices for every pt.
% The row index in these matrices is the recording channel number used in Step 01.
% ChannelMap2.mat is preferred when it exists because it includes newer-system
% channel mappings; otherwise ChannelMap.mat is used.
% Nill's update: I found some patients that ChannelMap2.mat existed but it
% seemed not to be correct! for example pt 202003 has only 132 channels,
% but in the ChannelMap2.mat there are channel numbers more than 132 which
% is till 137!! so further in the code I changed the ChannelMap2.mat to
% ChannelMap.mat


nPts = length(fileList);
Brainnetome_Atlas_All_tmp2 = cell(0, nPts);
allTrodeLabelsNeuro = cell(1, nPts);
all_isECoG = cell(1, nPts);
neurologistStructArray = struct([]);


%%
for p = 1:length(fileList)

    fileName = fileList(p).name;
    fileNameParts = strsplit(fileName, '.');

    ptID = fileNameParts{1};
    registeredDir = fullfile(electrodePath, ptID, 'Imaging', 'Registered');
    disp(registeredDir);
    channelMapFile2 = fullfile(registeredDir, 'ChannelMap.mat');
    channelMapFile1 = fullfile(registeredDir, 'ChannelMap2.mat');
    electrodesFile = fullfile(registeredDir, 'Electrodes.mat');

    if exist(channelMapFile2, 'file')
        channelMapData = load(channelMapFile2);
    elseif exist(channelMapFile1, 'file')
        channelMapData = load(channelMapFile1);
    else
        error('Could not find ChannelMap.mat or ChannelMap2.mat for %s.', ptID);
    end

    if exist(electrodesFile, 'file')
        electrodeData = load(electrodesFile, 'ElecMapRaw');
    else
        electrodeData = struct();
    end

    fields = fieldnames(channelMapData);
    for i = 1:numel(fields)
        neurologistStructArray(p).(fields{i}) = channelMapData.(fields{i});
    end
    if isfield(electrodeData, 'ElecMapRaw')
        neurologistStructArray(p).ElecMapRaw = electrodeData.ElecMapRaw;
    end

    if isfield(channelMapData, 'ChannelMap2') && any(isfinite(channelMapData.ChannelMap2(:)))
        channelMap = channelMapData.ChannelMap2;
    elseif isfield(channelMapData, 'ChannelMap1') && any(isfinite(channelMapData.ChannelMap1(:)))
        channelMap = channelMapData.ChannelMap1;
    else
        error('No finite channel numbers found in ChannelMap for %s.', ptID);
    end

    if ~isfield(channelMapData, 'LabelMap')
        error('ChannelMap file for %s does not contain LabelMap.', ptID);
    end
    if ~isfield(channelMapData, 'ElecAtlasProj')
        error('ChannelMap file for %s does not contain ElecAtlasProj.', ptID);
    end

    atlasNames = channelMapData.AtlasNames;
    brainnetomeCol = find(strcmp(atlasNames, 'Brainnetome'), 1, 'first');
    if isempty(brainnetomeCol)
        brainnetomeCol = 5;
    end

    channelNums = channelMap(:);
    labelMap = channelMapData.LabelMap(:);
    validMap = isfinite(channelNums) & channelNums > 0 & cellfun(@(x) ~strcmp(char(x), 'NaN'), labelMap);
    validChannels = channelNums(validMap);
    maxChan = max(validChannels);

    if size(Brainnetome_Atlas_All_tmp2, 1) < maxChan
        Brainnetome_Atlas_All_tmp2{maxChan, nPts} = '';
    end

    trodeLabels = repmat({''}, maxChan, 1);
    isECoG = false(maxChan, 1);

    validLinearIdx = find(validMap);
    for ii = 1:numel(validLinearIdx)
        mapIdx = validLinearIdx(ii);
        ch = channelNums(mapIdx);
        elecLabel = char(labelMap{mapIdx});

        trodeLabels{ch} = elecLabel;
        isECoG(ch) = ~startsWith(elecLabel, 'm', 'IgnoreCase', true);

        atlasLabel = '';
        if ch <= size(channelMapData.ElecAtlasProj, 1) && brainnetomeCol <= size(channelMapData.ElecAtlasProj, 2)
            atlasValue = channelMapData.ElecAtlasProj{ch, brainnetomeCol};
            if iscell(atlasValue) && ~isempty(atlasValue)
                atlasLabel = char(atlasValue{1});
            elseif ischar(atlasValue) || isstring(atlasValue)
                atlasLabel = char(atlasValue);
            end
        end
        Brainnetome_Atlas_All_tmp2{ch, p} = atlasLabel;
    end

    allTrodeLabelsNeuro{p} = trodeLabels;
    all_isECoG{p} = isECoG;
end

all_isECoG_logical_matrix = false(size(Brainnetome_Atlas_All_tmp2));
for p = 1:nPts
    current_logical_array = all_isECoG{p};
    all_isECoG_logical_matrix(1:numel(current_logical_array), p) = current_logical_array;
end
%% (5b) Replacing Unknowns with probabilistic location over 15%
probThreshold = 0.10; % updating labels if unknown but others are probable.
UnknownFilledCounts = zeros(1, length(neurologistStructArray));
UnknownRemainingCounts = zeros(1, length(neurologistStructArray));

for i = 1:length(neurologistStructArray) % loop over patients
    atlasNames = neurologistStructArray(i).AtlasNames;
    brainnetomeCol = find(strcmp(atlasNames, 'Brainnetome'), 1, 'first');
    if isempty(brainnetomeCol)
        brainnetomeCol = 5;
    end

    for k = 1:size(Brainnetome_Atlas_All_tmp2, 1) % loop over channel rows
        currentLabel = Brainnetome_Atlas_All_tmp2{k, i};
        if contains(char(currentLabel), 'Unkno') || contains(char(currentLabel), 'Unknown')
            candidateList = {};
            if k <= size(neurologistStructArray(i).ElecAtlasProbProj, 1) && ...
                    brainnetomeCol <= size(neurologistStructArray(i).ElecAtlasProbProj, 2)
                candidateList = neurologistStructArray(i).ElecAtlasProbProj{k, brainnetomeCol};
            end

            if ~isempty(candidateList) && size(candidateList, 2) >= 2
                candidateLabel = candidateList{1, 1};
                candidateProb  = candidateList{1, 2};
                if candidateProb > probThreshold
                    Brainnetome_Atlas_All_tmp2{k, i} = candidateLabel;
                    UnknownFilledCounts(i) = UnknownFilledCounts(i) + 1;
                else
                    UnknownRemainingCounts(i) = UnknownRemainingCounts(i) + 1;
                end
            else
                UnknownRemainingCounts(i) = UnknownRemainingCounts(i) + 1;
            end
        end
    end
end

UnknownFillTable = table((1:length(neurologistStructArray))', ...
    UnknownFilledCounts', UnknownRemainingCounts', ...
    'VariableNames', {'Patient', 'UnknownsFilled', 'UnknownsRemaining'});
disp(UnknownFillTable)

% (4) Cleaning electrode labels
Brainnetome_Atlas_All_wHemi = cell(size(Brainnetome_Atlas_All_tmp2));
Brainnetome_Atlas_All_woHemi = cell(size(Brainnetome_Atlas_All_tmp2));
for ii = 1:numel(Brainnetome_Atlas_All_tmp2)
    label = regexprep(char(Brainnetome_Atlas_All_tmp2{ii}), '[\d"]', '');
    if numel(label) >= 2
        wHemiLabel = label(1:end-2);
    else
        wHemiLabel = '';
    end
    if numel(wHemiLabel) >= 2
        woHemiLabel = wHemiLabel(1:end-2);
    else
        woHemiLabel = '';
    end
    Brainnetome_Atlas_All_wHemi{ii} = wHemiLabel;
    Brainnetome_Atlas_All_woHemi{ii} = woHemiLabel;
end



%% saving important data:
save(fullfile(outputFolderName, 'Brainnetome_Atlas_All_wHemi.mat'), ...
    'Brainnetome_Atlas_All_wHemi');

save(fullfile(outputFolderName, 'Brainnetome_Atlas_All_woHemi.mat'), ...
    'Brainnetome_Atlas_All_woHemi');

