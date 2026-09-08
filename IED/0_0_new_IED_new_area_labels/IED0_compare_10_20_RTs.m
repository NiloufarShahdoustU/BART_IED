

clear;
clc;
close all;

%% paths

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED_new_area_labels\IED1_find_number_of_IEDs\';



fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

[~, sortIndex] = sort({fileList.name});
fileList = fileList(sortIndex);


RTs_10_counts = zeros(1, length(fileList));
RTs_20_counts = zeros(1, length(fileList));

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
    RTs = LFPIED.RTs;

    RTs_10 = RTs(RTs>10);
    RTs_20 = RTs(RTs>20);

    RTs_10_counts(pt) = length(RTs_10);
    RTs_20_counts(pt) = length(RTs_20);
end

%% visualize

ptIDs = strings(1, length(fileList));

for pt = 1:length(fileList)
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptIDs(pt) = string(fileNameParts{1});
end

figure('Color', 'w', 'Position', [100, 100, 1400, 500]);

barData = [RTs_10_counts(:), RTs_20_counts(:)];
b = bar(barData, 'grouped');

b(1).FaceColor = [0.20, 0.55, 0.80];
b(2).FaceColor = [0.90, 0.40, 0.30];
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

xlabel('Participant ID');
ylabel('Number of trials');
title('seeing 10 and 20 sec trials');

xticks(1:length(ptIDs));
xticklabels(ptIDs);
xtickangle(45);

legend({'RT > 10 s', 'RT > 20 s'}, ...
       'Location', 'northwest', ...
       'Box', 'off');

set(gca, ...
    'FontSize', 10, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'LineWidth', 1);

