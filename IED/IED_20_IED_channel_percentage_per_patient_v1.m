% CCN2026
% Percentage of CHANNELS with IEDs (per patient, per epoch)
% Computed as:
%   mean( number of channels with >=1 IED per trial ) / total channels * 100
%
% AUTHOR: Nill

clear;
clc;
close all;
warning('off','all');

%% ===================== loading neural and event data =====================

inputFolderName_IEDtrials = ...
    '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';

fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));

outputFolderName = ...
    '\\155.100.91.44\d\Code\Nill\BART\IED\IED_20_IED_channel_percentage_per_patient_v1\';

PatientsNum = length(fileList);

% 4 epochs × participants
IED_channel_percentage_vector = nan(4, PatientsNum);

xTicks = strings(1, PatientsNum);

% -------- epoch indices --------
PostOnset    = 1;
PreResponse  = 2;
PostResponse = 3;
PreOutcome   = 4;

%% ===================== compute % CHANNELS with IEDs =====================

for pt = 1:PatientsNum

    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    xTicks(pt) = ptID;

    disp("patient: " + ptID);

    LFPIEDfile = fullfile(inputFolderName_IEDtrials, [ptID '.LFPIED.mat']);
    load(LFPIEDfile);   % loads struct LFPIED

    % =========================================================
    % PostOnset
    % =========================================================
    IEDmat = LFPIED.IEDtrialsPostOnset;   % channels × trials
    chans_with_IED_per_trial = sum(IEDmat > 0, 1);
    IED_channel_percentage_vector(PostOnset, pt) = ...
        mean(chans_with_IED_per_trial) / size(IEDmat,1) * 100;

    % =========================================================
    % PreResponse
    % =========================================================
    IEDmat = LFPIED.IEDtrialsPreResponse;
    chans_with_IED_per_trial = sum(IEDmat > 0, 1);
    IED_channel_percentage_vector(PreResponse, pt) = ...
        mean(chans_with_IED_per_trial) / size(IEDmat,1) * 100;

    % =========================================================
    % PostResponse
    % =========================================================
    IEDmat = LFPIED.IEDtrialsPostResponse;
    chans_with_IED_per_trial = sum(IEDmat > 0, 1);
    IED_channel_percentage_vector(PostResponse, pt) = ...
        mean(chans_with_IED_per_trial) / size(IEDmat,1) * 100;

    % =========================================================
    % PreOutcome
    % =========================================================
    IEDmat = LFPIED.IEDtrialsPreOutcome;
    chans_with_IED_per_trial = sum(IEDmat > 0, 1);
    IED_channel_percentage_vector(PreOutcome, pt) = ...
        mean(chans_with_IED_per_trial) / size(IEDmat,1) * 100;

end

%% ===================== SCATTER + BOXPLOT =====================

epoch_idx   = [PostOnset PreResponse PostResponse PreOutcome];
epoch_names = {'PostOnset','PreResponse','PostResponse','PreOutcome'};

epoch_colors = [
    0.20 0.45 0.70
    0.30 0.70 0.40
    0.85 0.33 0.10
    0.60 0.40 0.80
];

assert(length(epoch_idx) == size(epoch_colors,1), ...
    'Number of epoch colors must match number of epochs');

data  = [];
group = [];

for i = 1:length(epoch_idx)

    temp = squeeze(IED_channel_percentage_vector(epoch_idx(i), :));
    temp = temp(~isnan(temp));

    data  = [data temp];
    group = [group i * ones(1, length(temp))];

end

hFig = figure( ...
    'Units','normalized', ...
    'Position',[0.25 0.2 0.45 0.6], ...
    'Color','w');
hold on;

% ----- boxplot -----
boxplot(data, group, ...
    'Colors','k', ...
    'Symbol','', ...
    'Widths',0.5);

% ----- scatter overlay -----
for i = 1:length(epoch_idx)

    temp = squeeze(IED_channel_percentage_vector(epoch_idx(i), :));
    temp = temp(~isnan(temp));

    xj = i + 0.15 * (rand(size(temp)) - 0.5);

    scatter(xj, temp, ...
        60, ...
        epoch_colors(i,:), ...
        'filled', ...
        'MarkerFaceAlpha', 0.4, ...
        'MarkerEdgeColor', 'none');
end

set(gca, ...
    'XTick', 1:length(epoch_idx), ...
    'XTickLabel', epoch_names, ...
    'FontSize', 14, ...
    'TickDir', 'out', ...
    'PlotBoxAspectRatio', [1 1 1]);

ylabel('% of channels with IEDs','FontSize',16);
title('Percentage of channels exhibiting IEDs','FontSize',22);

box off;

outFile = fullfile(outputFolderName, ...
    'IED_channel_percentage_scatter_boxplot.pdf');

exportgraphics(hFig, outFile, 'ContentType','vector');


%% debug


