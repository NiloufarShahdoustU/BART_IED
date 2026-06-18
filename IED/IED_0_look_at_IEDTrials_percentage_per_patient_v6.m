% ===================== CCN 2026 =====================
% AUTHOR: Nill
% PURPOSE: Compute and plot IED trial percentages
% NOTE: PreOnset1 & PreOnset2 REMOVED COMPLETELY
% ADDED: Print ptIDs where IED% < 10 in ANY of the 4 epochs (per patient)

clear;
clc;
close all;
warning('off','all');

%% ===================== loading neural and event data =====================

inputFolderName_IEDtrials = ...
    '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';

fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));

outputFolderName = ...
    '\\155.100.91.44\d\Code\Nill\BART\IED\IED_0_outputs\v6\';

PatientsNum = length(fileList);

% ONLY 4 epochs (no PreOnsets)
IED_percentage_vector = nan(4, PatientsNum);

% Epoch indices
PostOnset    = 1;
PreResponse  = 2;
PostResponse = 3;
PreOutcome   = 4;

epoch_idx   = 1:4;
epoch_names = {'Post-Onset','Pre-Response','Post-Response','Pre-Outcome'};

epoch_colors = [
    0.20 0.45 0.70
    0.30 0.70 0.40
    0.85 0.33 0.10
    0.60 0.40 0.80
];

xTicks = strings(1, PatientsNum);

%% ===================== compute IED percentages =====================

for pt = 1:PatientsNum

    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    xTicks(pt) = ptID;

    disp("patient: " + ptID);

    LFPIEDfile = fullfile(inputFolderName_IEDtrials, [ptID '.LFPIED.mat']);
    load(LFPIEDfile);

    % -------- Post-Onset --------
    IEDtrialsAcross = any(LFPIED.IEDtrialsPostOnset, 1);
    IED_percentage_vector(PostOnset, pt) = ...
        (nansum(IEDtrialsAcross) / length(IEDtrialsAcross)) * 100;

    % -------- Pre-Response --------
    IEDtrialsAcross = any(LFPIED.IEDtrialsPreResponse, 1);
    IED_percentage_vector(PreResponse, pt) = ...
        (nansum(IEDtrialsAcross) / length(IEDtrialsAcross)) * 100;

    % -------- Post-Response --------
    IEDtrialsAcross = any(LFPIED.IEDtrialsPostResponse, 1);
    IED_percentage_vector(PostResponse, pt) = ...
        (nansum(IEDtrialsAcross) / length(IEDtrialsAcross)) * 100;

    % -------- Pre-Outcome --------
    IEDtrialsAcross = any(LFPIED.IEDtrialsPreOutcome, 1);
    IED_percentage_vector(PreOutcome, pt) = ...
        (nansum(IEDtrialsAcross) / length(IEDtrialsAcross)) * 100;
end

%% ===================== print ptIDs with IED% < 10 in ANY epoch =====================

thresholdPct = 10;

% mask: epoch x patient
belowMask = IED_percentage_vector < thresholdPct;

% patients where ANY of the 4 epochs is below threshold
patientsBelowAny = any(belowMask, 1);

ptIDs_belowAny = xTicks(patientsBelowAny);

disp("==============================================================");
disp("Patients with IED trials (%) < " + thresholdPct + " in ANY epoch:");
disp("==============================================================");

if isempty(ptIDs_belowAny)
    disp("None");
else
    % Print per patient + which epochs are below threshold + values
    for k = 1:numel(ptIDs_belowAny)
        thisID = ptIDs_belowAny(k);

        % find the patient column index
        ptCol = find(xTicks == thisID, 1, 'first');

        % epochs below threshold for this patient
        epochHits = find(belowMask(:, ptCol));

        % build a readable string
        s = "  " + thisID + "  |  ";

        for e = 1:numel(epochHits)
            ep = epochHits(e);
            val = IED_percentage_vector(ep, ptCol);
            s = s + epoch_names{ep} + ": " + string(val) + "%";
            if e < numel(epochHits)
                s = s + ", ";
            end
        end

        disp(s);
    end
end

disp("==============================================================");

%% ===================== prepare boxplot data =====================

data  = [];
group = [];

for i = 1:4
    temp = IED_percentage_vector(i, :);
    temp = temp(~isnan(temp));

    data  = [data, temp];
    group = [group, i * ones(1, length(temp))];
end

hFig = figure( ...
    'Units','normalized', ...
    'Position',[0.25 0.2 0.45 0.6], ...
    'Color','w');
hold on;

% ---- Boxplot ----
boxplot(data, group, ...
    'Colors','k', ...
    'Symbol','', ...
    'Widths',0.5);

% ---- Scatter overlay ----
for i = 1:length(epoch_idx)

    temp = IED_percentage_vector(epoch_idx(i), :);
    temp = temp(~isnan(temp));

    xj = i + 0.15 * (rand(size(temp)) - 0.5);

    scatter(xj, temp, ...
        60, ...
        epoch_colors(i,:), ...
        'filled', ...
        'MarkerFaceAlpha', 0.4, ...
        'MarkerEdgeColor','none');
end

% ---- Formatting ----
set(gca,'PlotBoxAspectRatio',[1 1 1]);
set(gca,'XTick',1:4, ...
        'XTickLabel',epoch_names, ...
        'FontSize',14);

ylabel('IED trials (%)','FontSize',16);
title('percentage of IED trials','FontSize',22);

box off;
set(gca,'TickDir','out');

outFile = fullfile(outputFolderName, ...
    'IED_percentage_scatter_boxplot_selected_epochs.pdf');

exportgraphics(hFig, outFile, 'ContentType','vector');
