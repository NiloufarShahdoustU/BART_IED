% Scatter-boxplot of:
% 1) Mean RT per participant
% 2) Mean IT per participant
% 3) Mean BR (bank rate) per participant
% 4) Mean PR (pop rate) per participant

% Only non-control trials

% Author: Nill

clear;
clc;
close all;

inputFolderName_LFPIED = 'D:\Nill\data\BART\0_0_new_IED_new_area_labels\IED1_find_number_of_IEDs\';
outputFolderName = 'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED5_meanRT_meanIT_BR_PR\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

% Colors
colorRT  = [0.204 0.459 0.702];   % blue
colorIT  = [0.847 0.333 0.153];   % orange
colorBR  = [0.250 0.600 0.250];   % green
colorPop = [0.494 0.184 0.556];   % purple

ptIDs = {};

meanRT_control = [];
meanIT_control = [];
meanBR_control = [];
meanPR_control = [];

nTrialsUsed_RT = [];
nTrialsUsed_IT = [];
nTrialsUsed_BR = [];
nTrialsUsed_PR = [];

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    ptID = fileNameParts{1};

    disp(' ');
    disp(['Processing patient ID: ' ptID]);

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName));
    LFPIED = loadedData.LFPIED;

    if ~isfield(LFPIED, 'BankedTrials')
        warning(['BankedTrials not found for patient ' ptID ...
            '. Skipping this participant.']);
        continue;
    end

    nTrials = LFPIED.nTrials;

    RTs = LFPIED.RTs(:);
    ITs = LFPIED.ITs(:);
    isControl = LFPIED.isControl(:);
    BRs = double(LFPIED.BankedTrials(:));

    % PoppedTrials is the opposite of BankedTrials
    PRs = 1 - BRs;

    minLen = min([ ...
        nTrials, ...
        length(RTs), ...
        length(ITs), ...
        length(isControl), ...
        length(BRs)]);

    nTrials   = minLen;
    RTs       = RTs(1:nTrials);
    ITs       = ITs(1:nTrials);
    isControl = isControl(1:nTrials);
    BRs       = BRs(1:nTrials);
    PRs       = PRs(1:nTrials);

    % Only non-control trials
    nonControlTrials = isControl == 0;

    % Remove trials with invalid RT or RT > 20 seconds
    validRT_20sec = isfinite(RTs) & RTs > 0 & RTs <= 20;

    % Valid trials
    validTrials_RT = ...
        nonControlTrials & validRT_20sec;

    validTrials_IT = ...
        nonControlTrials & validRT_20sec & ...
        isfinite(ITs) & ITs > 0;

    validTrials_BR = ...
        nonControlTrials & validRT_20sec & ...
        isfinite(BRs) & (BRs == 0 | BRs == 1);

    validTrials_PR = ...
        nonControlTrials & validRT_20sec & ...
        isfinite(PRs) & (PRs == 0 | PRs == 1);

    % Means per participant
    thisMeanRT = mean(RTs(validTrials_RT), 'omitnan');
    thisMeanIT = mean(ITs(validTrials_IT), 'omitnan');
    thisMeanBR = mean(BRs(validTrials_BR), 'omitnan');
    thisMeanPR = mean(PRs(validTrials_PR), 'omitnan');

    ptIDs{end+1,1} = ptID;

    meanRT_control(end+1,1) = thisMeanRT;
    meanIT_control(end+1,1) = thisMeanIT;
    meanBR_control(end+1,1) = thisMeanBR;
    meanPR_control(end+1,1) = thisMeanPR;

    nTrialsUsed_RT(end+1,1) = sum(validTrials_RT);
    nTrialsUsed_IT(end+1,1) = sum(validTrials_IT);
    nTrialsUsed_BR(end+1,1) = sum(validTrials_BR);
    nTrialsUsed_PR(end+1,1) = sum(validTrials_PR);

end

%% Save summary table
summaryTable = table( ...
    ptIDs, ...
    meanRT_control, ...
    meanIT_control, ...
    meanBR_control, ...
    meanPR_control, ...
    nTrialsUsed_RT, ...
    nTrialsUsed_IT, ...
    nTrialsUsed_BR, ...
    nTrialsUsed_PR, ...
    'VariableNames', { ...
    'PatientID', ...
    'MeanRT_NonControlTrials', ...
    'MeanIT_NonControlTrials', ...
    'MeanBR_NonControlTrials', ...
    'MeanPR_NonControlTrials', ...
    'NTrialsUsed_RT', ...
    'NTrialsUsed_IT', ...
    'NTrialsUsed_BR', ...
    'NTrialsUsed_PR'});

outputCSV = fullfile( ...
    outputFolderName, ...
    'mean_RT_IT_BR_PR_noncontrol_trials_summary.csv');

writetable(summaryTable, outputCSV);

%% Plot
fig = figure('Visible', 'off', 'Color', 'w');
set(fig, 'Units', 'pixels', 'Position', [100 100 1050 800]);

% Create four equal-size axes
leftMargin   = 0.07;
rightMargin  = 0.03;
bottomMargin = 0.18;
topMargin    = 0.10;
gap          = 0.045;

axWidth = ...
    (1 - leftMargin - rightMargin - 3*gap) / 4;

axHeight = ...
    1 - bottomMargin - topMargin;

ax1 = axes( ...
    'Parent', fig, ...
    'Position', ...
    [leftMargin, bottomMargin, axWidth, axHeight]);

ax2 = axes( ...
    'Parent', fig, ...
    'Position', ...
    [leftMargin + axWidth + gap, ...
    bottomMargin, axWidth, axHeight]);

ax3 = axes( ...
    'Parent', fig, ...
    'Position', ...
    [leftMargin + 2*(axWidth + gap), ...
    bottomMargin, axWidth, axHeight]);

ax4 = axes( ...
    'Parent', fig, ...
    'Position', ...
    [leftMargin + 3*(axWidth + gap), ...
    bottomMargin, axWidth, axHeight]);

% Subplot 1: RT
axes(ax1);

plotSingleScatterBox( ...
    meanRT_control, ...
    colorRT, ...
    'RT', ...
    'time (s)', ...
    'mean RT');

% Subplot 2: IT
axes(ax2);

plotSingleScatterBox( ...
    meanIT_control, ...
    colorIT, ...
    'IT', ...
    'time (s)', ...
    'mean IT');

% Subplot 3: BR
axes(ax3);

plotSingleScatterBox( ...
    meanBR_control, ...
    colorBR, ...
    'BR', ...
    'bank rate', ...
    'mean BR');

% Subplot 4: PR
axes(ax4);

plotSingleScatterBox( ...
    meanPR_control, ...
    colorPop, ...
    'PR', ...
    'pop rate', ...
    'mean PR');

outputPDF = fullfile( ...
    outputFolderName, ...
    'mean_RT_IT_BR_PR_noncontrol_trials_scatter_boxplot.pdf');

exportgraphics(fig, outputPDF, 'ContentType', 'vector');

close(fig);

disp(' ');
disp('Done.');
disp(['Saved PDF: ' outputPDF]);
disp(['Saved CSV: ' outputCSV]);

%% Function
function plotSingleScatterBox( ...
    dataVals, dataColor, xLabelText, yLabelText, titleText)

    dataVals = dataVals(:);
    dataVals = dataVals(isfinite(dataVals));

    hold on;

    boxplot(dataVals, ones(size(dataVals)), ...
        'Labels', {xLabelText}, ...
        'Colors', 'k', ...
        'Widths', 0.22, ...
        'Symbol', '');

    set(findobj(gca, 'Tag', 'Box'), ...
        'LineWidth', 0.75);

    set(findobj(gca, 'Tag', 'Median'), ...
        'LineWidth', 1.2);

    set(findobj(gca, 'Tag', 'Whisker'), ...
        'LineWidth', 0.75);

    set(findobj(gca, 'Tag', 'Upper Whisker'), ...
        'LineWidth', 0.75);

    set(findobj(gca, 'Tag', 'Lower Whisker'), ...
        'LineWidth', 0.75);

    set(findobj(gca, 'Tag', 'Upper Adjacent Value'), ...
        'LineWidth', 0.75);

    set(findobj(gca, 'Tag', 'Lower Adjacent Value'), ...
        'LineWidth', 0.75);

    % Jittered scatter points
    rng(1);
    jitterAmount = 0.055;
    xVals = 1 + jitterAmount * randn(size(dataVals));

    scatter(xVals, dataVals, 45, ...
        'MarkerFaceColor', dataColor, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.5, ...
        'LineWidth', 0.4);

    ylabel(yLabelText);
    title(titleText);

    xlim([0.65 1.35]);

    % Do not force y-axis to start at zero
    if ~isempty(dataVals)

        yMin = min(dataVals);
        yMax = max(dataVals);
        yRange = yMax - yMin;

        if yRange == 0
            yRange = max(abs(yMax), 1) * 0.1;
        end

        ylim([ ...
            yMin - 0.12*yRange, ...
            yMax + 0.12*yRange]);
    end

    grid off;
    box off;

    hold off;

end