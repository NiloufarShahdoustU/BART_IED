% Scatter-boxplot of:
% 1) Percent of IED trials per participant
% 2) Percent of channels showing IEDs per participant
%
% RT = blue
% IT = orange
%
% Only non-control trials
% Trials with RT > 10 seconds are removed from BOTH RT and IT calculations
% No log transform
% No outlier removal
% Each dot = one participant
% Author: Nill

clear;
clc;
close all;

inputFolderName_LFPIED = 'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';
outputFolderName = 'D:\Nill\code\BART\IED\0_0_new_IED\IED4_trials_chans_percentage\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

% Colors
colorRT = [0.204   0.459   0.702];   % blue RT
colorIT = [0.847   0.333   0.153];   % orange IT

ptIDs = {};
pctIEDtrials_RT = [];
pctIEDtrials_IT = [];
pctIEDchans_RT = [];
pctIEDchans_IT = [];

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    ptID = fileNameParts{1};

    disp(' ');
    disp(['Processing patient ID: ' ptID]);

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName));
    LFPIED = loadedData.LFPIED;

    nTrials = LFPIED.nTrials;

    RTs = LFPIED.RTs(:);
    ITs = LFPIED.ITs(:);
    isControl = LFPIED.isControl(:);
    nChans = length(LFPIED.selectedChans);

    minLen = min([nTrials, length(RTs), length(ITs), length(isControl)]);

    nTrials = minLen;
    RTs = RTs(1:nTrials);
    ITs = ITs(1:nTrials);
    isControl = isControl(1:nTrials);

    nonControlTrials = isControl == 0;
    validRT_10sec = isfinite(RTs) & RTs <= 10;

    % Count IEDs per trial
    nIED_RT = countIEDsPerTrial(LFPIED, 'IED_occurance_RT', nTrials);
    nIED_IT = countIEDsPerTrial(LFPIED, 'IED_occurance_IT', nTrials);

    % Denominators for trial percentage
    validTrials_RT = nonControlTrials & validRT_10sec & isfinite(RTs) & RTs > 0;
    validTrials_IT = nonControlTrials & validRT_10sec & isfinite(ITs) & ITs > 0;

    % Percent of trials that had at least one IED
    thisPctIEDtrials_RT = getPercentIEDTrials(nIED_RT, validTrials_RT);
    thisPctIEDtrials_IT = getPercentIEDTrials(nIED_IT, validTrials_IT);

    % Number of unique channels that showed at least one IED
    nIEDchans_RT = countUniqueChannelsWithIED(LFPIED, 'IED_occurance_RT', nTrials, validTrials_RT);
    nIEDchans_IT = countUniqueChannelsWithIED(LFPIED, 'IED_occurance_IT', nTrials, validTrials_IT);

    % Percent of selected channels that showed at least one IED
    if nChans > 0
        thisPctIEDchans_RT = 100 * nIEDchans_RT / nChans;
        thisPctIEDchans_IT = 100 * nIEDchans_IT / nChans;
    else
        thisPctIEDchans_RT = NaN;
        thisPctIEDchans_IT = NaN;
    end

    ptIDs{end+1, 1} = ptID;

    pctIEDtrials_RT(end+1, 1) = thisPctIEDtrials_RT;
    pctIEDtrials_IT(end+1, 1) = thisPctIEDtrials_IT;

    pctIEDchans_RT(end+1, 1) = thisPctIEDchans_RT;
    pctIEDchans_IT(end+1, 1) = thisPctIEDchans_IT;

end

% Save summary table
summaryTable = table( ...
    ptIDs, ...
    pctIEDtrials_RT, ...
    pctIEDtrials_IT, ...
    pctIEDchans_RT, ...
    pctIEDchans_IT, ...
    'VariableNames', { ...
        'PatientID', ...
        'PctIEDTrials_RT', ...
        'PctIEDTrials_IT', ...
        'PctIEDChannels_RT', ...
        'PctIEDChannels_IT'});

outputCSV = fullfile(outputFolderName, 'IED_trial_and_channel_percent_summary.csv');
writetable(summaryTable, outputCSV);

% Plot
fig = figure('Visible', 'off');
set(fig, 'Position', [100 100 1100 500]);

tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plotScatterBoxPair( ...
    pctIEDtrials_RT, ...
    pctIEDtrials_IT, ...
    colorRT, ...
    colorIT, ...
    '% of IED trials', ...
    '% trials with at least one IED');

nexttile;
plotScatterBoxPair( ...
    pctIEDchans_RT, ...
    pctIEDchans_IT, ...
    colorRT, ...
    colorIT, ...
    '% of channels showing IEDs', ...
    '% selected channels with at least one IED');

sgtitle('IED trial percentage and IED channel percentage across participants');

outputPDF = fullfile(outputFolderName, 'IED_trial_and_channel_percent_scatter_boxplot.pdf');
exportgraphics(fig, outputPDF, 'ContentType', 'vector');

close(fig);

disp(' ');
disp('Done.');
disp(['Saved PDF: ' outputPDF]);
disp(['Saved CSV: ' outputCSV]);

function nIED = countIEDsPerTrial(LFPIED, fieldName, nTrials)

    nIED = zeros(nTrials, 1);

    if ~isfield(LFPIED, fieldName)
        return;
    end

    IEDocc = LFPIED.(fieldName);

    if isempty(IEDocc)
        return;
    end

    trialIndices = IEDocc(:, 1);

    validRows = ...
        isfinite(trialIndices) & ...
        trialIndices >= 1 & ...
        trialIndices <= nTrials;

    trialIndices = round(trialIndices(validRows));

    if ~isempty(trialIndices)
        nIED = accumarray(trialIndices, 1, [nTrials 1], @sum, 0);
    end

end

function pctIEDtrials = getPercentIEDTrials(nIED, validTrials)

    nValidTrials = sum(validTrials);

    if nValidTrials == 0
        pctIEDtrials = NaN;
        return;
    end

    nTrialsWithIED = sum(validTrials & nIED > 0);

    pctIEDtrials = 100 * nTrialsWithIED / nValidTrials;

end

function nUniqueChans = countUniqueChannelsWithIED(LFPIED, fieldName, nTrials, validTrials)

    nUniqueChans = 0;

    if ~isfield(LFPIED, fieldName)
        return;
    end

    IEDocc = LFPIED.(fieldName);

    if isempty(IEDocc) || size(IEDocc, 2) < 2
        return;
    end

    trialIndices = IEDocc(:, 1);
    channelIndices = IEDocc(:, 2);

    validRows = ...
        isfinite(trialIndices) & ...
        isfinite(channelIndices) & ...
        trialIndices >= 1 & ...
        trialIndices <= nTrials & ...
        channelIndices >= 1;

    trialIndices = round(trialIndices(validRows));
    channelIndices = round(channelIndices(validRows));

    if isempty(trialIndices)
        return;
    end

    % Keep only IEDs from trials included in the analysis
    keepRows = validTrials(trialIndices);

    channelIndices = channelIndices(keepRows);

    nUniqueChans = length(unique(channelIndices));

end

function plotScatterBoxPair(dataRT, dataIT, colorRT, colorIT, yLabelText, titleText)

    dataRT = dataRT(:);
    dataIT = dataIT(:);

    dataRT = dataRT(isfinite(dataRT));
    dataIT = dataIT(isfinite(dataIT));

    hold on;

    boxData = [dataRT; dataIT];
    groupData = [ones(size(dataRT)); 2 * ones(size(dataIT))];

    boxplot(boxData, groupData, ...
        'Labels', {'RT', 'IT'}, ...
        'Colors', 'k', ...
        'Widths', 0.35, ...
        'Symbol', '');

    set(findobj(gca, 'Tag', 'Box'), 'LineWidth', 0.75);
    set(findobj(gca, 'Tag', 'Median'), 'LineWidth', 1.2);
    set(findobj(gca, 'Tag', 'Whisker'), 'LineWidth', 0.75);
    set(findobj(gca, 'Tag', 'Upper Whisker'), 'LineWidth', 0.75);
    set(findobj(gca, 'Tag', 'Lower Whisker'), 'LineWidth', 0.75);
    set(findobj(gca, 'Tag', 'Upper Adjacent Value'), 'LineWidth', 0.75);
    set(findobj(gca, 'Tag', 'Lower Adjacent Value'), 'LineWidth', 0.75);

    % Jittered scatter points
    rng(1);

    jitterAmount = 0.08;

    xRT = 1 + jitterAmount * randn(size(dataRT));
    xIT = 2 + jitterAmount * randn(size(dataIT));

    scatter(xRT, dataRT, 45, ...
        'MarkerFaceColor', colorRT, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.5, ...
        'LineWidth', 0.4);

    scatter(xIT, dataIT, 45, ...
        'MarkerFaceColor', colorIT, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.5, ...
        'LineWidth', 0.4);

    ylabel(yLabelText);
    title(titleText);

    ylim([0 100]);
    xlim([0.5 2.5]);

    grid off;
    box off;

    pbaspect([1 1 1]);

    hold off;

end