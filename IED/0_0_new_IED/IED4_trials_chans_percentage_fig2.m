% Scatter-boxplot of:
% 1) Mean number of IEDs per valid trial
% 2) Mean number of unique IED channels per valid trial


clear;
clc;
close all;

inputFolderName_LFPIED = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

outputFolderName = ...
    'D:\Nill\code\BART\IED\0_0_new_IED\IED4_trials_chans_percentage\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

% Colors
colorRT = [0.204 0.459 0.702];   % blue
colorIT = [0.847 0.333 0.153];   % orange

ptIDs = {};

nChannelsPerParticipant = [];
nValidTrials_RT = [];
nValidTrials_IT = [];

meanIEDs_RT = [];
meanIEDs_IT = [];

meanUniqueChans_RT = [];
meanUniqueChans_IT = [];

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    ptID = fileNameParts{1};

    disp(' ');
    disp(['Processing patient ID: ' ptID]);

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName));
    LFPIED = loadedData.LFPIED;

    % -------------------------------------------------------------
    % Number of channels for this participant
    % -------------------------------------------------------------
    nChans = length(LFPIED.selectedChans);

    nTrials = LFPIED.nTrials;

    RTs = LFPIED.RTs(:);
    ITs = LFPIED.ITs(:);
    isControl = LFPIED.isControl(:);

    % Make all trial-level variables the same length
    minLen = min([ ...
        nTrials, ...
        length(RTs), ...
        length(ITs), ...
        length(isControl)]);

    nTrials = minLen;

    RTs = RTs(1:nTrials);
    ITs = ITs(1:nTrials);
    isControl = isControl(1:nTrials);

    nonControlTrials = isControl == 0;

    % RT > 10 seconds is excluded from both RT and IT analyses
    validRT_10sec = isfinite(RTs) & RTs <= 20;

    % Valid trials for RT analysis
    validTrials_RT = ...
        nonControlTrials & ...
        validRT_10sec & ...
        isfinite(RTs) & ...
        RTs > 0;

    % Valid trials for IT analysis
    validTrials_IT = ...
        nonControlTrials & ...
        validRT_10sec & ...
        isfinite(ITs) & ...
        ITs > 0;

    thisNValidTrials_RT = sum(validTrials_RT);
    thisNValidTrials_IT = sum(validTrials_IT);

    % -------------------------------------------------------------
    % Number of IEDs in each trial
    % -------------------------------------------------------------
    nIEDsPerTrial_RT = countIEDsPerTrial( ...
        LFPIED, ...
        'IED_occurance_RT', ...
        nTrials);

    nIEDsPerTrial_IT = countIEDsPerTrial( ...
        LFPIED, ...
        'IED_occurance_IT', ...
        nTrials);

    % -------------------------------------------------------------
    % Number of unique IED channels in each trial
    % -------------------------------------------------------------
    nUniqueChansPerTrial_RT = countUniqueChannelsPerTrial( ...
        LFPIED, ...
        'IED_occurance_RT', ...
        nTrials);

    nUniqueChansPerTrial_IT = countUniqueChannelsPerTrial( ...
        LFPIED, ...
        'IED_occurance_IT', ...
        nTrials);

    % -------------------------------------------------------------
    % Calculate means across valid trials
    %
    % Trials with zero IEDs or zero unique channels are included.
    % -------------------------------------------------------------
    thisMeanIEDs_RT = getMeanAcrossValidTrials( ...
        nIEDsPerTrial_RT, ...
        validTrials_RT);

    thisMeanIEDs_IT = getMeanAcrossValidTrials( ...
        nIEDsPerTrial_IT, ...
        validTrials_IT);

    thisMeanUniqueChans_RT = getMeanAcrossValidTrials( ...
        nUniqueChansPerTrial_RT, ...
        validTrials_RT);

    thisMeanUniqueChans_IT = getMeanAcrossValidTrials( ...
        nUniqueChansPerTrial_IT, ...
        validTrials_IT);

    % -------------------------------------------------------------
    % Store participant results
    % -------------------------------------------------------------
    ptIDs{end+1, 1} = ptID;

    nChannelsPerParticipant(end+1, 1) = nChans;

    nValidTrials_RT(end+1, 1) = thisNValidTrials_RT;
    nValidTrials_IT(end+1, 1) = thisNValidTrials_IT;

    meanIEDs_RT(end+1, 1) = thisMeanIEDs_RT;
    meanIEDs_IT(end+1, 1) = thisMeanIEDs_IT;

    meanUniqueChans_RT(end+1, 1) = thisMeanUniqueChans_RT;
    meanUniqueChans_IT(end+1, 1) = thisMeanUniqueChans_IT;

    % Display participant results
    fprintf('Number of channels: %d\n', nChans);

    fprintf('Number of valid RT trials: %d\n', ...
        thisNValidTrials_RT);

    fprintf('Number of valid IT trials: %d\n', ...
        thisNValidTrials_IT);

    fprintf('Mean IEDs per RT trial: %.4f\n', ...
        thisMeanIEDs_RT);

    fprintf('Mean IEDs per IT trial: %.4f\n', ...
        thisMeanIEDs_IT);


    fprintf('Mean unique IED channels per RT trial: %.4f\n', ...
        thisMeanUniqueChans_RT);

    fprintf('Mean unique IED channels per IT trial: %.4f\n', ...
        thisMeanUniqueChans_IT);


end

%% Save summary table

summaryTable = table( ...
    ptIDs, ...
    nChannelsPerParticipant, ...
    nValidTrials_RT, ...
    nValidTrials_IT, ...
    meanIEDs_RT, ...
    meanIEDs_IT, ...
    meanUniqueChans_RT, ...
    meanUniqueChans_IT, ...
    'VariableNames', { ...
        'PatientID', ...
        'NumberOfChannels', ...
        'NumberOfValidTrials_RT', ...
        'NumberOfValidTrials_IT', ...
        'MeanIEDsPerTrial_RT', ...
        'MeanIEDsPerTrial_IT', ...
        'MeanUniqueIEDChannelsPerTrial_RT', ...
        'MeanUniqueIEDChannelsPerTrial_IT'});

outputCSV = fullfile( ...
    outputFolderName, ...
    'mean_IED_and_unique_channel_counts_per_trial.csv');

writetable(summaryTable, outputCSV);

%% Plot

fig = figure('Visible', 'off');

set(fig, ...
    'Position', [100 100 650 400], ...
    'Color', 'w');

tiledlayout(1, 2, ...
    'Padding', 'compact', ...
    'TileSpacing', 'compact');

% -------------------------------------------------------------
% Mean number of IEDs per trial
% -------------------------------------------------------------
nexttile;

plotScatterBoxPair( ...
    meanIEDs_RT, ...
    meanIEDs_IT, ...
    colorRT, ...
    colorIT, ...
    'Mean number of IEDs per trial', ...
    'Mean IEDs per trial');

% -------------------------------------------------------------
% Mean number of unique IED channels per trial
% -------------------------------------------------------------
nexttile;

plotScatterBoxPair( ...
    meanUniqueChans_RT, ...
    meanUniqueChans_IT, ...
    colorRT, ...
    colorIT, ...
    'Mean number of unique IED channels per trial', ...
    'Mean unique IED channels per trial');

sgtitle( ...
    'Mean IED and unique-channel counts across trials', ...
    'Color', 'k');

% Make all figure fonts and axis text black
set(findall(fig, 'Type', 'text'), 'Color', 'k');
set(findall(fig, 'Type', 'axes'), ...
    'XColor', 'k', ...
    'YColor', 'k');

outputPDF = fullfile( ...
    outputFolderName, ...
    'mean_IED_and_unique_channel_counts_per_trial_scatter_boxplot.pdf');

exportgraphics(fig, outputPDF, 'ContentType', 'vector');

close(fig);

disp(' ');
disp('Done.');
disp(['Saved PDF: ' outputPDF]);
disp(['Saved CSV: ' outputCSV]);


%% Local functions


function nIEDsPerTrial = countIEDsPerTrial( ...
    LFPIED, fieldName, nTrials)

    % Initialize every trial with zero IEDs
    nIEDsPerTrial = zeros(nTrials, 1);

    if ~isfield(LFPIED, fieldName)
        return;
    end

    IEDocc = LFPIED.(fieldName);

    if isempty(IEDocc) || size(IEDocc, 2) < 1
        return;
    end

    trialIndices = IEDocc(:, 1);

    validRows = ...
        isfinite(trialIndices) & ...
        trialIndices >= 1 & ...
        trialIndices <= nTrials;

    trialIndices = round(trialIndices(validRows));

    if ~isempty(trialIndices)

        nIEDsPerTrial = accumarray( ...
            trialIndices, ...
            1, ...
            [nTrials 1], ...
            @sum, ...
            0);

    end

end


function nUniqueChansPerTrial = countUniqueChannelsPerTrial( ...
    LFPIED, fieldName, nTrials)

    % Initialize every trial with zero unique IED channels
    nUniqueChansPerTrial = zeros(nTrials, 1);

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

    % Count different channels containing at least one IED
    % separately for each trial
    for trialNumber = 1:nTrials

        channelsThisTrial = ...
            channelIndices(trialIndices == trialNumber);

        if ~isempty(channelsThisTrial)

            nUniqueChansPerTrial(trialNumber) = ...
                length(unique(channelsThisTrial));

        end

    end

end


function meanValue = getMeanAcrossValidTrials( ...
    trialValues, validTrials)

    validTrials = logical(validTrials(:));
    trialValues = trialValues(:);

    if ~any(validTrials)

        meanValue = NaN;
        return;

    end

    % Trials containing zero IEDs or zero channels are included
    meanValue = mean(trialValues(validTrials));

end


function plotScatterBoxPair( ...
    dataRT, dataIT, colorRT, colorIT, yLabelText, titleText)

    dataRT = dataRT(:);
    dataIT = dataIT(:);

    dataRT = dataRT(isfinite(dataRT));
    dataIT = dataIT(isfinite(dataIT));

    hold on;

    % Bring the two groups closer together
    positionRT = 1.00;
    positionIT = 1.4;

    boxData = [dataRT; dataIT];

    groupData = [ ...
        ones(size(dataRT)); ...
        2 * ones(size(dataIT))];

    boxplot( ...
        boxData, ...
        groupData, ...
        'Positions', [positionRT positionIT], ...
        'Labels', {'RT', 'IT'}, ...
        'Colors', 'k', ...
        'Widths', 0.2, ...
        'Symbol', '');

    set(findobj(gca, 'Tag', 'Box'), ...
        'LineWidth', 0.5);

    set(findobj(gca, 'Tag', 'Median'), ...
        'LineWidth', 0.5);

    set(findobj(gca, 'Tag', 'Whisker'), ...
        'LineWidth', 0.5);

    set(findobj(gca, 'Tag', 'Upper Whisker'), ...
        'LineWidth', 0.5);

    set(findobj(gca, 'Tag', 'Lower Whisker'), ...
        'LineWidth', 0.5);

    set(findobj(gca, 'Tag', 'Upper Adjacent Value'), ...
        'LineWidth', 0.5);

    set(findobj(gca, 'Tag', 'Lower Adjacent Value'), ...
        'LineWidth', 0.5);

    % Reproducible jitter
    rng(1);

    jitterAmount = 0.03;

    xRT = positionRT + jitterAmount * randn(size(dataRT));
    xIT = positionIT + jitterAmount * randn(size(dataIT));

    % Smaller scatter points
    scatter( ...
        xRT, ...
        dataRT, ...
        18, ...
        'MarkerFaceColor', colorRT, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.5);

    scatter( ...
        xIT, ...
        dataIT, ...
        18, ...
        'MarkerFaceColor', colorIT, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.5);

    ylabel(yLabelText, ...
        'FontSize', 12, ...
        'Color', 'k');

    title(titleText, ...
        'FontSize', 12, ...
        'Color', 'k');

    % Narrower x-axis because groups are closer
    xlim([0.65 1.90]);

    % Automatically set an appropriate y-axis range
    allData = [dataRT; dataIT];

    if ~isempty(allData)

        dataMin = min(allData);
        dataMax = max(allData);
        dataRange = dataMax - dataMin;

        if dataRange == 0
            padding = max(0.01, 0.15 * abs(dataMax));
        else
            padding = 0.10 * dataRange;
        end

        lowerLimit = max(0, dataMin - padding);
        upperLimit = dataMax + padding;

        if upperLimit <= lowerLimit
            upperLimit = lowerLimit + 1;
        end

        ylim([lowerLimit upperLimit]);

    end

    % Set all axes text to 12 and point tick marks outward
    ax = gca;

    set(ax, ...
        'FontSize', 12, ...
        'TickDir', 'out', ...
        'TickLength', [0.02 0.02], ...
        'LineWidth', 0.5, ...
        'XColor', 'k', ...
        'YColor', 'k', ...
        'XTick', [positionRT positionIT], ...
        'XTickLabel', {'RT', 'IT'});

    grid off;
    box off;

    pbaspect([1 1 1]);

    hold off;

end