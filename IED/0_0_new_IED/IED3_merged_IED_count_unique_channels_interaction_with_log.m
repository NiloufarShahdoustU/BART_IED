% Test interaction between IED count and unique IED channels
% RT and IT outcomes are z-scored
% Unique channels and IED counts are z-scored
% BR stays as 0 and 1
% No log and no division by duration
% Author: Nill

clear;
clc;
close all;

inputFolder = 'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';
outputFolder = 'D:\Nill\code\BART\IED\0_0_new_IED\IED3_merged_IED_count_unique_channels_interaction_with_log\';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

files = dir(fullfile(inputFolder, '*.LFPIED.mat'));
allTrials = table();

%% 

for pt = 1:length(files)

    fileName = files(pt).name;
    parts = strsplit(fileName, '.');
    ptID = string(parts{1});

    disp("Processing patient: " + ptID)

    data = load(fullfile(inputFolder, fileName));
    LFPIED = data.LFPIED;

    nTrials = LFPIED.nTrials;

    RT = LFPIED.RTs(1:nTrials);
    IT = LFPIED.ITs(1:nTrials);
    BR = LFPIED.BankedTrials(1:nTrials);
    control = LFPIED.isControl(1:nTrials);

    RT = RT(:);
    IT = IT(:);
    BR = BR(:);
    control = control(:);

    [nIED_RT, nChannels_RT] = ...
        countIEDsAndChannels(LFPIED.IED_occurance_RT, nTrials);

    [nIED_IT, nChannels_IT] = ...
        countIEDsAndChannels(LFPIED.IED_occurance_IT, nTrials);

    % RT
    keep = control == 0 & ...
        isfinite(RT) & RT > 0 & RT <= 20;

    allTrials = addTrials( ...
        allTrials, ptID, "RT", RT, ...
        nIED_RT, nChannels_RT, keep);

    % IT
    keep = control == 0 & ...
        isfinite(RT) & RT > 0 & RT <= 20 & ...
        isfinite(IT) & IT > 0;

    allTrials = addTrials( ...
        allTrials, ptID, "IT", IT, ...
        nIED_IT, nChannels_IT, keep);

    % BR
    keep = control == 0 & ...
        isfinite(RT) & RT > 0 & RT <= 20 & ...
        isfinite(IT) & IT > 0 & ...
        (BR == 0 | BR == 1);

    allTrials = addTrials( ...
        allTrials, ptID, "BR", BR, ...
        nIED_IT, nChannels_IT, keep);
end

%% z scoring

allTrials.zY = NaN(height(allTrials), 1);
allTrials.zIEDs = NaN(height(allTrials), 1);
allTrials.zChannels = NaN(height(allTrials), 1);

outcomes = ["RT", "IT", "BR"];

for k = 1:length(outcomes)

    outcome = outcomes(k);
    rows = allTrials.outcome == outcome;

    % Z-score IED counts and unique channels
    allTrials.zIEDs(rows) = zscore(allTrials.nIED(rows));
    allTrials.zChannels(rows) = zscore(allTrials.nChannels(rows));

    % Z-score RT and IT, but not BR
    if outcome == "BR"
        allTrials.zY(rows) = allTrials.y(rows);
    else
        allTrials.zY(rows) = zscore(allTrials.y(rows));
    end
end

writetable(allTrials, ...
    fullfile(outputFolder, 'trial_level_data.csv'));

%% mixed-effects models

results = table();

for k = 1:length(outcomes)

    outcome = outcomes(k);

    T = allTrials(allTrials.outcome == outcome, :);

    T.participant = categorical(T.participant);
    T.channels = T.zChannels;
    T.IEDs = T.zIEDs;
    T.modelY = T.zY;

    if outcome == "BR"

        model = fitglme(T, ...
            'modelY ~ channels*IEDs + (1 | participant)', ...
            'Distribution', 'Binomial', ...
            'Link', 'logit');

        % different model for BR
        modelType = "Logistic mixed-effects";

    else

        model = fitlme(T, ...
            'modelY ~ channels*IEDs + (1 | participant)');

        modelType = "Linear mixed-effects";
    end

    coefficientTable = model.Coefficients;
    CI = coefCI(model);

    wantedTerms = ["channels", "IEDs", "channels:IEDs"];

    for j = 1:length(wantedTerms)

        % intercept is row 1
        % Channels, IEDs, and interaction are rows 2, 3, and 4
        row = j + 1;

        beta = coefficientTable.Estimate(row);
        SE = coefficientTable.SE(row);
        tStat = coefficientTable.tStat(row);
        pValue = coefficientTable.pValue(row);

        CILow = CI(row, 1);
        CIHigh = CI(row, 2);

        oddsRatio = NaN;
        oddsRatioLow = NaN;
        oddsRatioHigh = NaN;

        if outcome == "BR"
            oddsRatio = exp(beta);
            oddsRatioLow = exp(CILow);
            oddsRatioHigh = exp(CIHigh);
        end

        newRow = table( ...
            outcome, ...
            wantedTerms(j), ...
            length(unique(T.participant)), ...
            height(T), ...
            beta, ...
            SE, ...
            tStat, ...
            pValue, ...
            CILow, ...
            CIHigh, ...
            oddsRatio, ...
            oddsRatioLow, ...
            oddsRatioHigh, ...
            modelType, ...
            'VariableNames', ...
            {'outcome', 'term', ...
            'nParticipants', 'nTrials', ...
            'beta', 'SE', 'tStat', 'pValue', ...
            'CILow', 'CIHigh', ...
            'oddsRatio', 'oddsRatioCILow', ...
            'oddsRatioCIHigh', 'modelType'});

        results = [results; newRow];
    end
end

%% FDR correction 

results.interactionP_FDR = NaN(height(results), 1);

interactionRows = results.term == "channels:IEDs";
p = results.pValue(interactionRows);

[sortedP, order] = sort(p);

fdr = sortedP .* length(sortedP) ./ ...
    (1:length(sortedP))';

for k = length(fdr)-1:-1:1
    fdr(k) = min(fdr(k), fdr(k + 1));
end

fdr = min(fdr, 1);

correctedP = zeros(size(p));
correctedP(order) = fdr;

results.interactionP_FDR(interactionRows) = correctedP;

writetable(results, ...
    fullfile(outputFolder, 'mixed_effects_results.csv'));

%% visualization

colors = [
    0.204 0.459 0.702
    0.847 0.333 0.153
    0.250 0.600 0.250
];

termLabels = {
    'Unique channels'
    'Number of IEDs'
    'Interaction'
};

figure( ...
    'Color', 'w', ...
    'Position', [100 100 1200 430]);

tiledlayout(1, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for k = 1:3

    nexttile;

    rows = results(results.outcome == outcomes(k), :);
    x = 1:3;

    lowerError = rows.beta - rows.CILow;
    upperError = rows.CIHigh - rows.beta;

    hold on

    errorbar(x, rows.beta, lowerError, upperError, ...
        'o', ...
        'Color', colors(k, :), ...
        'MarkerFaceColor', colors(k, :), ...
        'MarkerEdgeColor', colors(k, :), ...
        'MarkerSize', 6, ...
        'LineWidth', 1, ...
        'CapSize', 10);

    yline(0, '--', ...
        'Color', [0.45 0.45 0.45], ...
        'LineWidth', 1);

    pFDR = rows.interactionP_FDR(3);

    title(outcomes(k) + ...
        ", interaction p = " + formatP(pFDR), ...
        'FontSize', 12, ...
        'FontWeight', 'bold');

    xticks(1:3)
    xticklabels(termLabels)
    xtickangle(25)
    xlim([0.5 3.5])

    if outcomes(k) == "BR"

        ylabel('Coefficient (log odds)', ...
            'FontWeight', 'bold');

    else

        ylabel('Standardized coefficient', ...
            'FontWeight', 'bold');
    end

    set(gca, ...
        'FontSize', 11, ...
        'LineWidth', 1, ...
        'Box', 'off', ...
        'TickDir', 'out');

    hold off
end

exportgraphics(gcf, ...
    fullfile(outputFolder, ...
    'merged_IED_count_channels_interaction.pdf'), ...
    'ContentType', 'vector');

close(gcf)

disp("Saved results in: " + outputFolder);

%% Functions

function [nIED, nChannels] = ...
    countIEDsAndChannels(IEDocc, nTrials)

    trials = round(IEDocc(:, 1));
    channels = round(IEDocc(:, 2));

    keep = isfinite(trials) & ...
        isfinite(channels) & ...
        trials >= 1 & ...
        trials <= nTrials & ...
        channels >= 1;

    trials = trials(keep);
    channels = channels(keep);

    nIED = accumarray(trials, 1, [nTrials 1]);

    trialChannelPairs = unique([trials channels], 'rows');

    nChannels = accumarray( ...
        trialChannelPairs(:, 1), 1, [nTrials 1]);
end

function allTrials = addTrials( ...
    allTrials, ptID, outcome, y, nIED, nChannels, keep)

    y = y(keep);
    nIED = nIED(keep);
    nChannels = nChannels(keep);

    n = length(y);

    newTrials = table( ...
        repmat(ptID, n, 1), ...
        repmat(outcome, n, 1), ...
        y, ...
        nIED, ...
        nChannels, ...
        'VariableNames', ...
        {'participant', 'outcome', 'y', ...
        'nIED', 'nChannels'});

    allTrials = [allTrials; newTrials];
end

function textP = formatP(p)

    if p < 0.001
        textP = "<0.001";
    else
        textP = string(sprintf('%.3f', p));
    end
end