% Plot number of IEDs and unique IED channels per trial vs RT and IT
% One 2x2 PDF per participant
% Only non-control trials
% Trials with 0 IEDs / 0 IED channels are NOT shown
% Trials with RT > 10 seconds are removed from BOTH RT and IT plots
% No outlier removal
% Only x datapoints are log10-transformed before plotting and fitting
% Linear fit using log10(x) and raw y
% Fitted line with SEM band
% Each subplot is square
% Author: Nill

clear;
clc;
close all;

inputFolderName_LFPIED = 'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';
outputFolderName = 'D:\Nill\code\BART\IED\0_0_new_IED\IED2_IED_chans_vs_IT_RT\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));

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
    validRT_10sec = isfinite(RTs) & RTs <= 20;

    nIED_RT = countIEDsPerTrial(LFPIED, 'IED_occurance_RT', nTrials);
    nIED_IT = countIEDsPerTrial(LFPIED, 'IED_occurance_IT', nTrials);

    nIEDchans_RT = countUniqueIEDchansPerTrial(LFPIED, 'IED_occurance_RT', nTrials);
    nIEDchans_IT = countUniqueIEDchansPerTrial(LFPIED, 'IED_occurance_IT', nTrials);

    keep_IED_RT = nonControlTrials & validRT_10sec & isfinite(RTs) & RTs > 0 & nIED_RT > 0;
    keep_IED_IT = nonControlTrials & validRT_10sec & isfinite(ITs) & ITs > 0 & nIED_IT > 0;

    keep_chans_RT = nonControlTrials & validRT_10sec & isfinite(RTs) & RTs > 0 & nIEDchans_RT > 0;
    keep_chans_IT = nonControlTrials & validRT_10sec & isfinite(ITs) & ITs > 0 & nIEDchans_IT > 0;

    x_IED_RT = nIED_RT(keep_IED_RT);
    y_IED_RT = RTs(keep_IED_RT);

    x_IED_IT = nIED_IT(keep_IED_IT);
    y_IED_IT = ITs(keep_IED_IT);

    x_chans_RT = nIEDchans_RT(keep_chans_RT);
    y_chans_RT = RTs(keep_chans_RT);

    x_chans_IT = nIEDchans_IT(keep_chans_IT);
    y_chans_IT = ITs(keep_chans_IT);

    fig = figure('Visible', 'off');
    set(fig, 'Position', [100 100 1100 1000]);

    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plotOnePanel( ...
        x_IED_RT, ...
        y_IED_RT, ...
        'log10(Number of IEDs in RT window)', ...
        'RT (s)', ...
        ['log10(IED count) vs RT - Patient ' ptID]);

    nexttile;
    plotOnePanel( ...
        x_IED_IT, ...
        y_IED_IT, ...
        'log10(Number of IEDs in IT window)', ...
        'IT (s)', ...
        ['log10(IED count) vs IT - Patient ' ptID]);

    nexttile;
    plotOnePanel( ...
        x_chans_RT, ...
        y_chans_RT, ...
        'log10(Number of unique IED channels in RT window)', ...
        'RT (s)', ...
        ['log10(Unique IED channels) vs RT - Patient ' ptID]);

    nexttile;
    plotOnePanel( ...
        x_chans_IT, ...
        y_chans_IT, ...
        'log10(Number of unique IED channels in IT window)', ...
        'IT (s)', ...
        ['log10(Unique IED channels) vs IT - Patient ' ptID]);

    sgtitle(['Patient ' ptID ...
        ' - nChans = ' num2str(nChans)], ...
        'Interpreter', 'none');

    outputPDF = fullfile(outputFolderName, [ptID '_IED_chans_vs_IT_RT.pdf']);

    exportgraphics(fig, outputPDF, 'ContentType', 'vector');

    close(fig);

end

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

function nIEDchans = countUniqueIEDchansPerTrial(LFPIED, fieldName, nTrials)

    nIEDchans = zeros(nTrials, 1);

    if ~isfield(LFPIED, fieldName)
        return;
    end

    IEDocc = LFPIED.(fieldName);

    if isempty(IEDocc) || size(IEDocc, 2) < 2
        return;
    end

    trialChan = IEDocc(:, 1:2);

    validRows = ...
        isfinite(trialChan(:, 1)) & ...
        isfinite(trialChan(:, 2)) & ...
        trialChan(:, 1) >= 1 & ...
        trialChan(:, 1) <= nTrials & ...
        trialChan(:, 2) >= 1;

    trialIndices = round(trialChan(validRows, 1));
    channelIndices = round(trialChan(validRows, 2));

    uniqueTrialChan = unique([trialIndices channelIndices], 'rows');

    if ~isempty(uniqueTrialChan)
        nIEDchans = accumarray(uniqueTrialChan(:, 1), 1, [nTrials 1], @sum, 0);
    end

end

function plotOnePanel(x, y, xLabelText, yLabelText, titleText)

    x = x(:);
    y = y(:);

    validIdx = isfinite(x) & isfinite(y) & x > 0 & y > 0;
    x = x(validIdx);
    y = y(validIdx);

    logX = log10(x);
    logY = y;

    hold on;

    scatter(logX, logY, 18, ...
        'MarkerFaceColor', [0.00 0.70 0.75], ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.5);

    xlabel(xLabelText);
    ylabel(yLabelText);
    title(titleText, 'Interpreter', 'none');

    grid off;
    box off;

    if length(logX) > 2 && length(unique(logX)) > 1

        [rhoVal, pVal] = corr(logX, logY, ...
            'Type', 'Spearman', ...
            'Rows', 'complete');

        subtitle(sprintf('Spearman rho = %.3f, p = %.4f', rhoVal, pVal));

        plotFitWithSEM(logX, logY);

    else

        subtitle('Not enough data for fit');

    end

    pbaspect([1 1 1]);

    hold off;

end

function plotFitWithSEM(logX, logY)

    logX = logX(:);
    logY = logY(:);

    validIdx = isfinite(logX) & isfinite(logY);
    logX = logX(validIdx);
    logY = logY(validIdx);

    if length(logX) < 3 || length(unique(logX)) < 2
        return;
    end

    mdl = fitlm(logX, logY);

    xFit = linspace(min(logX), max(logX), 100)';
    Xnew = [ones(length(xFit), 1), xFit];

    beta = mdl.Coefficients.Estimate;
    yFit = Xnew * beta;

    covBeta = mdl.CoefficientCovariance;
    semFit = sqrt(diag(Xnew * covBeta * Xnew'));

    upperSEM = yFit + semFit;
    lowerSEM = yFit - semFit;

    fill([xFit; flipud(xFit)], ...
        [upperSEM; flipud(lowerSEM)], ...
        [0.7 0.7 0.7], ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none');

    plot(xFit, yFit, '-', ...
        'Color', [0.35 0.35 0.35], ...
        'LineWidth', 1);

end