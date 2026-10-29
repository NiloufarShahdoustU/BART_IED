% Plot trial-level IED counts and unique IED channels versus RT and IT
% across all participants. Each dot represents one eligible trial.
% Author: Nill
%
% This is a descriptive pooled visualization. Participants with more
% eligible trials contribute more dots. Pooled Spearman correlations and
% fitted lines mix within- and between-participant variation; they are not
% participant-adjusted effects. No pooled p-values or SEM bands are shown.
% Use a repeated-measures analysis separately for population inference.
%
% Requires Statistics and Machine Learning Toolbox for corr (Spearman).

clear;
clc;
close all;

inputFolderName_LFPIED = 'D:\Nill\data\BART\0_0_new_IED_new_area_labels\IED1_find_number_of_IEDs\';
outputFolderName = 'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED2_IED_chans_vs_IT_RT_across_pts\';

maxRT = 20; % Seconds; same RT cutoff as the original script.
markerSize = 10; % Area in points squared.
markerAlpha = 0.25;
pointColor = [0.00 0.70 0.75];
showDescriptiveFit = true;

if ~exist(inputFolderName_LFPIED, 'dir')
    error('Input folder does not exist: %s', inputFolderName_LFPIED);
end
if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_LFPIED, '*.LFPIED.mat'));
if isempty(fileList)
    error('No *.LFPIED.mat files found in: %s', inputFolderName_LFPIED);
end

% Store each participant separately, then concatenate only after loading
% everyone. A participant may contribute to some panels but not others.
xByParticipant = cell(numel(fileList), 4);
yByParticipant = cell(numel(fileList), 4);
idByParticipant = cell(numel(fileList), 4);
participantIDs = strings(numel(fileList), 1);

for pt = 1:numel(fileList)
    fileName = fileList(pt).name;
    fileNameParts = strsplit(fileName, '.');
    ptID = string(fileNameParts{1});
    participantIDs(pt) = ptID;
    fprintf('Processing participant %s (%d/%d)\n', ...
        char(ptID), pt, numel(fileList));

    loadedData = load(fullfile(inputFolderName_LFPIED, fileName), 'LFPIED');
    if ~isfield(loadedData, 'LFPIED')
        error('Missing LFPIED structure in %s.', fileName);
    end
    LFPIED = loadedData.LFPIED;
    requiredFields = {'nTrials', 'RTs', 'ITs', 'isControl', ...
        'IED_occurance_RT', 'IED_occurance_IT'};
    for k = 1:numel(requiredFields)
        if ~isfield(LFPIED, requiredFields{k})
            error('Missing LFPIED.%s in %s.', requiredFields{k}, fileName);
        end
    end

    RTs = LFPIED.RTs(:);
    ITs = LFPIED.ITs(:);
    isControl = LFPIED.isControl(:);
    nTrials = min([LFPIED.nTrials, numel(RTs), numel(ITs), ...
        numel(isControl)]);
    RTs = RTs(1:nTrials);
    ITs = ITs(1:nTrials);
    isControl = isControl(1:nTrials);

    nonControlTrials = isControl == 0;
    validRT = isfinite(RTs) & RTs <= maxRT;

    nIED_RT = countIEDsPerTrial(LFPIED, 'IED_occurance_RT', nTrials);
    nIED_IT = countIEDsPerTrial(LFPIED, 'IED_occurance_IT', nTrials);
    nIEDchans_RT = countUniqueIEDchansPerTrial( ...
        LFPIED, 'IED_occurance_RT', nTrials);
    nIEDchans_IT = countUniqueIEDchansPerTrial( ...
        LFPIED, 'IED_occurance_IT', nTrials);

    % Preserve the original panel-specific filters, including exclusion of
    % zero-IED trials and application of the RT cutoff to IT panels too.
    keep = cell(1, 4);
    keep{1} = nonControlTrials & validRT & RTs > 0 & nIED_RT > 0;
    keep{2} = nonControlTrials & validRT & isfinite(ITs) & ...
        ITs > 0 & nIED_IT > 0;
    keep{3} = nonControlTrials & validRT & RTs > 0 & nIEDchans_RT > 0;
    keep{4} = nonControlTrials & validRT & isfinite(ITs) & ...
        ITs > 0 & nIEDchans_IT > 0;

    trialX = {nIED_RT, nIED_IT, nIEDchans_RT, nIEDchans_IT};
    trialY = {RTs, ITs, RTs, ITs};
    for panel = 1:4
        xByParticipant{pt, panel} = trialX{panel}(keep{panel});
        yByParticipant{pt, panel} = trialY{panel}(keep{panel});
        idByParticipant{pt, panel} = repmat(ptID, sum(keep{panel}), 1);
    end
end

%% One combined figure for all participants
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1300 1050]);
layout = tiledlayout(fig, 2, 2, ...
    'Padding', 'compact', 'TileSpacing', 'compact');

xLabels = { ...
    'log10(Number of IEDs in RT window)', ...
    'log10(Number of IEDs in IT window)', ...
    'log10(Number of unique IED channels in RT window)', ...
    'log10(Number of unique IED channels in IT window)'};
yLabels = {'RT (s)', 'IT (s)', 'RT (s)', 'IT (s)'};
panelTitles = { ...
    'a  IED count versus RT', ...
    'b  IED count versus IT', ...
    'c  Unique IED channels versus RT', ...
    'd  Unique IED channels versus IT'};

for panel = 1:4
    ax = nexttile(layout, panel);
    allX = vertcat(xByParticipant{:, panel});
    allY = vertcat(yByParticipant{:, panel});
    allIDs = vertcat(idByParticipant{:, panel});
    plotOnePanel(ax, allX, allY, allIDs, ...
        xLabels{panel}, yLabels{panel}, panelTitles{panel}, ...
        markerSize, markerAlpha, pointColor, showDescriptiveFit);
end

title(layout, sprintf( ...
    'IED counts, channel counts, and behavior across participants (%d loaded)', ...
    numel(unique(participantIDs))), ...
    'FontWeight', 'bold', 'Interpreter', 'none');

outputPDF = fullfile(outputFolderName, ...
    'IED2_IED_chans_vs_IT_RT_across_pts.pdf');
exportgraphics(fig, outputPDF, 'ContentType', 'vector');
close(fig);
fprintf('Saved visualization: %s\n', outputPDF);

%% Local functions
function nIED = countIEDsPerTrial(LFPIED, fieldName, nTrials)
    nIED = zeros(nTrials, 1);
    IEDocc = LFPIED.(fieldName);
    if isempty(IEDocc)
        return;
    end
    trialIndices = IEDocc(:, 1);
    validRows = isfinite(trialIndices) & trialIndices >= 1 & ...
        trialIndices <= nTrials;
    trialIndices = round(trialIndices(validRows));
    if ~isempty(trialIndices)
        nIED = accumarray(trialIndices, 1, [nTrials 1], @sum, 0);
    end
end

function nIEDchans = countUniqueIEDchansPerTrial(LFPIED, fieldName, nTrials)
    nIEDchans = zeros(nTrials, 1);
    IEDocc = LFPIED.(fieldName);
    if isempty(IEDocc)
        return;
    end
    if size(IEDocc, 2) < 2
        error('LFPIED.%s must contain trial and channel columns.', fieldName);
    end
    trialChan = IEDocc(:, 1:2);
    validRows = isfinite(trialChan(:, 1)) & isfinite(trialChan(:, 2)) & ...
        trialChan(:, 1) >= 1 & trialChan(:, 1) <= nTrials & ...
        trialChan(:, 2) >= 1;
    trialIndices = round(trialChan(validRows, 1));
    channelIndices = round(trialChan(validRows, 2));
    uniqueTrialChan = unique([trialIndices channelIndices], 'rows');
    if ~isempty(uniqueTrialChan)
        nIEDchans = accumarray(uniqueTrialChan(:, 1), 1, ...
            [nTrials 1], @sum, 0);
    end
end

function plotOnePanel(ax, x, y, participantID, ...
        xLabelText, yLabelText, titleText, ...
        markerSize, markerAlpha, pointColor, showDescriptiveFit)
    x = x(:);
    y = y(:);
    participantID = participantID(:);
    valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;
    logX = log10(x(valid));
    y = y(valid);
    participantID = participantID(valid);
    numberOfParticipants = numel(unique(participantID));

    hold(ax, 'on');
    scatter(ax, logX, y, markerSize, pointColor, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', markerAlpha);
    xlabel(ax, xLabelText, 'Interpreter', 'none');
    ylabel(ax, yLabelText);
    title(ax, titleText, 'Interpreter', 'none');
    grid(ax, 'off');
    box(ax, 'off');

    subtitleText = sprintf('%d participants; %d trials', ...
        numberOfParticipants, numel(y));
    enoughData = numel(y) > 2 && numel(unique(logX)) > 1;
    if enoughData
        if numel(unique(y)) > 1
            rho = corr(logX, y, 'Type', 'Spearman', 'Rows', 'complete');
            subtitleText = sprintf( ...
                '%s | pooled rho = %.3f (descriptive)', subtitleText, rho);
        end
        if showDescriptiveFit
            % Descriptive least-squares line, without an independence-based
            % confidence/SEM band. No pooled significance test is performed.
            coefficients = [ones(numel(logX), 1), logX] \ y;
            xFit = linspace(min(logX), max(logX), 100)';
            yFit = [ones(numel(xFit), 1), xFit] * coefficients;
            plot(ax, xFit, yFit, '-', 'Color', [0.30 0.30 0.30], ...
                'LineWidth', 1.5);
        end
    elseif isempty(y)
        text(ax, 0.5, 0.5, 'No eligible trials', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
    else
        subtitleText = [subtitleText ' | insufficient variation for fit'];
    end
    subtitle(ax, subtitleText, 'Interpreter', 'none');
    pbaspect(ax, [1 1 1]);
    hold(ax, 'off');
    fprintf('%s: %d participants, %d trials\n', ...
        titleText, numberOfParticipants, numel(y));
end
