% ied vs non-ied spectrograms for every significant ied7 brain area

clear;
clc;
close all;

% folders

lfpiedFolder = ...
    'D:\Nill\data\BART\0_0_new_IED_new_area_labels\IED1_find_number_of_IEDs\';

rawDataFolder = ...
    'D:\Nill\data\BART_preprocessed\';

ied7ResultsFile ='D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED7_Cox_IT_RT_BR_postIED_by_brain_area\IT_RT_BR_brain_area_cox_all_results.mat';

outputFolder = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\' ...
     'IED7_significant_areas_IED_vs_nonIED_spectrogram\'];

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

outputPDF = fullfile(outputFolder, ...
    'IED7_significant_areas_IED_vs_nonIED_spectrogram.pdf');

if exist(outputPDF, 'file')
    delete(outputPDF);
end

% spectrogram settings

preIEDSeconds = 1;
postIEDSeconds = 1;
extraSeconds = 0.5;

baselineStartSeconds = -0.8;
baselineEndSeconds = -0.3;

spectrogramFrequencies = logspace(log10(5), log10(200), 60);
spectrogramTime = -0.9:0.02:0.9;

amplitudeThreshold = 5000;
minimumDistanceBetweenIEDsSeconds = 0.5;
maximumRTSeconds = 20;

% cluster permutation settings

numberOfPermutations = 1000;
alphaBin = 0.05;
alphaCluster = 0.05;
minimumClusterSize = 5;
minimumClusterDurationMs = 50;
nonSignificantAlpha = 0.20;

rng(1);

if exist('basewaveERP', 'file') ~= 2
    error('basewaveERP.m was not found');
end

if exist('bwconncomp', 'file') ~= 2
    error('bwconncomp was not found');
end

% get every significant area from rt, it, and br

ied7 = load(ied7ResultsFile);
outcomeNames = ["RT", "IT", "BR"];

allResults = struct( ...
    'outcome', {}, ...
    'area', {}, ...
    'iedTF', {}, ...
    'nonIEDTF', {}, ...
    'participantIDs', {}, ...
    'numberOfIEDEvents', {}, ...
    'numberOfNonIEDEvents', {});

for oo = 1:length(outcomeNames)

    outcome = outcomeNames(oo);
    resultTable = ied7.(char(outcome + "Results"));

    useRows = ...
        string(resultTable.status) == "fitted" & ...
        logical(resultTable.significantUncorrected) & ...
        isfinite(double(resultTable.beta_logHazard));

    significantAreas = string(resultTable.anatomicalArea(useRows));
    significantAreas = clean_area_names(significantAreas);
    significantAreas = unique(significantAreas, 'stable');

    significantAreas( ...
        ismissing(significantAreas) | ...
        strlength(significantAreas) == 0 | ...
        lower(significantAreas) == "unknown") = [];

    fprintf('\n%s significant areas:\n', outcome);
    disp(significantAreas);

    for aa = 1:length(significantAreas)

        rr = length(allResults) + 1;

        allResults(rr).outcome = outcome;
        allResults(rr).area = significantAreas(aa);

        allResults(rr).iedTF = zeros( ...
            length(spectrogramFrequencies), ...
            length(spectrogramTime), 0);

        allResults(rr).nonIEDTF = zeros( ...
            length(spectrogramFrequencies), ...
            length(spectrogramTime), 0);

        allResults(rr).participantIDs = strings(0, 1);
        allResults(rr).numberOfIEDEvents = [];
        allResults(rr).numberOfNonIEDEvents = [];

    end

end

if isempty(allResults)
    error('no significant ied7 brain areas were found');
end

% load each participant only once

fileList = dir(fullfile(lfpiedFolder, '*.LFPIED.mat'));

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileParts = strsplit(fileName, '.');
    participantID = string(fileParts{1});

    fprintf('\nparticipant %d/%d: %s\n', ...
        pt, length(fileList), participantID);

    loadedData = load(fullfile(lfpiedFolder, fileName));

    if ~isfield(loadedData, 'LFPIED')
        fprintf('skipped: LFPIED was missing\n');
        continue;
    end

    LFPIED = loadedData.LFPIED;

    patientRawFolder = fullfile( ...
        rawDataFolder, char(participantID), 'Data');

    ns2List = dir(fullfile(patientRawFolder, '*.ns2'));
    nevList = dir(fullfile(patientRawFolder, '*.nev'));

    if length(ns2List) ~= 1 || length(nevList) ~= 1
        fprintf('skipped: need one ns2 and one nev file\n');
        continue;
    end

    % get the response times from the triggers

    NEV = openNEV( ...
        fullfile(nevList.folder, nevList.name), 'overwrite');

    trigs = NEV.Data.SerialDigitalIO.UnparsedData;
    trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
    responseTimes = trigTimes(trigs == 23);

    clear NEV;

    % load the lfp

    NSX = openNSx(fullfile(ns2List.folder, ns2List.name));
    Fs = double(NSX.MetaTags.SamplingFreq);
    rawLFP = NSX.Data;

    if iscell(rawLFP)
        rawLFP = [rawLFP{:}];
    end

    if Fs <= 400
        fprintf('skipped: sampling rate is too low\n');
        clear NSX rawLFP;
        continue;
    end

    selectedChans = round(double(LFPIED.selectedChans(:)));
    allAreaLabels = make_string_labels(LFPIED.anatomicalLocs);

    if isempty(selectedChans) || isempty(allAreaLabels)
        fprintf('skipped: channel labels were missing\n');
        clear NSX rawLFP;
        continue;
    end

    if length(allAreaLabels) >= max(selectedChans)
        selectedAreaLabels = allAreaLabels(selectedChans);
    elseif length(allAreaLabels) == length(selectedChans)
        selectedAreaLabels = allAreaLabels;
    else
        fprintf('skipped: channel labels did not match selectedChans\n');
        clear NSX rawLFP;
        continue;
    end

    selectedAreaLabels = clean_area_names(selectedAreaLabels);

    % do every significant outcome-area pair

    for rr = 1:length(allResults)

        outcome = allResults(rr).outcome;
        areaName = allResults(rr).area;

        areaChannels = find(strcmpi(selectedAreaLabels, areaName));

        if isempty(areaChannels)
            continue;
        end

        RTs = double(LFPIED.RTs(:));
        ITs = double(LFPIED.ITs(:));
        isControl = double(LFPIED.isControl(:));
        balloonType = double(LFPIED.balloonType(:));
        bankedTrials = double(LFPIED.BankedTrials(:));
        balloonTimes = double(LFPIED.balloonTimes(:));

        if outcome == "RT"
            IEDrows = double(LFPIED.IED_occurance_RT);
            startTimes = balloonTimes;
            duration = RTs;
        else
            IEDrows = double(LFPIED.IED_occurance_IT);
            startTimes = double(responseTimes(:));
            duration = ITs;
        end

        nTrials = min([ ...
            length(RTs), length(ITs), length(isControl), ...
            length(balloonType), length(bankedTrials), ...
            length(startTimes), length(duration)]);

        if nTrials < 1 || isempty(IEDrows)
            continue;
        end

        RTs = RTs(1:nTrials);
        duration = duration(1:nTrials);
        isControl = isControl(1:nTrials);
        balloonType = balloonType(1:nTrials);
        bankedTrials = bankedTrials(1:nTrials);
        startTimes = startTimes(1:nTrials);

        colorCode = NaN(size(balloonType));
        goodColors = ismember(round(balloonType), [1 2 3 11 12 13]);

        colorCode(goodColors) = ...
            mod(round(balloonType(goodColors)) - 1, 10) + 1;

        validTrials = ...
            isControl == 0 & ...
            isfinite(RTs) & RTs > 0 & RTs <= maximumRTSeconds & ...
            isfinite(duration) & duration > 0 & ...
            isfinite(startTimes) & ...
            ismember(colorCode, [1 2 3]);

        % it only uses banked trials

        if outcome == "IT"
            validTrials = validTrials & bankedTrials == 1;
        end

        % br uses banked and lost trials

        if outcome == "BR"
            validTrials = ...
                validTrials & ismember(bankedTrials, [0 1]);
        end

        % keep ieds from this area and valid trials

        goodIEDrows = ...
            isfinite(IEDrows(:, 1)) & ...
            isfinite(IEDrows(:, 2)) & ...
            isfinite(IEDrows(:, 3)) & ...
            ismember(round(IEDrows(:, 2)), areaChannels);

        areaIEDrows = IEDrows(goodIEDrows, :);

        if isempty(areaIEDrows)
            continue;
        end

        trialNumbers = round(areaIEDrows(:, 1));

        keepIEDrows = ...
            trialNumbers >= 1 & trialNumbers <= nTrials;

        goodRowNumbers = find(keepIEDrows);

        if ~isempty(goodRowNumbers)
            keepIEDrows(goodRowNumbers) = ...
                validTrials(trialNumbers(goodRowNumbers));
        end

        areaIEDrows = areaIEDrows(keepIEDrows, :);

        if isempty(areaIEDrows)
            continue;
        end

        % make one list of ied times per trial
        % close detections across channels count as the same event

        eventTimesByTrial = cell(nTrials, 1);
        IEDtrials = unique(round(areaIEDrows(:, 1)), 'stable');

        for ii = 1:length(IEDtrials)

            thisTrial = IEDtrials(ii);

            thisTrialTimes = sort(round(areaIEDrows( ...
                round(areaIEDrows(:, 1)) == thisTrial, 3)));

            if isempty(thisTrialTimes)
                continue;
            end

            keepTimes = [true; ...
                diff(thisTrialTimes) >= ...
                round(minimumDistanceBetweenIEDsSeconds * Fs)];

            eventTimesByTrial{thisTrial} = ...
                thisTrialTimes(keepTimes);

        end

        IEDtrials = IEDtrials( ...
            ~cellfun(@isempty, eventTimesByTrial(IEDtrials)));

        % non-ied trials have no ied in this brain area

        nonIEDtrials = find( ...
            validTrials & ...
            ~ismember((1:nTrials)', IEDtrials));

        if isempty(IEDtrials) || isempty(nonIEDtrials)
            continue;
        end

        % match without replacement inside this participant
        % first match color and bank/pop, then use closest duration

        [matchedIEDtrials, matchedNonIEDtrials] = ...
            match_trials( ...
                IEDtrials, nonIEDtrials, ...
                colorCode, bankedTrials, duration);

        if isempty(matchedIEDtrials)
            continue;
        end

        participantIEDevents = zeros( ...
            length(spectrogramFrequencies), ...
            length(spectrogramTime), 0);

        participantNonIEDevents = zeros( ...
            length(spectrogramFrequencies), ...
            length(spectrogramTime), 0);

        for mm = 1:length(matchedIEDtrials)

            IEDtrial = matchedIEDtrials(mm);
            nonIEDtrial = matchedNonIEDtrials(mm);

            realTimes = eventTimesByTrial{IEDtrial};

            for ee = 1:length(realTimes)

                realFraction = ...
                    ((realTimes(ee) - 1) / Fs) / duration(IEDtrial);

                realFraction = min(max(realFraction, 0), 1);

                pseudoTime = ...
                    round(realFraction * duration(nonIEDtrial) * Fs) + 1;

                realAbsoluteSample = ...
                    floor(Fs * startTimes(IEDtrial)) + ...
                    realTimes(ee) - 1;

                pseudoAbsoluteSample = ...
                    floor(Fs * startTimes(nonIEDtrial)) + ...
                    pseudoTime - 1;

                realChannelTF = zeros( ...
                    length(spectrogramFrequencies), ...
                    length(spectrogramTime), 0);

                pseudoChannelTF = zeros( ...
                    length(spectrogramFrequencies), ...
                    length(spectrogramTime), 0);

                % use the same area channels for each real-pseudo pair

                for cc = 1:length(areaChannels)

                    localChannel = areaChannels(cc);
                    rawChannelNumber = selectedChans(localChannel);

                    if rawChannelNumber < 1 || ...
                            rawChannelNumber > size(rawLFP, 1)
                        continue;
                    end

                    signal = double(rawLFP(rawChannelNumber, :));

                    realTF = make_one_tf( ...
                        signal, realAbsoluteSample, Fs, ...
                        preIEDSeconds, postIEDSeconds, extraSeconds, ...
                        baselineStartSeconds, baselineEndSeconds, ...
                        spectrogramFrequencies, spectrogramTime, ...
                        amplitudeThreshold);

                    pseudoTF = make_one_tf( ...
                        signal, pseudoAbsoluteSample, Fs, ...
                        preIEDSeconds, postIEDSeconds, extraSeconds, ...
                        baselineStartSeconds, baselineEndSeconds, ...
                        spectrogramFrequencies, spectrogramTime, ...
                        amplitudeThreshold);

                    if isempty(realTF) || isempty(pseudoTF)
                        continue;
                    end

                    realChannelTF(:, :, end + 1) = realTF;
                    pseudoChannelTF(:, :, end + 1) = pseudoTF;

                end

                if isempty(realChannelTF)
                    continue;
                end

                participantIEDevents(:, :, end + 1) = ...
                    mean(realChannelTF, 3, 'omitnan');

                participantNonIEDevents(:, :, end + 1) = ...
                    mean(pseudoChannelTF, 3, 'omitnan');

            end

        end

        if isempty(participantIEDevents)
            continue;
        end

        % first average channels/events inside the participant

        participantIEDTF = ...
            mean(participantIEDevents, 3, 'omitnan');

        participantNonIEDTF = ...
            mean(participantNonIEDevents, 3, 'omitnan');

        allResults(rr).iedTF(:, :, end + 1) = participantIEDTF;
        allResults(rr).nonIEDTF(:, :, end + 1) = participantNonIEDTF;
        allResults(rr).participantIDs(end + 1, 1) = participantID;

        allResults(rr).numberOfIEDEvents(end + 1, 1) = ...
            size(participantIEDevents, 3);

        allResults(rr).numberOfNonIEDEvents(end + 1, 1) = ...
            size(participantNonIEDevents, 3);

        fprintf('%s %s: %d matched events\n', ...
            outcome, areaName, size(participantIEDevents, 3));

    end

    clear NSX rawLFP LFPIED loadedData;

end

% make one pdf page for every significant outcome-area pair

pageNumber = 0;

for rr = 1:length(allResults)

    outcome = allResults(rr).outcome;
    areaName = allResults(rr).area;
    IEDall = allResults(rr).iedTF;
    nonIEDall = allResults(rr).nonIEDTF;

    numberOfParticipants = size(IEDall, 3);

    if numberOfParticipants < 1
        fprintf('\n%s %s skipped: no matched participants\n', ...
            outcome, areaName);
        continue;
    end

    meanIED = mean(IEDall, 3, 'omitnan');
    meanNonIED = mean(nonIEDall, 3, 'omitnan');
    meanDifference = mean(IEDall - nonIEDall, 3, 'omitnan');

    significantMask = false(size(meanDifference));
    clusterPvalues = [];

    if numberOfParticipants >= 2

        fprintf('\ncbpt: %s %s, n = %d\n', ...
            outcome, areaName, numberOfParticipants);

        [significantMask, clusterPvalues] = ...
            run_paired_cbpt( ...
                IEDall, nonIEDall, ...
                spectrogramTime, ...
                numberOfPermutations, ...
                alphaBin, alphaCluster, ...
                minimumClusterSize, ...
                minimumClusterDurationMs);

    end

    groupValues = [meanIED(:); meanNonIED(:)];
    groupValues = abs(groupValues(isfinite(groupValues)));

    if isempty(groupValues)
        groupColorLimit = 1;
    else
        groupColorLimit = prctile(groupValues, 98);
    end

    if ~isfinite(groupColorLimit) || groupColorLimit <= 0
        groupColorLimit = 1;
    end

    differenceValues = abs(meanDifference(isfinite(meanDifference)));

    if isempty(differenceValues)
        differenceColorLimit = 1;
    else
        differenceColorLimit = prctile(differenceValues, 98);
    end

    if ~isfinite(differenceColorLimit) || differenceColorLimit <= 0
        differenceColorLimit = 1;
    end

    alphaMask = nonSignificantAlpha .* ones(size(significantMask));
    alphaMask(significantMask) = 1;

    fig = figure('Visible', 'off', 'Color', 'w', 'Renderer', 'painters');
    set(fig, 'Position', [40 40 1900 720]);

    layout = tiledlayout(fig, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(layout);
    plot_vector_heatmap(ax1, ...
        spectrogramTime, spectrogramFrequencies, meanIED, []);
    set(ax1, ...
        'YScale', 'log', ...
        'YDir', 'normal', ...
        'YLim', [5 200], ...
        'YTick', [10 100], ...
        'FontSize', 12, ...
        'TickDir', 'out');
    xline(ax1, 0, 'r', 'LineWidth', 1);
    xlabel(ax1, 'time from ied peak (s)');
    ylabel(ax1, 'frequency (hz)');
    title(ax1, sprintf('ied, n = %d participants', ...
        numberOfParticipants));
    colormap(ax1, parula);
    axes(ax1);
    caxis([-groupColorLimit groupColorLimit]);
    colorbar(ax1);
    box(ax1, 'on');

    ax2 = nexttile(layout);
    plot_vector_heatmap(ax2, ...
        spectrogramTime, spectrogramFrequencies, meanNonIED, []);
    set(ax2, ...
        'YScale', 'log', ...
        'YDir', 'normal', ...
        'YLim', [5 200], ...
        'YTick', [10 100], ...
        'FontSize', 12, ...
        'TickDir', 'out');
    xline(ax2, 0, 'r', 'LineWidth', 1);
    xlabel(ax2, 'time from pseudo-ied (s)');
    ylabel(ax2, 'frequency (hz)');
    title(ax2, sprintf('matched non-ied, n = %d participants', ...
        numberOfParticipants));
    colormap(ax2, parula);
    axes(ax2);
    caxis([-groupColorLimit groupColorLimit]);
    colorbar(ax2);
    box(ax2, 'on');

    ax3 = nexttile(layout);
    plot_vector_heatmap(ax3, ...
        spectrogramTime, spectrogramFrequencies, ...
        meanDifference, alphaMask);
    set(ax3, ...
        'YScale', 'log', ...
        'YDir', 'normal', ...
        'YLim', [5 200], ...
        'YTick', [10 100], ...
        'FontSize', 12, ...
        'TickDir', 'out');
    xline(ax3, 0, 'r', 'LineWidth', 1);
    xlabel(ax3, 'time from ied/pseudo-ied (s)');
    ylabel(ax3, 'frequency (hz)');
    title(ax3, 'ied - non-ied, cbpt');
    colormap(ax3, parula);
    axes(ax3);
    caxis([-differenceColorLimit differenceColorLimit]);
    colorbar(ax3);
    box(ax3, 'on');

    title(layout, sprintf([ ...
        'IED7 %s | %s | baseline %.1f to %.1f s | ' ...
        '%d ied and %d pseudo-ied events | %d significant pixels'], ...
        outcome, areaName, ...
        baselineStartSeconds, baselineEndSeconds, ...
        sum(allResults(rr).numberOfIEDEvents), ...
        sum(allResults(rr).numberOfNonIEDEvents), ...
        nnz(significantMask)), ...
        'Interpreter', 'none', ...
        'FontSize', 15, ...
        'FontWeight', 'bold');

    if ~isempty(clusterPvalues)
        fprintf('significant cluster p values: ');
        fprintf('%.4f ', clusterPvalues);
        fprintf('\n');
    else
        fprintf('no significant clusters\n');
    end

    pageNumber = pageNumber + 1;

    % save this brain area in its own pdf too

    safeAreaName = char(areaName);
    safeAreaName = strtrim(safeAreaName);
    safeAreaName = regexprep(safeAreaName, '\s+', '_');
    safeAreaName = regexprep(safeAreaName, '[<>:"/\\|?*]', '_');
    safeAreaName = regexprep(safeAreaName, '_+', '_');

    separatePDF = fullfile(outputFolder, sprintf( ...
        '%s_%s_spectrogram.pdf', ...
        safeAreaName, char(outcome)));

    separateSVG = fullfile(outputFolder, sprintf( ...
        '%s_%s_spectrogram.svg', ...
        safeAreaName, char(outcome)));

    separateFIG = fullfile(outputFolder, sprintf( ...
        '%s_%s_spectrogram.fig', ...
        safeAreaName, char(outcome)));

    exportgraphics(fig, separatePDF, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    print(fig, separateSVG, '-dsvg', '-painters');

    savefig(fig, separateFIG);

    fprintf('saved editable files:\n%s\n%s\n%s\n', ...
        separatePDF, separateSVG, separateFIG);

    % also add it to the combined pdf

    if pageNumber == 1
        exportgraphics(fig, outputPDF, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white');
    else
        exportgraphics(fig, outputPDF, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white', ...
            'Append', true);
    end

    close(fig);

end

if pageNumber == 0
    warning('no pdf was made because no area had matched participants');
else
    fprintf('\nsaved:\n%s\n', outputPDF);
end

function heatmapSurface = plot_vector_heatmap( ...
        ax, timeCenters, frequencyCenters, colorData, alphaData)

    % Draw every time-frequency bin as a vector surface face. Unlike
    % imagesc, this keeps the heatmap cells editable in PDF/SVG exports.

    timeCenters = double(timeCenters(:)');
    frequencyCenters = double(frequencyCenters(:));

    timeEdges = centers_to_edges(timeCenters, false);
    frequencyEdges = centers_to_edges(frequencyCenters, true);

    paddedColorData = [colorData, colorData(:, end)];
    paddedColorData = [paddedColorData; paddedColorData(end, :)];

    [timeGrid, frequencyGrid] = meshgrid(timeEdges, frequencyEdges);

    heatmapSurface = surf(ax, ...
        timeGrid, frequencyGrid, zeros(size(timeGrid)), ...
        paddedColorData, ...
        'FaceColor', 'flat', ...
        'EdgeColor', 'none');

    view(ax, 2);

    if ~isempty(alphaData)
        paddedAlphaData = [alphaData, alphaData(:, end)];
        paddedAlphaData = [paddedAlphaData; paddedAlphaData(end, :)];

        set(heatmapSurface, ...
            'FaceAlpha', 'flat', ...
            'AlphaData', paddedAlphaData, ...
            'AlphaDataMapping', 'none');
    end

end

function edges = centers_to_edges(centers, useLogSpacing)

    centers = double(centers(:)');

    if numel(centers) < 2
        if useLogSpacing
            edges = centers .* [1 / sqrt(2), sqrt(2)];
        else
            edges = centers + [-0.5, 0.5];
        end
        return;
    end

    if useLogSpacing
        logCenters = log(centers);
        logEdges = [ ...
            logCenters(1) - (logCenters(2) - logCenters(1)) / 2, ...
            (logCenters(1:end-1) + logCenters(2:end)) / 2, ...
            logCenters(end) + ...
                (logCenters(end) - logCenters(end-1)) / 2];
        edges = exp(logEdges);
    else
        edges = [ ...
            centers(1) - (centers(2) - centers(1)) / 2, ...
            (centers(1:end-1) + centers(2:end)) / 2, ...
            centers(end) + (centers(end) - centers(end-1)) / 2];
    end

end

function labels = make_string_labels(rawLabels)

    if iscell(rawLabels)

        labels = strings(numel(rawLabels), 1);

        for ii = 1:numel(rawLabels)

            oneLabel = rawLabels{ii};

            if isempty(oneLabel)
                labels(ii) = "Unknown";
            elseif iscell(oneLabel)
                labels(ii) = string(oneLabel{1});
            else
                labels(ii) = string(oneLabel);
            end

        end

    elseif ischar(rawLabels)
        labels = string(cellstr(rawLabels));
    else
        labels = string(rawLabels);
    end

    labels = labels(:);

end

function areaNames = clean_area_names(areaNames)

    areaNames = string(areaNames(:));
    areaNames = strip(areaNames);

    areaNames( ...
        ismissing(areaNames) | strlength(areaNames) == 0) = ...
        "Unknown";

    areaNames = regexprep(areaNames, ...
        '^\s*\(?\s*(Left|Right|LH|RH|L|R)\s*\)?[_\-\s]+', ...
        '', 'ignorecase');

    areaNames = regexprep(areaNames, ...
        '[_\-\s]+\(?\s*(Left|Right|LH|RH|L|R)\s*\)?\s*$', ...
        '', 'ignorecase');

    areaNames = regexprep(areaNames, '[_\-]+', ' ');
    areaNames = regexprep(areaNames, '\s+', ' ');
    areaNames = strip(areaNames);

end

function [matchedIED, matchedNonIED] = ...
        match_trials( ...
            IEDtrials, nonIEDtrials, ...
            colorCode, bankedTrials, duration)

    matchedIED = [];
    matchedNonIED = [];
    availableNonIED = nonIEDtrials(:);

    % do the harder-to-match trials first

    numberOfExactMatches = zeros(length(IEDtrials), 1);

    for ii = 1:length(IEDtrials)

        thisTrial = IEDtrials(ii);

        numberOfExactMatches(ii) = sum( ...
            colorCode(availableNonIED) == colorCode(thisTrial) & ...
            bankedTrials(availableNonIED) == bankedTrials(thisTrial));

    end

    [~, trialOrder] = sort(numberOfExactMatches, 'ascend');
    IEDtrials = IEDtrials(trialOrder);

    for ii = 1:length(IEDtrials)

        if isempty(availableNonIED)
            break;
        end

        thisTrial = IEDtrials(ii);

        exactPool = availableNonIED( ...
            colorCode(availableNonIED) == colorCode(thisTrial) & ...
            bankedTrials(availableNonIED) == bankedTrials(thisTrial));

        if isempty(exactPool)
            exactPool = availableNonIED( ...
                colorCode(availableNonIED) == colorCode(thisTrial));
        end

        if isempty(exactPool)
            exactPool = availableNonIED;
        end

        [~, closestIndex] = min( ...
            abs(duration(exactPool) - duration(thisTrial)));

        chosenNonIED = exactPool(closestIndex);

        matchedIED(end + 1, 1) = thisTrial;
        matchedNonIED(end + 1, 1) = chosenNonIED;

        availableNonIED(availableNonIED == chosenNonIED) = [];

    end

end

function oneTF = make_one_tf( ...
        signal, centerSample, Fs, ...
        preIEDSeconds, postIEDSeconds, extraSeconds, ...
        baselineStartSeconds, baselineEndSeconds, ...
        spectrogramFrequencies, spectrogramTime, ...
        amplitudeThreshold)

    oneTF = [];

    longPreSamples = round((preIEDSeconds + extraSeconds) * Fs);
    longPostSamples = round((postIEDSeconds + extraSeconds) * Fs);

    firstSample = centerSample - longPreSamples;
    lastSample = centerSample + longPostSamples;

    if firstSample < 1 || lastSample > length(signal)
        return;
    end

    oneEpoch = signal(firstSample:lastSample);
    oneEpoch = detrend(oneEpoch, 0);

    if any(~isfinite(oneEpoch)) || ...
            max(abs(oneEpoch)) > amplitudeThreshold
        return;
    end

    % this is the same basewaveerp idea as the older code

    [wave, period] = basewaveERP( ...
        oneEpoch(:)', Fs, ...
        min(spectrogramFrequencies), ...
        max(spectrogramFrequencies), 6, 0);

    [waveletFrequencies, frequencyOrder] = sort(1 ./ period);
    wave = wave(frequencyOrder, :);

    waveletTime = linspace( ...
        -(preIEDSeconds + extraSeconds), ...
        postIEDSeconds + extraSeconds, ...
        length(oneEpoch));

    power = abs(wave).^2;

    baselineRows = ...
        waveletTime >= baselineStartSeconds & ...
        waveletTime <= baselineEndSeconds;

    baselinePower = mean(power(:, baselineRows), 2, 'omitnan');

    if any(~isfinite(baselinePower)) || any(baselinePower <= 0)
        return;
    end

    powerDB = 10 .* log10(power ./ baselinePower);

    [oldTimeGrid, oldFrequencyGrid] = meshgrid( ...
        waveletTime, waveletFrequencies);

    [newTimeGrid, newFrequencyGrid] = meshgrid( ...
        spectrogramTime, spectrogramFrequencies);

    oneTF = interp2( ...
        oldTimeGrid, oldFrequencyGrid, powerDB, ...
        newTimeGrid, newFrequencyGrid, ...
        'linear', NaN);

end

function [significantMask, significantClusterPvalues] = ...
        run_paired_cbpt( ...
            IEDall, nonIEDall, timeVector, ...
            numberOfPermutations, alphaBin, alphaCluster, ...
            minimumClusterSize, minimumClusterDurationMs)

    differenceData = IEDall - nonIEDall;
    [realT, realP] = paired_t_map(differenceData);

    positiveMask = realP < alphaBin & realT > 0;
    negativeMask = realP < alphaBin & realT < 0;

    [positiveStats, positivePixels] = ...
        get_cluster_stats( ...
            positiveMask, realT, minimumClusterSize, "positive");

    [negativeStats, negativePixels] = ...
        get_cluster_stats( ...
            negativeMask, realT, minimumClusterSize, "negative");

    numberOfParticipants = size(differenceData, 3);
    maximumNullCluster = zeros(numberOfPermutations, 1);

    for pp = 1:numberOfPermutations

        randomSigns = ...
            2 .* randi([0 1], 1, 1, numberOfParticipants) - 1;

        permutedDifference = differenceData .* randomSigns;
        [permutedT, permutedP] = paired_t_map(permutedDifference);

        permutedPositive = ...
            permutedP < alphaBin & permutedT > 0;

        permutedNegative = ...
            permutedP < alphaBin & permutedT < 0;

        positiveNullStats = get_cluster_stats( ...
            permutedPositive, permutedT, ...
            minimumClusterSize, "positive");

        negativeNullStats = get_cluster_stats( ...
            permutedNegative, permutedT, ...
            minimumClusterSize, "negative");

        allNullStats = [positiveNullStats; negativeNullStats];

        if ~isempty(allNullStats)
            maximumNullCluster(pp) = max(allNullStats);
        end

        if mod(pp, 100) == 0
            fprintf('permutation %d/%d\n', ...
                pp, numberOfPermutations);
        end

    end

    significantMask = false(size(realT));
    significantClusterPvalues = [];
    timeStepMs = median(diff(timeVector)) * 1000;

    allRealStats = [positiveStats; negativeStats];
    allRealPixels = [positivePixels; negativePixels];

    for cc = 1:length(allRealStats)

        pixels = allRealPixels{cc};

        [~, timeIndices] = ind2sub(size(realT), pixels);

        clusterDurationMs = ...
            (max(timeVector(timeIndices)) - ...
            min(timeVector(timeIndices))) * 1000 + timeStepMs;

        clusterP = ...
            (1 + sum(maximumNullCluster >= allRealStats(cc))) / ...
            (numberOfPermutations + 1);

        if clusterP < alphaCluster && ...
                clusterDurationMs >= minimumClusterDurationMs

            significantMask(pixels) = true;
            significantClusterPvalues(end + 1, 1) = clusterP;

        end

    end

end

function [tMap, pMap] = paired_t_map(differenceData)

    numberOfPairs = sum(isfinite(differenceData), 3);
    meanDifference = mean(differenceData, 3, 'omitnan');
    standardDeviation = std(differenceData, 0, 3, 'omitnan');
    standardError = standardDeviation ./ sqrt(numberOfPairs);

    tMap = meanDifference ./ standardError;
    degreesFreedom = numberOfPairs - 1;
    pMap = 2 .* tcdf(-abs(tMap), degreesFreedom);

    badPixels = ...
        numberOfPairs < 2 | ...
        ~isfinite(tMap) | ...
        standardError <= 0;

    tMap(badPixels) = NaN;
    pMap(badPixels) = NaN;

end

function [clusterStats, clusterPixels] = ...
        get_cluster_stats( ...
            candidateMask, tMap, minimumClusterSize, direction)

    candidateMask(~isfinite(tMap)) = false;
    connectedClusters = bwconncomp(candidateMask, 8);

    clusterStats = [];
    clusterPixels = {};

    for cc = 1:connectedClusters.NumObjects

        pixels = connectedClusters.PixelIdxList{cc};

        if numel(pixels) < minimumClusterSize
            continue;
        end

        values = tMap(pixels);

        if direction == "positive"
            clusterMass = sum(values(values > 0), 'omitnan');
        else
            clusterMass = sum(abs(values(values < 0)), 'omitnan');
        end

        if isfinite(clusterMass) && clusterMass > 0
            clusterStats(end + 1, 1) = clusterMass;
            clusterPixels{end + 1, 1} = pixels;
        end

    end

end
