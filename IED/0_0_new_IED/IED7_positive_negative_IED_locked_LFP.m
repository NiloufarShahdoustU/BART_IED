% compare ied-locked spectrograms for positive and negative ied7 results

clear;
clc;
close all;

% folders

lfpiedFolder = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

rawDataFolder = ...
    'D:\Nill\data\BART_preprocessed\';

ied7ResultsFile = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED\' ...
     'IED7_Cox_IT_RT_BR_postIED_by_brain_area\' ...
     'IT_RT_BR_brain_area_cox_all_results.mat'];

outputFolder = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED\' ...
     'IED7_positive_negative_IED_locked_LFP\'];

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

outputPDF = fullfile(outputFolder, ...
    'IED7_positive_negative_IED_locked_LFP_simple.pdf');

if exist(outputPDF, 'file')
    delete(outputPDF);
end

% settings

preIEDSeconds = 1;
postIEDSeconds = 1;
extraSeconds = 0.5;

baselineStartSeconds = -0.8;
baselineEndSeconds = -0.3;

% log-spaced frequencies make the y axis look like the attached figure
spectrogramFrequencies = logspace(log10(5), log10(200), 60);
spectrogramTime = -0.9:0.02:0.9;

minimumDistanceBetweenIEDsSeconds = 0.5;
amplitudeThreshold = 5000;
minimumIEDsPerChannel = 3;
maximumRTSeconds = 20;

numberOfPermutations = 1000;
alphaBin = 0.10;
alphaCluster = 0.05;
minimumClusterSize = 5;
minimumClusterDurationMs = 50;
minimumSignificantPixelRunMs = 40;
nonSignificantAlpha = 0.20;

rng(1);

if exist('basewaveERP', 'file') ~= 2
    error('basewaveERP.m was not found');
end

if exist('bwconncomp', 'file') ~= 2
    error('bwconncomp was not found');
end

% find the significant positive and negative areas

ied7 = load(ied7ResultsFile);
outcomeNames = ["RT", "IT", "BR"];

positiveAreas = cell(1, 3);
negativeAreas = cell(1, 3);

for oo = 1:3

    outcome = outcomeNames(oo);
    resultTable = ied7.(char(outcome + "Results"));

    useRows = ...
        resultTable.status == "fitted" & ...
        resultTable.significantUncorrected & ...
        isfinite(resultTable.beta_logHazard);

    positiveAreas{oo} = string(resultTable.anatomicalArea( ...
        useRows & resultTable.beta_logHazard > 0));

    negativeAreas{oo} = string(resultTable.anatomicalArea( ...
        useRows & resultTable.beta_logHazard < 0));

    positiveAreas{oo} = strip(positiveAreas{oo});
    negativeAreas{oo} = strip(negativeAreas{oo});

    positiveAreas{oo}( ...
        ismissing(positiveAreas{oo}) | strlength(positiveAreas{oo}) == 0) = ...
        "Unknown";

    negativeAreas{oo}( ...
        ismissing(negativeAreas{oo}) | strlength(negativeAreas{oo}) == 0) = ...
        "Unknown";

    fprintf('\n%s positive areas:\n', outcome);
    disp(positiveAreas{oo});

    fprintf('%s negative areas:\n', outcome);
    disp(negativeAreas{oo});

    positiveAll{oo} = zeros( ...
        length(spectrogramFrequencies), length(spectrogramTime), 0);

    negativeAll{oo} = zeros( ...
        length(spectrogramFrequencies), length(spectrogramTime), 0);

    positiveIEDcount{oo} = [];
    negativeIEDcount{oo} = [];

end

% load the participants

fileList = dir(fullfile(lfpiedFolder, '*.LFPIED.mat'));

for pt = 1:length(fileList)

    fileName = fileList(pt).name;
    fileParts = strsplit(fileName, '.');
    patientID = string(fileParts{1});

    fprintf('\nparticipant %d/%d: %s\n', ...
        pt, length(fileList), patientID);

    loadedData = load(fullfile(lfpiedFolder, fileName));

    if ~isfield(loadedData, 'LFPIED')
        fprintf('skipped: LFPIED was missing\n');
        continue;
    end

    LFPIED = loadedData.LFPIED;

    patientRawFolder = fullfile(rawDataFolder, char(patientID), 'Data');
    ns2List = dir(fullfile(patientRawFolder, '*.ns2'));
    nevList = dir(fullfile(patientRawFolder, '*.nev'));

    if length(ns2List) ~= 1 || length(nevList) ~= 1
        fprintf('skipped: need one ns2 and one nev file\n');
        continue;
    end

    % response times are needed for the it and br windows

    NEV = openNEV(fullfile(nevList.folder, nevList.name), 'overwrite');
    trigs = NEV.Data.SerialDigitalIO.UnparsedData;
    trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
    responseTimes = trigTimes(trigs == 23);
    clear NEV;

    % load lfp

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

    % turn the area labels into strings

    rawLabels = LFPIED.anatomicalLocs;

    if iscell(rawLabels)

        allAreaLabels = strings(numel(rawLabels), 1);

        for ii = 1:numel(rawLabels)

            oneLabel = rawLabels{ii};

            if isempty(oneLabel)
                allAreaLabels(ii) = "Unknown";
            elseif iscell(oneLabel)
                allAreaLabels(ii) = string(oneLabel{1});
            else
                allAreaLabels(ii) = string(oneLabel);
            end

        end

    elseif ischar(rawLabels)
        allAreaLabels = string(cellstr(rawLabels));
    else
        allAreaLabels = string(rawLabels);
    end

    allAreaLabels = allAreaLabels(:);

    if length(allAreaLabels) >= max(selectedChans)
        selectedAreaLabels = allAreaLabels(selectedChans);
    else
        selectedAreaLabels = allAreaLabels;
    end

    selectedAreaLabels = strip(selectedAreaLabels);
    selectedAreaLabels( ...
        ismissing(selectedAreaLabels) | strlength(selectedAreaLabels) == 0) = ...
        "Unknown";

    % remove left and right from the area names

    selectedAreaLabels = regexprep(selectedAreaLabels, ...
        '^\s*\(?\s*(Left|Right|LH|RH|L|R)\s*\)?[_\-\s]+', ...
        '', 'ignorecase');

    selectedAreaLabels = regexprep(selectedAreaLabels, ...
        '[_\-\s]+\(?\s*(Left|Right|LH|RH|L|R)\s*\)?\s*$', ...
        '', 'ignorecase');

    selectedAreaLabels = regexprep(selectedAreaLabels, '[_\-]+', ' ');
    selectedAreaLabels = regexprep(selectedAreaLabels, '\s+', ' ');
    selectedAreaLabels = strip(selectedAreaLabels);

    % run rt, it, and br

    for oo = 1:3

        outcome = outcomeNames(oo);

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
            length(startTimes)]);

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

        % it uses banked trials only
        if outcome == "IT"
            validTrials = validTrials & bankedTrials == 1;
        end

        % br uses both banked and lost trials
        if outcome == "BR"
            validTrials = validTrials & ismember(bankedTrials, [0 1]);
        end

        if isempty(IEDrows) || ~any(validTrials)
            continue;
        end

        % group 1 is positive and group 2 is negative

        for gg = 1:2

            if gg == 1
                groupAreas = positiveAreas{oo};
            else
                groupAreas = negativeAreas{oo};
            end

            if isempty(groupAreas)
                continue;
            end

            participantAreaTF = zeros( ...
                length(spectrogramFrequencies), ...
                length(spectrogramTime), 0);

            participantIEDcount = 0;

            for aa = 1:length(groupAreas)

                areaName = groupAreas(aa);
                areaChannels = find(selectedAreaLabels == areaName);

                if isempty(areaChannels)
                    continue;
                end

                areaChannelTF = zeros( ...
                    length(spectrogramFrequencies), ...
                    length(spectrogramTime), 0);

                areaIEDcount = 0;

                for cc = 1:length(areaChannels)

                    localChannel = areaChannels(cc);

                    channelRows = ...
                        isfinite(IEDrows(:, 1)) & ...
                        isfinite(IEDrows(:, 2)) & ...
                        isfinite(IEDrows(:, 3)) & ...
                        round(IEDrows(:, 2)) == localChannel;

                    thisIED = IEDrows(channelRows, :);

                    if isempty(thisIED)
                        continue;
                    end

                    trialNumber = round(thisIED(:, 1));

                    keepRows = ...
                        trialNumber >= 1 & ...
                        trialNumber <= nTrials;

                    goodRowNumbers = find(keepRows);

                    if ~isempty(goodRowNumbers)
                        keepRows(goodRowNumbers) = ...
                            validTrials(trialNumber(goodRowNumbers));
                    end

                    thisIED = thisIED(keepRows, :);

                    if isempty(thisIED)
                        continue;
                    end

                    trialNumber = round(thisIED(:, 1));

                    absoluteSample = ...
                        floor(Fs .* startTimes(trialNumber)) + ...
                        round(thisIED(:, 3)) - 1;

                    % remove ieds that are too close to another ied

                    separatedRows = true(length(absoluteSample), 1);
                    minimumDistanceSamples = ...
                        minimumDistanceBetweenIEDsSeconds * Fs;

                    for ii = 1:length(absoluteSample)

                        sameTrial = trialNumber == trialNumber(ii);
                        distance = abs(absoluteSample - absoluteSample(ii));

                        if any(sameTrial & distance > 0 & ...
                                distance < minimumDistanceSamples)
                            separatedRows(ii) = false;
                        end

                    end

                    absoluteSample = absoluteSample(separatedRows);

                    if length(absoluteSample) < minimumIEDsPerChannel
                        continue;
                    end

                    rawChannelNumber = selectedChans(localChannel);

                    if rawChannelNumber < 1 || ...
                            rawChannelNumber > size(rawLFP, 1)
                        continue;
                    end

                    signal = double(rawLFP(rawChannelNumber, :));

                    channelEventTF = zeros( ...
                        length(spectrogramFrequencies), ...
                        length(spectrogramTime), 0);

                    for ee = 1:length(absoluteSample)

                        longPreSamples = round( ...
                            (preIEDSeconds + extraSeconds) * Fs);

                        longPostSamples = round( ...
                            (postIEDSeconds + extraSeconds) * Fs);

                        firstSample = absoluteSample(ee) - longPreSamples;
                        lastSample = absoluteSample(ee) + longPostSamples;

                        if firstSample < 1 || lastSample > length(signal)
                            continue;
                        end

                        oneEpoch = signal(firstSample:lastSample);
                        oneEpoch = detrend(oneEpoch, 0);

                        if any(~isfinite(oneEpoch)) || ...
                                max(abs(oneEpoch)) > amplitudeThreshold
                            continue;
                        end

                        % basewaveerp makes the spectrogram

                        [wave, period] = basewaveERP( ...
                            oneEpoch(:)', Fs, ...
                            min(spectrogramFrequencies), ...
                            max(spectrogramFrequencies), 6, 0);

                        [waveletFrequencies, frequencyOrder] = ...
                            sort(1 ./ period);

                        wave = wave(frequencyOrder, :);

                        waveletTime = linspace( ...
                            -(preIEDSeconds + extraSeconds), ...
                            postIEDSeconds + extraSeconds, ...
                            length(oneEpoch));

                        power = abs(wave).^2;

                        baselineRows = ...
                            waveletTime >= baselineStartSeconds & ...
                            waveletTime <= baselineEndSeconds;

                        baselinePower = ...
                            mean(power(:, baselineRows), 2, 'omitnan');

                        if any(~isfinite(baselinePower)) || ...
                                any(baselinePower <= 0)
                            continue;
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

                        channelEventTF(:, :, end + 1) = oneTF;

                    end

                    if size(channelEventTF, 3) < minimumIEDsPerChannel
                        continue;
                    end

                    areaChannelTF(:, :, end + 1) = ...
                        mean(channelEventTF, 3, 'omitnan');

                    areaIEDcount = ...
                        areaIEDcount + size(channelEventTF, 3);

                end

                if isempty(areaChannelTF)
                    continue;
                end

                participantAreaTF(:, :, end + 1) = ...
                    mean(areaChannelTF, 3, 'omitnan');

                participantIEDcount = ...
                    participantIEDcount + areaIEDcount;

            end

            if isempty(participantAreaTF)
                continue;
            end

            participantTF = mean(participantAreaTF, 3, 'omitnan');

            if gg == 1
                positiveAll{oo}(:, :, end + 1) = participantTF;
                positiveIEDcount{oo}(end + 1, 1) = participantIEDcount;
            else
                negativeAll{oo}(:, :, end + 1) = participantTF;
                negativeIEDcount{oo}(end + 1, 1) = participantIEDcount;
            end

        end

    end

    clear NSX rawLFP LFPIED loadedData;

end

%% make one pdf page for rt, it, and br

for oo = 1:3

    outcome = outcomeNames(oo);
    positiveTF = positiveAll{oo};
    negativeTF = negativeAll{oo};

    if isempty(positiveTF)
        meanPositiveTF = NaN( ...
            length(spectrogramFrequencies), ...
            length(spectrogramTime));
    else
        meanPositiveTF = mean(positiveTF, 3, 'omitnan');
    end

    if isempty(negativeTF)
        meanNegativeTF = NaN( ...
            length(spectrogramFrequencies), ...
            length(spectrogramTime));
    else
        meanNegativeTF = mean(negativeTF, 3, 'omitnan');
    end

    % cluster based permutation test between positive and negative

    numberOfPositiveParticipants = size(positiveTF, 3);
    numberOfNegativeParticipants = size(negativeTF, 3);

    significantMask = false(size(meanPositiveTF));

    if numberOfPositiveParticipants >= 2 && ...
            numberOfNegativeParticipants >= 2

        allParticipantTF = cat(3, positiveTF, negativeTF);

        positiveParticipantIdx = ...
            1:numberOfPositiveParticipants;

        negativeParticipantIdx = ...
            numberOfPositiveParticipants + ...
            (1:numberOfNegativeParticipants);

        [~, significantPositive, significantNegative] = ...
            run_ttest_cbpt_sign_separated( ...
            allParticipantTF, ...
            positiveParticipantIdx, ...
            negativeParticipantIdx, ...
            spectrogramTime, ...
            numberOfPermutations, ...
            alphaBin, ...
            alphaCluster, ...
            minimumClusterSize, ...
            minimumClusterDurationMs);

        significantPositive = remove_short_time_runs( ...
            significantPositive, spectrogramTime, ...
            minimumSignificantPixelRunMs);

        significantNegative = remove_short_time_runs( ...
            significantNegative, spectrogramTime, ...
            minimumSignificantPixelRunMs);

        significantMask = ...
            significantPositive | significantNegative;

    end

    fprintf('%s: %d significant time-frequency pixels\n', ...
        outcome, nnz(significantMask));

    rowValues = [meanPositiveTF(:); meanNegativeTF(:)];
    rowValues = rowValues(isfinite(rowValues));

    if isempty(rowValues)
        rowColorLimits = [-1 1];
    else
        rowColorLimits = [ ...
            prctile(rowValues, 5), ...
            prctile(rowValues, 95)];
    end

    if any(~isfinite(rowColorLimits)) || ...
            rowColorLimits(1) >= rowColorLimits(2)
        rowColorLimits = [-1 1];
    end

    alphaMask = nonSignificantAlpha .* ...
        ones(size(significantMask));

    alphaMask(significantMask) = 1;

    positiveAreaText = strjoin(positiveAreas{oo}, ', ');
    negativeAreaText = strjoin(negativeAreas{oo}, ', ');

    if strlength(positiveAreaText) == 0
        positiveAreaText = "none";
    end

    if strlength(negativeAreaText) == 0
        negativeAreaText = "none";
    end

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [50 50 1800 900]);

    layout = tiledlayout(fig, 1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(layout);
    positiveImage = imagesc(ax1, ...
        spectrogramTime, spectrogramFrequencies, ...
        meanPositiveTF);
    set(positiveImage, 'AlphaData', alphaMask);
    set(ax1, 'YScale', 'log', ...
        'YDir', 'normal', ...
        'YLim', [5 200], ...
        'YTick', [10 100], ...
        'FontSize', 12, ...
        'TickDir', 'out');
    xline(ax1, 0, 'r', 'LineWidth', 1);
    xlabel(ax1, 'time from ied peak (s)');
    ylabel(ax1, 'frequency (hz)');
    title(ax1, { ...
        sprintf('positive coefficient, n = %d', size(positiveTF, 3)), ...
        char(positiveAreaText)}, ...
        'Interpreter', 'none');
    colormap(ax1, parula);
    axes(ax1);
    caxis(rowColorLimits);
    colorbar(ax1);
    box(ax1, 'on');

    ax2 = nexttile(layout);
    negativeImage = imagesc(ax2, ...
        spectrogramTime, spectrogramFrequencies, ...
        meanNegativeTF);
    set(negativeImage, 'AlphaData', alphaMask);
    set(ax2, 'YScale', 'log', ...
        'YDir', 'normal', ...
        'YLim', [5 200], ...
        'YTick', [10 100], ...
        'FontSize', 12, ...
        'TickDir', 'out');
    xline(ax2, 0, 'r', 'LineWidth', 1);
    xlabel(ax2, 'time from ied peak (s)');
    ylabel(ax2, 'frequency (hz)');
    title(ax2, { ...
        sprintf('negative coefficient, n = %d', size(negativeTF, 3)), ...
        char(negativeAreaText)}, ...
        'Interpreter', 'none');
    colormap(ax2, parula);
    axes(ax2);
    caxis(rowColorLimits);
    colorbar(ax2);
    box(ax2, 'on');

    title(layout, sprintf([ ...
        'IED7 %s positive-negative spectrogram CBPT | ' ...
        'baseline %.1f to %.1f s | ' ...
        'positive %d ieds | negative %d ieds'], ...
        outcome, baselineStartSeconds, baselineEndSeconds, ...
        sum(positiveIEDcount{oo}), sum(negativeIEDcount{oo})), ...
        'FontSize', 16, ...
        'FontWeight', 'bold');

    if oo == 1
        exportgraphics(fig, outputPDF, ...
            'ContentType', 'image', ...
            'Resolution', 220);
    else
        exportgraphics(fig, outputPDF, ...
            'ContentType', 'image', ...
            'Resolution', 220, ...
            'Append', true);
    end

    close(fig);

end

function [significantMask, ...
          significantPositive, ...
          significantNegative] = ...
          run_ttest_cbpt_sign_separated( ...
            allData, ...
            condition1Idx, ...
            condition2Idx, ...
            timeVector, ...
            numberOfPermutations, ...
            alphaBin, ...
            alphaCluster, ...
            minimumClusterSize, ...
            minimumClusterDurationMs)

    numberOfFrequencies = size(allData, 1);
    numberOfTimes = size(allData, 2);
    numberOfObservations = size(allData, 3);

    [tMap, pMap] = compute_ttest_maps_fast( ...
        allData, condition1Idx, condition2Idx);

    candidatePositive = pMap < alphaBin & tMap > 0;
    candidateNegative = pMap < alphaBin & tMap < 0;

    candidatePositive(isnan(candidatePositive)) = false;
    candidateNegative(isnan(candidateNegative)) = false;

    [realStatsPositive, realPixelsPositive] = ...
        get_cluster_stats( ...
        candidatePositive, tMap, ...
        minimumClusterSize, "positive");

    [realStatsNegative, realPixelsNegative] = ...
        get_cluster_stats( ...
        candidateNegative, tMap, ...
        minimumClusterSize, "negative");

    labels = zeros(numberOfObservations, 1);
    labels(condition1Idx) = 1;
    labels(condition2Idx) = 2;

    validObservations = find(labels > 0);
    validLabels = labels(validObservations);

    maximumNullPositive = ...
        zeros(numberOfPermutations, 1);

    maximumNullNegative = ...
        zeros(numberOfPermutations, 1);

    for permutationNumber = 1:numberOfPermutations

        shuffledLabels = ...
            validLabels(randperm(numel(validLabels)));

        permutedCondition1Idx = ...
            validObservations(shuffledLabels == 1);

        permutedCondition2Idx = ...
            validObservations(shuffledLabels == 2);

        [permutedTMap, permutedPMap] = ...
            compute_ttest_maps_fast( ...
            allData, ...
            permutedCondition1Idx, ...
            permutedCondition2Idx);

        permutedPositive = ...
            permutedPMap < alphaBin & permutedTMap > 0;

        permutedNegative = ...
            permutedPMap < alphaBin & permutedTMap < 0;

        permutedPositive(isnan(permutedPositive)) = false;
        permutedNegative(isnan(permutedNegative)) = false;

        permutedStatsPositive = get_cluster_stats( ...
            permutedPositive, permutedTMap, ...
            minimumClusterSize, "positive");

        permutedStatsNegative = get_cluster_stats( ...
            permutedNegative, permutedTMap, ...
            minimumClusterSize, "negative");

        if ~isempty(permutedStatsPositive)
            maximumNullPositive(permutationNumber) = ...
                max(permutedStatsPositive);
        end

        if ~isempty(permutedStatsNegative)
            maximumNullNegative(permutationNumber) = ...
                max(permutedStatsNegative);
        end

        if mod(permutationNumber, 100) == 0
            fprintf('permutation %d / %d complete\n', ...
                permutationNumber, numberOfPermutations);
        end

    end

    positiveClusterThreshold = prctile( ...
        maximumNullPositive, 100 * (1 - alphaCluster));

    negativeClusterThreshold = prctile( ...
        maximumNullNegative, 100 * (1 - alphaCluster));

    significantPositive = ...
        false(numberOfFrequencies, numberOfTimes);

    significantNegative = ...
        false(numberOfFrequencies, numberOfTimes);

    timeStepMs = median(diff(timeVector)) * 1000;

    for clusterNumber = 1:numel(realStatsPositive)

        clusterPixels = ...
            realPixelsPositive{clusterNumber};

        [~, timeIndices] = ind2sub( ...
            [numberOfFrequencies, numberOfTimes], ...
            clusterPixels);

        clusterDurationMs = ...
            (max(timeVector(timeIndices)) - ...
            min(timeVector(timeIndices))) * 1000 + timeStepMs;

        if realStatsPositive(clusterNumber) > ...
                positiveClusterThreshold && ...
                clusterDurationMs >= minimumClusterDurationMs

            significantPositive(clusterPixels) = true;

        end

    end

    for clusterNumber = 1:numel(realStatsNegative)

        clusterPixels = ...
            realPixelsNegative{clusterNumber};

        [~, timeIndices] = ind2sub( ...
            [numberOfFrequencies, numberOfTimes], ...
            clusterPixels);

        clusterDurationMs = ...
            (max(timeVector(timeIndices)) - ...
            min(timeVector(timeIndices))) * 1000 + timeStepMs;

        if realStatsNegative(clusterNumber) > ...
                negativeClusterThreshold && ...
                clusterDurationMs >= minimumClusterDurationMs

            significantNegative(clusterPixels) = true;

        end

    end

    significantMask = ...
        significantPositive | significantNegative;

    fprintf('positive cluster threshold = %.4f\n', ...
        positiveClusterThreshold);

    fprintf('negative cluster threshold = %.4f\n', ...
        negativeClusterThreshold);

end

function [clusterStats, clusterPixels] = ...
        get_cluster_stats( ...
        candidateMask, tMap, minimumClusterSize, direction)

    connectedClusters = bwconncomp(candidateMask, 8);

    clusterStats = [];
    clusterPixels = {};

    for clusterNumber = 1:connectedClusters.NumObjects

        pixels = ...
            connectedClusters.PixelIdxList{clusterNumber};

        if numel(pixels) < minimumClusterSize
            continue;
        end

        values = tMap(pixels);

        if direction == "positive"
            clusterMass = ...
                sum(values(values > 0), 'omitnan');
        else
            clusterMass = ...
                sum(abs(values(values < 0)), 'omitnan');
        end

        if ~isnan(clusterMass) && clusterMass > 0
            clusterStats(end + 1, 1) = clusterMass;
            clusterPixels{end + 1, 1} = pixels;
        end

    end

end

function [tMap, pMap] = compute_ttest_maps_fast( ...
        allData, condition1Idx, condition2Idx)

    condition1Data = allData(:, :, condition1Idx);
    condition2Data = allData(:, :, condition2Idx);

    numberCondition1 = ...
        sum(~isnan(condition1Data), 3);

    numberCondition2 = ...
        sum(~isnan(condition2Data), 3);

    meanCondition1 = ...
        mean(condition1Data, 3, 'omitnan');

    meanCondition2 = ...
        mean(condition2Data, 3, 'omitnan');

    varianceCondition1 = ...
        var(condition1Data, 0, 3, 'omitnan');

    varianceCondition2 = ...
        var(condition2Data, 0, 3, 'omitnan');

    standardError = sqrt( ...
        varianceCondition1 ./ numberCondition1 + ...
        varianceCondition2 ./ numberCondition2);

    tMap = ...
        (meanCondition1 - meanCondition2) ./ standardError;

    degreesNumerator = ( ...
        varianceCondition1 ./ numberCondition1 + ...
        varianceCondition2 ./ numberCondition2) .^ 2;

    degreesDenominator = ...
        (varianceCondition1 ./ numberCondition1) .^ 2 ./ ...
        (numberCondition1 - 1) + ...
        (varianceCondition2 ./ numberCondition2) .^ 2 ./ ...
        (numberCondition2 - 1);

    degreesFreedom = ...
        degreesNumerator ./ degreesDenominator;

    pMap = ...
        2 .* tcdf(-abs(tMap), degreesFreedom);

    badValues = ...
        numberCondition1 < 2 | ...
        numberCondition2 < 2 | ...
        standardError == 0 | ...
        isnan(standardError) | ...
        isnan(degreesFreedom) | ...
        degreesFreedom <= 0;

    tMap(badValues) = NaN;
    pMap(badValues) = NaN;

end

function cleanMask = remove_short_time_runs( ...
        significantMask, timeVector, minimumRunMs)

    cleanMask = false(size(significantMask));
    timeStepMs = median(diff(timeVector)) * 1000;

    for frequencyNumber = 1:size(significantMask, 1)

        rowMask = significantMask(frequencyNumber, :);
        connectedRuns = bwconncomp(rowMask, 4);

        for runNumber = 1:connectedRuns.NumObjects

            pixels = ...
                connectedRuns.PixelIdxList{runNumber};

            runDurationMs = ...
                (max(timeVector(pixels)) - ...
                min(timeVector(pixels))) * 1000 + timeStepMs;

            if runDurationMs >= minimumRunMs
                cleanMask(frequencyNumber, pixels) = true;
            end

        end

    end

end
