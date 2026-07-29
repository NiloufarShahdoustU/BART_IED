% compare ied-locked spectrograms for positive and negative ied7 results

clear;
clc;
close all;

% folders

lfpiedFolder = ...
    'D:\Nill\data\BART\0_0_new_IED_last_1000_ms\IED1_find_number_of_IEDs\';

rawDataFolder = ...
    'D:\Nill\data\BART_preprocessed\';

ied7ResultsFile = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED_last_1000_ms\' ...
     'IED7_Cox_IT_RT_BR_postIED_by_brain_area\' ...
     'IT_RT_BR_brain_area_cox_all_results.mat'];

outputFolder = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED_last_1000_ms\' ...
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
% a slightly softer threshold for forming clusters
% the final cluster still needs p < 0.05
clusterStartingP = 0.5;
clusterSignificanceP = 0.05;
nonSignificantAlpha = 0.18;

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

    differenceTF = meanPositiveTF - meanNegativeTF;

    positiveAndNegativeValues = abs([ ...
        meanPositiveTF(isfinite(meanPositiveTF)); ...
        meanNegativeTF(isfinite(meanNegativeTF))]);

    if isempty(positiveAndNegativeValues)
        groupColorLimit = 1;
    else
        groupColorLimit = prctile(positiveAndNegativeValues, 98);
    end

    differenceValues = abs(differenceTF(isfinite(differenceTF)));

    if isempty(differenceValues)
        differenceColorLimit = 1;
    else
        differenceColorLimit = prctile(differenceValues, 98);
    end

    if ~isfinite(groupColorLimit) || groupColorLimit <= 0
        groupColorLimit = 1;
    end

    if ~isfinite(differenceColorLimit) || differenceColorLimit <= 0
        differenceColorLimit = 1;
    end

    % cluster based permutation test for the difference plot

    significantDifference = false(size(differenceTF));

    numberOfPositiveParticipants = size(positiveTF, 3);
    numberOfNegativeParticipants = size(negativeTF, 3);

    if numberOfPositiveParticipants >= 2 && ...
            numberOfNegativeParticipants >= 2

        positiveRows = reshape(permute(positiveTF, [3 1 2]), ...
            numberOfPositiveParticipants, []);

        negativeRows = reshape(permute(negativeTF, [3 1 2]), ...
            numberOfNegativeParticipants, []);

        [~, observedP, ~, observedStats] = ttest2( ...
            positiveRows, negativeRows, ...
            'Vartype', 'unequal');

        observedP = reshape(observedP, size(differenceTF));
        observedT = reshape(observedStats.tstat, size(differenceTF));

        observedClusters = cell(1, 2);
        observedClusters{1} = bwconncomp( ...
            observedP < clusterStartingP & observedT > 0, 4);

        observedClusters{2} = bwconncomp( ...
            observedP < clusterStartingP & observedT < 0, 4);

        allRows = [positiveRows; negativeRows];
        totalParticipants = size(allRows, 1);
        largestPermutedCluster = zeros(numberOfPermutations, 1);

        rng(1);

        for pp = 1:numberOfPermutations

            shuffledRows = randperm(totalParticipants);

            permutedPositive = allRows( ...
                shuffledRows(1:numberOfPositiveParticipants), :);

            permutedNegative = allRows( ...
                shuffledRows(numberOfPositiveParticipants + 1:end), :);

            [~, permutedP, ~, permutedStats] = ttest2( ...
                permutedPositive, permutedNegative, ...
                'Vartype', 'unequal');

            permutedP = reshape(permutedP, size(differenceTF));
            permutedT = reshape( ...
                permutedStats.tstat, size(differenceTF));

            for ss = 1:2

                if ss == 1
                    permutedMask = ...
                        permutedP < clusterStartingP & permutedT > 0;
                else
                    permutedMask = ...
                        permutedP < clusterStartingP & permutedT < 0;
                end

                permutedClusters = bwconncomp(permutedMask, 4);

                for cc = 1:permutedClusters.NumObjects

                    clusterMass = sum(abs(permutedT( ...
                        permutedClusters.PixelIdxList{cc})), ...
                        'omitnan');

                    largestPermutedCluster(pp) = max( ...
                        largestPermutedCluster(pp), clusterMass);

                end

            end

        end

        for ss = 1:2

            for cc = 1:observedClusters{ss}.NumObjects

                clusterPixels = ...
                    observedClusters{ss}.PixelIdxList{cc};

                observedClusterMass = ...
                    sum(abs(observedT(clusterPixels)), 'omitnan');

                clusterP = (1 + sum( ...
                    largestPermutedCluster >= observedClusterMass)) / ...
                    (numberOfPermutations + 1);

                if clusterP < clusterSignificanceP
                    significantDifference(clusterPixels) = true;
                end

            end

        end

    end

    fprintf('%s: %d significant time-frequency pixels\n', ...
        outcome, nnz(significantDifference));

    positiveAreaText = strjoin(positiveAreas{oo}, ', ');
    negativeAreaText = strjoin(negativeAreas{oo}, ', ');

    if strlength(positiveAreaText) == 0
        positiveAreaText = "none";
    end

    if strlength(negativeAreaText) == 0
        negativeAreaText = "none";
    end

    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [50 50 2100 720]);

    layout = tiledlayout(fig, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax1 = nexttile(layout);
    surf(ax1, spectrogramTime, spectrogramFrequencies, ...
        meanPositiveTF, 'EdgeColor', 'none');
    view(ax1, 2);
    set(ax1, 'YScale', 'log', ...
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
    caxis([-groupColorLimit groupColorLimit]);
    colorbar(ax1);
    box(ax1, 'off');

    ax2 = nexttile(layout);
    surf(ax2, spectrogramTime, spectrogramFrequencies, ...
        meanNegativeTF, 'EdgeColor', 'none');
    view(ax2, 2);
    set(ax2, 'YScale', 'log', ...
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
    caxis([-groupColorLimit groupColorLimit]);
    colorbar(ax2);
    box(ax2, 'off');

    ax3 = nexttile(layout);
    differenceSurface = surf(ax3, ...
        spectrogramTime, spectrogramFrequencies, ...
        differenceTF, 'EdgeColor', 'none');

    set(differenceSurface, ...
        'AlphaData', nonSignificantAlpha + ...
        (1 - nonSignificantAlpha) .* ...
        double(significantDifference), ...
        'AlphaDataMapping', 'none', ...
        'FaceAlpha', 'flat');

    view(ax3, 2);
    set(ax3, 'YScale', 'log', ...
        'YLim', [5 200], ...
        'YTick', [10 100], ...
        'FontSize', 12, ...
        'TickDir', 'out');
    xline(ax3, 0, 'r', 'LineWidth', 1);
    xlabel(ax3, 'time from ied peak (s)');
    ylabel(ax3, 'frequency (hz)');
    title(ax3, { ...
        'positive minus negative', ...
        'opaque areas: cbpt p < 0.05'});
    colormap(ax3, parula);
    axes(ax3);
    caxis([-differenceColorLimit differenceColorLimit]);
    colorbar(ax3);
    box(ax3, 'off');

    title(layout, sprintf([ ...
        'IED7 %s ied-locked lfp | baseline %.1f to %.1f s | ' ...
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
