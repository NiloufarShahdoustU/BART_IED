% positive and negative ied8 interaction frequency band curves

clear;
clc;
close all;

% folders

lfpiedFolder = ...
    'D:\Nill\data\BART\0_0_new_IED\IED1_find_number_of_IEDs\';

rawDataFolder = ...
    'D:\Nill\data\BART_preprocessed\';

modelingFolder = ...
    ['D:\Nill\data\BART\0_0_new_IED\' ...
     'context_modeling\param_recovery_1_modeling\'];

ied8ResultsFile = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED\' ...
     'IED8_cox_expected_reward_by_brain_area\' ...
     'IT_RT_BR_mechanistic_IED_x_expected_reward_all_results.mat'];

outputFolder = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED\' ...
     'IED8_positive_negative_interaction_frequency_band_curves\'];

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

outputPDF = fullfile(outputFolder, ...
    'IED8_positive_negative_interaction_frequency_band_curves.pdf');

if exist(outputPDF, 'file')
    delete(outputPDF);
end

% settings

timeVector = -1:0.01:1;
baselineStart = -0.8;
baselineEnd = -0.3;

bandNames = ["Delta", "Theta", "Alpha", "Beta", "Gamma"];
bandRanges = [1 4; 4 8; 8 13; 13 30; 30 100];

outcomeNames = ["IT", "RT", "BR", "PR"];

minimumDistanceBetweenIEDsSeconds = 0.5;
minimumIEDsPerChannel = 3;
amplitudeThreshold = 5000;
maximumRTSeconds = 20;

numberOfPermutations = 1000;
sampleAlpha = 0.05;
clusterAlpha = 0.05;

positiveColor = [0.10 0.35 0.85];
negativeColor = [0.85 0.15 0.15];

rng(1);

if exist('butter', 'file') ~= 2 || ...
        exist('filtfilt', 'file') ~= 2 || ...
        exist('hilbert', 'file') ~= 2 || ...
        exist('zp2sos', 'file') ~= 2
    error('butter, filtfilt, hilbert, and zp2sos are needed');
end

% get positive and negative interaction areas

ied8 = load(ied8ResultsFile);

positiveAreas = cell(1, 4);
negativeAreas = cell(1, 4);

for oo = 1:4

    outcome = outcomeNames(oo);

    if outcome == "PR"
        resultTable = ied8.BRResults;
        beta = resultTable.betaPop_IEDxExpectedReward;
        significant = ...
            resultTable.significantPopUncorrected_IEDxExpectedReward;
    else
        resultTable = ied8.(char(outcome + "Results"));
        beta = resultTable.beta_IEDxExpectedReward;
        significant = ...
            resultTable.significantUncorrected_IEDxExpectedReward;
    end

    useRows = ...
        resultTable.status == "fitted" & ...
        significant & ...
        isfinite(beta);

    allAreas = clean_area_names(resultTable.anatomicalArea);

    positiveAreas{oo} = unique( ...
        allAreas(useRows & beta > 0), 'stable');

    negativeAreas{oo} = unique( ...
        allAreas(useRows & beta < 0), 'stable');

    fprintf('\n%s positive areas:\n', outcome);
    disp(positiveAreas{oo});

    fprintf('%s negative areas:\n', outcome);
    disp(negativeAreas{oo});

end

% each array is band x time x participant

positiveAll = cell(1, 4);
negativeAll = cell(1, 4);

for oo = 1:4
    positiveAll{oo} = zeros( ...
        length(bandNames), length(timeVector), 0);

    negativeAll{oo} = zeros( ...
        length(bandNames), length(timeVector), 0);
end

% load participants

fileList = dir(fullfile(lfpiedFolder, '*.LFPIED.mat'));
modelFileList = dir(fullfile(modelingFolder, '*.mat'));

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

    expectedRewardAll = load_expected_reward( ...
        modelFileList, modelingFolder, patientID);

    if isempty(expectedRewardAll)
        fprintf('skipped: expected reward data were not found\n');
        continue;
    end

    neededFields = { ...
        'selectedChans', 'anatomicalLocs', ...
        'RTs', 'ITs', 'isControl', 'balloonType', ...
        'BankedTrials', 'balloonTimes', ...
        'IED_occurance_RT', 'IED_occurance_IT'};

    missingField = false;

    for ff = 1:length(neededFields)
        if ~isfield(LFPIED, neededFields{ff})
            missingField = true;
        end
    end

    if missingField
        fprintf('skipped: a needed LFPIED field was missing\n');
        continue;
    end

    patientRawFolder = fullfile( ...
        rawDataFolder, char(patientID), 'Data');

    ns2List = dir(fullfile(patientRawFolder, '*.ns2'));
    nevList = dir(fullfile(patientRawFolder, '*.nev'));

    if length(ns2List) ~= 1 || length(nevList) ~= 1
        fprintf('skipped: need one ns2 and one nev file\n');
        continue;
    end

    % get response times for it, br, and pr

    NEV = openNEV( ...
        fullfile(nevList.folder, nevList.name), 'overwrite');

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

    if Fs <= 2.2 * max(bandRanges(:, 2))
        fprintf('skipped: sampling rate is too low\n');
        clear NSX rawLFP LFPIED loadedData;
        continue;
    end

    selectedChans = round(double(LFPIED.selectedChans(:)));

    selectedAreaLabels = get_selected_area_labels( ...
        LFPIED.anatomicalLocs, selectedChans);

    % run it, rt, br, and pr

    for oo = 1:4

        outcome = outcomeNames(oo);
        expectedReward = expectedRewardAll;

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
            length(RTs), length(ITs), ...
            length(isControl), length(balloonType), ...
            length(bankedTrials), length(startTimes), ...
            length(expectedReward)]);

        if nTrials < 1
            continue;
        end

        RTs = RTs(1:nTrials);
        duration = duration(1:nTrials);
        isControl = isControl(1:nTrials);
        balloonType = balloonType(1:nTrials);
        bankedTrials = bankedTrials(1:nTrials);
        startTimes = startTimes(1:nTrials);
        expectedReward = expectedReward(1:nTrials);

        colorCode = NaN(size(balloonType));
        goodColors = ismember( ...
            round(balloonType), [1 2 3 11 12 13]);

        colorCode(goodColors) = ...
            mod(round(balloonType(goodColors)) - 1, 10) + 1;

        validTrials = ...
            isControl == 0 & ...
            isfinite(RTs) & RTs > 0 & RTs <= maximumRTSeconds & ...
            isfinite(duration) & duration > 0 & ...
            isfinite(startTimes) & ...
            isfinite(expectedReward) & ...
            ismember(colorCode, [1 2 3]);

        % it and br use bank trials

        if outcome == "IT" || outcome == "BR"
            validTrials = validTrials & bankedTrials == 1;
        end

        % pr uses pop/loss trials

        if outcome == "PR"
            validTrials = validTrials & bankedTrials == 0;
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

            % band x time x area

            participantAreaCurves = zeros( ...
                length(bandNames), length(timeVector), 0);

            for aa = 1:length(groupAreas)

                areaName = groupAreas(aa);
                areaChannels = find( ...
                    strcmpi(selectedAreaLabels, areaName));

                if isempty(areaChannels)
                    continue;
                end

                % band x time x channel

                areaChannelCurves = zeros( ...
                    length(bandNames), length(timeVector), 0);

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

                    % remove close ieds in the same trial

                    separatedRows = true(length(absoluteSample), 1);
                    minimumDistanceSamples = ...
                        minimumDistanceBetweenIEDsSeconds * Fs;

                    for ii = 1:length(absoluteSample)

                        sameTrial = trialNumber == trialNumber(ii);
                        distance = abs( ...
                            absoluteSample - absoluteSample(ii));

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

                    % remove events outside the data or with big artifacts

                    goodEvents = false(length(absoluteSample), 1);
                    firstOffset = round(timeVector(1) * Fs);
                    lastOffset = round(timeVector(end) * Fs);

                    for ee = 1:length(absoluteSample)

                        firstSample = ...
                            absoluteSample(ee) + firstOffset;

                        lastSample = ...
                            absoluteSample(ee) + lastOffset;

                        if firstSample < 1 || ...
                                lastSample > length(signal)
                            continue;
                        end

                        oneEpoch = signal(firstSample:lastSample);

                        if any(~isfinite(oneEpoch)) || ...
                                max(abs(oneEpoch)) > amplitudeThreshold
                            continue;
                        end

                        goodEvents(ee) = true;

                    end

                    absoluteSample = absoluteSample(goodEvents);

                    if length(absoluteSample) < minimumIEDsPerChannel
                        continue;
                    end

                    channelCurves = get_band_curves( ...
                        signal, absoluteSample, Fs, ...
                        timeVector, baselineStart, baselineEnd, ...
                        bandRanges, minimumIEDsPerChannel);

                    if ~any(isfinite(channelCurves(:)))
                        continue;
                    end

                    areaChannelCurves(:, :, end + 1) = ...
                        channelCurves;

                end

                if isempty(areaChannelCurves)
                    continue;
                end

                participantAreaCurves(:, :, end + 1) = ...
                    mean(areaChannelCurves, 3, 'omitnan');

            end

            if isempty(participantAreaCurves)
                continue;
            end

            % each area has equal weight inside this participant

            participantCurve = ...
                mean(participantAreaCurves, 3, 'omitnan');

            if gg == 1
                positiveAll{oo}(:, :, end + 1) = ...
                    participantCurve;
            else
                negativeAll{oo}(:, :, end + 1) = ...
                    participantCurve;
            end

        end

    end

    clear NSX rawLFP LFPIED loadedData expectedRewardAll;

end

% make the 5 by 4 figure

fig = figure('Visible', 'off', 'Color', 'w');
set(fig, 'Position', [20 20 2100 1500]);

layout = tiledlayout(fig, 5, 4, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for bb = 1:length(bandNames)

    for oo = 1:length(outcomeNames)

        ax = nexttile(layout, ...
            (bb - 1) * length(outcomeNames) + oo);

        positiveData = reshape( ...
            permute(positiveAll{oo}(bb, :, :), [3 2 1]), ...
            size(positiveAll{oo}, 3), length(timeVector));

        negativeData = reshape( ...
            permute(negativeAll{oo}(bb, :, :), [3 2 1]), ...
            size(negativeAll{oo}, 3), length(timeVector));

        meanPositive = mean(positiveData, 1, 'omitnan');
        meanNegative = mean(negativeData, 1, 'omitnan');

        semPositive = get_sem(positiveData);
        semNegative = get_sem(negativeData);

        hold(ax, 'on');

        % sem shading

        if any(isfinite(meanPositive))
            fill(ax, ...
                [timeVector fliplr(timeVector)], ...
                [meanPositive - semPositive, ...
                 fliplr(meanPositive + semPositive)], ...
                positiveColor, ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end

        if any(isfinite(meanNegative))
            fill(ax, ...
                [timeVector fliplr(timeVector)], ...
                [meanNegative - semNegative, ...
                 fliplr(meanNegative + semNegative)], ...
                negativeColor, ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end

        plot(ax, timeVector, meanPositive, ...
            'Color', positiveColor, ...
            'LineWidth', 1.5, ...
            'DisplayName', 'positive');

        plot(ax, timeVector, meanNegative, ...
            'Color', negativeColor, ...
            'LineWidth', 1.5, ...
            'DisplayName', 'negative');

        xline(ax, 0, ':k', ...
            'LineWidth', 0.8, ...
            'HandleVisibility', 'off');

        yline(ax, 0, ':', ...
            'Color', [0.6 0.6 0.6], ...
            'LineWidth', 0.7, ...
            'HandleVisibility', 'off');

        % cbpt for this band and outcome
        % there is no fdr or bonferroni across the 20 subplots

        significantTime = false(size(timeVector));

        if size(positiveData, 1) >= 2 && ...
                size(negativeData, 1) >= 2

            significantTime = run_cluster_permutation( ...
                positiveData, negativeData, ...
                numberOfPermutations, ...
                sampleAlpha, clusterAlpha);
        end

        xlim(ax, [timeVector(1), timeVector(end)]);
        set(ax, ...
            'FontSize', 9, ...
            'TickDir', 'out', ...
            'Box', 'off');

        drawnow;
        yLimits = ylim(ax);

        if any(significantTime)
            yBar = yLimits(1) + 0.06 * diff(yLimits);
            plot_significance_bar( ...
                ax, timeVector, significantTime, yBar);
        end

        if oo == 1
            ylabel(ax, sprintf('%s\n%d-%d Hz', ...
                char(bandNames(bb)), ...
                bandRanges(bb, 1), bandRanges(bb, 2)));
        end

        if bb == length(bandNames)
            xlabel(ax, 'time from IED (s)');
        else
            set(ax, 'XTickLabel', []);
        end

        if bb == 1

            titleLines = make_area_title( ...
                outcomeNames(oo), ...
                positiveAreas{oo}, negativeAreas{oo});

            title(ax, titleLines, ...
                'Interpreter', 'none', ...
                'FontSize', 8, ...
                'FontWeight', 'normal');
        end

        if bb == 1 && oo == 1
            legend(ax, 'Location', 'best', ...
                'Box', 'off', 'FontSize', 8);
        end

        fprintf('%s %s: positive n = %d, negative n = %d\n', ...
            outcomeNames(oo), bandNames(bb), ...
            size(positiveData, 1), size(negativeData, 1));

    end

end

exportgraphics(fig, outputPDF, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

close(fig);

fprintf('\nfigure saved:\n%s\n', outputPDF);

%% local functions

function expectedReward = load_expected_reward( ...
        modelFileList, modelingFolder, patientID)

    expectedReward = [];
    modelNames = string({modelFileList.name})';

    if isempty(modelNames)
        return;
    end

    patientIDlower = lower(char(patientID));
    escapedID = regexptranslate('escape', patientIDlower);
    tokenPattern = ...
        ['(^|[^a-z0-9])' escapedID '([^a-z0-9]|$)'];

    matchedRows = false(length(modelNames), 1);

    for ii = 1:length(modelNames)
        matchedRows(ii) = ~isempty(regexp( ...
            lower(char(modelNames(ii))), ...
            tokenPattern, 'once'));
    end

    matchedRows = find(matchedRows);

    if isempty(matchedRows)
        matchedRows = find(contains( ...
            lower(modelNames), lower(string(patientID))));
    end

    if length(matchedRows) ~= 1
        return;
    end

    modelData = load(fullfile( ...
        modelingFolder, char(modelNames(matchedRows))));

    if ~isfield(modelData, 'TDdataParamRecovery')
        return;
    end

    TD = modelData.TDdataParamRecovery;

    if ~isfield(TD, 'bestApIdx') || ...
            ~isfield(TD, 'bestAnIdx') || ...
            ~isfield(TD, 'expectedReward')
        return;
    end

    expectedReward = squeeze(TD.expectedReward( ...
        TD.bestApIdx, TD.bestAnIdx, :))';

    expectedReward = double(expectedReward(:));

end

function areaNames = clean_area_names(areaNames)

    areaNames = string(areaNames(:));
    areaNames = strip(areaNames);

    areaNames(ismissing(areaNames) | ...
        strlength(areaNames) == 0) = "Unknown";

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

function selectedAreaLabels = get_selected_area_labels( ...
        rawLabels, selectedChans)

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

    if ~isempty(selectedChans) && ...
            length(allAreaLabels) >= max(selectedChans)
        selectedAreaLabels = allAreaLabels(selectedChans);
    else
        selectedAreaLabels = allAreaLabels;
    end

    selectedAreaLabels = clean_area_names(selectedAreaLabels);

end

function channelCurves = get_band_curves( ...
        signal, eventSamples, Fs, ...
        timeVector, baselineStart, baselineEnd, ...
        bandRanges, minimumIEDsPerChannel)

    numberOfBands = size(bandRanges, 1);
    channelCurves = NaN(numberOfBands, length(timeVector));

    firstOffset = round(timeVector(1) * Fs);
    lastOffset = round(timeVector(end) * Fs);
    sampleOffsets = firstOffset:lastOffset;
    sampleTimes = sampleOffsets ./ Fs;

    baselineRows = ...
        sampleTimes >= baselineStart & ...
        sampleTimes <= baselineEnd;

    for bb = 1:numberOfBands

        lowFrequency = bandRanges(bb, 1);
        highFrequency = bandRanges(bb, 2);

        [filterZ, filterP, filterK] = butter( ...
            4, [lowFrequency highFrequency] ./ (Fs / 2), ...
            'bandpass');

        [filterSOS, filterG] = zp2sos( ...
            filterZ, filterP, filterK);

        filteredSignal = filtfilt( ...
            filterSOS, filterG, signal);
        powerSignal = abs(hilbert(filteredSignal)).^2;

        eventCurves = NaN( ...
            length(eventSamples), length(timeVector));

        for ee = 1:length(eventSamples)

            indices = eventSamples(ee) + sampleOffsets;
            onePower = powerSignal(indices);

            baselinePower = mean( ...
                onePower(baselineRows), 'omitnan');

            if ~isfinite(baselinePower) || baselinePower <= 0
                continue;
            end

            oneCurve = 10 .* log10(onePower ./ baselinePower);

            eventCurves(ee, :) = interp1( ...
                sampleTimes, oneCurve, timeVector, ...
                'linear', NaN);

        end

        goodRows = sum(isfinite(eventCurves), 2) == ...
            length(timeVector);

        if sum(goodRows) >= minimumIEDsPerChannel
            channelCurves(bb, :) = ...
                mean(eventCurves(goodRows, :), 1, 'omitnan');
        end

        clear filteredSignal powerSignal;

    end

end

function semValues = get_sem(data)

    if isempty(data)
        semValues = NaN(1, size(data, 2));
        return;
    end

    numberValues = sum(isfinite(data), 1);
    semValues = std(data, 0, 1, 'omitnan') ./ sqrt(numberValues);
    semValues(numberValues < 2) = 0;

end

function significantTime = run_cluster_permutation( ...
        positiveData, negativeData, ...
        numberOfPermutations, sampleAlpha, clusterAlpha)

    [realT, realP] = ...
        get_welch_t(positiveData, negativeData);

    positiveClusters = get_time_clusters( ...
        realP < sampleAlpha & realT > 0);

    negativeClusters = get_time_clusters( ...
        realP < sampleAlpha & realT < 0);

    realClusters = [positiveClusters; negativeClusters];
    realMass = zeros(length(realClusters), 1);

    for cc = 1:length(realClusters)
        realMass(cc) = sum( ...
            abs(realT(realClusters{cc})), 'omitnan');
    end

    allData = [positiveData; negativeData];
    numberPositive = size(positiveData, 1);
    numberAll = size(allData, 1);

    maximumNullMass = zeros(numberOfPermutations, 1);

    for pp = 1:numberOfPermutations

        shuffledRows = randperm(numberAll);

        permutedPositive = ...
            allData(shuffledRows(1:numberPositive), :);

        permutedNegative = ...
            allData(shuffledRows(numberPositive + 1:end), :);

        [permutedT, permutedP] = ...
            get_welch_t(permutedPositive, permutedNegative);

        permutedPositiveClusters = get_time_clusters( ...
            permutedP < sampleAlpha & permutedT > 0);

        permutedNegativeClusters = get_time_clusters( ...
            permutedP < sampleAlpha & permutedT < 0);

        permutedClusters = [ ...
            permutedPositiveClusters; ...
            permutedNegativeClusters];

        for cc = 1:length(permutedClusters)

            thisMass = sum( ...
                abs(permutedT(permutedClusters{cc})), ...
                'omitnan');

            maximumNullMass(pp) = max( ...
                maximumNullMass(pp), thisMass);
        end

    end

    significantTime = false(1, size(allData, 2));

    for cc = 1:length(realClusters)

        clusterP = ...
            (1 + sum(maximumNullMass >= realMass(cc))) ./ ...
            (numberOfPermutations + 1);

        if clusterP < clusterAlpha
            significantTime(realClusters{cc}) = true;
        end

    end

end

function [tValues, pValues] = ...
        get_welch_t(positiveData, negativeData)

    numberPositive = sum(isfinite(positiveData), 1);
    numberNegative = sum(isfinite(negativeData), 1);

    meanPositive = mean(positiveData, 1, 'omitnan');
    meanNegative = mean(negativeData, 1, 'omitnan');

    variancePositive = var( ...
        positiveData, 0, 1, 'omitnan');

    varianceNegative = var( ...
        negativeData, 0, 1, 'omitnan');

    positivePart = variancePositive ./ numberPositive;
    negativePart = varianceNegative ./ numberNegative;

    standardError = sqrt(positivePart + negativePart);

    tValues = ...
        (meanPositive - meanNegative) ./ standardError;

    degreesFreedom = ...
        (positivePart + negativePart).^2 ./ ...
        (positivePart.^2 ./ (numberPositive - 1) + ...
         negativePart.^2 ./ (numberNegative - 1));

    pValues = 2 .* tcdf(-abs(tValues), degreesFreedom);

    badValues = ...
        numberPositive < 2 | ...
        numberNegative < 2 | ...
        ~isfinite(standardError) | ...
        standardError <= 0 | ...
        ~isfinite(degreesFreedom) | ...
        degreesFreedom <= 0;

    tValues(badValues) = NaN;
    pValues(badValues) = NaN;

end

function clusters = get_time_clusters(mask)

    mask = logical(mask(:)');
    changes = diff([false mask false]);

    clusterStarts = find(changes == 1);
    clusterEnds = find(changes == -1) - 1;

    clusters = cell(length(clusterStarts), 1);

    for cc = 1:length(clusterStarts)
        clusters{cc} = clusterStarts(cc):clusterEnds(cc);
    end

end

function plot_significance_bar( ...
        ax, timeVector, significantTime, yBar)

    clusters = get_time_clusters(significantTime);
    timeStep = median(diff(timeVector));

    for cc = 1:length(clusters)

        timeRows = clusters{cc};

        if length(timeRows) == 1
            xValues = [ ...
                timeVector(timeRows) - timeStep / 2, ...
                timeVector(timeRows) + timeStep / 2];
        else
            xValues = timeVector( ...
                [timeRows(1), timeRows(end)]);
        end

        plot(ax, xValues, [yBar yBar], ...
            'k', 'LineWidth', 2.5, ...
            'HandleVisibility', 'off');

    end

end

function titleLines = make_area_title( ...
        outcome, positiveAreas, negativeAreas)

    titleLines = {char(outcome)};

    if isempty(positiveAreas)
        positiveAreas = "none";
    end

    if isempty(negativeAreas)
        negativeAreas = "none";
    end

    areasPerLine = 3;

    for firstArea = 1:areasPerLine:length(positiveAreas)

        lastArea = min( ...
            firstArea + areasPerLine - 1, ...
            length(positiveAreas));

        areaText = strjoin( ...
            positiveAreas(firstArea:lastArea), ', ');

        if firstArea == 1
            titleLines{end + 1} = ...
                ['blue: ' char(areaText)];
        else
            titleLines{end + 1} = char(areaText);
        end

    end

    for firstArea = 1:areasPerLine:length(negativeAreas)

        lastArea = min( ...
            firstArea + areasPerLine - 1, ...
            length(negativeAreas));

        areaText = strjoin( ...
            negativeAreas(firstArea:lastArea), ', ');

        if firstArea == 1
            titleLines{end + 1} = ...
                ['red: ' char(areaText)];
        else
            titleLines{end + 1} = char(areaText);
        end

    end

end
