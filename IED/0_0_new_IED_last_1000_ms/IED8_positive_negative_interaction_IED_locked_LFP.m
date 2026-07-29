% compare ied-locked spectrograms for positive and negative
% ied x expected reward cox coefficients from ied8
%
% this code makes one pdf with 3 pages:
% page 1 = rt
% page 2 = it
% page 3 = br
%
% positive and negative groups only include significant interaction areas
% from ied8
% no csv or mat files are saved

clear;
clc;
close all;

rng(1);

% folders

lfpiedFolder = ...
    'D:\Nill\data\BART\0_0_new_IED_last_1000_ms\IED1_find_number_of_IEDs\';

rawDataFolder = ...
    'D:\Nill\data\BART_preprocessed\';

modelingFolder = ...
    ['D:\Nill\data\BART\0_0_new_IED_last_1000_ms\' ...
     'context_modeling\param_recovery_1_modeling\'];

ied8ResultsFile = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED_last_1000_ms\' ...
     'IED8_cox_expected_reward_by_brain_area\' ...
     'IT_RT_BR_mechanistic_IED_x_expected_reward_all_results.mat'];

outputFolder = ...
    ['D:\Nill\code\BART\IED\0_0_new_IED_last_1000_ms\' ...
     'IED8_positive_negative_interaction_IED_locked_LFP\'];

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

outputPDF = fullfile(outputFolder, ...
    'IED8_positive_negative_interaction_IED_locked_LFP.pdf');

if exist(outputPDF, 'file')
    delete(outputPDF);
end

% analysis settings

preIEDSeconds = 1;
postIEDSeconds = 1;
extraSeconds = 0.5;

baselineStartSeconds = -0.8;
baselineEndSeconds = -0.3;

% log-spaced frequencies make the y axis logarithmic
spectrogramFrequencies = logspace(log10(5), log10(200), 60);
spectrogramTime = -0.9:0.02:0.9;

minimumDistanceBetweenIEDsSeconds = 0.5;
amplitudeThreshold = 5000;
minimumIEDsPerChannel = 3;

nPermutations = 1000;
clusterFormingP = 0.5;
clusterP = 0.05;
nonSignificantAlpha = 0.18;



maximumRTSeconds = 20;

if exist('basewaveERP', 'file') ~= 2
    error('basewaveERP.m was not found');
end

if exist('bwconncomp', 'file') ~= 2
    error('bwconncomp was not found');
end

% get significant positive and negative interaction areas from ied8

ied8 = load(ied8ResultsFile);

outcomeNames = ["RT", "IT", "BR"];

for oo = 1:length(outcomeNames)

    outcome = outcomeNames(oo);
    tableName = char(outcome + "Results");
    resultTable = ied8.(tableName);

    useRows = ...
        resultTable.status == "fitted" & ...
        resultTable.significantUncorrected_IEDxExpectedReward & ...
        isfinite(resultTable.beta_IEDxExpectedReward);

    positiveAreas{oo} = string(resultTable.anatomicalArea( ...
        useRows & resultTable.beta_IEDxExpectedReward > 0));

    negativeAreas{oo} = string(resultTable.anatomicalArea( ...
        useRows & resultTable.beta_IEDxExpectedReward < 0));

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

    groupData(oo).positiveTF = zeros( ...
        length(spectrogramFrequencies), length(spectrogramTime), 0);
    groupData(oo).negativeTF = zeros( ...
        length(spectrogramFrequencies), length(spectrogramTime), 0);
    groupData(oo).positiveIEDCount = [];
    groupData(oo).negativeIEDCount = [];

end

% load every participant

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

    expectedReward = loadExpectedReward( ...
        modelFileList, modelingFolder, patientID);

    if isempty(expectedReward)
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
            fprintf('skipped: %s was missing\n', neededFields{ff});
            missingField = true;
        end
    end

    if missingField
        continue;
    end

    patientRawFolder = fullfile(rawDataFolder, char(patientID), 'Data');
    ns2List = dir(fullfile(patientRawFolder, '*.ns2'));
    nevList = dir(fullfile(patientRawFolder, '*.nev'));

    if length(ns2List) ~= 1 || length(nevList) ~= 1
        fprintf('skipped: need exactly one ns2 and one nev file\n');
        continue;
    end

    % get the exact response trigger times for the it and br events

    NEV = openNEV(fullfile(nevList.folder, nevList.name), 'overwrite');
    trigs = NEV.Data.SerialDigitalIO.UnparsedData;
    trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
    responseTimes = trigTimes(trigs == 23);
    clear NEV;

    % load the raw lfp

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
    selectedAreaLabels = getSelectedChannelLabels( ...
        LFPIED.anatomicalLocs, selectedChans);
    selectedAreaLabels = cleanAreaLabels(selectedAreaLabels);

    % run rt, it, and br separately

    for oo = 1:length(outcomeNames)

        outcome = outcomeNames(oo);

        if isempty(positiveAreas{oo}) && isempty(negativeAreas{oo})
            continue;
        end

        [IEDrows, startTimes, validTrials] = getOutcomeIEDdata( ...
            LFPIED, responseTimes, expectedReward, ...
            outcome, maximumRTSeconds);

        if isempty(IEDrows) || ~any(validTrials)
            continue;
        end

        settings.Fs = Fs;
        settings.preIEDSeconds = preIEDSeconds;
        settings.postIEDSeconds = postIEDSeconds;
        settings.extraSeconds = extraSeconds;
        settings.baselineStartSeconds = baselineStartSeconds;
        settings.baselineEndSeconds = baselineEndSeconds;
        settings.spectrogramFrequencies = spectrogramFrequencies;
        settings.spectrogramTime = spectrogramTime;
        settings.minimumDistanceBetweenIEDsSeconds = ...
            minimumDistanceBetweenIEDsSeconds;
        settings.amplitudeThreshold = amplitudeThreshold;
        settings.minimumIEDsPerChannel = minimumIEDsPerChannel;

        [positiveTF, positiveCounts] = ...
            getParticipantGroupAverage( ...
                rawLFP, selectedChans, selectedAreaLabels, ...
                IEDrows, startTimes, validTrials, ...
                positiveAreas{oo}, settings);

        [negativeTF, negativeCounts] = ...
            getParticipantGroupAverage( ...
                rawLFP, selectedChans, selectedAreaLabels, ...
                IEDrows, startTimes, validTrials, ...
                negativeAreas{oo}, settings);

        if ~isempty(positiveTF)
            groupData(oo).positiveTF(:, :, end + 1) = positiveTF;
            groupData(oo).positiveIEDCount(end + 1, 1) = ...
                positiveCounts.nIEDs;
        end

        if ~isempty(negativeTF)
            groupData(oo).negativeTF(:, :, end + 1) = negativeTF;
            groupData(oo).negativeIEDCount(end + 1, 1) = ...
                negativeCounts.nIEDs;
        end

    end

    clear NSX rawLFP LFPIED loadedData;

end

% make one pdf page for each outcome

for oo = 1:length(outcomeNames)

    outcome = outcomeNames(oo);

    positiveTF = groupData(oo).positiveTF;
    negativeTF = groupData(oo).negativeTF;

    meanPositiveTF = meanOrNaN( ...
        positiveTF, length(spectrogramFrequencies), ...
        length(spectrogramTime));

    meanNegativeTF = meanOrNaN( ...
        negativeTF, length(spectrogramFrequencies), ...
        length(spectrogramTime));

    differenceTF = meanPositiveTF - meanNegativeTF;

    positiveLimit = getColorLimit([meanPositiveTF(:); meanNegativeTF(:)]);
    differenceLimit = getColorLimit(differenceTF(:));

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
            observedP < clusterFormingP & observedT > 0, 4);

        observedClusters{2} = bwconncomp( ...
            observedP < clusterFormingP & observedT < 0, 4);

        allRows = [positiveRows; negativeRows];
        totalParticipants = size(allRows, 1);
        largestPermutedCluster = zeros(nPermutations, 1);

        rng(1);

        for pp = 1:nPermutations

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
                        permutedP < clusterFormingP & permutedT > 0;
                else
                    permutedMask = ...
                        permutedP < clusterFormingP & permutedT < 0;
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

                thisClusterP = (1 + sum( ...
                    largestPermutedCluster >= observedClusterMass)) / ...
                    (nPermutations + 1);

                if thisClusterP < clusterP
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
    plotTimeFrequency(ax1, spectrogramTime, ...
        spectrogramFrequencies, meanPositiveTF, positiveLimit);
    title(ax1, { ...
        sprintf('positive coefficient, n = %d', size(positiveTF, 3)), ...
        char(positiveAreaText)}, ...
        'Interpreter', 'none');

    ax2 = nexttile(layout);
    plotTimeFrequency(ax2, spectrogramTime, ...
        spectrogramFrequencies, meanNegativeTF, positiveLimit);
    title(ax2, { ...
        sprintf('negative coefficient, n = %d', size(negativeTF, 3)), ...
        char(negativeAreaText)}, ...
        'Interpreter', 'none');

    ax3 = nexttile(layout);
    differenceSurface = plotTimeFrequency(ax3, spectrogramTime, ...
        spectrogramFrequencies, differenceTF, differenceLimit);

    set(differenceSurface, ...
        'AlphaData', nonSignificantAlpha + ...
        (1 - nonSignificantAlpha) .* ...
        double(significantDifference), ...
        'AlphaDataMapping', 'none', ...
        'FaceAlpha', 'flat');

    title(ax3, { ...
        'positive minus negative', ...
        'opaque areas: cbpt p < 0.05'});

    title(layout, sprintf([ ...
        'IED8 %s ied-locked lfp | baseline %.1f to %.1f s | ' ...
        'positive %d ieds | negative %d ieds'], ...
        outcome, baselineStartSeconds, baselineEndSeconds, ...
        sum(groupData(oo).positiveIEDCount), ...
        sum(groupData(oo).negativeIEDCount)), ...
        'FontSize', 16, ...
        'FontWeight', 'bold');

    if oo == 1
        exportgraphics(fig, outputPDF, ...
            'ContentType', 'image', 'Resolution', 220);
    else
        exportgraphics(fig, outputPDF, ...
            'ContentType', 'image', 'Resolution', 220, ...
            'Append', true);
    end

    close(fig);

end

fprintf('\nfinished. pdf saved here:\n%s\n', outputPDF);


function expectedReward = loadExpectedReward( ...
    modelFileList, modelingFolder, patientID)

    expectedReward = [];
    modelNames = string({modelFileList.name})';

    if isempty(modelNames)
        return;
    end

    patientIDlower = lower(char(patientID));
    escapedID = regexptranslate('escape', patientIDlower);
    tokenPattern = ['(^|[^a-z0-9])' escapedID '([^a-z0-9]|$)'];

    matchedRows = false(length(modelNames), 1);

    for ii = 1:length(modelNames)
        matchedRows(ii) = ~isempty(regexp( ...
            lower(char(modelNames(ii))), tokenPattern, 'once'));
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


function [IEDrows, startTimes, validTrials] = ...
    getOutcomeIEDdata(LFPIED, responseTimes, expectedReward, ...
    outcome, maximumRTSeconds)

    RTs = double(LFPIED.RTs(:));
    ITs = double(LFPIED.ITs(:));
    isControl = double(LFPIED.isControl(:));
    balloonType = double(LFPIED.balloonType(:));
    BankedTrials = double(LFPIED.BankedTrials(:));
    balloonTimes = double(LFPIED.balloonTimes(:));
    responseTimes = double(responseTimes(:));
    expectedReward = double(expectedReward(:));

    if outcome == "RT"
        IEDrows = double(LFPIED.IED_occurance_RT);
        startTimes = balloonTimes;
        duration = RTs;
        lengths = [length(RTs), length(isControl), ...
            length(balloonType), length(startTimes), ...
            length(expectedReward)];
    else
        IEDrows = double(LFPIED.IED_occurance_IT);
        startTimes = responseTimes;
        duration = ITs;
        lengths = [length(RTs), length(ITs), length(isControl), ...
            length(balloonType), length(startTimes), ...
            length(expectedReward)];
    end

    if outcome == "IT" || outcome == "BR"
        lengths(end + 1) = length(BankedTrials);
    end

    if isfield(LFPIED, 'nTrials')
        lengths(end + 1) = double(LFPIED.nTrials);
    end

    nTrials = floor(min(lengths));

    RTs = RTs(1:nTrials);
    duration = duration(1:nTrials);
    isControl = isControl(1:nTrials);
    balloonType = balloonType(1:nTrials);
    startTimes = startTimes(1:nTrials);
    expectedReward = expectedReward(1:nTrials);

    colorCode = NaN(size(balloonType));
    goodColors = ismember(round(balloonType), [1 2 3 11 12 13]);
    colorCode(goodColors) = ...
        mod(round(balloonType(goodColors)) - 1, 10) + 1;

    validTrials = ...
        isControl == 0 & ...
        isfinite(RTs) & RTs > 0 & RTs <= maximumRTSeconds & ...
        isfinite(duration) & duration > 0 & ...
        isfinite(startTimes) & ...
        isfinite(expectedReward) & ...
        ismember(colorCode, [1 2 3]);

    if outcome == "IT" || outcome == "BR"
        BankedTrials = BankedTrials(1:nTrials);
    end

    % it uses banked trials only
    if outcome == "IT"
        validTrials = validTrials & ...
            isfinite(BankedTrials) & BankedTrials == 1;
    end

    % br uses both banked and lost trials
    if outcome == "BR"
        validTrials = validTrials & ...
            isfinite(BankedTrials) & ...
            ismember(BankedTrials, [0 1]);
    end

end


function [groupTF, counts] = getParticipantGroupAverage( ...
    rawLFP, selectedChans, selectedAreaLabels, ...
    IEDrows, startTimes, validTrials, groupAreas, settings)

    groupTF = [];

    counts.nIEDs = 0;
    counts.nChannels = 0;
    counts.nAreas = 0;

    if isempty(groupAreas)
        return;
    end

    areaTF = zeros( ...
        length(settings.spectrogramFrequencies), ...
        length(settings.spectrogramTime), 0);

    for aa = 1:length(groupAreas)

        areaName = groupAreas(aa);
        areaChannels = find(selectedAreaLabels == areaName);

        if isempty(areaChannels)
            continue;
        end

        channelTF = zeros( ...
            length(settings.spectrogramFrequencies), ...
            length(settings.spectrogramTime), 0);
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
            goodRows = ...
                trialNumber >= 1 & ...
                trialNumber <= length(validTrials);

            goodIndices = find(goodRows);

            if ~isempty(goodIndices)
                goodRows(goodIndices) = ...
                    validTrials(trialNumber(goodIndices));
            end

            thisIED = thisIED(goodRows, :);

            if isempty(thisIED)
                continue;
            end

            trialNumber = round(thisIED(:, 1));
            absoluteSample = ...
                floor(settings.Fs .* startTimes(trialNumber)) + ...
                round(thisIED(:, 3)) - 1;

            [absoluteSample, keepRows] = keepSeparatedIEDs( ...
                absoluteSample, trialNumber, ...
                settings.minimumDistanceBetweenIEDsSeconds * settings.Fs);

            thisIED = thisIED(keepRows, :);

            if length(absoluteSample) < settings.minimumIEDsPerChannel
                continue;
            end

            rawChannelNumber = selectedChans(localChannel);

            if rawChannelNumber < 1 || ...
                    rawChannelNumber > size(rawLFP, 1)
                continue;
            end

            signal = double(rawLFP(rawChannelNumber, :));

            eventTF = zeros( ...
                length(settings.spectrogramFrequencies), ...
                length(settings.spectrogramTime), 0);

            for ee = 1:length(absoluteSample)

                oneTF = getOneIEDfeatures( ...
                    signal, absoluteSample(ee), settings);

                if isempty(oneTF)
                    continue;
                end

                eventTF(:, :, end + 1) = oneTF;

            end

            if size(eventTF, 3) < settings.minimumIEDsPerChannel
                continue;
            end

            channelTF(:, :, end + 1) = mean(eventTF, 3, 'omitnan');

            areaIEDcount = areaIEDcount + size(eventTF, 3);
            counts.nChannels = counts.nChannels + 1;

        end

        if isempty(channelTF)
            continue;
        end

        areaTF(:, :, end + 1) = mean(channelTF, 3, 'omitnan');

        counts.nIEDs = counts.nIEDs + areaIEDcount;
        counts.nAreas = counts.nAreas + 1;

    end

    if isempty(areaTF)
        return;
    end

    groupTF = mean(areaTF, 3, 'omitnan');

end


function oneTF = getOneIEDfeatures( ...
    signal, peakSample, settings)

    oneTF = [];

    Fs = settings.Fs;

    longPreSamples = round( ...
        (settings.preIEDSeconds + settings.extraSeconds) * Fs);

    longPostSamples = round( ...
        (settings.postIEDSeconds + settings.extraSeconds) * Fs);

    firstSample = peakSample - longPreSamples;
    lastSample = peakSample + longPostSamples;

    if firstSample < 1 || lastSample > length(signal)
        return;
    end

    oneEpoch = signal(firstSample:lastSample);
    oneEpoch = detrend(oneEpoch, 0);

    if any(~isfinite(oneEpoch)) || ...
            max(abs(oneEpoch)) > settings.amplitudeThreshold
        return;
    end

    % basewaveerp makes the spectrogram

    [wave, period] = basewaveERP( ...
        oneEpoch(:)', Fs, ...
        min(settings.spectrogramFrequencies), ...
        max(settings.spectrogramFrequencies), 6, 0);

    [waveletFrequencies, frequencyOrder] = sort(1 ./ period);
    wave = wave(frequencyOrder, :);

    waveletTime = linspace( ...
        -(settings.preIEDSeconds + settings.extraSeconds), ...
        settings.postIEDSeconds + settings.extraSeconds, ...
        length(oneEpoch));

    power = abs(wave).^2;

    baselineRows = ...
        waveletTime >= settings.baselineStartSeconds & ...
        waveletTime <= settings.baselineEndSeconds;

    baselinePower = mean(power(:, baselineRows), 2, 'omitnan');

    if any(~isfinite(baselinePower)) || any(baselinePower <= 0)
        return;
    end

    powerDB = 10 .* log10(power ./ baselinePower);

    [oldTimeGrid, oldFrequencyGrid] = meshgrid( ...
        waveletTime, waveletFrequencies);

    [newTimeGrid, newFrequencyGrid] = meshgrid( ...
        settings.spectrogramTime, ...
        settings.spectrogramFrequencies);

    oneTF = interp2( ...
        oldTimeGrid, oldFrequencyGrid, powerDB, ...
        newTimeGrid, newFrequencyGrid, ...
        'linear', NaN);

end


function [keptSamples, keepRows] = keepSeparatedIEDs( ...
    absoluteSample, trialNumber, minimumDistanceSamples)

    absoluteSample = round(absoluteSample(:));
    trialNumber = round(trialNumber(:));

    keepRows = true(length(absoluteSample), 1);

    for ii = 1:length(absoluteSample)
        sameTrial = trialNumber == trialNumber(ii);
        distance = abs(absoluteSample - absoluteSample(ii));

        hasCloseIED = any( ...
            sameTrial & distance > 0 & ...
            distance < minimumDistanceSamples);

        if hasCloseIED
            keepRows(ii) = false;
        end
    end

    keptSamples = absoluteSample(keepRows);

end


function labels = getSelectedChannelLabels(anatomicalLocs, selectedChans)

    labels = convertLabelsToString(anatomicalLocs);
    labels = labels(:);

    if length(labels) >= max(selectedChans)
        labels = labels(selectedChans);
    elseif length(labels) == length(selectedChans)
        % labels are already only for the selected channels
    else
        error('cannot map anatomical locations to selected channels');
    end

end


function labels = convertLabelsToString(rawLabels)

    if isstring(rawLabels)
        labels = rawLabels;
    elseif iscell(rawLabels)
        labels = strings(numel(rawLabels), 1);

        for ii = 1:numel(rawLabels)
            value = rawLabels{ii};

            if isempty(value)
                labels(ii) = "";
            elseif iscell(value) && numel(value) == 1
                labels(ii) = string(value{1});
            else
                labels(ii) = string(value);
            end
        end

    elseif ischar(rawLabels)
        labels = string(cellstr(rawLabels));
    elseif iscategorical(rawLabels)
        labels = string(rawLabels);
    else
        labels = string(rawLabels);
    end

    labels = labels(:);

end


function labels = cleanAreaLabels(labels)

    labels = strip(string(labels));
    labels(ismissing(labels) | strlength(labels) == 0) = "Unknown";

    labels = regexprep(labels, ...
        '^\s*\(?\s*(Left|Right|LH|RH|L|R)\s*\)?[_\-\s]+', ...
        '', 'ignorecase');

    labels = regexprep(labels, ...
        '[_\-\s]+\(?\s*(Left|Right|LH|RH|L|R)\s*\)?\s*$', ...
        '', 'ignorecase');

    labels = regexprep(labels, ...
        '^\s*(Left|Right)\s+hemisphere[_\-\s]+', ...
        '', 'ignorecase');

    labels = regexprep(labels, ...
        '[_\-\s]+(Left|Right)\s+hemisphere\s*$', ...
        '', 'ignorecase');

    labels = regexprep(labels, '[_\-]+', ' ');
    labels = regexprep(labels, '\s+', ' ');
    labels = strip(labels);
    labels(strlength(labels) == 0) = "Unknown";

end


function meanMatrix = meanOrNaN(data, nFrequencies, nTimes)

    if isempty(data)
        meanMatrix = NaN(nFrequencies, nTimes);
    else
        meanMatrix = mean(data, 3, 'omitnan');
    end

end


function colorLimit = getColorLimit(values)

    values = abs(values(isfinite(values)));

    if isempty(values)
        colorLimit = 1;
    else
        colorLimit = prctile(values, 98);

        if ~isfinite(colorLimit) || colorLimit <= 0
            colorLimit = max(values);
        end

        if ~isfinite(colorLimit) || colorLimit <= 0
            colorLimit = 1;
        end
    end

end


function surfaceHandle = plotTimeFrequency( ...
    ax, time, frequencies, data, colorLimit)

    surfaceHandle = surf(ax, time, frequencies, ...
        data, 'EdgeColor', 'none');

    view(ax, 2);
    set(ax, ...
        'YScale', 'log', ...
        'YLim', [5 200], ...
        'YTick', [10 100], ...
        'FontSize', 12, ...
        'TickDir', 'out');

    xline(ax, 0, 'r', 'LineWidth', 1);
    xlabel(ax, 'time from ied peak (s)');
    ylabel(ax, 'frequency (hz)');
    colormap(ax, parula);
    axes(ax);
    caxis([-colorLimit colorLimit]);
    colorbar(ax);
    box(ax, 'off');

end
