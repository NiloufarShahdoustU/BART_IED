%% IED detection sanity-check figures
% By Niloufar Shahdoust (niloufar.shahdoust@utah.edu)
%
% Finds accepted IEDs during the final 1000 ms of RT and IT.
% Plots mean +/- SEM from -500 to +500 ms around each IED peak.
% Saves only PDF figures; no LFPIED.mat files are saved.

clear;
clc;
close all;


%% Paths

behaviorFolder = ...
    'D:\Nill\data\BART\bhvStruct_Nill_made';

anatomicalFolder = ...
    'D:\Nill\data\BART\0_0_new_IED_new_area_labels\IED1_find_number_of_IEDs';

figureOutputFolder = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_new_area_labels\IED_0_IED_detection_sanity_check';

if ~exist(figureOutputFolder, 'dir')
    mkdir(figureOutputFolder);
end

fileList = dir(fullfile(behaviorFolder, '*.bhvStruct.mat'));

[~, sortIndex] = sort({fileList.name});
fileList = fileList(sortIndex);


%% Process participants

for pt = 1:length(fileList)

    data = load(fullfile( ...
        behaviorFolder, ...
        fileList(pt).name));

    fileNameParts = strsplit( ...
        fileList(pt).name, ...
        '.');

    ptID = fileNameParts{1};

    disp(' ');
    disp(['Processing patient ID: ' ptID]);


    %% Load existing LFPIED anatomical locations

    anatomicalFile = fullfile( ...
        anatomicalFolder, ...
        [ptID '.LFPIED.mat']);

    if ~exist(anatomicalFile, 'file')

        warning( ...
            ['No LFPIED anatomical-location file was found for ' ...
             '%s. Skipping this patient.'], ...
            ptID);

        continue;

    end

    anatomicalData = load( ...
        anatomicalFile, ...
        'LFPIED');

    if ~isfield(anatomicalData, 'LFPIED') || ...
            ~isfield(anatomicalData.LFPIED, 'anatomicalLocs')

        warning( ...
            ['LFPIED.anatomicalLocs was not found for %s. ' ...
             'Skipping this patient.'], ...
            ptID);

        continue;

    end

    referenceLFPIED = anatomicalData.LFPIED;


    %% Find participant NEV file

    nevList = dir(fullfile( ...
        'D:\Nill\data\BART_preprocessed', ...
        ptID, ...
        'Data', ...
        '*.nev'));

    if length(nevList) > 1

        error( ...
            'Many NEV files are available for patient %s.', ...
            ptID);

    elseif isempty(nevList)

        error( ...
            'No NEV file was found for patient %s.', ...
            ptID);

    else

        nevFile = fullfile( ...
            nevList.folder, ...
            nevList.name);

    end

    [~, isECoG, ~, ~, ~] = ptTrodesBART_2(ptID);


    %% Load triggers

    NEV = openNEV( ...
        nevFile, ...
        'overwrite');

    trigs = ...
        NEV.Data.SerialDigitalIO.UnparsedData;

    trigTimes = ...
        NEV.Data.SerialDigitalIO.TimeStampSec;


    %% Load neural data

    [nevPath, nevName, ~] = fileparts(nevFile);

    NSX = openNSx( ...
        fullfile(nevPath, [nevName '.ns2']));

    selectedChans = find(isECoG);

    nChans = length(selectedChans);

    Fs = NSX.MetaTags.SamplingFreq;


    %% Notch-filter ECoG channels

    notchFilter = true;

    data2K = zeros( ...
        nChans, ...
        size(NSX.Data, 2));

    for ch = 1:nChans

        tempSignal = resample( ...
            double(NSX.Data(selectedChans(ch), :)), ...
            Fs, ...
            Fs);

        if notchFilter

            [b1, a1] = iirnotch( ...
                60 / (Fs / 2), ...
                (60 / (Fs / 2)) / 50);

            tempSignal = filtfilt( ...
                b1, ...
                a1, ...
                tempSignal);

            [b2, a2] = iirnotch( ...
                120 / (Fs / 2), ...
                (120 / (Fs / 2)) / 50);

            tempSignal = filtfilt( ...
                b2, ...
                a2, ...
                tempSignal);

        end

        data2K(ch, :) = tempSignal;

    end

    clear NSX NEV tempSignal a1 a2 b1 b2


    %% Task times

    balloonTimes = trigTimes( ...
        trigs == 1  | ...
        trigs == 2  | ...
        trigs == 3  | ...
        trigs == 4  | ...
        trigs == 11 | ...
        trigs == 12 | ...
        trigs == 13 | ...
        trigs == 14);

    inflateTimes = trigTimes(trigs == 23);

    if length(balloonTimes) > length(inflateTimes)
        balloonTimes(end) = [];
    end

    respTimes = trigTimes(trigs == 23);

    outcomeType = trigs( ...
        sort([ ...
            find(trigs == 25); ...
            find(trigs == 26)])) - 24;

    nTrials = length(outcomeType);


    %% Find all IEDs in neural signal

    IED_timepoints = zeros(size(data2K));

    for ch = 1:nChans

        mySig = data2K(ch, :);

        IEDStruct = detectIEDs_single_array_v8( ...
            mySig, ...
            Fs);

        IED_indices = ...
            IEDStruct.foundPeaks.locs;

        IED_timepoints(ch, IED_indices) = 1;

    end


    %% Keep IEDs from RT - 1000 ms to RT
    %  and from IT - 1000 ms to IT

    RTs = data.bhvStruct.allRTs;
    ITs = data.bhvStruct.allITs;

    detectionWindowSeconds = 1;

    detectionWindowSamples = floor( ...
        Fs * detectionWindowSeconds);

    IED_timepointsRT = cell(nTrials, 1);
    IED_timepointsIT = cell(nTrials, 1);

    nSamplesRT = zeros(nTrials, 1);
    nSamplesIT = zeros(nTrials, 1);

    for trial = 1:nTrials

        % If the interval is shorter than 1000 ms,
        % use its complete duration.

        nSamplesRT(trial) = min( ...
            floor(Fs * RTs(trial)), ...
            detectionWindowSamples);

        nSamplesIT(trial) = min( ...
            floor(Fs * ITs(trial)), ...
            detectionWindowSamples);

        IED_timepointsRT{trial} = ...
            nan(nChans, nSamplesRT(trial));

        IED_timepointsIT{trial} = ...
            nan(nChans, nSamplesIT(trial));

    end


    %% Remove noisy and large-amplitude segments

    threshold = 7;
    ampThreshold = 5000;

    for ch = 1:nChans

        for trial = 1:nTrials

            %% RT - 1000 ms to RT

            RT_end = ...
                floor(Fs * balloonTimes(trial)) + ...
                floor(Fs * RTs(trial)) - 1;

            RT_start = ...
                RT_end - nSamplesRT(trial) + 1;

            LFP_temp = data2K( ...
                ch, ...
                RT_start:RT_end);

            IED_temp = IED_timepoints( ...
                ch, ...
                RT_start:RT_end);

            isNoisy = classifyNoisySignal( ...
                LFP_temp, ...
                threshold);

            if ~isNoisy

                isOutofRange = detectLargeAmplitude( ...
                    LFP_temp, ...
                    ampThreshold);

                if ~isOutofRange

                    IED_timepointsRT{trial}(ch, :) = ...
                        IED_temp;

                end

            end


            %% IT - 1000 ms to IT

            IT_end = ...
                floor(Fs * respTimes(trial)) + ...
                floor(Fs * ITs(trial)) - 1;

            IT_start = ...
                IT_end - nSamplesIT(trial) + 1;

            LFP_temp = data2K( ...
                ch, ...
                IT_start:IT_end);

            IED_temp = IED_timepoints( ...
                ch, ...
                IT_start:IT_end);

            isNoisy = classifyNoisySignal( ...
                LFP_temp, ...
                threshold);

            if ~isNoisy

                isOutofRange = detectLargeAmplitude( ...
                    LFP_temp, ...
                    ampThreshold);

                if ~isOutofRange

                    IED_timepointsIT{trial}(ch, :) = ...
                        IED_temp;

                end

            end

        end

    end


    %% Plot accepted IEDs from -500 to +500 ms

    plotWindowSeconds = 0.5;

    plotWindowSamples = floor( ...
        Fs * plotWindowSeconds);

    plotTime = ...
        (-plotWindowSamples:plotWindowSamples) ...
        / Fs * 1000;


    for ch = 1:nChans

        channelNumber = selectedChans(ch);

        IED_absolute_indices = [];

        for trial = 1:nTrials

            %% Absolute indices of RT IEDs

            RT_end = ...
                floor(Fs * balloonTimes(trial)) + ...
                floor(Fs * RTs(trial)) - 1;

            RT_start = ...
                RT_end - nSamplesRT(trial) + 1;

            IED_idx_RT = find( ...
                IED_timepointsRT{trial}(ch, :) == 1);

            RT_absolute_indices = ...
                RT_start + IED_idx_RT - 1;

            IED_absolute_indices = [ ...
                IED_absolute_indices, ...
                RT_absolute_indices];


            %% Absolute indices of IT IEDs

            IT_end = ...
                floor(Fs * respTimes(trial)) + ...
                floor(Fs * ITs(trial)) - 1;

            IT_start = ...
                IT_end - nSamplesIT(trial) + 1;

            IED_idx_IT = find( ...
                IED_timepointsIT{trial}(ch, :) == 1);

            IT_absolute_indices = ...
                IT_start + IED_idx_IT - 1;

            IED_absolute_indices = [ ...
                IED_absolute_indices, ...
                IT_absolute_indices];

        end


        % Do not plot an IED twice if it appears
        % in both RT and IT windows.

        IED_absolute_indices = unique( ...
            IED_absolute_indices);

        IED_waveforms = zeros( ...
            0, ...
            2 * plotWindowSamples + 1);


        for ied = 1:length(IED_absolute_indices)

            waveformStart = ...
                IED_absolute_indices(ied) - ...
                plotWindowSamples;

            waveformEnd = ...
                IED_absolute_indices(ied) + ...
                plotWindowSamples;

            if waveformStart >= 1 && ...
                    waveformEnd <= size(data2K, 2)

                waveform = data2K( ...
                    ch, ...
                    waveformStart:waveformEnd);

                IED_waveforms(end + 1, :) = waveform;

            end

        end


        %% Make one figure for this participant and channel

        if ~isempty(IED_waveforms)

            anatomicalLocation = getAnatomicalLocation( ...
                referenceLFPIED, ...
                channelNumber, ...
                ch);


            %% Calculate mean and SEM

            meanWaveform = mean( ...
                IED_waveforms, ...
                1);

            semWaveform = std( ...
                IED_waveforms, ...
                0, ...
                1) / sqrt(size(IED_waveforms, 1));

            upperSEM = meanWaveform + semWaveform;
            lowerSEM = meanWaveform - semWaveform;


            %% Create figure

            fig = figure( ...
                'Visible', 'off', ...
                'Color', 'white', ...
                'Position', [100 100 1000 650]);

            hold on;


            %% Plot mean +/- SEM

            fill( ...
                [plotTime, fliplr(plotTime)], ...
                [upperSEM, fliplr(lowerSEM)], ...
                [0.75 0.85 1.00], ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.45);


            %% Plot mean waveform

            plot( ...
                plotTime, ...
                meanWaveform, ...
                'Color', [0.10 0.35 0.75], ...
                'LineWidth', 2);


            %% Mark IED peak

            plot( ...
                0, ...
                meanWaveform(plotWindowSamples + 1), ...
                'r.', ...
                'MarkerSize', 14);

            xline( ...
                0, ...
                '--r', ...
                'IED peak', ...
                'LineWidth', 1);


            %% Figure labels

            xlabel('Time relative to IED peak (ms)');
            ylabel('Mean LFP amplitude');

            title([ ...
                ptID ...
                ', chan ' num2str(channelNumber) ...
                ', ' anatomicalLocation ...
                ', number of IEDs = ' ...
                num2str(size(IED_waveforms, 1))], ...
                'Interpreter', 'none');

            xlim([-500 500]);

            box off;


            %% Save only PDF

            pdfName = [ ...
                ptID ...
                '_chan' ...
                num2str(channelNumber) ...
                '.pdf'];

            exportgraphics( ...
                fig, ...
                fullfile(figureOutputFolder, pdfName), ...
                'ContentType', 'vector');

            close(fig);

        end

    end


    clear data data2K IED_timepoints ...
          IED_timepointsRT IED_timepointsIT ...
          anatomicalData referenceLFPIED

end


%% Match displayed channel to anatomical location

function anatomicalLocation = getAnatomicalLocation( ...
    LFPIED, channelNumber, channelPosition)

    anatomicalLocation = 'Unknown';

    anatomicalLocs = LFPIED.anatomicalLocs;

    locationIndex = [];


    % Case 1:
    % anatomicalLocs contains one entry
    % for every original channel.

    if numel(anatomicalLocs) >= channelNumber && ...
            (~isfield(LFPIED, 'selectedChans') || ...
             numel(anatomicalLocs) ~= ...
             numel(LFPIED.selectedChans))

        locationIndex = channelNumber;


    % Case 2:
    % anatomicalLocs contains only selected channels.

    elseif isfield(LFPIED, 'selectedChans')

        locationIndex = find( ...
            LFPIED.selectedChans == channelNumber, ...
            1, ...
            'first');

    end


    % Fallback if selectedChans is unavailable.

    if isempty(locationIndex) && ...
            numel(anatomicalLocs) >= channelPosition

        locationIndex = channelPosition;

    end


    if isempty(locationIndex) || ...
            locationIndex > numel(anatomicalLocs)

        return;

    end


    %% Extract anatomical-location text

    if iscell(anatomicalLocs)

        locationValue = ...
            anatomicalLocs{locationIndex};

    else

        locationValue = ...
            anatomicalLocs(locationIndex);

    end


    if iscell(locationValue)
        locationValue = locationValue{1};
    end


    if isstring(locationValue) || ...
            iscategorical(locationValue)

        anatomicalLocation = char(locationValue);

    elseif ischar(locationValue)

        anatomicalLocation = locationValue;

    elseif isnumeric(locationValue)

        anatomicalLocation = num2str(locationValue);

    end


    anatomicalLocation = strtrim( ...
        anatomicalLocation);

    if isempty(anatomicalLocation)
        anatomicalLocation = 'Unknown';
    end

end