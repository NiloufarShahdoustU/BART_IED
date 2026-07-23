% Find IEDs in the last 500 ms of RT and the last 500 ms of IT
% Visualize the original LFP from -1000 to 1000 ms around each IED peak

clear;
clc;
close all;


% Load behavioral files

inputFolderName = 'D:\Nill\data\BART\bhvStruct_Nill_made';

fileList = dir(fullfile(inputFolderName, '*.bhvStruct.mat'));

outputFolderName = ...
    'D:\Nill\data\BART\0_0_new_IED_last_500_ms\IED1_find_number_of_IEDs\';

figureOutputFolder = ...
    'D:\Nill\code\BART\IED\0_0_new_IED_last_500_ms\IED1_find_number_of_IEDs\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

if ~exist(figureOutputFolder, 'dir')
    mkdir(figureOutputFolder);
end


for pt = 1:length(fileList)

    data = load(fullfile(inputFolderName, fileList(pt).name));

    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1};

    disp(' ');
    disp(['Processing patient ID: ' patientID]);

    ptID = patientID;

    nevList = dir(['D:\Nill\data\BART_preprocessed\' ptID '\Data\*.nev']);

    if length(nevList) > 1
        error('Many NEV files are available for this patient.')
    elseif length(nevList) < 1
        error('No NEV file was found for this patient.')
    else
        nevFile = fullfile(nevList.folder, nevList.name);
    end

    [trodeLabels, isECoG, ~, ~, anatomicalLocs] = ptTrodesBART_2(ptID);


    % Load triggers

    NEV = openNEV(nevFile, 'overwrite');

    trigs = NEV.Data.SerialDigitalIO.UnparsedData;
    trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;


    % Load neural data

    [nevPath, nevName, ~] = fileparts(nevFile);

    NSX = openNSx(fullfile(nevPath, [nevName '.ns2']));

    selectedChans = find(isECoG);
    nChans = length(selectedChans);
    Fs = NSX.MetaTags.SamplingFreq;

    poppedTrials = data.bhvStruct.poppedTrials;
    BankedTrials = ~poppedTrials;


    % Notch filter

    notchFilter = true;

    data2K = zeros(nChans, size(NSX.Data, 2));

    for ch = 1:nChans

        if notchFilter

            [b1, a1] = iirnotch(60/(Fs/2), (60/(Fs/2))/50);

            tempSignal = resample( ...
                double(NSX.Data(selectedChans(ch), :)), Fs, Fs);

            tempSignal = filtfilt(b1, a1, tempSignal);

            [b2, a2] = iirnotch(120/(Fs/2), (120/(Fs/2))/50);

            data2K(ch, :) = filtfilt(b2, a2, tempSignal);

        else

            data2K(ch, :) = resample( ...
                double(NSX.Data(selectedChans(ch), :)), Fs, Fs);

        end

    end

    clear NSX NEV tempSignal a1 a2 b1 b2


    % Task times

    balloonTimes = trigTimes( ...
        trigs == 1 | trigs == 2 | trigs == 3 | trigs == 4 | ...
        trigs == 11 | trigs == 12 | trigs == 13 | trigs == 14);

    inflateTimes = trigTimes(trigs == 23);

    balloonType = trigs( ...
        trigs == 1 | trigs == 2 | trigs == 3 | trigs == 4 | ...
        trigs == 11 | trigs == 12 | trigs == 13 | trigs == 14);

    if length(balloonTimes) > length(inflateTimes)
        balloonTimes(end) = [];
    end

    respTimes = trigTimes(trigs == 23);

    outcomeTimes = trigTimes(trigs == 25 | trigs == 26);

    outcomeType = trigs( ...
        sort([find(trigs == 25); find(trigs == 26)])) - 24;

    nTrials = length(outcomeType);

    if length(balloonType) > nTrials
        balloonType = balloonType(1:end-1);
    end


    % Find all IEDs in the neural signal

    IED_timepoints = zeros(size(data2K));

    for ch = 1:nChans

        mySig = data2K(ch, :);

        IEDStruct = detectIEDs_single_array_v8(mySig, Fs);

        IED_indices = IEDStruct.foundPeaks.locs;

        IED_timepoints(ch, IED_indices) = 1;

    end


    % Keep only the last 500 ms of RT and IT

    RTs = data.bhvStruct.allRTs;
    ITs = data.bhvStruct.allITs;

    windowSeconds = 0.5;
    windowSamples = floor(Fs * windowSeconds);

    IED_timepointsRT = cell(nTrials, 1);
    IED_timepointsIT = cell(nTrials, 1);

    nSamplesRT = zeros(nTrials, 1);
    nSamplesIT = zeros(nTrials, 1);

    for trial = 1:nTrials

        % If the interval is shorter than 500 ms, use its whole duration

        nSamplesRT(trial) = min( ...
            floor(Fs * RTs(trial)), windowSamples);

        nSamplesIT(trial) = min( ...
            floor(Fs * ITs(trial)), windowSamples);

        IED_timepointsRT{trial} = ...
            nan(nChans, nSamplesRT(trial));

        IED_timepointsIT{trial} = ...
            nan(nChans, nSamplesIT(trial));

    end


    % Remove noisy channels and save IED time points

    threshold = 7;
    ampThreshold = 5000;

    for ch = 1:nChans

        for trial = 1:nTrials

            % Last 500 ms of RT

            RT_end = floor(Fs * balloonTimes(trial)) + ...
                     floor(Fs * RTs(trial)) - 1;

            RT_start = RT_end - nSamplesRT(trial) + 1;

            LFP_temp = data2K(ch, RT_start:RT_end);
            IED_temp = IED_timepoints(ch, RT_start:RT_end);

            isNoisy = classifyNoisySignal(LFP_temp, threshold);

            if ~isNoisy

                isOutofRange = ...
                    detectLargeAmplitude(LFP_temp, ampThreshold);

                if ~isOutofRange
                    IED_timepointsRT{trial}(ch, :) = IED_temp;
                end

            end


            % Last 500 ms of IT

            IT_end = floor(Fs * respTimes(trial)) + ...
                     floor(Fs * ITs(trial)) - 1;

            IT_start = IT_end - nSamplesIT(trial) + 1;

            LFP_temp = data2K(ch, IT_start:IT_end);
            IED_temp = IED_timepoints(ch, IT_start:IT_end);

            isNoisy = classifyNoisySignal(LFP_temp, threshold);

            if ~isNoisy

                isOutofRange = ...
                    detectLargeAmplitude(LFP_temp, ampThreshold);

                if ~isOutofRange
                    IED_timepointsIT{trial}(ch, :) = IED_temp;
                end

            end

        end

    end


    % Save every IED occurrence
    % Column 1 = trial
    % Column 2 = channel
    % Column 3 = time index inside the 500 ms window

    IED_occurance_RT = zeros(0, 3);
    IED_occurance_IT = zeros(0, 3);

    for ch = 1:nChans

        for trial = 1:nTrials

            IED_idx_RT = ...
                find(IED_timepointsRT{trial}(ch, :) == 1);

            newRowsRT = [ ...
                trial * ones(length(IED_idx_RT), 1), ...
                ch * ones(length(IED_idx_RT), 1), ...
                IED_idx_RT(:)];

            IED_occurance_RT = ...
                [IED_occurance_RT; newRowsRT];


            IED_idx_IT = ...
                find(IED_timepointsIT{trial}(ch, :) == 1);

            newRowsIT = [ ...
                trial * ones(length(IED_idx_IT), 1), ...
                ch * ones(length(IED_idx_IT), 1), ...
                IED_idx_IT(:)];

            IED_occurance_IT = ...
                [IED_occurance_IT; newRowsIT];

        end

    end

    IED_occurance_RT = sortrows(IED_occurance_RT, 1);
    IED_occurance_IT = sortrows(IED_occurance_IT, 1);


    % Visualize all accepted IEDs from each channel
    % The plotted signal is from -1000 to 1000 ms around each IED.
    % The plotting window is independent of the saved 500 ms RT and IT data.

    plotWindowSeconds = 1;
    plotWindowSamples = floor(Fs * plotWindowSeconds);

    plotTime = ...
        (-plotWindowSamples:plotWindowSamples) / Fs * 1000;

    for ch = 1:nChans

        IED_absolute_indices = [];

        for trial = 1:nTrials

            % Absolute indices of RT IEDs

            RT_end = floor(Fs * balloonTimes(trial)) + ...
                     floor(Fs * RTs(trial)) - 1;

            RT_start = RT_end - nSamplesRT(trial) + 1;

            IED_idx_RT = ...
                find(IED_timepointsRT{trial}(ch, :) == 1);

            RT_absolute_indices = ...
                RT_start + IED_idx_RT - 1;

            IED_absolute_indices = ...
                [IED_absolute_indices, RT_absolute_indices];


            % Absolute indices of IT IEDs

            IT_end = floor(Fs * respTimes(trial)) + ...
                     floor(Fs * ITs(trial)) - 1;

            IT_start = IT_end - nSamplesIT(trial) + 1;

            IED_idx_IT = ...
                find(IED_timepointsIT{trial}(ch, :) == 1);

            IT_absolute_indices = ...
                IT_start + IED_idx_IT - 1;

            IED_absolute_indices = ...
                [IED_absolute_indices, IT_absolute_indices];

        end


        % Do not plot the same IED twice if it appears in both windows

        IED_absolute_indices = unique(IED_absolute_indices);

        IED_waveforms = [];

        for ied = 1:length(IED_absolute_indices)

            waveformStart = ...
                IED_absolute_indices(ied) - plotWindowSamples;

            waveformEnd = ...
                IED_absolute_indices(ied) + plotWindowSamples;

            if waveformStart >= 1 && ...
                    waveformEnd <= size(data2K, 2)

                waveform = ...
                    data2K(ch, waveformStart:waveformEnd);

                IED_waveforms = ...
                    [IED_waveforms; waveform];

            end

        end


        % Make one plot for this patient and channel

        if ~isempty(IED_waveforms)

            fig = figure( ...
                'Visible', 'off', ...
                'Color', 'white', ...
                'Position', [100 100 1000 650]);

            hold on;

            plot(plotTime, IED_waveforms', ...
                'Color', [0.35 0.55 0.85], ...
                'LineWidth', 0.5);

            plot(zeros(size(IED_waveforms, 1), 1), ...
                IED_waveforms(:, plotWindowSamples + 1), ...
                'r.', ...
                'MarkerSize', 10);

            xline(0, '--r', 'IED peak', 'LineWidth', 1);

            xlabel('Time relative to IED peak (ms)');
            ylabel('LFP amplitude');

            title([ptID ...
                ', channel ' num2str(selectedChans(ch)) ...
                ', number of IEDs = ' ...
                num2str(size(IED_waveforms, 1))]);

            xlim([-1000 1000]);
            box off;


            % Save PDF

            pdfName = [ ...
                ptID ...
                '_chan' ...
                num2str(selectedChans(ch)) ...
                '.pdf'];

            exportgraphics( ...
                fig, ...
                fullfile(figureOutputFolder, pdfName), ...
                'ContentType', 'vector');

            close(fig);

        end

    end


    % Save the same LFPIED data

    LFPIED.IED_occurance_RT = IED_occurance_RT;
    LFPIED.IED_occurance_IT = IED_occurance_IT;

    LFPIED.IED_occurance_RT_columns = ...
        {'trial', 'channel', 'time_index_within_RT'};

    LFPIED.IED_occurance_IT_columns = ...
        {'trial', 'channel', 'time_index_within_IT'};

    LFPIED.outcomeType = outcomeType;
    LFPIED.anatomicalLocs = anatomicalLocs;
    LFPIED.balloonType = balloonType;
    LFPIED.balloonTimes = balloonTimes;
    LFPIED.nTrials = nTrials;
    LFPIED.selectedChans = selectedChans;
    LFPIED.trodeLabels = trodeLabels;

    LFPIED.RTs = RTs;
    LFPIED.ITs = ITs;
    LFPIED.isControl = data.bhvStruct.isCtrl;
    LFPIED.poppedTrials = poppedTrials;
    LFPIED.BankedTrials = BankedTrials;

    save( ...
        [outputFolderName ptID '.LFPIED.mat'], ...
        'LFPIED');


    clear data2K IED_timepoints ...
          IED_timepointsRT IED_timepointsIT LFPIED

end