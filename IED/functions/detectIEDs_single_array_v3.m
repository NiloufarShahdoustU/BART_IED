%% this is with hilbert transfrom

function [IEDdata] = detectIEDs_single_array_v3(data, Fo)
% DETECTIEDS finds interictal epileptiform discharges (IEDs) in a single EEG channel
% data: A vector representing the EEG signal
% Fo: Original sampling rate of the data

Fs = Fo; % Sampling rate for processing, kept the same as original in this case

% Define the desired frequency range
freqRange = [5 25]; % Frequency range for bandpass filter
smoothing_var = 10;

% Design the bandpass filter
[b, a] = butter(4, freqRange / (Fs / 2), 'bandpass');

% Apply bandpass filtering
filteredData = filtfilt(b, a, data);

% Apply Hilbert transform
hilbertTransform = abs(hilbert(filteredData));

% Downsampling to Fs
dataFs = resample(hilbertTransform, Fs, Fo);
% disp(size(dataFs));
% Finding peaks
DeviateFromMax = 0.9; % 90 percents deviate from max
maxValue = max(smooth(abs(dataFs), Fs / smoothing_var));
minPeakHeight = DeviateFromMax * maxValue;
[amp, locs, ~, ~] = findpeaks(smooth(abs(dataFs), Fs / smoothing_var), 'MinPeakHeight', minPeakHeight);

% Initialize variables for filtered peaks
filteredLocs = [];
lastAcceptedLoc = 0;

minSamplesApart = 0.8 * Fs; % 400 ms in terms of samples (it doubles)

for i = 1:length(locs)
    if (locs(i) - lastAcceptedLoc) >= minSamplesApart
        filteredLocs(end+1) = locs(i); % Accept the current location
        lastAcceptedLoc = locs(i); % Update the last accepted location
    end
end

% Update IEDdata structure with filtered locations
IEDdata.foundPeaks(1).locs = filteredLocs;
IEDdata.detections(1).times = filteredLocs / Fs; % Convert locations to time in seconds

end
