%% this is Nill's crappy version that happens to work well! :D


function [IEDdata] = detectIEDs_single_array_v2(data, Fo)
% DETECTIEDS finds interictal epileptiform discharges (IEDs) in a single EEG channel
% data: A vector representing the EEG signal
% Fo: Original sampling rate of the data


Fs = Fo; % Sampling rate for processing, kept the same as original in this case

% Filtering and downsampling
[b, a] = butter(4, [5 25] / (Fs / 2), 'bandpass');

% Filtering and thresholding
tmpSig = resample(data, Fs, Fo) - mean(resample(data, Fs, Fo));
data2070 = filtfilt(b, a, tmpSig);

% Finding peaks
DeviateFromMax = 0.5;
maxValue = max(smooth(abs(data2070), Fs / 50));
minPeakHeight = DeviateFromMax * maxValue;
[~, locs, ~, ~] = ...
    findpeaks(smooth(abs(data2070), Fs / 50), 'MinPeakHeight', minPeakHeight);

% Initialize variables for filtered peaks
filteredLocs = [];
lastAcceptedLoc = 0;
minSamplesApart = 0.15 * Fs; % 150 ms in terms of samples

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
