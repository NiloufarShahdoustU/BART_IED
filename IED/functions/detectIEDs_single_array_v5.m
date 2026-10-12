function [IEDdata] = detectIEDs_single_array_v5(data, Fo)
% DETECTIEDS_SINGLE_ARRAY_V5
% Detect IEDs using Hilbert envelope for detection
% BUT save ONLY true LFP waveform peaks
%
% data : raw LFP (1D)
% Fo   : sampling rate (Hz)

Fs = Fo;

%% ---------------- parameters ----------------
freqRange    = [5 25];
smoothingVar = 10;
minPeakFrac  = 0.6;          % 60% of max envelope
minISI       = 0.250;        % 250 ms
searchWinSec = 0.050;        % ±50 ms for true peak

%% ---------------- bandpass ----------------
[b, a] = butter(4, freqRange/(Fs/2), 'bandpass');
filteredData = filtfilt(b, a, data);

%% ---------------- envelope detection (INTERNAL ONLY) ----------------
env = abs(hilbert(filteredData));
envSmooth = smooth(env, Fs/smoothingVar);

minPeakProminence = minPeakFrac * max(envSmooth);
[~, locs] = findpeaks(envSmooth, ...
    'MinPeakProminence', minPeakProminence);

%% ---------------- enforce refractory ----------------
minSamplesApart = round(minISI * Fs);

keep = true(size(locs));
lastLoc = -inf;

for i = 1:length(locs)
    if locs(i) - lastLoc < minSamplesApart
        keep(i) = false;
    else
        lastLoc = locs(i);
    end
end

locs = locs(keep);

%% ---------------- TRUE LFP PEAK REFINEMENT ----------------
searchWin = round(searchWinSec * Fs);
truePeaks = zeros(size(locs));

for i = 1:length(locs)

    c = locs(i);
    i1 = max(1, c - searchWin);
    i2 = min(length(filteredData), c + searchWin);

    seg = filteredData(i1:i2);

    % TRUE waveform peak (absolute)
    [~, idx] = max(abs(seg));
    truePeaks(i) = i1 + idx - 1;

end

%% ---------------- OUTPUT (ONLY TRUE PEAKS) ----------------
IEDdata.foundPeaks(1).locs  = truePeaks;
IEDdata.detections(1).times = truePeaks / Fs;

end
