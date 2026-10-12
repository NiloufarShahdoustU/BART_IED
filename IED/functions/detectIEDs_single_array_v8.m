function [IEDdata] = detectIEDs_single_array_v8(data, Fo)
% DETECTIEDS_SINGLE_ARRAY_V5
% Detect IEDs using Hilbert envelope
% Save true LFP waveform peaks after refinement

Fs = Fo;

%---------------- parameters ----------------
freqRange    = [5 25];
smoothingVar = 10;
minPeakFrac  = 0.6;          % 60% of max envelope height
minISI       = 0.250;        % 250 ms
searchWinSec = 0.050;        % +/- 50 ms for true peak

% ---------------- make data row vector ----------------
data = double(data(:)');

% ---------------- bandpass ----------------
[b, a] = butter(4, freqRange/(Fs/2), 'bandpass');
filteredData = filtfilt(b, a, data);

% ---------------- envelope detection ----------------
env = abs(hilbert(filteredData));

smoothWin = round(Fs / smoothingVar);
envSmooth = smooth(env, smoothWin);

minPeakHeight = minPeakFrac * max(envSmooth);

[~, locs] = findpeaks(envSmooth, ...
    'MinPeakHeight', minPeakHeight, ...
    'MinPeakDistance', round(minISI * Fs));

% ---------------- TRUE LFP PEAK REFINEMENT ----------------
searchWin = round(searchWinSec * Fs);
truePeaks = zeros(size(locs));

for i = 1:length(locs)

    c = locs(i);
    i1 = max(1, c - searchWin);
    i2 = min(length(filteredData), c + searchWin);

    seg = filteredData(i1:i2);

    % true absolute waveform peak in filtered LFP
    [~, idx] = max(abs(seg));
    truePeaks(i) = i1 + idx - 1;

end

% ---------------- enforce refractory again after refinement ----------------
truePeaks = sort(truePeaks);

minSamplesApart = round(minISI * Fs);
keep = true(size(truePeaks));
lastLoc = -inf;

for i = 1:length(truePeaks)
    if truePeaks(i) - lastLoc < minSamplesApart
        keep(i) = false;
    else
        lastLoc = truePeaks(i);
    end
end

truePeaks = truePeaks(keep);

% ---------------- OUTPUT ----------------
IEDdata.foundPeaks(1).locs  = truePeaks;
IEDdata.detections(1).times = truePeaks / Fs;

end