%% this is with hilbert transfrom
% in this version I am using minPeakProminent instead of minPeakHeight


function [IEDdata] = detectIEDs_single_array_v4(data, Fo)
    % DETECTIEDS finds interictal epileptiform discharges (IEDs) in a single EEG channel
    % data: A vector representing the EEG signal
    % Fo: Original sampling rate of the data
    
    Fs = Fo; % Sampling rate for processing
    
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
    
    % Compute the smoothed signal for peak detection
    smoothedData = smooth(abs(dataFs), Fs / smoothing_var);
    
    % Determine the maximum amplitude for reference
    maxValue = max(smoothedData);
    
    % Set MinPeakProminence based on a percentage of the maximum value
    minPeakProminence = 0.05 * maxValue; % For example, 5% of the max
    
    % Find peaks based on prominence
    [amp, locs] = findpeaks(smoothedData, 'MinPeakProminence', minPeakProminence);
    
    % Initialize variables for filtered peaks
    filteredLocs = [];
    lastAcceptedLoc = 0;
    
    minSamplesApart = 0.8 * Fs; % Minimum distance between peaks
    
    % Filter peaks based on the distance criterion
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

