function isOutOfRange = detectLargeAmplitude(lfpData, threshold)
    % Detect if any channels in the LFP data have an amplitude outside the threshold.

    % threshold is the amplitude limit (e.g., 1000).
    % Output isOutOfRange is 1 if there are outliers, 0 otherwise.
    
    % Check if any element in the array exceeds the threshold in magnitude
    if any(abs(lfpData(:)) > threshold)
        isOutOfRange = 1;  % Return 1 if outliers are found
    else
        isOutOfRange = 0;  % Return 0 if no outliers
    end
end
