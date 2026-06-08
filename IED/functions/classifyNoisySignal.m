function isNoisy = classifyNoisySignal(signal, threshold)
    % Calculate the power of the signal
    signalPower = var(signal);

    % Estimate the noise power
    noiseEstimate = medfilt1(signal, 10) - signal;  % Use median filtering to estimate noise
    noisePower = var(noiseEstimate);

    % Compute the Signal-to-Noise Ratio (SNR)
    snrValue = 10 * log10(signalPower / noisePower);

    % Classify as noisy or not based on the threshold
    if snrValue < threshold
        isNoisy = 1;  % Signal is super noisy
    else
        isNoisy = 0;  % Signal is not super noisy
    end
end
