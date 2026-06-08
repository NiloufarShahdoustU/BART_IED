function plotWaveletSpectrogram(signal, samplingRate, loF, hiF)
    % Perform Continuous Wavelet Transform
    [cwtCoeffs, frequencies, coi] = cwt(signal, 'amor', samplingRate, 'FrequencyLimits', [loF hiF]);
    
    % Create time vector
    t = (0:length(signal)-1) / samplingRate;
    
    % Plot the spectrogram (scalogram)
    figure;
    tms = t * 1000; % Convert to milliseconds for the x-axis
    surface(tms, frequencies, abs(cwtCoeffs));
    axis tight;
    shading flat;
    colormap(jet);
    colorbar;
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    title('Continuous Wavelet Transform (Spectrogram)');
    set(gca, 'YScale', 'log');
    
    % Plot cone of influence
    hold on;
    plot(tms, coi, 'w--', 'LineWidth', 2);
    hold off;
end
