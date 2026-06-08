% find out the spectogram of IED and non-IED trials:


clear;
clc;
close all;
%% reading data
inputFolderName_LFPmat = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\LFPmat_bad_chans_removed';
inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\IEDdata_bad_chans_removed';

outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_17_spectogram_v1\';
fileList = dir(fullfile(inputFolderName_LFPmat, '*.LFPmat.mat'));
PatientsNum = length(fileList);
%% parameters
Fs = 500;               % sampling rate
loF = 1;                % lower frequency of the range for which the transform is to be done
hiF = 150;              % higher frequency of the range for which the transform is to be done
MotherWaveParam = 6;    % the mother wavelet parameter (wavenumber); constant that has to do with converting from scales to fourier space
waitc = 0;              % a handle to the qiqi waitbar (usually 0)

%%
for pt = 1:PatientsNum
% for pt = 1:1
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    %IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.IEDdata.mat'];
    load(IEDdata);
    IEDtrials = IEDdata.IEDtrials;
    
    % LFP read
    LFPmat = [inputFolderName_LFPmat '\' ptID '.LFPmat.mat'];
    load(LFPmat);
    LFPmat = LFPmatStruct_new.LFPmat;

    nTrials = size(LFPmat,3);
    nSignal = size(LFPmat,2);
    nChans = size(LFPmat,1)-1;
    
    LFPmat_IED = nan(size(LFPmat));
    LFPmat_nonIED = nan(size(LFPmat));


    % TODO: TOOOOOOOOOOOOOODDDDDDOOOOOOOOOOOOOOOO: separate IED and nonIED lfpmat
    for chan=1:nChans
        for trial=1:nTrials
            if IEDtrials(chan, trial) == 1
                LFPmat_IED(chan,:, trial) = LFPmat(chan,:, trial);
            else
                LFPmat_nonIED(chan,:, trial) = LFPmat(chan,:, trial);
            end

        end
    end
     
    % we need to normalize these means! We count the number of IEDs and
    % nonIEDs and when finding AverageLFP_IED and AverageLFP_nonIED we need
    % to normalize that amount. 

    ChanTrialNum = size(IEDtrials,1)*size(IEDtrials,1);
    IEDnum = sum(IEDtrials(:));
    nonIEDnum = ChanTrialNum - IEDnum;

    %average over channels and trials and create a single signal
    AverageLFP_IED = (mean(LFPmat_IED, [1 3], 'omitnan') * IEDnum)/ChanTrialNum;
    AverageLFP_nonIED = (mean(LFPmat_nonIED, [1 3], 'omitnan')*nonIEDnum)/ChanTrialNum;


    % analysis
    % in the code below, the main input which is the first input is a
    % vector and the outputs are wave_IED (the size is freqs X time), period_IED and scale_IED
    % the size is the number of freqs and also coi_IED (the size is time)

    [wave_IED,period_IED,scale_IED,coi_IED] = basewaveERP(AverageLFP_IED,Fs, loF, hiF, MotherWaveParam,waitc);
    [wave_nonIED,period_nonIED,scale_nonIED,coi_nonIED] = basewaveERP(AverageLFP_nonIED,Fs, loF, hiF, MotherWaveParam,waitc);



    % Visualization
    t = 1:nSignal;
    % Define your time range and corresponding labels
    new_ticks = linspace(1, nSignal, 7); % Adjust '7' based on desired number of ticks
    new_labels = linspace(-2, 4, length(new_ticks)); % Corresponding labels from -2 to 4
    
    frequency_IED = 1 ./ period_IED;
    frequency_nonIED = 1 ./ period_nonIED;
    
    % Compute power in decibels
    PowerDecibel_IED = abs(10 * log10(wave_IED));
    PowerDecibel_nonIED = abs(10 * log10(wave_nonIED));
    
    % Compute the color axis limits based on both datasets
    cmin = min([PowerDecibel_IED(:); PowerDecibel_nonIED(:)]);
    cmax = max([PowerDecibel_IED(:); PowerDecibel_nonIED(:)]);
    
    % Visualization
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.5], 'Visible', 'off'); 
    
    % IED subplot
    subplot(1,2,1)
    pcolor(t, 1 ./ period_IED, PowerDecibel_IED);
    set(gca, 'YScale', 'log', 'YDir', 'normal');
    shading flat;
    colorbar;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('IED');
    caxis([cmin cmax]); % Set the color axis limits
    
    % Set x-axis ticks and labels
    xticks(new_ticks);
    xticklabels(arrayfun(@num2str, new_labels, 'UniformOutput', false));
    
    % Add cone of influence
    hold on;
    plot(t, 1 ./ coi_IED, 'w--', 'LineWidth', 2);
    
    % Add vertical dashed line at x = 1001
    xline(1001, '--k', 'LineWidth', 3, 'Color', 'r');
    
    hold off;
    
    % Non-IED subplot
    subplot(1,2,2)
    pcolor(t, 1 ./ period_nonIED, PowerDecibel_nonIED);
    set(gca, 'YScale', 'log', 'YDir', 'normal');
    shading flat;
    colorbar;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Non IED');
    caxis([cmin cmax]); % Set the color axis limits
    
    % Set x-axis ticks and labels
    xticks(new_ticks);
    xticklabels(arrayfun(@num2str, new_labels, 'UniformOutput', false));
    
    % Add cone of influence
    hold on;
    plot(t, 1 ./ coi_nonIED, 'w--', 'LineWidth', 2);
    
    % Add vertical dashed line at x = 1001
    xline(1001, '--k', 'LineWidth', 3, 'Color', 'r');
    
    hold off;
    
    % Set the overall title for the figure
    sgtitle(['patient = ' char(string(ptID))], 'FontSize', 18, 'FontWeight', 'bold');



    

    
    % Save figure
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = string(ptID) + '_scalogram';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

    clear AverageLFP_IED AverageLFP_IED IEDtrials LFPmat AverageLFP_IED AverageLFP_nonIED
    clear wave_IED period_IED scale_IED coi_IED wave_nonIED period_nonIED scale_nonIED coi_nonIED

end



%% debug:
