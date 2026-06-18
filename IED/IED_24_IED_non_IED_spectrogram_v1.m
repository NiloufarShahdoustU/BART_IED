%this is only a helper code, do not use this 
% this is bullshit
%use version v2
% find out the spectogram of IED and non-IED trials:
% I am using baseWaveERP for finding the spectrograms

clear;
clc;
close all;
%% reading data
inputFolderName_LFPmat = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\LFPmat_bad_chans_removed';
inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\IEDdata_bad_chans_removed';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_24_IED_non_IED_spectrogram_v1\';
fileList = dir(fullfile(inputFolderName_LFPmat, '*.LFPmat.mat'));
PatientsNum = length(fileList);

%% parameters
Fs = 1000;               % sampling rate
loF = 1;                % lower frequency of the range for which the transform is to be done
hiF = 150;              % higher frequency of the range for which the transform is to be done
MotherWaveParam = 6;    % the mother wavelet parameter (wavenumber); constant that has to do with converting from scales to fourier space
waitc = 0;              % a handle to the qiqi waitbar (usually 0)

%%
% for pt = 1:PatientsNum
for pt = 1:3
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


% ===================== VISUALIZATION: IED vs non-IED =====================
    
    freq_IED    = 1 ./ period_IED;
    freq_nonIED = 1 ./ period_nonIED;
    
    Spec_IED    = 10*log10(abs(wave_IED).^2);
    Spec_nonIED = 10*log10(abs(wave_nonIED).^2);
    
    timeVec = 1:size(Spec_IED,2);
    
    % shared color limits (robust, like Starling code)
    allVals = [Spec_IED(:); Spec_nonIED(:)];
    climVals = [prctile(allVals,25), prctile(allVals,100)];
    
    fig = figure('Visible','off','Position',[100 100 1200 500]);
    sgtitle(sprintf('%s | IED vs non-IED', ptID), ...
        'FontWeight','bold','Interpreter','none', 'Visible', 'off');
    
    
    subplot(1,2,1)
    imagesc(timeVec, freq_IED, Spec_IED);
    set(gca,'YDir','normal','YScale','log')
    ylim([loF hiF])
    caxis(climVals)
    axis square tight
    title('IED trials')
    xlabel('Time (samples)')
    ylabel('Frequency (Hz)')
    set(gca,'FontSize',10)
    
    
    subplot(1,2,2)
    imagesc(timeVec, freq_nonIED, Spec_nonIED);
    set(gca,'YDir','normal','YScale','log')
    ylim([loF hiF])
    caxis(climVals)
    axis square tight
    title('Non-IED trials')
    xlabel('Time (samples)')
    set(gca,'FontSize',10)
    
    
    h = colorbar('southoutside');
    h.Position = [0.25 0.08 0.5 0.03];
    h.Label.String = 'Power (dB)';
    h.FontSize = 10;
    
    set(fig,'Renderer','painters');
    
    exportgraphics(fig, ...
        fullfile(outputFolderName, sprintf('%s_IED_vs_nonIED_scalogram.pdf', ptID)), ...
        'ContentType','vector', ...
        'BackgroundColor','none', ...
        'Resolution',600);
    



    clear AverageLFP_IED AverageLFP_IED IEDtrials LFPmat AverageLFP_IED AverageLFP_nonIED
    clear wave_IED period_IED scale_IED coi_IED wave_nonIED period_nonIED scale_nonIED coi_nonIED

end



%% debug:

% aaaa = squeeze(IEDdata.IED_timepoints(5,:,:));