
% this is the visulization function for AES 2024 
% we created the last visualization using this code
% author: Nill

clc;
clear;
close all;


ptID = '202314a';


% nevList = dir(['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);
nevList = dir(['\\155.100.91.44\d\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);

if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesBART(ptID);

%%

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));

% data parameters
selectedChans = find(isECoG);
nChans = length(selectedChans)-1;
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;

% resampling LFP at Fnew sampling frequency
notchFilter = true;
Fnew = 500;
for ch = 1:nChans
    if notchFilter
        [b1,a1] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        tmp(ch,:) = filtfilt(b1,a1,resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs));
        
        [b2,a2] = iirnotch(120/(Fnew/2),(120/(Fnew/2))/50);
        data2K(ch,:) = filtfilt(b2,a2,tmp(ch,:));
    else
        data2K(ch,:) = resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs);
    end
end
clear NSX NEV
Fs = Fnew;

% common average re-referencing.
%   data2K = remove1stPC(data2K);

% timing parameters.
pre = 2;
post = 4;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% TF parameters
params.fpass = [1 250]; % pick a value that's 50 higher for wavelets
params.Fs = Fs;
params.dBconversion = false;
params.normalized = false;
params.theoreticalNorm = false;
params.baseline = true;

% picking pre-bird baseline period
if params.baseline
    baselineType = 'preCue';
    switch baselineType
        case {'preCue'}
            bP = [-2.5 -1.5];
        case {'preTask'}
            % normalizing based on the mean spectrum from a pre-task baseline epoch.
            secsPreTask = 50;
            if (ppData.Event.trigTimes(1)/3e4)<secsPreTask
                fprintf('only %d seconds before first trigger.',ppData.Event.trigTimes(1)./3e4);
                % TODO:: then do spectrum
            else
                % TODO:: do spectrum.
            end
    end
end

% task parameters in chronological order..
% There aren't any trigs that == 4
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
if length(balloonTimes)>length(inflateTimes)
    balloonTimes(end) = [];
end


% balloon color colormap - [yellow, orange, read, gray]
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% task parameters in chronological order..
respTimes = trigTimes(trigs==23);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
nTrials = length(outcomeType);


if length(balloonType)>nTrials; balloonType=balloonType(1:end-1); end


%% TODO:: don't analyz 'NaC' trodes...
 
% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
% plotting bank/pop responses
for ch2 = 1:nChans
    % epoch the spectral data for each channel.
    for tt = 1:nTrials
        updateUser('finished spectral calculations',tt,50,nTrials);
        % epoch the data here [channels X samples X trials]
        LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
    end
end
clear ch2 tt

%% with this code I found the best channels of 4 different areas with IEDs
% AES 2024 figures for each channel


plotRawTrials = true;
start_channel = 1;
end_channel = 65;

if plotRawTrials
    for chz = start_channel:end_channel
        figure('Position', [100, 100, 800, 600]);  % Adjust the position and size as needed

        hold on 

        IEDtrials = outliers(range(squeeze(LFPmat(chz,:,:))));
        portion_of_trials = 1;
        for trl = 1:size(LFPmat,3)/portion_of_trials
            if IEDtrials(trl)
                plot(tSec, LFPmat(chz,:,trl) + (500 * trl), 'r')
            else
                plot(tSec, LFPmat(chz,:,trl) + (500 * trl), 'k')
            end
        end

        hold off
        xlabel('time (s)')
        ylabel('Amplitude')
        title(['Electrode in ' anatomicalLocs(selectedChans(chz))])


    end
end


%%


plotRawTrials = true;
chosen_channels = [18, 57, 17, 16];

if plotRawTrials
  figure('Units', 'normalized', 'Position', [0.5, 0.5, 0.25, 0.28]);

    numRows = 2;
    numCols = 2;

    for chz = 1:length(chosen_channels)
        subplot(numRows, numCols, chz);
        hold on;

        IEDtrials = outliers(range(squeeze(LFPmat(chosen_channels(chz), :, :))));
        portion_of_trials = 1;
        start_from_trial = 1;

        for trl = start_from_trial:size(LFPmat, 3) / portion_of_trials
            if IEDtrials(trl)
                plot(tSec, LFPmat(chosen_channels(chz), :, trl) + (500 * trl), 'r');
            else
                plot(tSec, LFPmat(chosen_channels(chz), :, trl) + (500 * trl), 'k');
            end
        end

        hold off;

        % Remove y labels from the right subplots
        if mod(chz, numCols) ~= 1
            ylabel('');
            ylim([30000, 50000]);
            yticklabels({});
        else
            ylim([30000,50000]);
            yticklabels({});
        end

        title([anatomicalLocs(selectedChans(chosen_channels(chz)))], 'FontSize', 12);
        
        % Set x-axis tick labels to bold and font size 12
        set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontWeight', 'bold', 'FontSize', 12);
    end

    % Common x-axis label for the entire figure
    xlabel('time to stimulus onset (s)', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized', 'Position', [-0.15, -0.15, 0]);
    annotation('textbox', [0.09, 0.07, 0.5, 0.06], 'String', 'brain signal in different trials', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14, 'Rotation', 90);


    title_m = suptitle('brain areas with IED occurrence');
    set(title_m, 'FontSize', 22, 'FontWeight', 'bold');

    set(gca, 'box', 'off', 'tickdir', 'out');
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    saveas(gcf, '\\155.100.91.44\d\Code\Nill\BART\AES2024\IED6_8_10_output\IEDs.pdf');

end
