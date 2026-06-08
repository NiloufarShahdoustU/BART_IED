% in this code I'm going to read the lfpmat created in the previous version
% of this code and then recreate another version of lfpmat that does not
% have noisy and bad channels and trials.
% this is only for the older version of analysis, the one for ccn 2024.
% from SfN 2024 forward I have IED_4_create_lfpmat-v2_6chunks


% IMPORTANT:
% you need to run IED_11_IED_nonIED_saving_IED_trials_v3 after this for
% have an updated dataset.

% author: Nill



clear;
clc;
close all;
warning('off','all');

%% loading data

inputFolderName_LFPmat = '\\155.100.91.44\d\Data\Nill\BART\LFPmat';
fileList = dir(fullfile(inputFolderName_LFPmat, '*.LFPmat.mat'));
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\LFPmat_bad_chans_removed\'; 
PatientsNum = length(fileList);

%%

for pt = 1:PatientsNum
% for pt = 57:61
    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1}; 
    disp(' ');
    disp(['Processing patient ID: ' patientID ' (' int2str(pt) '/' int2str(PatientsNum) ')']);
    ptID = patientID;


    matFile_bhvStruct = [inputFolderName_LFPmat '\' ptID '.LFPmat.mat'];
    load(matFile_bhvStruct);

    nTrials = LFPmatStruct.nTrials;
    nChans = length(LFPmatStruct.selectedChans);
    LFPmat = LFPmatStruct.LFPmat;



    LFPmat_new = nan(size(LFPmat)); % some of the channels and trials will be nan in here
    
    for chan = 1:nChans
        for trial = 1:nTrials
            LFPmat_temp_vec = squeeze(LFPmat(chan,:,trial));
            % now I need to check if this specific trial for each specific
            % channel is not a good signal. Then I won't even save it. 
            threshold = 7; % in dB

            % Classify signal
            isNoisy = classifyNoisySignal(LFPmat_temp_vec, threshold);
    
            if ~isNoisy
                LFPmat_new(chan,:,trial) = LFPmat_temp_vec;
            end
        end
        clear LFPmat_temp_vec
    end

    LFPmatStruct_new.LFPmat = LFPmat_new;
    LFPmatStruct_new.anatomicalLocs = LFPmatStruct.anatomicalLocs;
    LFPmatStruct_new.balloonType = LFPmatStruct.balloonType;
    LFPmatStruct_new.balloonTimes = LFPmatStruct.balloonTimes;
    LFPmatStruct_new.nTrials = LFPmatStruct.nTrials;
    LFPmatStruct_new.selectedChans = LFPmatStruct.selectedChans;
    LFPmatStruct_new.trodeLabels = LFPmatStruct.trodeLabels;

    save([outputFolderName ptID '.LFPmat.mat'],'LFPmatStruct_new')

    clear LFPmat_new nTrials nChans LFPmat
end 



%% debug


% 
% cleanSignal = squeeze(LFPmat(94,:,178));
% 
% noisySignal = squeeze(LFPmat(103,:,108));
% 
% 
% % Define SNR threshold
% threshold = 7; % in dB
% 
% % Classify both signals
% isNoisyClean = classifyNoisySignal(cleanSignal, threshold);
% isNoisyNoisy = classifyNoisySignal(noisySignal, threshold);
% 
% % Display results
% fprintf('Clean Signal is classified as noisy: %d\n', isNoisyClean);
% fprintf('Noisy Signal is classified as noisy: %d\n', isNoisyNoisy);
% 
% % Visualization
% figure;
% subplot(2,1,1);
% plot(cleanSignal);
% title('Clean Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
% 
% subplot(2,1,2);
% plot(noisySignal);
% title('Noisy Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
% 
% % Mark the detected noisy signal
% if isNoisyClean
%     subplot(2,1,1);
%     title('Clean Signal (Detected as Noisy)');
% else
%     subplot(2,1,1);
%     title('Clean Signal (Detected as Not Noisy)');
% end
% 
% if isNoisyNoisy
%     subplot(2,1,2);
%     title('Noisy Signal (Detected as Noisy)');
% else
%     subplot(2,1,2);
%     title('Noisy Signal (Detected as Not Noisy)');
% end
% 
% 
% 
% aaaaa = squeeze(LFPmat_new(103,:,108));