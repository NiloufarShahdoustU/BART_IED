% find out the spectogram of IED and non-IED trials:
% you can find explanation in notes

% TODO: you need to change this code in the ranksum part and use cluster
% based permuation test. 

clear;
clc;
close all;
%% reading data
inputFolderName_LFPmat = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\LFPmat_bad_chans_removed';
inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\IEDdata_bad_chans_removed';

outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\IED_nonIED_freq_time_pValues\';
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
% for pt = 50:50
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    %IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.IEDdata.mat'];
    load(IEDdata);
    IEDtrials = IEDdata.IEDtrials;
    SelectedChans = IEDdata.selectedChans;
    anatomicalLocs = IEDdata.anatomicalLocs;
    SelectedAnatomicalLoc = anatomicalLocs(SelectedChans(1:end-1));

    
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
    
   %the second dimension of this should be 38, cause we do have 38 sFreqs
   pValueMatrix = nan(nChans,38, nSignal);

    % analysis

    for ch = 1:nChans
    % for ch = 1:1
        for tt = 1: nTrials
            fprintf('\nchannel %d of %d',ch,nChans)
            fprintf('\ntrial %d of %d',tt,nTrials)

    % IED analysis:
        % calculate scalograms
        % if(LFPmat_IED(ch,:,tt))
            [wave_IED,period_IED,scale_IED,coi_IED] = basewaveERP(LFPmat_IED(ch,:,tt),Fs, loF, hiF, MotherWaveParam,waitc);   
            sFreqs_IED = period_IED.^-1; % first find frequencies
            % convert to power in decibels
            SdB_IED(:,:,tt) = abs(10*log10(wave_IED));
            % just getting components
            % spectrum [freqs,samples,trials]
            S_IED(:,:,tt) = abs(wave_IED);
            % phase
            phi_IED(:,:,tt) = angle(wave_IED);

            
        % end


    % nonIED analysis:
            [wave_nonIED,period_nonIED,scale_nonIED,coi_nonIED] = basewaveERP(squeeze(LFPmat_nonIED(ch,:,tt)),Fs, loF, hiF, MotherWaveParam,waitc);
            sFreqs_nonIED = period_nonIED.^-1; % first find frequencies
            % convert to power in decibels
            SdB_nonIED(:,:,tt) = abs(10*log10(wave_nonIED));
            % just getting components
            % spectrum [scales,samples,trials]
            S_nonIED(:,:,tt) = abs(wave_nonIED);
            % phase
            phi_nonIED(:,:,tt) = angle(wave_nonIED);
    
        end
    
    
        % doing baseline normalization
        % pick baseline window
        bWin = 1:nSignal;
        Sbl_IED = S_IED./repmat(nanmean(nanmean(S_IED(:,bWin,:),2),3),1,size(S_IED,2),size(S_IED,3));
        Sbl_nonIED = S_nonIED./repmat(nanmean(nanmean(S_nonIED(:,bWin,:),2),3),1,size(S_nonIED,2),size(S_nonIED,3));
        
        % now I need to delete the Nan trials
        nan_slices_IED = all(all(isnan(Sbl_IED), 1), 2);
        Sbl_IED(:,:,nan_slices_IED) = [];

        nan_slices_nonIED = all(all(isnan(Sbl_nonIED), 1), 2);
        Sbl_nonIED(:,:,nan_slices_nonIED) = [];
        

        if(~all(isnan(Sbl_IED), 'all')) % here I am checking whether we have IED trials or not :D
            for freq=1:size(Sbl_IED,1)
                for time=1:size(Sbl_IED,2)
                    FreqTime_IED = squeeze(Sbl_IED(freq, time, :));
                    FreqTime_nonIED = squeeze(Sbl_nonIED(freq, time, :));
    
                    %test
                    %TODO: this ranksum test is not correct you need to
                    %do cluster based permutation test
                    pValue = ranksum(FreqTime_IED, FreqTime_nonIED);
                    pValueMatrix(ch, freq,time) = pValue;
                end
            end    
        end

        clear wave_IED period_IED scale_IED coi_IED Sbl_IED
        clear wave_nonIED period_nonIED scale_nonIED coi_nonIED Sbl_nonIED
    end 

    % I need to save pValue matrix
    pValueMatrix_Data.pValueMatrix = pValueMatrix;
    pValueMatrix_Data.sFreqs = sFreqs_IED;
    pValueMatrix_Data.period = period_IED;
    pValueMatrix_Data.SelectedAnatomicalLoc = SelectedAnatomicalLoc;
    save([outputFolderName ptID '.pValue.mat'],'pValueMatrix_Data')

    clear LFPmat_IED LFPmat_nonIED ;
end



%% debug:
% bbbb = squeeze(LFPmat_IED(1,:,1));
% aaaa =  squeeze(SdB_IED(24,:,120));


% for ch = 11:-1:1
%     disp(ch);
% end

% 
% % Example 3D matrix
% A = rand(4, 5, 3); % Creating a random 4x5x3 matrix
% A(:,:,2) = NaN;    % Setting all elements of the second "slice" to NaN
% %
% % Identify slices where all elements are NaN
% nan_slices = all(all(isnan(A), 1), 2);
% 
% % Remove those slices from the third dimension
% A(:,:,nan_slices) = [];
