% ccn 2026
% I am going to have spectrograms (WITH epoch dimension)
% UPDATED: store ONLY MEAN spectrogram over TRIALS (no per-trial cells)
% For each epoch, for each channel group, we compute mean POWER spectrogram across trials.
% AUTHOR: Nill (edited)

clear; clc; close all;
warning('off','all');

SigChanFolder_IT = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_ITs_2_chunks';
matFile_IT  = dir(fullfile(SigChanFolder_IT, 'sigChanNums_ITs_PostResponse_PreOutcome.mat'));
filePath_IT = fullfile(SigChanFolder_IT, matFile_IT(1).name);
SigChans_IT = load(filePath_IT);

SigChanFolder_RT = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_RTs_2_chunks';
matFile_RT  = dir(fullfile(SigChanFolder_RT, 'sigChanNums_RTs_PostOnset_PreResponse.mat'));
filePath_RT = fullfile(SigChanFolder_RT, matFile_RT(1).name);
SigChans_RT = load(filePath_RT);

SigChanFolder_ACC = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_acc_2_chunks';
matFile_ACC  = dir(fullfile(SigChanFolder_ACC, 'sigChanNums_acc_PostResponse_PreOutcome.mat'));
filePath_ACC = fullfile(SigChanFolder_ACC, matFile_ACC(1).name);
SigChans_ACC = load(filePath_ACC);

inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_24_IED_non_IED_spectrogram_v2\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));
PatientsNum = length(fileList);

Fs = 1000;
loF = 1;
hiF = 150;
MotherWaveParam = 6;
waitc = 0;

% Keep same middle window (500 ms if 1000 samples): 250:749
midIdx = 250:749;

for pt = 1:PatientsNum
% for pt = 1:1

    ptID = erase(fileList(pt).name, '.LFPIED.mat');
    disp("Processing patient: " + ptID);

    load(fullfile(inputFolderName_IEDtrials, fileList(pt).name));

    epochList = {
        LFPIED.IEDtrialsPostOnset;
        LFPIED.IEDtrialsPreResponse;
        LFPIED.IEDtrialsPostResponse;
        LFPIED.IEDtrialsPreOutcome;
    };

    LFPepochList = {
        LFPIED.LFPmatPostOnset;
        LFPIED.LFPmatPreResponse;
        LFPIED.LFPmatPostResponse;
        LFPIED.LFPmatPreOutcome;
    };

    epochNames = {'PostOnset','PreResponse','PostResponse','PreOutcome'};

    isControl = LFPIED.isControl;

    pre_IT  = SigChans_IT.sigChanNums_PreOutcome{pt,1};
    post_IT = SigChans_IT.sigChanNums_PostResponse{pt,1};

    pre_ACC  = SigChans_ACC.sigChanNums_PreOutcome{pt,1};
    post_ACC = SigChans_ACC.sigChanNums_PostResponse{pt,1};

    pre_RT  = SigChans_RT.sigChanNums_PostOnset{pt,1};
    post_RT = SigChans_RT.sigChanNums_PreResponse{pt,1};

    significantChannelIndices_temp = unique([pre_IT(:); post_IT(:); pre_ACC(:); post_ACC(:); pre_RT(:); post_RT(:)], 'stable');

    [~, IEDsignificantChannelIndices] = ismember(significantChannelIndices_temp, LFPIED.selectedChans);
    IEDsignificantChannelIndices = sort(IEDsignificantChannelIndices);

    IEDchannels_all = [];

    for e = 1:length(epochList)
        IEDtrials = epochList{e};
        chansWithIED = find(any(IEDtrials(:, ~isControl) ~= 0, 2));
        IEDchannels_all = [IEDchannels_all; chansWithIED(:)];
    end

    IEDchannels_all = unique(IEDchannels_all, 'stable');

    IEDNonSignificantChannelIndices = setdiff(IEDchannels_all, IEDsignificantChannelIndices, 'stable');
    IEDNonSignificantChannelIndices = sort(IEDNonSignificantChannelIndices);

    allChannelIndices = 1:length(LFPIED.selectedChans);
    nonIEDChannelIndices = setdiff(allChannelIndices, [IEDNonSignificantChannelIndices(:); IEDsignificantChannelIndices(:)], 'stable');

    % ==========================
    % STORE MEAN SPECS PER EPOCH
    % ==========================
    % Each epoch stores a 1D cell:
    %   SpecMean_..._epoch{e}{c} = [freq x time] mean POWER over trials
    SpecMean_IED_sig_epoch    = cell(length(epochList), 1);
    SpecMean_IED_nonsig_epoch = cell(length(epochList), 1);
    SpecMean_nonIED_epoch     = cell(length(epochList), 1);

    % keep refs per epoch
    period_ref_epoch = cell(length(epochList), 1);
    scale_ref_epoch  = cell(length(epochList), 1);
    coi_ref_epoch    = cell(length(epochList), 1);

    % keep trial indices per epoch
    IEDtrialIdx_epoch    = cell(length(epochList), 1);
    nonIEDtrialIdx_epoch = cell(length(epochList), 1);

    for e = 1:length(epochList)

        IEDtrials_epoch = epochList{e};
        LFP_epoch       = LFPepochList{e};
        epochName       = epochNames{e};

        % ---- trial definitions (same logic as your code) ----
        IEDtrialIdx    = find(~isControl & any(IEDtrials_epoch(IEDsignificantChannelIndices, :) ~= 0, 1));
        nonIEDtrialIdx = find(~isControl & all(IEDtrials_epoch(:, :) == 0, 1));

        IEDtrialIdx_epoch{e}    = IEDtrialIdx;
        nonIEDtrialIdx_epoch{e} = nonIEDtrialIdx;

        fprintf('%s | epoch: %d/%d\n', epochName, e, length(epochList));

        period_ref = [];
        scale_ref  = [];
        coi_ref    = [];

        % allocate epoch cells (no trial dim)
        SpecMean_IED_sig_epoch{e}    = cell(length(IEDsignificantChannelIndices), 1);
        SpecMean_IED_nonsig_epoch{e} = cell(length(IEDNonSignificantChannelIndices), 1);
        SpecMean_nonIED_epoch{e}     = cell(length(nonIEDChannelIndices), 1);

        % =======================
        % IED SIGNIFICANT CHANNELS
        % =======================
        for c = 1:length(IEDsignificantChannelIndices)

            ch = IEDsignificantChannelIndices(c);

            if isempty(IEDtrialIdx)
                SpecMean_IED_sig_epoch{e}{c} = [];
                continue
            end

            sumPow = [];
            nUsed  = 0;

            for t = 1:length(IEDtrialIdx)
                tr = IEDtrialIdx(t);

                signal = squeeze(LFP_epoch(ch, :, tr));
                signal = double(signal(:));

                targetN = 1000;
                if length(signal) < targetN
                    padN = targetN - length(signal);
                    leftPad  = floor(padN/2);
                    rightPad = padN - leftPad;
                    signal = [zeros(leftPad,1); signal; zeros(rightPad,1)];
                end

                [wave, period, scale, coi] = basewaveERP(signal, Fs, loF, hiF, MotherWaveParam, waitc);

                % Use POWER for spectrogram
                powMid = abs(wave(:, midIdx)).^2;

                if isempty(sumPow)
                    sumPow = zeros(size(powMid));
                end

                sumPow = sumPow + powMid;
                nUsed  = nUsed + 1;

                if isempty(period_ref)
                    period_ref = period;
                    scale_ref  = scale;
                    coi_ref    = coi;
                end
            end

            SpecMean_IED_sig_epoch{e}{c} = sumPow ./ max(nUsed, 1);

        end

        period_ref_epoch{e} = period_ref;
        scale_ref_epoch{e}  = scale_ref;
        coi_ref_epoch{e}    = coi_ref;

        % ==========================
        % IED NON-SIGNIFICANT CHANNELS
        % ==========================
        for c = 1:length(IEDNonSignificantChannelIndices)

            ch = IEDNonSignificantChannelIndices(c);

            if isempty(IEDtrialIdx)
                SpecMean_IED_nonsig_epoch{e}{c} = [];
                continue
            end

            sumPow = [];
            nUsed  = 0;

            for t = 1:length(IEDtrialIdx)
                tr = IEDtrialIdx(t);

                signal = squeeze(LFP_epoch(ch, :, tr));
                signal = double(signal(:));

                targetN = 1000;
                if length(signal) < targetN
                    padN = targetN - length(signal);
                    leftPad  = floor(padN/2);
                    rightPad = padN - leftPad;
                    signal = [zeros(leftPad,1); signal; zeros(rightPad,1)];
                end

                [wave, ~, ~, ~] = basewaveERP(signal, Fs, loF, hiF, MotherWaveParam, waitc);

                powMid = abs(wave(:, midIdx)).^2;

                if isempty(sumPow)
                    sumPow = zeros(size(powMid));
                end

                sumPow = sumPow + powMid;
                nUsed  = nUsed + 1;
            end

            SpecMean_IED_nonsig_epoch{e}{c} = sumPow ./ max(nUsed, 1);

        end

        % ==============
        % non-IED CHANNELS
        % ==============
        for c = 1:length(nonIEDChannelIndices)

            ch = nonIEDChannelIndices(c);

            if isempty(nonIEDtrialIdx)
                SpecMean_nonIED_epoch{e}{c} = [];
                continue
            end

            sumPow = [];
            nUsed  = 0;

            for t = 1:length(nonIEDtrialIdx)
                tr = nonIEDtrialIdx(t);

                signal = squeeze(LFP_epoch(ch, :, tr));
                signal = double(signal(:));

                targetN = 1000;
                if length(signal) < targetN
                    padN = targetN - length(signal);
                    leftPad  = floor(padN/2);
                    rightPad = padN - leftPad;
                    signal = [zeros(leftPad,1); signal; zeros(rightPad,1)];
                end

                [wave, ~, ~, ~] = basewaveERP(signal, Fs, loF, hiF, MotherWaveParam, waitc);

                powMid = abs(wave(:, midIdx)).^2;

                if isempty(sumPow)
                    sumPow = zeros(size(powMid));
                end

                sumPow = sumPow + powMid;
                nUsed  = nUsed + 1;
            end

            SpecMean_nonIED_epoch{e}{c} = sumPow ./ max(nUsed, 1);

        end
    end
        outFile = fullfile(outputFolderName, [ptID '_SpecMeans_epoch.mat']);

    % Save important refs + indices + means
    save(outFile, ...
        'ptID', 'epochNames', ...
        'IEDsignificantChannelIndices', 'IEDNonSignificantChannelIndices', 'nonIEDChannelIndices', ...
        'IEDtrialIdx_epoch', 'nonIEDtrialIdx_epoch', ...
        'SpecMean_IED_sig_epoch', 'SpecMean_IED_nonsig_epoch', 'SpecMean_nonIED_epoch', ...
        'period_ref_epoch', 'scale_ref_epoch', 'coi_ref_epoch', ...
        'Fs', 'loF', 'hiF', 'MotherWaveParam', 'midIdx', ...
        '-v7.3');


end




%% stuff
