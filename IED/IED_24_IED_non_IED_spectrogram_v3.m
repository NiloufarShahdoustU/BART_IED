% ccn 2026
% Spectrograms WITH epoch dimension
% UPDATED: keep per-trial spectrograms + also store mean over trials per channel
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
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_24_IED_non_IED_spectrogram_v3\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));
PatientsNum = length(fileList);

Fs = 1000;
loF = 1;
hiF = 200;
MotherWaveParam = 6;
waitc = 0;

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


    %+++++++++++++++++++++++++++++++++++++++++++++++++++
    % dearling with finding the common brain areas:

    AnatomicalLocs = LFPIED.anatomicalLocs;
    selectedChannels = LFPIED.selectedChans;

    %IED significant
    BrainAreaIEDsignificant = AnatomicalLocs(selectedChannels(IEDsignificantChannelIndices));
    BrainAreaIEDsignificant = erase( BrainAreaIEDsignificant, ["Left ", "Right ", " Left", " Right"]);

    %IED Non-significant
    BrainAreaIEDNonSignificant = AnatomicalLocs(selectedChannels(IEDNonSignificantChannelIndices));
    BrainAreaIEDNonSignificant = erase( BrainAreaIEDNonSignificant, ["Left ", "Right ", " Left", " Right"]);

    %nonIED
    BrainAreaNonIED = AnatomicalLocs(selectedChannels(nonIEDChannelIndices));
    BrainAreaNonIED = erase( BrainAreaNonIED, ["Left ", "Right ", " Left", " Right"]);

    BrainAreasCommon = intersect( ...
        intersect(BrainAreaIEDsignificant, BrainAreaIEDNonSignificant), ...
        BrainAreaNonIED);

    % IMPORTANT condition:
    % OLNLY if we have any common brain areas between our 3 conditions, we can
    % continue on our analysis!!

    if(~isempty(BrainAreasCommon))
        % --- IED significant ---
        isCommon_IEDsig = ismember(BrainAreaIEDsignificant, BrainAreasCommon);
        CommonChansIndicesIEDSignificant = IEDsignificantChannelIndices(isCommon_IEDsig);

        CommonBrainAreaIEDsignificant = BrainAreaIEDsignificant(isCommon_IEDsig);

        % --- IED non-significant ---
        isCommon_IEDnonsig = ismember(BrainAreaIEDNonSignificant, BrainAreasCommon);
        CommonChansIndicesIEDNonSignificant = IEDNonSignificantChannelIndices(isCommon_IEDnonsig);

        CommonBrainAreaIEDNonSignificant = BrainAreaIEDNonSignificant(isCommon_IEDnonsig);

        % --- non-IED ---
        isCommon_nonIED = ismember(BrainAreaNonIED, BrainAreasCommon);
        CommonChansIndicesNonIED = nonIEDChannelIndices(isCommon_nonIED);

        CommonBrainAreaNonIED = BrainAreaNonIED(isCommon_nonIED);

        %+++++++++++++++++++++++++++++++++++++++++++++++++++

        % ==========================
        % STORE SPECS PER EPOCH HERE
        % ==========================
        % Each epoch stores a 2D cell: {channel, trial}
        % IMPORTANT: store LINEAR power (abs(wave).^2) and convert to dB only for plotting
        Spec_IED_sig_epoch    = cell(length(epochList), 1);
        Spec_IED_nonsig_epoch = cell(length(epochList), 1);
        Spec_nonIED_epoch     = cell(length(epochList), 1);

        % keep refs per epoch
        period_ref_epoch = cell(length(epochList), 1);
        scale_ref_epoch  = cell(length(epochList), 1);
        coi_ref_epoch    = cell(length(epochList), 1);

        % keep trial indices per epoch so you know what you computed
        IEDtrialIdx_epoch    = cell(length(epochList), 1);
        nonIEDtrialIdx_epoch = cell(length(epochList), 1);

        for e = 1:length(epochList)

            IEDtrials_epoch = epochList{e};
            LFP_epoch       = LFPepochList{e};
            epochName       = epochNames{e};

            IEDtrialIdx    = find(~isControl & any(IEDtrials_epoch(IEDsignificantChannelIndices, :) ~= 0, 1));
            nonIEDtrialIdx = find(~isControl & all(IEDtrials_epoch(:, :) == 0, 1));

            IEDtrialIdx_epoch{e}    = IEDtrialIdx;
            nonIEDtrialIdx_epoch{e} = nonIEDtrialIdx;

            fprintf('%s , epoch: %d/%d\n', epochName, e, length(epochList));

            period_ref = [];
            scale_ref  = [];
            coi_ref    = [];

            % Allocate per-epoch containers
            Spec_IED_sig_epoch{e}    = cell(length(CommonChansIndicesIEDSignificant), length(IEDtrialIdx));
            Spec_IED_nonsig_epoch{e} = cell(length(CommonChansIndicesIEDNonSignificant), length(IEDtrialIdx));
            Spec_nonIED_epoch{e}     = cell(length(CommonBrainAreaNonIED), length(nonIEDtrialIdx));

            % ------------------ IED significant ------------------
            for c = 1:length(CommonChansIndicesIEDSignificant)

                ch = CommonChansIndicesIEDSignificant(c);

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

                    powtemp = abs(wave).^2;
                    Spec_IED_sig_epoch{e}{c,t} = powtemp(:, midIdx);

                    if isempty(period_ref)
                        period_ref = period;
                        scale_ref  = scale;
                        coi_ref    = coi;
                    end
                end

            end

            period_ref_epoch{e} = period_ref;
            scale_ref_epoch{e}  = scale_ref;
            coi_ref_epoch{e}    = coi_ref;

            % ---------------- IED non-significant ----------------
            for c = 1:length(CommonChansIndicesIEDNonSignificant)

                ch = CommonChansIndicesIEDNonSignificant(c);

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
                    powtemp = abs(wave).^2;
                    Spec_IED_nonsig_epoch{e}{c,t} = powtemp(:, midIdx);
                end

            end

            % ---------------------- non-IED ----------------------
            for c = 1:length(CommonChansIndicesNonIED)

                ch = CommonChansIndicesNonIED(c);

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
                    powtemp = abs(wave).^2;
                    Spec_nonIED_epoch{e}{c,t} = powtemp(:, midIdx);

                end

            end

        end
    end

    % visualization per participants:

    % finding which data to look at:
    idxIEDsignificant    = cell(numel(BrainAreasCommon),1);
    idxIEDNonSignificant = cell(numel(BrainAreasCommon),1);
    idxNonIED            = cell(numel(BrainAreasCommon),1);
    for i = 1:numel(BrainAreasCommon)
        area = BrainAreasCommon{i};
        idxIEDsignificant{i} = find(strcmp(CommonBrainAreaIEDsignificant, area));
        idxIEDNonSignificant{i} = find(strcmp(CommonBrainAreaIEDNonSignificant, area));
        idxNonIED{i} = find(strcmp(CommonBrainAreaNonIED, area));
    end

    t_ms = (0:(numel(midIdx)-1)) / Fs * 1000;

    % frequency axis (from period_ref_epoch). basewave typically returns period in seconds.
    % Use first epoch as reference (they should match across epochs).
    if ~isempty(period_ref_epoch{1})
        f_Hz = 1 ./ period_ref_epoch{1};
    else
        f_Hz = 1:size(Spec_nonIED_epoch{1}{1,1},1); % fallback
    end

    % helper: average spectrogram over selected channels + all trials, then average across epochs
    avgSpecAcrossEpochs = @(SpecEpoch, chanIdxPerArea) local_avg_spec_across_epochs(SpecEpoch, chanIdxPerArea);

    for iA = 1:numel(BrainAreasCommon)

        areaName = BrainAreasCommon{iA};

        % indices into the *common* channel lists (i.e., rows of Spec_*_epoch{e})
        chIdx_sig    = idxIEDsignificant{iA};
        chIdx_nonsig = idxIEDNonSignificant{iA};
        chIdx_nonIED = idxNonIED{iA};

        % compute mean spectrograms (freq x time) in LINEAR power
        S_sig_lin    = avgSpecAcrossEpochs(Spec_IED_sig_epoch,    chIdx_sig);
        S_nonsig_lin = avgSpecAcrossEpochs(Spec_IED_nonsig_epoch, chIdx_nonsig);
        S_nonIED_lin = avgSpecAcrossEpochs(Spec_nonIED_epoch,     chIdx_nonIED);

        % convert to dB only for plotting
        S_sig    = 10*log10(S_sig_lin    + eps);
        S_nonsig = 10*log10(S_nonsig_lin + eps);
        S_nonIED = 10*log10(S_nonIED_lin + eps);

        % ---- plot ----
        fig = figure('Color','w','Position',[100 100 1200 400], 'Visible','off');
        tl = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

        % ---- robust shared CLim across the 3 conditions ----
        allVals = [S_sig(:); S_nonsig(:); S_nonIED(:)];
        allVals = allVals(isfinite(allVals));
        clim = [prctile(allVals, 5), prctile(allVals, 95)];

        % ---- y-axis (log) settings (show as 10^k) ----
        fPlot = f_Hz(:);
        fWin  = [min(fPlot) max(fPlot)];   % or e.g. [1 150]

        kmin = ceil(log10(fWin(1)));
        kmax = floor(log10(fWin(2)));
        yt   = 10.^(kmin:kmax);
        ytLbl = arrayfun(@(k) sprintf('10^{%d}', k), kmin:kmax, 'UniformOutput', false);

        % time axis in ms should match your data columns
        x = t_ms;
        xlimWanted = [x(1) x(end)];

        % --------- panel 1 ---------
        ax1 = nexttile;
        imagesc(ax1, x, fPlot, S_sig);
        set(ax1,'YDir','normal','YScale','log');
        xlim(ax1, xlimWanted);
        ylim(ax1, fWin);
        yticks(ax1, yt);
        yticklabels(ax1, ytLbl);
        ax1.CLim = clim;
        xlabel(ax1,'Time (ms)'); ylabel(ax1,'Frequency (Hz)');
        title(ax1,'IED significant');

        % --------- panel 2 ---------
        ax2 = nexttile;
        imagesc(ax2, x, fPlot, S_nonsig);
        set(ax2,'YDir','normal','YScale','log');
        xlim(ax2, xlimWanted);
        ylim(ax2, fWin);
        yticks(ax2, yt);
        yticklabels(ax2, ytLbl);
        ax2.CLim = clim;
        xlabel(ax2,'Time (ms)'); ylabel(ax2,'Frequency (Hz)');
        title(ax2,'IED non-significant');

        % --------- panel 3 ---------
        ax3 = nexttile;
        imagesc(ax3, x, fPlot, S_nonIED);
        set(ax3,'YDir','normal','YScale','log');
        xlim(ax3, xlimWanted);
        ylim(ax3, fWin);
        yticks(ax3, yt);
        yticklabels(ax3, ytLbl);
        ax3.CLim = clim;
        xlabel(ax3,'Time (ms)'); ylabel(ax3,'Frequency (Hz)');
        title(ax3,'non-IED');

        % ---- one shared colorbar on the far right (vertical) ----
        cb = colorbar(ax3, 'eastoutside');
        cb.Layout.Tile = 'east';
        cb.Label.String = 'power (dB)';

        sgtitle(sprintf('%s | %s', ptID, areaName), 'Interpreter','none');

        % safe filename parts
        areaSafe = regexprep(areaName, '[^\w\-]+', '_');
        fnameBase = sprintf('%s_%s_spectrogram', ptID, areaSafe);

        set(fig,'Renderer','painters');
        exportgraphics(fig, fullfile(outputFolderName, [fnameBase '.pdf']), ...
            'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);

    end

end
%%
function Smean = local_avg_spec_across_epochs(SpecEpoch, chanIdxPerArea)
% SpecEpoch{e} is a cell matrix: (channels_in_common_list x trials_in_that_condition)
% Each entry is (freq x time) LINEAR power (abs(wave).^2)

    Smean = [];
    if isempty(chanIdxPerArea)
        return;
    end

    epochMeans = {};

    for e = 1:numel(SpecEpoch)
        C = SpecEpoch{e};
        if isempty(C)
            continue;
        end

        mats = {};
        for ii = 1:numel(chanIdxPerArea)
            r = chanIdxPerArea(ii);
            if r < 1 || r > size(C,1)
                continue;
            end
            for t = 1:size(C,2)
                if ~isempty(C{r,t})
                    mats{end+1} = C{r,t};
                end
            end
        end

        if ~isempty(mats)
            S3 = cat(3, mats{:});
            epochMeans{end+1} = mean(S3, 3, 'omitnan');
        end
    end

    if isempty(epochMeans)
        return;
    end

    Sepoch3 = cat(3, epochMeans{:});
    Smean = mean(Sepoch3, 3, 'omitnan');
end
