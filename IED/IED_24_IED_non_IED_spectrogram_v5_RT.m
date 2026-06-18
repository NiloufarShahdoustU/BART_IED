% ccn 2026
% in this code, I'm taking a look at different spectrograms seperately for
% each epoch.
% UPDATED: only for RTs
% UPDATED2: put different epochs in the SAME FIGURE (subplots) per brain area
% AUTHOR: Nill (edited)

clear; clc; close all;
warning('off','all');

SigChanFolder_RT = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_RTs_2_chunks';
matFile_RT  = dir(fullfile(SigChanFolder_RT, 'sigChanNums_RTs_PostOnset_PreResponse.mat'));
filePath_RT = fullfile(SigChanFolder_RT, matFile_RT(1).name);
SigChans_RT = load(filePath_RT);

inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_24_IED_non_IED_spectrogram_v5_RT\';

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
    };

    LFPepochList = {
        LFPIED.LFPmatPostOnset;
        LFPIED.LFPmatPreResponse;
    };

    epochNames = {'PostOnset','PreResponse'};

    isControl = LFPIED.isControl;

    pre_RT  = SigChans_RT.sigChanNums_PostOnset{pt,1};
    post_RT = SigChans_RT.sigChanNums_PreResponse{pt,1};

    significantChannelIndices_temp = unique([pre_RT(:); post_RT(:)], 'stable');

    [~, IEDsignificantChannelIndices] = ismember(significantChannelIndices_temp, LFPIED.selectedChans);
    IEDsignificantChannelIndices = sort(IEDsignificantChannelIndices);

    AnatomicalLocs = LFPIED.anatomicalLocs;
    selectedChannels = LFPIED.selectedChans;
    allChannelIndices = 1:length(selectedChannels);

    f_Hz_ref = [];
    t_ms = (0:(numel(midIdx)-1)) / Fs * 1000;

    % ---- store per-area, per-epoch averages so we can plot in one figure later ----
    areaData = struct(); % dynamic fields: areaSafe -> struct with fields per epoch

    for e = 1:length(epochList)

        IEDtrials_epoch = epochList{e};
        LFP_epoch       = LFPepochList{e};
        epochName       = epochNames{e};

        fprintf('%s , epoch: %d/%d\n', epochName, e, length(epochList));

        chansWithIED_epoch = find(any(IEDtrials_epoch(:, ~isControl) ~= 0, 2));
        chansWithIED_epoch = unique(chansWithIED_epoch, 'stable');

        IEDsig_epoch_idx = intersect(IEDsignificantChannelIndices(:), chansWithIED_epoch(:), 'stable');

        IEDnonsig_epoch_idx = setdiff(chansWithIED_epoch(:), IEDsignificantChannelIndices(:), 'stable');
        IEDnonsig_epoch_idx = sort(IEDnonsig_epoch_idx);

        nonIED_epoch_idx = setdiff(allChannelIndices(:), chansWithIED_epoch(:), 'stable');

        BrainAreaIEDsig_epoch = AnatomicalLocs(selectedChannels(IEDsig_epoch_idx));
        BrainAreaIEDsig_epoch = erase(BrainAreaIEDsig_epoch, ["Left ", "Right ", " Left", " Right"]);

        BrainAreaIEDnonsig_epoch = AnatomicalLocs(selectedChannels(IEDnonsig_epoch_idx));
        BrainAreaIEDnonsig_epoch = erase(BrainAreaIEDnonsig_epoch, ["Left ", "Right ", " Left", " Right"]);

        BrainAreaNonIED_epoch = AnatomicalLocs(selectedChannels(nonIED_epoch_idx));
        BrainAreaNonIED_epoch = erase(BrainAreaNonIED_epoch, ["Left ", "Right ", " Left", " Right"]);

        BrainAreasCommon_epoch = intersect( ...
            intersect(BrainAreaIEDsig_epoch, BrainAreaIEDnonsig_epoch), ...
            BrainAreaNonIED_epoch);

        if isempty(BrainAreasCommon_epoch)
            disp("  No common brain areas for this epoch. Skipping epoch: " + epochName);
            continue;
        end

        idxIEDsig_epoch    = cell(numel(BrainAreasCommon_epoch),1);
        idxIEDnonsig_epoch = cell(numel(BrainAreasCommon_epoch),1);
        idxNonIED_epoch    = cell(numel(BrainAreasCommon_epoch),1);
        for iA = 1:numel(BrainAreasCommon_epoch)
            area = BrainAreasCommon_epoch{iA};
            idxIEDsig_epoch{iA}    = find(strcmp(BrainAreaIEDsig_epoch, area));
            idxIEDnonsig_epoch{iA} = find(strcmp(BrainAreaIEDnonsig_epoch, area));
            idxNonIED_epoch{iA}    = find(strcmp(BrainAreaNonIED_epoch, area));
        end

        IEDtrialIdx_sig = find(~isControl & any(IEDtrials_epoch(IEDsig_epoch_idx, :) ~= 0, 1));
        IEDtrialIdx_nonsig = find(~isControl & any(IEDtrials_epoch(IEDnonsig_epoch_idx, :) ~= 0, 1));
        nonIEDtrialIdx = find(~isControl & all(IEDtrials_epoch(:, :) == 0, 1));

        if isempty(IEDtrialIdx_sig) || isempty(IEDtrialIdx_nonsig) || isempty(nonIEDtrialIdx)
            disp("  Not enough trials for one or more conditions in epoch. Skipping epoch: " + epochName);
            continue;
        end

        Spec_IED_sig    = cell(length(IEDsig_epoch_idx),    length(IEDtrialIdx_sig));
        Spec_IED_nonsig = cell(length(IEDnonsig_epoch_idx), length(IEDtrialIdx_nonsig));
        Spec_nonIED     = cell(length(nonIED_epoch_idx),    length(nonIEDtrialIdx));

        period_ref = [];
        scale_ref  = [];
        coi_ref    = [];

        for c = 1:length(IEDsig_epoch_idx)
            ch = IEDsig_epoch_idx(c);
            for t = 1:length(IEDtrialIdx_sig)
                tr = IEDtrialIdx_sig(t);

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
                Spec_IED_sig{c,t} = powtemp(:, midIdx);

                if isempty(period_ref)
                    period_ref = period;
                    scale_ref  = scale;
                    coi_ref    = coi;
                end
            end
        end

        for c = 1:length(IEDnonsig_epoch_idx)
            ch = IEDnonsig_epoch_idx(c);
            for t = 1:length(IEDtrialIdx_nonsig)
                tr = IEDtrialIdx_nonsig(t);

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
                Spec_IED_nonsig{c,t} = powtemp(:, midIdx);
            end
        end

        for c = 1:length(nonIED_epoch_idx)
            ch = nonIED_epoch_idx(c);
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
                Spec_nonIED{c,t} = powtemp(:, midIdx);
            end
        end

        if isempty(f_Hz_ref)
            if ~isempty(period_ref)
                f_Hz_ref = 1 ./ period_ref;
            else
                f_Hz_ref = 1:size(Spec_nonIED{1,1},1);
            end
        end

        for iA = 1:numel(BrainAreasCommon_epoch)

            areaName = BrainAreasCommon_epoch{iA};
            areaSafe = regexprep(areaName, '[^\w\-]+', '_');

            chIdx_sig_local    = idxIEDsig_epoch{iA};
            chIdx_nonsig_local = idxIEDnonsig_epoch{iA};
            chIdx_nonIED_local = idxNonIED_epoch{iA};

            S_sig_lin    = local_avg_spec_one_epoch(Spec_IED_sig,    chIdx_sig_local);
            S_nonsig_lin = local_avg_spec_one_epoch(Spec_IED_nonsig, chIdx_nonsig_local);
            S_nonIED_lin = local_avg_spec_one_epoch(Spec_nonIED,     chIdx_nonIED_local);

            if isempty(S_sig_lin) || isempty(S_nonsig_lin) || isempty(S_nonIED_lin)
                continue;
            end

            S_sig    = 10*log10(S_sig_lin    + eps);
            S_nonsig = 10*log10(S_nonsig_lin + eps);
            S_nonIED = 10*log10(S_nonIED_lin + eps);

            if ~isfield(areaData, areaSafe)
                areaData.(areaSafe).areaName = areaName;
            end

            areaData.(areaSafe).(epochName).S_sig    = S_sig;
            areaData.(areaSafe).(epochName).S_nonsig = S_nonsig;
            areaData.(areaSafe).(epochName).S_nonIED = S_nonIED;

        end

    end

    % ---- after all epochs: plot per area with epochs as subplot rows ----
    areaFields = fieldnames(areaData);
    for a = 1:numel(areaFields)

        areaSafe = areaFields{a};
        areaName = areaData.(areaSafe).areaName;

        haveEpoch = false(numel(epochNames),1);
        for e = 1:numel(epochNames)
            haveEpoch(e) = isfield(areaData.(areaSafe), epochNames{e});
        end
        if ~any(haveEpoch)
            continue;
        end

        % compute ONE shared clim across all epochs + conditions for this area
        allVals = [];
        for e = 1:numel(epochNames)
            ep = epochNames{e};
            if isfield(areaData.(areaSafe), ep)
                allVals = [allVals; ...
                    areaData.(areaSafe).(ep).S_sig(:); ...
                    areaData.(areaSafe).(ep).S_nonsig(:); ...
                    areaData.(areaSafe).(ep).S_nonIED(:)];
            end
        end
        allVals = allVals(isfinite(allVals));
        if isempty(allVals)
            continue;
        end
        clim = [prctile(allVals, 5), prctile(allVals, 99)];

        fPlot = f_Hz_ref(:);
        fWin  = [min(fPlot) max(fPlot)];

        kmin = ceil(log10(fWin(1)));
        kmax = floor(log10(fWin(2)));
        yt   = 10.^(kmin:kmax);
        ytLbl = arrayfun(@(k) sprintf('10^{%d}', k), kmin:kmax, 'UniformOutput', false);

        x = t_ms;
        xlimWanted = [x(1) x(end)];

        nEpochPlot = sum(haveEpoch);
        fig = figure('Color','w','Position',[100 100 1200 380*nEpochPlot], 'Visible','off');

        tiledlayout(nEpochPlot, 3, 'Padding','compact', 'TileSpacing','compact');

        row = 0;
        for e = 1:numel(epochNames)
            epochName = epochNames{e};
            if ~isfield(areaData.(areaSafe), epochName)
                continue;
            end
            row = row + 1;

            S_sig    = areaData.(areaSafe).(epochName).S_sig;
            S_nonsig = areaData.(areaSafe).(epochName).S_nonsig;
            S_nonIED = areaData.(areaSafe).(epochName).S_nonIED;

            ax1 = nexttile;
            imagesc(ax1, x, fPlot, S_sig);
            set(ax1,'YDir','normal','YScale','log');
            xlim(ax1, xlimWanted); ylim(ax1, fWin);
            yticks(ax1, yt); yticklabels(ax1, ytLbl);
            ax1.CLim = clim;
            xlabel(ax1,'Time (ms)'); ylabel(ax1,'Frequency (Hz)');
            title(ax1, sprintf('%s | IED significant', epochName), 'Interpreter','none');

            ax2 = nexttile;
            imagesc(ax2, x, fPlot, S_nonsig);
            set(ax2,'YDir','normal','YScale','log');
            xlim(ax2, xlimWanted); ylim(ax2, fWin);
            yticks(ax2, yt); yticklabels(ax2, ytLbl);
            ax2.CLim = clim;
            xlabel(ax2,'Time (ms)'); ylabel(ax2,'Frequency (Hz)');
            title(ax2, sprintf('%s | IED non-significant', epochName), 'Interpreter','none');

            ax3 = nexttile;
            imagesc(ax3, x, fPlot, S_nonIED);
            set(ax3,'YDir','normal','YScale','log');
            xlim(ax3, xlimWanted); ylim(ax3, fWin);
            yticks(ax3, yt); yticklabels(ax3, ytLbl);
            ax3.CLim = clim;
            xlabel(ax3,'Time (ms)'); ylabel(ax3,'Frequency (Hz)');
            title(ax3, sprintf('%s | non-IED', epochName), 'Interpreter','none');

            if row == 1
                cb = colorbar(ax3, 'eastoutside');
                cb.Layout.Tile = 'east';
                cb.Label.String = 'power (dB)';
            end
        end

        sgtitle(sprintf('%s | %s', ptID, areaName), 'Interpreter','none');

        fnameBase = sprintf('%s_%s_spectrograms_ALLepochs', ptID, areaSafe);

        set(fig,'Renderer','painters');
        exportgraphics(fig, fullfile(outputFolderName, [fnameBase '.pdf']), ...
            'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);

        close(fig);

    end

end

function Smean = local_avg_spec_one_epoch(SpecCell, chanIdxPerArea)
% SpecCell is a cell matrix: (channels x trials)
% Each entry is (freq x time) LINEAR power (abs(wave).^2)

    Smean = [];
    if isempty(chanIdxPerArea) || isempty(SpecCell)
        return;
    end

    mats = {};
    for ii = 1:numel(chanIdxPerArea)
        r = chanIdxPerArea(ii);
        if r < 1 || r > size(SpecCell,1)
            continue;
        end
        for t = 1:size(SpecCell,2)
            if ~isempty(SpecCell{r,t})
                mats{end+1} = SpecCell{r,t};
            end
        end
    end

    if isempty(mats)
        return;
    end

    S3 = cat(3, mats{:});
    Smean = mean(S3, 3, 'omitnan');
end
