% ccn 2026
% power box plot (BOXPLOT + SCATTER per trial, band 8-80 Hz)
% 2 rows = epochs, each column = brain area
% each subplot = 3 boxes (IED sig / IED non-sig / non-IED) + scatters
%
% - y-axis top fixed to 100
% - pairwise Wilcoxon rank-sum (Mann–Whitney), exact when small
% - Holm correction within each subplot (3 pairwise tests)
% - significance bars starting at y=80
% - boxplot lines BLACK, linewidth 0.5, no mean lines
% - NO top/right subplot borders (Box OFF enforced after boxplot)
%
% AUTHOR: Nill (edited)

clear; clc; close all;
warning('off','all');

SigChanFolder_ACC = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_acc_2_chunks';
matFile_ACC  = dir(fullfile(SigChanFolder_ACC, 'sigChanNums_acc_PostResponse_PreOutcome.mat'));
filePath_ACC = fullfile(SigChanFolder_ACC, matFile_ACC(1).name);
SigChans_ACC = load(filePath_ACC);

inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_26_IED_non_IED_power_boxplot_v1_acc\';

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

bandLo = 8;
bandHi = 80;

yTopFixed = 100;

boxLW = 0.5;
sigLW = 0.5;

for pt = 1:PatientsNum
% for pt = 1:1

    ptID = erase(fileList(pt).name, '.LFPIED.mat');
    disp("Processing patient: " + ptID);

    load(fullfile(inputFolderName_IEDtrials, fileList(pt).name));

    epochList = {
        LFPIED.IEDtrialsPostResponse;
        LFPIED.IEDtrialsPreOutcome;
    };

    LFPepochList = {
        LFPIED.LFPmatPostResponse;
        LFPIED.LFPmatPreOutcome;
    };

    epochNames = {'PostResponse','PreOutcome'};

    isControl = LFPIED.isControl;

    pre_ACC  = SigChans_ACC.sigChanNums_PreOutcome{pt,1};
    post_ACC = SigChans_ACC.sigChanNums_PostResponse{pt,1};

    significantChannelIndices_temp = unique([pre_ACC(:); post_ACC(:)], 'stable');

    [~, IEDsignificantChannelIndices] = ismember(significantChannelIndices_temp, LFPIED.selectedChans);
    IEDsignificantChannelIndices = sort(IEDsignificantChannelIndices);

    AnatomicalLocs = LFPIED.anatomicalLocs;
    selectedChannels = LFPIED.selectedChans;
    allChannelIndices = 1:length(selectedChannels);

    f_Hz_ref = [];
    t_ms = (0:(numel(midIdx)-1)) / Fs * 1000; %#ok<NASGU>

    areaData = struct();

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
        scale_ref  = []; %#ok<NASGU>
        coi_ref    = []; %#ok<NASGU>

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

                [wave, period, scale, coi] = basewaveERP(signal, Fs, loF, hiF, MotherWaveParam, waitc); %#ok<ASGLU>
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
                f_Hz_ref = (1:size(Spec_nonIED{1,1},1)).';
            end
        end

        fPlot = f_Hz_ref(:);
        freqMask = (fPlot >= bandLo) & (fPlot <= bandHi);

        for iA = 1:numel(BrainAreasCommon_epoch)

            areaName = BrainAreasCommon_epoch{iA};
            areaSafe = regexprep(areaName, '[^\w\-]+', '_');

            chIdx_sig_local    = idxIEDsig_epoch{iA};
            chIdx_nonsig_local = idxIEDnonsig_epoch{iA};
            chIdx_nonIED_local = idxNonIED_epoch{iA};

            pts_sig    = local_trial_bandpower_points(Spec_IED_sig,    chIdx_sig_local,    freqMask);
            pts_nonsig = local_trial_bandpower_points(Spec_IED_nonsig, chIdx_nonsig_local, freqMask);
            pts_nonIED = local_trial_bandpower_points(Spec_nonIED,     chIdx_nonIED_local, freqMask);

            if isempty(pts_sig) || isempty(pts_nonsig) || isempty(pts_nonIED)
                continue;
            end

            if ~isfield(areaData, areaSafe)
                areaData.(areaSafe).areaName = areaName;
            end

            areaData.(areaSafe).(epochName).pts_sig    = pts_sig;
            areaData.(areaSafe).(epochName).pts_nonsig = pts_nonsig;
            areaData.(areaSafe).(epochName).pts_nonIED = pts_nonIED;

        end

    end

    areaFields = fieldnames(areaData);
    if isempty(areaFields)
        continue;
    end

    keep = false(numel(areaFields),1);
    for a = 1:numel(areaFields)
        areaSafe = areaFields{a};
        for e = 1:numel(epochNames)
            if isfield(areaData.(areaSafe), epochNames{e})
                keep(a) = true;
                break;
            end
        end
    end
    areaFields = areaFields(keep);
    nAreas = numel(areaFields);
    if nAreas == 0
        continue;
    end

    tilePx = 340;
    figW = max(1200, tilePx*nAreas);
    figH = max(760,  tilePx*2 + 140);
    fig = figure('Color','w','Position',[60 60 figW figH], 'Visible','off');

    tiledlayout(2, nAreas, 'Padding','compact', 'TileSpacing','compact');

    groupNames = {'IED sig','IED non-sig','non-IED'};

    for e = 1:numel(epochNames)
        epochName = epochNames{e};

        for a = 1:nAreas
            areaSafe = areaFields{a};
            areaName = areaData.(areaSafe).areaName;

            ax = nexttile;
            hold(ax,'on');
            set(ax,'TickDir','out');
            set(ax,'Box','off'); % (will be forced again AFTER boxplot)

            if isfield(areaData.(areaSafe), epochName)

                D = areaData.(areaSafe).(epochName);

                pts1 = D.pts_sig(:);
                pts2 = D.pts_nonsig(:);
                pts3 = D.pts_nonIED(:);

                y = [pts1; pts2; pts3];
                g = [ones(size(pts1)); 2*ones(size(pts2)); 3*ones(size(pts3))];

                % --- boxplot (boxplot turns Box ON, so we turn it OFF right after) ---
                boxplot(ax, y, g, ...
                    'Labels', groupNames, ...
                    'Symbol','', ...
                    'Whisker', 1.5, ...
                    'Colors', 'k');

                local_force_boxplot_black(ax, boxLW);

                % FORCE: remove top/right borders AFTER boxplot
                set(ax,'Box','off');

                jitter = 0.18;
                x1 = 1 + (rand(size(pts1))-0.5)*2*jitter;
                x2 = 2 + (rand(size(pts2))-0.5)*2*jitter;
                x3 = 3 + (rand(size(pts3))-0.5)*2*jitter;

                scatter(ax, x1, pts1, 12, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.0);
                scatter(ax, x2, pts2, 12, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.0);
                scatter(ax, x3, pts3, 12, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.0);

                xlim(ax, [0.4 3.6]);
                set(ax,'XTick',1:3,'XTickLabel',groupNames);
                ax.XTickLabelRotation = 25;

                yAll = y(isfinite(y));
                if isempty(yAll)
                    yLow = 0;
                else
                    yLow = min(yAll);
                    pad = 0.08 * (yTopFixed - yLow);
                    if ~isfinite(pad) || pad<=0, pad = 1; end
                    yLow = yLow - pad;
                end
                ylim(ax, [yLow, yTopFixed]);

                [p12, ok12] = local_ranksum_smallN(pts1, pts2);
                [p13, ok13] = local_ranksum_smallN(pts1, pts3);
                [p23, ok23] = local_ranksum_smallN(pts2, pts3);

                pRaw = [p12; p13; p23];
                ok   = [ok12; ok13; ok23];

                pAdj = nan(3,1);
                if any(ok)
                    pAdj(ok) = local_holm_bonferroni(pRaw(ok));
                end

                pairs = [1 2; 1 3; 2 3];
                yBase = 80;
                yStep = 4;
                yLevels = [yBase; yBase + yStep; yBase + 2*yStep];

                for k = 1:3
                    if ok(k) && isfinite(pAdj(k))
                        lab = local_p_to_stars_or_ns(pAdj(k));
                    else
                        lab = 'n/a';
                    end
                    local_sig_bar(ax, pairs(k,1), pairs(k,2), yLevels(k), lab, sigLW);
                end

                % FORCE AGAIN (in case anything toggled it)
                set(ax,'Box','off');

            else
                text(ax, 0.5, 0.5, 'No data', 'Units','normalized', ...
                    'HorizontalAlignment','center', 'Color',[0.4 0.4 0.4]);
                xlim(ax, [0 1]); ylim(ax, [0 1]);
                set(ax,'XTick',[],'YTick',[]);
                set(ax,'Box','off');
            end

            axis(ax,'square');
            set(ax,'FontSize',10);

            if a == 1
                ylabel(ax, sprintf('%s\nMean band power (%d–%d Hz) [dB]', epochName, bandLo, bandHi), 'Interpreter','none');
            end
            if e == 1
                title(ax, areaName, 'Interpreter','none');
            end

            % FINAL FORCE OFF (guarantee no top/right)
            set(ax,'Box','off');
        end
    end

    sgtitle(sprintf('%s | Band power boxplots', ptID), 'Interpreter','none');

    fnameBase = sprintf('%s_BoxplotBandPower_%dto%dHz_2epochs_ALLareas_yTop%d_stats', ...
        ptID, bandLo, bandHi, round(yTopFixed));

    set(fig,'Renderer','painters');
    exportgraphics(fig, fullfile(outputFolderName, [fnameBase '.pdf']), ...
        'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);

    close(fig);

end

%% ===================== HELPERS =====================

function pts_dB = local_trial_bandpower_points(SpecCell, chanIdxPerArea, freqMask)
    pts_dB = [];
    if isempty(SpecCell) || isempty(chanIdxPerArea)
        return;
    end

    nTrials = size(SpecCell, 2);
    trialVals = nan(nTrials, 1);

    for t = 1:nTrials
        valsThisTrial = [];
        for ii = 1:numel(chanIdxPerArea)
            r = chanIdxPerArea(ii);
            if r < 1 || r > size(SpecCell,1), continue; end

            M = SpecCell{r,t};
            if isempty(M), continue; end
            if size(M,1) ~= numel(freqMask), continue; end

            Mb = M(freqMask, :);
            v_lin = mean(Mb(:), 'omitnan');
            v_dB  = 10*log10(v_lin + eps);

            if isfinite(v_dB)
                valsThisTrial(end+1,1) = v_dB; %#ok<AGROW>
            end
        end
        if ~isempty(valsThisTrial)
            trialVals(t) = mean(valsThisTrial, 'omitnan');
        end
    end

    pts_dB = trialVals(isfinite(trialVals));
end

function local_force_boxplot_black(ax, lineW)
    if nargin < 2 || isempty(lineW), lineW = 0.5; end
    tags = {'Box','Median','Whisker','Cap','Outliers'};
    for i = 1:numel(tags)
        h = findobj(ax, 'Tag', tags{i});
        for k = 1:numel(h)
            if isgraphics(h(k))
                if isprop(h(k),'Color'), h(k).Color = 'k'; end
                if isprop(h(k),'MarkerEdgeColor'), h(k).MarkerEdgeColor = 'k'; end
                if isprop(h(k),'MarkerFaceColor'), h(k).MarkerFaceColor = 'k'; end
                if isprop(h(k),'LineWidth'), h(k).LineWidth = lineW; end
            end
        end
    end
end

function [p, ok] = local_ranksum_smallN(x, y)
    x = x(isfinite(x)); y = y(isfinite(y));
    ok = (numel(x) >= 2) && (numel(y) >= 2);
    p = nan; if ~ok, return; end
    try
        if numel(x) <= 10 && numel(y) <= 10
            p = ranksum(x, y, 'method', 'exact');
        else
            p = ranksum(x, y, 'method', 'approximate');
        end
    catch
        p = ranksum(x, y);
    end
end

function pAdj = local_holm_bonferroni(pVals)
    pVals = pVals(:);
    m = numel(pVals);
    [pSort, idx] = sort(pVals, 'ascend');
    pHolm = nan(m,1);
    for i = 1:m
        pHolm(i) = min(1, (m - i + 1) * pSort(i));
    end
    for i = 2:m
        pHolm(i) = max(pHolm(i), pHolm(i-1));
    end
    pAdj = nan(m,1);
    pAdj(idx) = pHolm;
end

function lab = local_p_to_stars_or_ns(p)
    if ~isfinite(p), lab = 'n/a'; return; end
    if p < 0.001
        lab = '***';
    elseif p < 0.01
        lab = '**';
    elseif p < 0.05
        lab = '*';
    else
        lab = 'n.s.';
    end
end

function local_sig_bar(ax, x1, x2, y, label, lineW)
    if nargin < 6 || isempty(lineW), lineW = 0.5; end
    if ~isfinite(y), return; end

    yl = ylim(ax);
    ySpan = yl(2) - yl(1);
    if ~isfinite(ySpan) || ySpan<=0, ySpan = 1; end
    h = 0.015 * ySpan;

    plot(ax, [x1 x1 x2 x2], [y-h y y y-h], 'k-', 'LineWidth', lineW);

    text(ax, mean([x1 x2]), y + 0.01*ySpan, label, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'Color','k', 'FontSize',10, 'FontWeight','normal');

    % extra safety
    set(ax,'Box','off');
end
