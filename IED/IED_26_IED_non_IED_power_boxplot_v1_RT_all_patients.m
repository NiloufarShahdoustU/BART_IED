% ccn 2026
% RT: power box plot (ALL PATIENTS, BASELINE-CORRECTED)
% 2 rows = epochs, each column = brain area
%
% Visualization:
%   - 2 boxplots per subplot:
%       (1) IED sig  - nonIED  (within-patient delta)
%       (2) IED non-sig - nonIED (within-patient delta)
%   - each scatter dot = one patient's mean delta
%
% Stats per subplot (paired across patients), then Holm within subplot (3 tests):
%   (1) deltaSig    vs 0
%   (2) deltaNonSig vs 0
%   (3) deltaSig    vs deltaNonSig
%
% - y-axis top fixed to 100
% - boxplot lines BLACK, linewidth 0.5, no mean lines
% - NO top/right subplot borders (Box OFF enforced after boxplot)
%
% AUTHOR: Nill (edited)

clear; clc; close all;
warning('off','all');

SigChanFolder_RT = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_RTs_2_chunks';
matFile_RT  = dir(fullfile(SigChanFolder_RT, 'sigChanNums_RTs_PostOnset_PreResponse.mat'));
filePath_RT = fullfile(SigChanFolder_RT, matFile_RT(1).name);
SigChans_RT = load(filePath_RT);

inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_26_IED_non_IED_power_boxplot_v1_RT_all_patients\';

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

% ---- band for boxplots ----
bandLo = 1;
bandHi = 150;

% ---- plotting/stats params ----
yTopFixed = 100;
boxLW = 0.5;
sigLW = 0.5;

epochNames = {'PostOnset','PreResponse'};
groupNames2 = {'IED sig - nonIED','IED non-sig - nonIED'};

% ===================== COLLECT ACROSS PATIENTS =====================
areaDataAll = struct();  % areaDataAll.(areaSafe).(epochName) has ptIdx, sigMeans, nonsigMeans, nonIEDMeans
f_Hz_ref = [];

for pt = 1:PatientsNum
% for pt = 1:5

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

    isControl = LFPIED.isControl;

    pre_RT  = SigChans_RT.sigChanNums_PostOnset{pt,1};
    post_RT = SigChans_RT.sigChanNums_PreResponse{pt,1};

    significantChannelIndices_temp = unique([pre_RT(:); post_RT(:)], 'stable');

    [~, IEDsignificantChannelIndices] = ismember(significantChannelIndices_temp, LFPIED.selectedChans);
    IEDsignificantChannelIndices = sort(IEDsignificantChannelIndices);

    AnatomicalLocs = LFPIED.anatomicalLocs;
    selectedChannels = LFPIED.selectedChans;
    allChannelIndices = 1:length(selectedChannels);

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

        IEDtrialIdx_sig    = find(~isControl & any(IEDtrials_epoch(IEDsig_epoch_idx, :) ~= 0, 1));
        IEDtrialIdx_nonsig = find(~isControl & any(IEDtrials_epoch(IEDnonsig_epoch_idx, :) ~= 0, 1));
        nonIEDtrialIdx     = find(~isControl & all(IEDtrials_epoch(:, :) == 0, 1));

        if isempty(IEDtrialIdx_sig) || isempty(IEDtrialIdx_nonsig) || isempty(nonIEDtrialIdx)
            disp("  Not enough trials for one or more conditions in epoch. Skipping epoch: " + epochName);
            continue;
        end

        Spec_IED_sig    = cell(length(IEDsig_epoch_idx),    length(IEDtrialIdx_sig));
        Spec_IED_nonsig = cell(length(IEDnonsig_epoch_idx), length(IEDtrialIdx_nonsig));
        Spec_nonIED     = cell(length(nonIED_epoch_idx),    length(nonIEDtrialIdx));

        period_ref_local = [];

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

                [wave, period] = basewaveERP(signal, Fs, loF, hiF, MotherWaveParam, waitc); %#ok<ASGLU>
                powtemp = abs(wave).^2;
                Spec_IED_sig{c,t} = powtemp(:, midIdx);

                if isempty(period_ref_local)
                    period_ref_local = period;
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

                [wave, ~] = basewaveERP(signal, Fs, loF, hiF, MotherWaveParam, waitc);
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

                [wave, ~] = basewaveERP(signal, Fs, loF, hiF, MotherWaveParam, waitc);
                powtemp = abs(wave).^2;
                Spec_nonIED{c,t} = powtemp(:, midIdx);
            end
        end

        if isempty(f_Hz_ref)
            if ~isempty(period_ref_local)
                f_Hz_ref = 1 ./ period_ref_local;
            else
                f_Hz_ref = (1:size(Spec_nonIED{1,1},1)).';
            end
        end

        fPlot = f_Hz_ref(:);
        freqMask = (fPlot >= bandLo) & (fPlot <= bandHi);

        % ---- for each common area: compute patient mean per condition (mean of trial points) ----
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

            % patient mean (each patient contributes ONE dot per condition)
            m_sig    = mean(pts_sig,    'omitnan');
            m_nonsig = mean(pts_nonsig, 'omitnan');
            m_nonIED = mean(pts_nonIED, 'omitnan');

            if ~isfinite(m_sig) || ~isfinite(m_nonsig) || ~isfinite(m_nonIED)
                continue;
            end

            if ~isfield(areaDataAll, areaSafe)
                areaDataAll.(areaSafe).areaName = areaName;
            end
            if ~isfield(areaDataAll.(areaSafe), epochName)
                areaDataAll.(areaSafe).(epochName).ptIdx       = [];
                areaDataAll.(areaSafe).(epochName).sigMeans    = [];
                areaDataAll.(areaSafe).(epochName).nonsigMeans = [];
                areaDataAll.(areaSafe).(epochName).nonIEDMeans = [];
            end

            areaDataAll.(areaSafe).(epochName).ptIdx(end+1,1)       = pt;       %#ok<AGROW>
            areaDataAll.(areaSafe).(epochName).sigMeans(end+1,1)    = m_sig;    %#ok<AGROW>
            areaDataAll.(areaSafe).(epochName).nonsigMeans(end+1,1) = m_nonsig; %#ok<AGROW>
            areaDataAll.(areaSafe).(epochName).nonIEDMeans(end+1,1) = m_nonIED; %#ok<AGROW>
        end
    end
end



%%
% ===================== PLOT ONE FIGURE (ALL PATIENTS) =====================
areaFields = fieldnames(areaDataAll);
if isempty(areaFields)
    disp("No areas collected across patients.");
    return;
end

% keep only areas that have at least one epoch with any data
keep = false(numel(areaFields),1);
for a = 1:numel(areaFields)
    areaSafe = areaFields{a};
    for e = 1:numel(epochNames)
        if isfield(areaDataAll.(areaSafe), epochNames{e})
            keep(a) = true;
            break;
        end
    end
end
areaFields = areaFields(keep);
nAreas = numel(areaFields);
if nAreas == 0
    disp("No areas with data after filtering.");
    return;
end

tilePx = 340;
figW = max(1200, tilePx*nAreas);
figH = max(760,  tilePx*2 + 140);
fig = figure('Color','w','Position',[60 60 figW figH], 'Visible','off');

tiledlayout(2, nAreas, 'Padding','compact', 'TileSpacing','compact');

for e = 1:numel(epochNames)
    epochName = epochNames{e};

    for a = 1:nAreas
        areaSafe = areaFields{a};
        areaName = areaDataAll.(areaSafe).areaName;

        ax = nexttile;
        hold(ax,'on');
        set(ax,'TickDir','out');
        set(ax,'Box','off');

        if isfield(areaDataAll.(areaSafe), epochName)

            D = areaDataAll.(areaSafe).(epochName);

            % ---- within-patient deltas (paired) ----
            % FIX: force column vectors (avoids size mismatch in scatter & signrank)
            [dSig, dNonSig] = local_within_patient_deltas( ...
                D.ptIdx(:), D.sigMeans(:), D.nonsigMeans(:), D.nonIEDMeans(:));

            okStats = (numel(dSig) >= 2) && (numel(dNonSig) >= 2);

            % ---- 2-group boxplot ----
            y = [dSig(:); dNonSig(:)];
            g = [ones(size(dSig(:))); 2*ones(size(dNonSig(:)))];


            boxplot(ax, y, g, ...
                'Labels', groupNames2, ...
                'Symbol','', ...
                'Whisker', 1.5, ...
                'Colors', 'k');

            local_force_boxplot_black(ax, boxLW);
            set(ax,'Box','off');

            % FIX: jitter vectors match dSig(:)/dNonSig(:)
            jitter = 0.18;
            dSigV    = dSig(:);
            dNonSigV = dNonSig(:);
            x1 = 1 + (rand(size(dSigV))-0.5)*2*jitter;
            x2 = 2 + (rand(size(dNonSigV))-0.5)*2*jitter;

            scatter(ax, x1, dSigV,    12, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.0);
            scatter(ax, x2, dNonSigV, 12, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.0);

            xlim(ax, [0.4 2.6]);
            set(ax,'XTick',1:2,'XTickLabel',groupNames2);
            ax.XTickLabelRotation = 25;

            % ---- CHANGED: FIXED Y LIMITS ----
            ylim(ax, [-20 40]);

            % ---- tests (paired), Holm within subplot ----
            % FIX: use vector versions for stats
            pRaw = nan(3,1);
            ok   = false(3,1);

            if okStats
                try
                    pRaw(1) = signrank(dSigV, 0);            % sig vs nonIED
                    pRaw(2) = signrank(dNonSigV, 0);         % nonsig vs nonIED
                    pRaw(3) = signrank(dSigV, dNonSigV);     % sig vs nonsig
                    ok = isfinite(pRaw);
                catch
                    ok = false(3,1);
                end
            end

            pAdj = nan(3,1);
            if any(ok)
                pAdj(ok) = local_holm_bonferroni(pRaw(ok));
            end

            % ---- labels above boxes for vs nonIED ----
            yl = ylim(ax);
            ySpan = yl(2) - yl(1);
            if ~isfinite(ySpan) || ySpan<=0, ySpan = 1; end

            yTop1 = local_box_top_y(dSigV);
            yTop2 = local_box_top_y(dNonSigV);

            yTextPad = 0.03 * ySpan;
            if ~isfinite(yTextPad) || yTextPad<=0, yTextPad = 1; end

            lab1 = 'n/a'; lab2 = 'n/a';
            if ok(1) && isfinite(pAdj(1)), lab1 = local_p_to_stars_or_ns(pAdj(1)); end
            if ok(2) && isfinite(pAdj(2)), lab2 = local_p_to_stars_or_ns(pAdj(2)); end

            text(ax, 1, yTop1 + yTextPad, lab1, 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', 'Color','k', 'FontSize',10);
            text(ax, 2, yTop2 + yTextPad, lab2, 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', 'Color','k', 'FontSize',10);

            % ---- significance bar between boxes for sig vs nonsig ----
            yBar = max([yTop1, yTop2]) + 0.10*ySpan;
            if ~isfinite(yBar), yBar = 80; end
            if yBar > yl(2) - 0.02*ySpan
                yBar = yl(2) - 0.06*ySpan;
            end

            labBar = 'n/a';
            if ok(3) && isfinite(pAdj(3)), labBar = local_p_to_stars_or_ns(pAdj(3)); end
            local_sig_bar(ax, 1, 2, yBar, labBar, sigLW);

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
            ylabel(ax, sprintf('%s\nRT Band power Δ from non-IED\n(%d–%d Hz) [dB]', ...
                epochName, bandLo, bandHi), 'Interpreter','none');
        end
        if e == 1
            title(ax, areaName, 'Interpreter','none');
        end

        set(ax,'Box','off');
    end
end

sgtitle(sprintf('ALL PATIENTS | RT Band power (baseline-corrected within patient)'), 'Interpreter','none');

fnameBase = sprintf('ALLPATIENTS_RT_BoxplotBandPower_DELTAfromNonIED_%dto%dHz_2epochs_ALLareas_yTop%d_signrank_Holm', ...
    bandLo, bandHi, round(yTopFixed));

set(fig,'Renderer','painters');
exportgraphics(fig, fullfile(outputFolderName, [fnameBase '.pdf']), ...
    'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);

close(fig);


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

function [dSig, dNonSig] = local_within_patient_deltas(ptIdx, mSig, mNonSig, mNonIED)
    ptIdx = ptIdx(:);
    mSig = mSig(:); mNonSig = mNonSig(:); mNonIED = mNonIED(:);

    ok = isfinite(ptIdx) & isfinite(mSig) & isfinite(mNonSig) & isfinite(mNonIED);
    ptIdx = ptIdx(ok);
    mSig = mSig(ok); mNonSig = mNonSig(ok); mNonIED = mNonIED(ok);

    u = unique(ptIdx, 'stable');

    dSig = nan(numel(u),1);
    dNonSig = nan(numel(u),1);

    for i = 1:numel(u)
        p = u(i);
        r = (ptIdx == p);

        s  = mean(mSig(r),    'omitnan');
        ns = mean(mNonSig(r), 'omitnan');
        b  = mean(mNonIED(r), 'omitnan');

        if ~isfinite(s) || ~isfinite(ns) || ~isfinite(b)
            continue;
        end

        dSig(i)    = s  - b;
        dNonSig(i) = ns - b;
    end

    keep = isfinite(dSig) & isfinite(dNonSig);
    dSig = dSig(keep);
    dNonSig = dNonSig(keep);
end

function yTop = local_box_top_y(x)
    x = x(isfinite(x));
    if isempty(x)
        yTop = 0;
        return;
    end
    q1 = prctile(x,25);
    q3 = prctile(x,75);
    iqrV = q3 - q1;
    if ~isfinite(iqrV), iqrV = 0; end
    yTop = q3 + 1.0*iqrV;
    yTop = max(yTop, max(x));
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

    set(ax,'Box','off');
end

