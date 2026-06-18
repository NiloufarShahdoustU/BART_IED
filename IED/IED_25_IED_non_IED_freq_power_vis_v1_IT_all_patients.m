% ccn 2026
% group plot: ONE FIGURE + ONE PDF for ALL patients (IT)
% per patient: compute power(dB) vs frequency(Hz) curves per epoch/area/condition
% group: mean across patients of patient means, SEM across patients (NOT sem-of-sems)
% AUTHOR: Nill (edited)

clear; clc; close all;
warning('off','all');

SigChanFolder_IT = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_ITs_2_chunks';
matFile_IT  = dir(fullfile(SigChanFolder_IT, 'sigChanNums_ITs_PostResponse_PreOutcome.mat'));
filePath_IT = fullfile(SigChanFolder_IT, matFile_IT(1).name);
SigChans_IT = load(filePath_IT);

inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_25_IED_non_IED_freq_power_vis_v1_IT_all_patients\';

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

xTickStep = 50;   % Hz
xPad      = 5;    % Hz

epochNames = {'PostResponse','PreOutcome'};

f_ref = [];

% group container:
% groupData.(areaSafe).areaName
% groupData.(areaSafe).(epochName).sig.muMat    (freq x PatientsNum)
% groupData.(areaSafe).(epochName).nonsig.muMat (freq x PatientsNum)
% groupData.(areaSafe).(epochName).nonIED.muMat (freq x PatientsNum)
groupData = struct();

for pt = 1:PatientsNum
    

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

    isControl = LFPIED.isControl;

    pre_IT  = SigChans_IT.sigChanNums_PreOutcome{pt,1};
    post_IT = SigChans_IT.sigChanNums_PostResponse{pt,1};

    significantChannelIndices_temp = unique([pre_IT(:); post_IT(:)], 'stable');

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

        if isempty(IEDsig_epoch_idx) || isempty(IEDnonsig_epoch_idx) || isempty(nonIED_epoch_idx)
            disp("  Missing channels for one or more conditions. Skipping epoch: " + epochName);
            continue;
        end

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
            disp("  No common brain areas for this epoch in this patient. Skipping epoch: " + epochName);
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
            disp("  Not enough trials for one or more conditions. Skipping epoch: " + epochName);
            continue;
        end

        Spec_IED_sig    = cell(length(IEDsig_epoch_idx),    length(IEDtrialIdx_sig));
        Spec_IED_nonsig = cell(length(IEDnonsig_epoch_idx), length(IEDtrialIdx_nonsig));
        Spec_nonIED     = cell(length(nonIED_epoch_idx),    length(nonIEDtrialIdx));

        period_ref = [];

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

                [wave, period, ~, ~] = basewaveERP(signal, Fs, loF, hiF, MotherWaveParam, waitc);
                powtemp = abs(wave).^2;
                Spec_IED_sig{c,t} = powtemp(:, midIdx);

                if isempty(period_ref)
                    period_ref = period;
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

        if isempty(f_ref)
            if ~isempty(period_ref)
                f_ref = (1 ./ period_ref);
                f_ref = f_ref(:);
            else
                f_ref = (1:size(Spec_nonIED{1,1},1)).';
            end
        end

        f_ref = f_ref(:);

        for iA = 1:numel(BrainAreasCommon_epoch)

            areaName = BrainAreasCommon_epoch{iA};
            areaSafe = regexprep(areaName, '[^\w\-]+', '_');

            chIdx_sig_local    = idxIEDsig_epoch{iA};
            chIdx_nonsig_local = idxIEDnonsig_epoch{iA};
            chIdx_nonIED_local = idxNonIED_epoch{iA};

            [mu_sig, ~]    = local_power_freq_stats(Spec_IED_sig,    chIdx_sig_local);
            [mu_nonsig, ~] = local_power_freq_stats(Spec_IED_nonsig, chIdx_nonsig_local);
            [mu_nonIED, ~] = local_power_freq_stats(Spec_nonIED,     chIdx_nonIED_local);

            if isempty(mu_sig) || isempty(mu_nonsig) || isempty(mu_nonIED)
                continue;
            end

            mu_sig    = mu_sig(:);
            mu_nonsig = mu_nonsig(:);
            mu_nonIED = mu_nonIED(:);

            if ~isfield(groupData, areaSafe)
                groupData.(areaSafe).areaName = areaName;
            end

            groupData = local_store_patient_curve(groupData, areaSafe, epochName, 'sig',    mu_sig,    f_ref, pt, PatientsNum);
            groupData = local_store_patient_curve(groupData, areaSafe, epochName, 'nonsig', mu_nonsig, f_ref, pt, PatientsNum);
            groupData = local_store_patient_curve(groupData, areaSafe, epochName, 'nonIED', mu_nonIED, f_ref, pt, PatientsNum);

        end
    end
end

areaFields = fieldnames(groupData);
if isempty(areaFields)
    error('No group data found.');
end

nAreas = numel(areaFields);
fPlot = f_ref(:);

col_sig    = [0.0000 0.4470 0.7410];
col_nonsig = [0.8500 0.3250 0.0980];
col_nonIED = [0.4660 0.6740 0.1880];
shadeAlpha = 0.4;

tilePx = 300;
figW = max(1200, tilePx*nAreas);
figH = max(700,  tilePx*2 + 120);
fig = figure('Color','w','Position',[60 60 figW figH], 'Visible','off');

tiledlayout(2, nAreas, 'Padding','compact', 'TileSpacing','compact');

hLeg = gobjects(3,1);
legSet = false;

for e = 1:numel(epochNames)
    epochName = epochNames{e};

    for a = 1:nAreas
        areaSafe = areaFields{a};
        areaName = groupData.(areaSafe).areaName;

        ax = nexttile;
        hold(ax,'on'); box(ax,'off');

        if isfield(groupData.(areaSafe), epochName)

            D = groupData.(areaSafe).(epochName);

            [mu_sig, sem_sig]       = local_group_mean_sem(D.sig.muMat);
            [mu_nonsig, sem_nonsig] = local_group_mean_sem(D.nonsig.muMat);
            [mu_nonIED, sem_nonIED] = local_group_mean_sem(D.nonIED.muMat);

            local_plot_mean_sem(ax, fPlot, mu_sig,    sem_sig,    col_sig,    shadeAlpha);
            local_plot_mean_sem(ax, fPlot, mu_nonsig, sem_nonsig, col_nonsig, shadeAlpha);
            local_plot_mean_sem(ax, fPlot, mu_nonIED, sem_nonIED, col_nonIED, shadeAlpha);

            xt = 0:xTickStep:ceil(max(fPlot)/xTickStep)*xTickStep;
            set(ax,'XTick', xt);
            xlim(ax, [-xPad, max(fPlot)]);

            if ~legSet
                htmp = findobj(ax, 'Type','line');
                htmp = flipud(htmp(:));
                if numel(htmp) >= 3
                    hLeg(1) = htmp(1);
                    hLeg(2) = htmp(2);
                    hLeg(3) = htmp(3);
                    legSet = true;
                end
            end
        else
            text(ax, 0.5, 0.5, 'No data', 'Units','normalized', ...
                'HorizontalAlignment','center', 'Color',[0.4 0.4 0.4]);

            xt = 0:xTickStep:ceil(max(fPlot)/xTickStep)*xTickStep;
            set(ax,'XTick', xt);
            xlim(ax, [-xPad, max(fPlot)]);
        end

        axis(ax,'square');

        if a == 1
            ylabel(ax, sprintf('%s\npower (dB)', epochName), 'Interpreter','none');
        else
            ylabel(ax, '');
        end

        if e == 2
            xlabel(ax, 'frequency (Hz)');
        else
            xlabel(ax, '');
        end

        if e == 1
            title(ax, areaName, 'Interpreter','none');
        end

        set(ax,'FontSize',10);
    end
end

if legSet
    lg = legend(hLeg, {'IED significant','IED non-significant','non-IED'}, ...
        'Location','southoutside', 'Orientation','horizontal');
    lg.Box = 'off';

    txtObjs = findobj(lg, 'Type', 'Text');
    txtObjs = flipud(txtObjs(:));
    for k = 1:min(numel(txtObjs), numel(hLeg))
        if isgraphics(hLeg(k)) && isprop(hLeg(k), 'Color')
            txtObjs(k).Color = hLeg(k).Color;
        end
    end
end

sgtitle('all patients | power vs frequency', 'Interpreter','none');

fnameBase = 'ALLpatients_PowerFreq_2epochs_ALLareas_IT';

set(fig,'Renderer','painters');
exportgraphics(fig, fullfile(outputFolderName, [fnameBase '.pdf']), ...
    'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);

close(fig);

function groupData = local_store_patient_curve(groupData, areaSafe, epochName, condName, muVec, f_ref, pt, PatientsNum)

    if ~isfield(groupData.(areaSafe), epochName)
        groupData.(areaSafe).(epochName) = struct();
    end
    if ~isfield(groupData.(areaSafe).(epochName), condName)
        groupData.(areaSafe).(epochName).(condName).muMat = nan(numel(f_ref), PatientsNum);
    end

    M = groupData.(areaSafe).(epochName).(condName).muMat;

    muVec = muVec(:);
    nmin = min(numel(f_ref), numel(muVec));
    tmp = nan(numel(f_ref),1);
    tmp(1:nmin) = muVec(1:nmin);

    M(:,pt) = tmp;

    groupData.(areaSafe).(epochName).(condName).muMat = M;
end

function [mu_dB, sem_dB] = local_power_freq_stats(SpecCell, chanIdxPerArea)

    mu_dB  = [];
    sem_dB = [];

    if isempty(SpecCell) || isempty(chanIdxPerArea)
        return;
    end

    vecs = {};
    for ii = 1:numel(chanIdxPerArea)
        r = chanIdxPerArea(ii);
        if r < 1 || r > size(SpecCell,1)
            continue;
        end
        for t = 1:size(SpecCell,2)
            M = SpecCell{r,t};
            if ~isempty(M)
                v_lin = mean(M, 2, 'omitnan');
                v_dB  = 10*log10(v_lin + eps);
                if all(isfinite(v_dB))
                    vecs{end+1} = v_dB; %#ok<AGROW>
                end
            end
        end
    end

    if isempty(vecs)
        return;
    end

    V = cat(2, vecs{:});
    mu_dB = mean(V, 2, 'omitnan');

    n = sum(isfinite(V), 2);
    sd = std(V, 0, 2, 'omitnan');
    sem_dB = sd ./ sqrt(max(n,1));
end

function [mu, sem] = local_group_mean_sem(muMat)
    muMat = double(muMat);
    mu = mean(muMat, 2, 'omitnan');
    n = sum(isfinite(muMat), 2);
    sd = std(muMat, 0, 2, 'omitnan');
    sem = sd ./ sqrt(max(n,1));
end

function local_plot_mean_sem(ax, x, mu, sem, col, alphaVal)

    x = x(:); mu = mu(:); sem = sem(:);

    good = isfinite(x) & isfinite(mu) & isfinite(sem);
    x = x(good); mu = mu(good); sem = sem(good);

    if isempty(x)
        return;
    end

    y1 = mu - sem;
    y2 = mu + sem;

    patch(ax, [x; flipud(x)], [y1; flipud(y2)], col, ...
        'FaceAlpha', alphaVal, 'EdgeColor','none');

    plot(ax, x, mu, 'LineWidth', 0.5, 'Color', col);
end
