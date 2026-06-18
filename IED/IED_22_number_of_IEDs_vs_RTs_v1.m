% CCN2026
% # of IEDs vs RTs (per trial, all patients)
% NOTE: EXCLUDE isControl trials
% note: this code uses data that are the output of codes 23!!! keep that in
% mind. 
% VIS: bins = 1 vs 2+
% AUTHOR: Nill

clear; clc; close all;
warning('off','all');

SigChanFolder = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_RTs_2_chunks';

matFile = dir(fullfile(SigChanFolder, 'sigChanNums_RTs_PostOnset_PreResponse.mat'));
filePath = fullfile(SigChanFolder, matFile(1).name);
SigChansNumPatients = load(filePath);

inputFolderName_IEDtrials = ...
    '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';

outputFolderName = ...
    '\\155.100.91.44\d\Code\Nill\BART\IED\IED_22_number_of_IEDs_vs_RTs_v1\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));
PatientsNum = length(fileList);

IEDsNumbersAndRTs = struct();

for pt = 1:PatientsNum

    ptID = erase(fileList(pt).name, '.LFPIED.mat');
    disp("Processing patient: " + ptID);

    load(fullfile(inputFolderName_IEDtrials, fileList(pt).name)); % loads LFPIED

    epochList = {
        LFPIED.IEDtrialsPostOnset;
        LFPIED.IEDtrialsPreResponse;
    };

    [~, nTrials] = size(epochList{1});

    pre  = SigChansNumPatients.sigChanNums_PostOnset{pt,1};
    post = SigChansNumPatients.sigChanNums_PreResponse{pt,1};
    significantChannelIndices_temp = unique([pre(:); post(:)], 'stable');

    [~, significantChannelIndices] = ismember( ...
        significantChannelIndices_temp, ...
        LFPIED.selectedChans ...
    );

    IED_count = zeros(nTrials,1);

    for e = 1:length(epochList)
        IEDtrials = epochList{e};
        IEDtrials_sig = IEDtrials(significantChannelIndices, :);
        IED_count = IED_count + sum(IEDtrials_sig, 1)';
    end

    RTs = LFPIED.RTs(:);

    if isfield(LFPIED,'isControl')
        isControl = logical(LFPIED.isControl(:));
    elseif isfield(LFPIED,'isControlTrial')
        isControl = logical(LFPIED.isControlTrial(:));
    elseif isfield(LFPIED,'isControlTrials')
        isControl = logical(LFPIED.isControlTrials(:));
    else
        error('No isControl field for patient %s', ptID);
    end

    if numel(isControl) ~= nTrials
        error('isControl length mismatch for patient %s', ptID);
    end

    % --------- RT outlier removal (>10s) ---------
    RTsThreshold   = 10;
    OutlierIndices = RTs > RTsThreshold;

    RTs       = RTs(~OutlierIndices);
    IED_count = IED_count(~OutlierIndices);
    isControl = isControl(~OutlierIndices);
    % ---------------------------------------------

    IEDsNumbersAndRTs(pt).ptID      = ptID;
    IEDsNumbersAndRTs(pt).IED_count = IED_count;
    IEDsNumbersAndRTs(pt).RTs       = RTs;
    IEDsNumbersAndRTs(pt).isControl = isControl;
end


%%
minTrialsPerPoint = 5;

outFile = fullfile(outputFolderName, ...
    'ALLPATIENTS_IEDcount_vs_RTs_1_vs_2plus_nonControl.pdf');

xVals     = [1 2];
xTickLbls = {'1','2+'};

allBin = [];
allRT  = [];

A = nan(PatientsNum,2);
nPerBin = zeros(PatientsNum,2);

for pt = 1:PatientsNum

    IED_count = IEDsNumbersAndRTs(pt).IED_count;
    RTs       = IEDsNumbersAndRTs(pt).RTs;
    isControl = IEDsNumbersAndRTs(pt).isControl;

    good = ~isnan(RTs) & ~isControl & (IED_count >= 1);
    IED_count = IED_count(good);
    RTs       = RTs(good);

    IED_bin = IED_count;
    IED_bin(IED_bin >= 2) = 2;

    for b = 1:2
        idx = (IED_bin == b);
        n   = sum(idx);

        if n < minTrialsPerPoint
            continue
        end

        meanRT = mean(RTs(idx));

        allBin = [allBin; b];
        allRT  = [allRT; meanRT];

        A(pt,b) = meanRT;
        nPerBin(pt,b) = n;
    end
end

rowsComplete = all(~isnan(A),2);
A2 = A(rowsComplete,:);

p_perm = nan;

if size(A2,1) >= 3

    nPerm   = 10000;
    obsDiff = mean(A2(:,2) - A2(:,1));

    permDiff = zeros(nPerm,1);
    nPts = size(A2,1);

    for k = 1:nPerm
        flip = rand(nPts,1) > 0.5;
        Aperm = A2;
        Aperm(flip,:) = Aperm(flip,[2 1]);
        permDiff(k) = mean(Aperm(:,2) - Aperm(:,1));
    end

    p_perm = mean(abs(permDiff) >= abs(obsDiff));

    fprintf('\nPermutation RTs 1 vs 2+: p = %.4g (N=%d patients)\n', ...
        p_perm, nPts);
end

pooledN = nansum(nPerBin,1);

hFig = figure('Color','w','Units','inches','Position',[2 2 5.2 4.2]);
hold on;

xJ = allBin + 0.14*(rand(size(allBin))-0.5);
scatter(xJ, allRT, 45, 'filled', ...
    'MarkerFaceAlpha',0.35, ...
    'MarkerEdgeColor','none');

meanRTs = nanmean(A,1);
semRTs  = nanstd(A,[],1) ./ sqrt(sum(~isnan(A),1));

errorbar(xVals, meanRTs, semRTs, '-o', ...
    'LineWidth',2,'MarkerSize',7,...
    'MarkerFaceColor','k','Color','k','CapSize',10);

for i = 1:2
    text(xVals(i), meanRTs(i) + semRTs(i) + 0.03*max(allRT), ...
        sprintf('Trials=%d', pooledN(i)), ...
        'HorizontalAlignment','center','FontSize',10);
end

xlabel('# of IEDs across channels in trial','FontSize',14);
ylabel('RT (s)','FontSize',14);
title('RT vs #IEDs (1 vs 2+)', 'FontSize',14);

xlim([0.5 2.5]);
xticks(xVals);
xticklabels(xTickLbls);

set(gca,'FontSize',12,'TickDir','out');
box off;

if ~isnan(p_perm)
    y = max(allRT) * 1.15;
    plot([1 1 2 2],[y y*1.03 y*1.03 y],'k-','LineWidth',1.5);
    text(1.5, y*1.06, ...
        sprintf('perm p = %.3g', p_perm), ...
        'HorizontalAlignment','center', ...
        'FontSize',12,'FontWeight','bold');
end

exportgraphics(hFig, outFile, 'ContentType','vector');
disp("Saved group plot to: " + outFile);
