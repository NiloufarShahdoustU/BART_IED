% CCN2026
% # of IEDs vs accuracy (per trial, all patients)
% NOTE: EXCLUDE isControl trials
% note: this code uses data that are the output of codes 23!!! keep that in
% mind. 
% VIS: bins = 0, 1, 2, 3+
% AUTHOR: Nill

clear; clc; close all;
warning('off','all');

%%
SigChanFolder = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_acc_2_chunks';

matFile = dir(fullfile(SigChanFolder, 'sigChanNums_acc_PostResponse_PreOutcome.mat'));
filePath = fullfile(SigChanFolder, matFile(1).name);  
SigChansNumPatients = load(filePath);


%%
inputFolderName_IEDtrials = ...
    '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';

outputFolderName = ...
    '\\155.100.91.44\d\Code\Nill\BART\IED\IED_22_number_of_IEDs_vs_acc_v1\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));
PatientsNum = length(fileList);


IEDsNumbersAndOutcomes = struct();

for pt = 1:PatientsNum
% for pt = 1:1

    ptID = erase(fileList(pt).name, '.LFPIED.mat');
    disp("Processing patient: " + ptID);

    load(fullfile(inputFolderName_IEDtrials, fileList(pt).name)); % loads LFPIED

    % ---------------- epoch list ----------------
    epochList = {
        % LFPIED.IEDtrialsPostOnset;
        % LFPIED.IEDtrialsPreResponse;
        LFPIED.IEDtrialsPostResponse;
        LFPIED.IEDtrialsPreOutcome;
    };

    % number of trials
    [~, nTrials] = size(epochList{1});

    pre  = SigChansNumPatients.sigChanNums_PreOutcome{pt,1};
    post = SigChansNumPatients.sigChanNums_PostResponse{pt,1};

    significantChannelIndices_temp = unique([pre(:); post(:)], 'stable');


    [~, significantChannelIndices] = ismember( ...
        significantChannelIndices_temp, ...
        LFPIED.selectedChans ...
    );





% ---------------- IED count (SIGNIFICANT CHANNELS ONLY) ----------------
    IED_count = zeros(nTrials,1);   % trials × 1
    
    for e = 1:length(epochList)
        IEDtrials = epochList{e};   % chans × trials
    
    
        % keep only significant channels
        IEDtrials_sig = IEDtrials(significantChannelIndices, :);
    
        % sum across significant channels
        IED_count = IED_count + sum(IEDtrials_sig, 1)';
    end



    % ---------------- outcome ----------------
    outcome = LFPIED.outcomeType(:);   % trials × 1

    % ---------------- isControl ----------------
    % Expecting LFPIED.isControl as trials×1 logical (or 0/1).
    % If your field name differs, add it here.
    if isfield(LFPIED, 'isControl')
        isControl = logical(LFPIED.isControl(:));
    elseif isfield(LFPIED, 'isControlTrial')
        isControl = logical(LFPIED.isControlTrial(:));
    elseif isfield(LFPIED, 'isControlTrials')
        isControl = logical(LFPIED.isControlTrials(:));
    else
        error('No isControl field found in LFPIED for patient %s. Add the correct field name.', ptID);
    end

    % Safety check
    if numel(isControl) ~= nTrials
        error('isControl length (%d) does not match nTrials (%d) for patient %s.', numel(isControl), nTrials, ptID);
    end

    % ---------------- store ----------------
    IEDsNumbersAndOutcomes(pt).ptID      = ptID;
    IEDsNumbersAndOutcomes(pt).IED_count = IED_count;
    IEDsNumbersAndOutcomes(pt).outcome   = outcome;
    IEDsNumbersAndOutcomes(pt).isControl = isControl;

end



%% ===================== ALL PATIENTS: 1 vs 2+ =====================
minTrialsPerPoint = 5;

outFile = fullfile(outputFolderName, ...
    'ALLPATIENTS_IEDcount_vs_accuracy_1_vs_2plus_nonControl.pdf');

xVals     = [1 2];
xTickLbls = {'1','2+'};

allBin = [];
allAcc = [];

pooledWins = zeros(1,2);
pooledN    = zeros(1,2);

A = nan(PatientsNum,2);
nPerBin = zeros(PatientsNum,2);

for pt = 1:PatientsNum

    IED_count = IEDsNumbersAndOutcomes(pt).IED_count;
    outcome   = IEDsNumbersAndOutcomes(pt).outcome;
    isControl = IEDsNumbersAndOutcomes(pt).isControl;

    good = ~isnan(outcome) & ~isControl & (IED_count >= 1);
    IED_count = IED_count(good);
    outcome   = outcome(good);

    % ---- binning: 1 vs 2+ ----
    IED_bin = IED_count;
    IED_bin(IED_bin >= 2) = 2;

    for b = 1:2
        idx = (IED_bin == b);
        n   = sum(idx);

        if n == 0
            continue
        end

        wins = sum(outcome(idx) == 1);

        pooledWins(b) = pooledWins(b) + wins;
        pooledN(b)    = pooledN(b) + n;
        nPerBin(pt,b) = n;

        if n >= minTrialsPerPoint
            acc = wins / n;
            allBin = [allBin; b];
            allAcc = [allAcc; acc];
            A(pt,b) = acc;
        end
    end
end

% ===================== POOLED CI =====================
pooledAcc  = nan(1,2);
pooledLoCI = nan(1,2);
pooledHiCI = nan(1,2);

for b = 1:2
    if pooledN(b) > 0
        [p, pci] = binofit(pooledWins(b), pooledN(b), 0.05);
        pooledAcc(b)  = p;
        pooledLoCI(b) = pci(1);
        pooledHiCI(b) = pci(2);
    end
end

% ===================== STATS =====================
rowsComplete = all(~isnan(A),2);
A2 = A(rowsComplete,:);

p_wilcoxon = nan;
p_perm     = nan;

if size(A2,1) >= 3

    % ----- Wilcoxon signed-rank -----
    p_wilcoxon = signrank(A2(:,1), A2(:,2));

    % ----- Paired permutation test -----
    nPerm = 10000;
    obsDiff = mean(A2(:,2) - A2(:,1));

    permDiff = zeros(nPerm,1);
    nPts = size(A2,1);

    for k = 1:nPerm
        flip = rand(nPts,1) > 0.5;
        Aperm = A2;
        Aperm(flip,:) = Aperm(flip,[2 1]);   % swap within patient
        permDiff(k) = mean(Aperm(:,2) - Aperm(:,1));
    end

    p_perm = mean(abs(permDiff) >= abs(obsDiff));

    fprintf('\n1 vs 2+: p = %.4g (N=%d patients)\n', ...
        p_wilcoxon, nPts);
    fprintf('Permutation test (10,000 perms): p = %.4g\n', p_perm);

else
    fprintf('\nNot enough complete patients (N=%d)\n', size(A2,1));
end

% ===================== PLOT =====================
hFig = figure('Color','w','Units','inches','Position',[2 2 5.2 4.2]);
hold on;

xJ = allBin + 0.14*(rand(size(allBin))-0.5);
scatter(xJ, allAcc, 45, 'filled', ...
    'MarkerFaceAlpha',0.35, ...
    'MarkerEdgeColor','none');

errLow  = pooledAcc - pooledLoCI;
errHigh = pooledHiCI - pooledAcc;

errorbar(xVals, pooledAcc, errLow, errHigh, '-o', ...
    'LineWidth',2,'MarkerSize',7,...
    'MarkerFaceColor','k','Color','k','CapSize',10);

for i = 1:2
    text(xVals(i), pooledAcc(i)+0.06, ...
        sprintf('Trials=%d', pooledN(i)), ...
        'HorizontalAlignment','center','FontSize',10);
end

xlabel('# of IEDs across channels in trial','FontSize',14);
ylabel('accuracy','FontSize',14);
title('accuracy vs #IEDs (1 vs 2+)', 'FontSize',14);

xlim([0.5 2.5]);
xticks(xVals);
xticklabels(xTickLbls);
ylim([0 1.25]);

set(gca,'FontSize',12,'TickDir','out');
box off;

if ~isnan(p_perm)
    y = 1.12;
    plot([1 1 2 2],[y y+0.02 y+0.02 y],'k-','LineWidth',1.5);
    text(1.5, y+0.035, ...
        sprintf('perm test p = %.3g', p_perm), ...
        'HorizontalAlignment','center', ...
        'FontSize',12,'FontWeight','bold');
end

exportgraphics(hFig, outFile, 'ContentType','vector');
disp("Saved group plot to: " + outFile);

% ===================== UTIL =====================
function out = ternary(cond,a,b)
if cond
    out = a;
else
    out = b;
end
end
