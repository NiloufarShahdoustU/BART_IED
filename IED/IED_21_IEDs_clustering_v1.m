% CCN2026
% IED clustering using t-SNE (channel × trial based)
% AUTHOR: Nill

clear; clc; close all;
warning('off','all');

preSamp  = 150;
postSamp = 150;
winLen   = preSamp + postSamp + 1;
nPCs = 10;               % PCA components
nNeighbors = 20;         % UMAP neighborhood size
minDist = 0.1;           % UMAP cluster tightness
minClusterSize = 30;     % HDBSCAN parameter
maxClusters = 2;
epsilon = 0.8;   

inputFolderName_IEDtrials = ...
    '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';

outputFolderName = ...
    '\\155.100.91.44\d\Code\Nill\BART\IED\IED_21_IEDs_clustering_v1\';

fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));
PatientsNum = length(fileList);

IED_snippets = [];   % [nIEDs × winLen]

for pt = 1:PatientsNum
% for pt = 1:2   % debug

    ptID = erase(fileList(pt).name, '.LFPIED.mat');
    disp("Processing patient: " + ptID);

    load(fullfile(inputFolderName_IEDtrials, fileList(pt).name)); % LFPIED

    epochList = {
        LFPIED.IEDtrialsPostOnset,    LFPIED.IED_timepointsPostOnset,    LFPIED.LFPmatPostOnset;
        LFPIED.IEDtrialsPreResponse,  LFPIED.IED_timepointsPreResponse,  LFPIED.LFPmatPreResponse;
        LFPIED.IEDtrialsPostResponse, LFPIED.IED_timepointsPostResponse, LFPIED.LFPmatPostResponse;
        LFPIED.IEDtrialsPreOutcome,   LFPIED.IED_timepointsPreOutcome,   LFPIED.LFPmatPreOutcome;
    };

    for e = 1:size(epochList,1)

        IEDtrials    = epochList{e,1};   % chans × trials
        IED_tp       = epochList{e,2};   % chans × time × trials
        LFPmat       = epochList{e,3};   % chans × time × trials

        if isempty(IEDtrials)
            continue
        end

        [nCh, nTime, nTrials] = size(LFPmat);

        [chIdx, trIdx] = find(IEDtrials == 1);

        if isempty(chIdx)
            continue
        end

        for k = 1:length(chIdx)

            ch = chIdx(k);
            tr = trIdx(k);

            % find IED timepoint for THIS channel & trial
            t0 = find(IED_tp(ch,:,tr), 1, 'first');

            if isempty(t0)
                continue
            end

            if t0-preSamp < 1 || t0+postSamp > nTime
                continue
            end

            snippet = squeeze( ...
                LFPmat(ch, t0-preSamp:t0+postSamp, tr));

            % normalize per waveform
            snippet = zscore(snippet);

            IED_snippets = [IED_snippets; snippet(:)'];
        end
    end
end

fprintf('Total IED waveforms collected: %d\n', size(IED_snippets,1));

%% ===================== t-SNE + k-means clustering =====================
rng(1);

Y = tsne(IED_snippets, ...
    'NumDimensions', 2, ...
    'Perplexity', 40, ...
    'Standardize', false);

k = 2;   % HARD CAP: cannot exceed 2 clusters

idx = kmeans(IED_snippets, k, ...
    'Replicates', 20, ...
    'MaxIter', 1000);

cluster_colors = [
    0.1216 0.4667 0.7059;   % Cluster 1 - blue
    1.0000 0.4980 0.0549;   % Cluster 2 - orange
];


hFig1 = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[2 2 5.5 5.5]);

hold on;

centroids = nan(k,2);

for c = 1:k
    % --- scatter points ---
    scatter( ...
        Y(idx==c,1), ...
        Y(idx==c,2), ...
        12, ...
        cluster_colors(c,:), ...
        'filled', ...
        'MarkerFaceAlpha',0.5);

    % --- centroid ---
    centroids(c,:) = mean(Y(idx==c,:), 1, 'omitnan');
end

% --- plot centroids ---
scatter( ...
    centroids(:,1), ...
    centroids(:,2), ...
    140, ...
    cluster_colors, ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.2);

xlabel('t-SNE 1');
ylabel('t-SNE 2');
title('2D t-SNE of IED LFP waveforms');

axis equal;
box off;
set(gca,'TickDir','out','FontSize',12);


% ---- save ----
outFile1 = fullfile(outputFolderName, 'IED_tSNE_2D_clusters.pdf');
exportgraphics(hFig1, outFile1, 'ContentType','vector');


hFig2 = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[2 2 5 5]);

hold on;

for c = 1:k

    % --- data for this cluster ---
    data_c = IED_snippets(idx == c, :);   % trials × time

    % --- mean & SEM ---
    mu  = mean(data_c, 1, 'omitnan');
    sem = std(data_c, 0, 1, 'omitnan');

    t = 1:length(mu);

    % --- SEM shaded area ---
    fill([t fliplr(t)], ...
         [mu+sem fliplr(mu-sem)], ...
         cluster_colors(c,:), ...
         'FaceAlpha',0.25, ...
         'EdgeColor','none');

    % --- mean curve ---
    plot(t, mu, ...
        'LineWidth', 2, ...
        'Color', cluster_colors(c,:));
end

xlabel('Time (samples)');
ylabel('Amplitude (z)');
title('Mean IED waveform per cluster');

axis square;
box off;
set(gca,'TickDir','out','FontSize',12);

% ---- save ----
outFile2 = fullfile(outputFolderName, 'IED_tSNE_2D_clusters_mean_waveforms.pdf');
exportgraphics(hFig2, outFile2, 'ContentType','vector');


%% ===================== UMAP + DBSCAN =====================

[umapXYZ, ~, ~] = run_umap( ...
    IED_snippets, ...
    'n_neighbors', nNeighbors, ...
    'min_dist', minDist, ...
    'n_components', 2, ...
    'metric', 'euclidean', ...
    'randomize', false, ...
    'verbose', 'none');

labels_raw = dbscan( ...
    umapXYZ, ...
    epsilon, ...
    minClusterSize);

labels = -1 * ones(size(labels_raw));   % initialize as noise

validClusters = unique(labels_raw(labels_raw > 0));
clusterSizes  = arrayfun(@(c) sum(labels_raw == c), validClusters);

[~, order] = sort(clusterSizes, 'descend');

nKeep = min(2, length(order));   % HARD CAP: max 2 clusters
keptClusters = validClusters(order(1:nKeep));

for i = 1:length(keptClusters)
    labels(labels_raw == keptClusters(i)) = i;
end

K = length(keptClusters);
disp("Clusters kept: " + K);

cluster_colors = lines(K);

hFig1 = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[2 2 6 6]);

hold on;

for c = 1:K
    scatter( ...
        umapXYZ(labels==c,1), ...
        umapXYZ(labels==c,2), ...
        12, ...
        cluster_colors(c,:), ...
        'filled', ...
        'MarkerFaceAlpha',0.5);
end

xlabel('UMAP 1');
ylabel('UMAP 2');
title('clustering of IED waveforms');

axis equal;
box off;
set(gca,'TickDir','out','FontSize',12);

outFile1 = fullfile(outputFolderName,'IED_UMAP_2D_clusters.pdf');
exportgraphics(hFig1,outFile1,'ContentType','vector');

hFig2 = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[2 2 5 5]);

hold on;

t = -preSamp:postSamp;

for c = 1:K

    data_c = IED_snippets(labels==c,:);

    mu  = mean(data_c,1,'omitnan');
    sem = std(data_c,0,1,'omitnan');

    fill([t fliplr(t)], ...
         [mu+sem fliplr(mu-sem)], ...
         cluster_colors(c,:), ...
         'FaceAlpha',0.25, ...
         'EdgeColor','none');

    plot(t, mu, ...
        'LineWidth',2, ...
        'Color',cluster_colors(c,:));
end

xlabel('time (samples)');
ylabel('amplitude (z)');
title('mean IED waveform per cluster');

axis square;
box off;
set(gca,'TickDir','out','FontSize',12);

outFile2 = fullfile(outputFolderName,'IED_UMAP_2D_clusters_mean_waveforms.pdf');
exportgraphics(hFig2,outFile2,'ContentType','vector');
