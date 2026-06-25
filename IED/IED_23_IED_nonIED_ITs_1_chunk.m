% the final version, mix of 5 and 2
% Author: Nill

clear;
clc;
close all;
warning('off','all');

inputFolderName_IEDdata = 'D:\Nill\data\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = 'D:\Nill\code\BART\IED\IED_23_IED_nonIED_ITs_1_chunk\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

fileList = dir(fullfile(inputFolderName_IEDdata, '*.LFPIED.mat'));
PatientsNum = length(fileList);

PostResponse = 1;
alpha = 0.05;

epoch_names = {'Post-Response'};
epoch_colors = [
    0.85 0.33 0.10
];

nonIEDtrials_bhvr_measure_mean = nan(PatientsNum,1);
IEDtrials_bhvr_measure_mean    = nan(PatientsNum,1);

for pt = 1:PatientsNum

    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    disp("patient: " + ptID);

    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    nChans = length(LFPIED.selectedChans)-1;

    RTs = LFPIED.RTs;
    ITs = LFPIED.ITs;
    RTsThreshold = 10;
    OutlierIndices = RTs >= RTsThreshold;

    isControl = LFPIED.isControl;      % 1 = control, 0 = not control
    isControl = isControl(:)';        
    keepIdx   = (~OutlierIndices) & (~isControl);

    RTs = RTs(keepIdx);
    ITs = ITs(keepIdx);

    IEDtrialsPostResponse = LFPIED.IEDtrialsPostResponse(:, keepIdx);

    allTimePoints = IEDtrialsPostResponse;
    nonIEDIndices = find(all(allTimePoints == 0, 1));

    nonIEDtrials_bhvr_measure = ITs(nonIEDIndices);
    nonIEDtrials_bhvr_measure = nonIEDtrials_bhvr_measure(~isnan(nonIEDtrials_bhvr_measure));
    nonIEDtrials_bhvr_measure_mean(pt) = mean(nonIEDtrials_bhvr_measure);

    IEDTrials_bhvr_measure_MeanPerChan_PostResponse = nan(nChans,1);

    IEDtrialsPostResponse_only = IEDtrialsPostResponse;

    pVal_PostResponse = nan(1,nChans);

    NumberofPermutations = 10000;

    for chz = 1:nChans

        IEDTrials_bhvr_measure = IEDtrialsPostResponse_only(chz,:);
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure .* ITs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(~isnan(IEDTrials_bhvr_measure));

        if (numel(IEDTrials_bhvr_measure) > 0)
            pVal_PostResponse(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PostResponse(chz) < alpha
                IEDTrials_bhvr_measure_MeanPerChan_PostResponse(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PostResponse(chz) = NaN;
        end

        clear IEDTrials_bhvr_measure

    end

    IEDtrials_bhvr_measure_mean(pt,PostResponse) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PostResponse);

end

%%

name = 'IEDtrials_bhvr_measure_mean_ITs_PostResponse';
save([outputFolderName name '.mat'],'IEDtrials_bhvr_measure_mean');

% =========================
% Full participant vectors
% =========================
vec1_full = nonIEDtrials_bhvr_measure_mean(:);
vec2_full = IEDtrials_bhvr_measure_mean(:,PostResponse);

% =========================
% Filtered vectors for boxplot
% ONLY keep non-IED datapoints for participants
% who also have valid Post-Response IED IT datapoints
% =========================
validPair = ~isnan(vec1_full) & vec1_full ~= 0 & ...
            ~isnan(vec2_full) & vec2_full ~= 0;

vec1 = vec1_full(validPair);
vec2 = vec2_full(validPair);

allVecs = [vec1; vec2];
group = [ones(length(vec1),1); 2*ones(length(vec2),1)];

figSize = 4;

figure( ...
    'Units','inches', ...
    'Position',[1 1 figSize figSize], ...   
    'Visible','on');

boxplot(allVecs, group, ...
    'Labels', {'non-IED','Post-Response'}, ...
    'Color', 'k', ...
    'Symbol', '');

xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('normalized mean IT', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

set(findobj(gca, 'Type', 'line'), 'LineWidth', 0.5);

hold on;
jitterAmount = 0.1;
markerSize = 18;

% =========================
% Jittered x positions
% =========================
x1_full = nan(size(vec1_full));
x2_full = nan(size(vec2_full));

x1_full(validPair) = 1 + (rand(sum(validPair),1) - 0.5) * jitterAmount;
x2_full(validPair) = 2 + (rand(sum(validPair),1) - 0.5) * jitterAmount;

% =========================
% Connect same participants
% =========================
nParticipants = length(vec1_full);

for i = 1:nParticipants
    xline = [x1_full(i), x2_full(i)];
    yline = [vec1_full(i), vec2_full(i)];

    validLine = ~isnan(xline) & ~isnan(yline) & yline ~= 0;

    if sum(validLine) >= 2
        plot(xline(validLine), yline(validLine), ...
            '-', ...
            'Color', [0.82 0.82 0.82], ...
            'LineWidth', 0.25);
    end
end

scatter(x1_full(validPair), vec1_full(validPair), markerSize, ...
    'filled', ...
    'MarkerFaceColor', [0.5 0.5 0.5], ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'none');

scatter(x2_full(validPair), vec2_full(validPair), markerSize, ...
    'filled', ...
    'MarkerFaceColor', epoch_colors(1,:), ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'none');

hold off;

set(gca, 'box', 'off', 'tickdir', 'out');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 figSize figSize], ...
         'PaperSize', [figSize figSize]);

filename = 'ITs_boxplot_PostResponse';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

pValuePostResponse_ITs_pt = ranksum(vec1,vec2);

outputFileName = fullfile(outputFolderName, 'ITs_pValues_PostResponse.txt');
fileID = fopen(outputFileName, 'w');

fprintf(fileID, 'P-values from Rank Sum Tests:\n');
fprintf(fileID, 'P-value (non-IED vs Post-Response): %.4f\n', pValuePostResponse_ITs_pt);

fclose(fileID);