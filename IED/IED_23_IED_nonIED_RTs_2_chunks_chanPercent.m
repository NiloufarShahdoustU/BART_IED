% ===================== CCN 2026 =====================
% Author: Nill
% PURPOSE: IED vs non-IED RT analysis (ONLY PostOnset + PreResponse)
% NOTE: PreOnset1 & PreOnset2 REMOVED COMPLETELY

clear;
clc;
close all;
warning('off','all');

%% ===================== paths =====================
inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName        = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_RTs_2_chunks\';

fileList    = dir(fullfile(inputFolderName_IEDdata, '*.LFPIED.mat'));
PatientsNum = length(fileList);

%% ===================== epochs (ONLY 2) =====================
PostOnset   = 1;
PreResponse = 2;

epoch_names  = {'post-onset','pre-response'};
epoch_colors = [
    0.20 0.45 0.70
    0.30 0.70 0.40
];

alpha = 0.05;

%% ===================== outputs =====================
% means are over trials (per patient)
nonIEDtrials_bhvr_measure_mean = nan(PatientsNum,1);
IEDtrials_bhvr_measure_mean    = nan(PatientsNum, 2);  % [PostOnset, PreResponse]

%% ===================== anatomical bookkeeping =====================
AnatoicalLocsNums = 150; % arbitrary number

AnatomicalLocsPatientsPostOnset   = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsPreResponse = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsAll         = zeros(AnatoicalLocsNums, PatientsNum);

nan_cell_array = repmat({'nan'},AnatoicalLocsNums, 1);
AnatomicalLocsVecPostOnset   = string(nan_cell_array);
AnatomicalLocsVecPreResponse = string(nan_cell_array);
AnatomicalLocsVecAll         = string(nan_cell_array);

% -------------------- ADDED: store significant channel numbers per patient (per epoch) --------------------
sigChanNums_PostOnset   = cell(PatientsNum,1);
sigChanNums_PreResponse = cell(PatientsNum,1);
ptID_list               = strings(PatientsNum,1);
% ---------------------------------------------------------------------------------------------------------

%% ===================== main loop =====================
for pt = 1:PatientsNum

    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    disp("patient: " + ptID);

    % -------------------- ADDED --------------------
    ptID_list(pt) = string(ptID);
    % -----------------------------------------------

    % load
    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    nChans = length(LFPIED.selectedChans)-1; % removing the last selected chan

    anatomicalLocs = LFPIED.anatomicalLocs;

    % ===================== COMBINE LEFT/RIGHT =====================
    anatomicalLocs = regexprep(string(anatomicalLocs), '^(Left|Right)\s+', '');
    % =============================================================

    selectedChans = LFPIED.selectedChans;
    selectedChans = selectedChans(1:end-1);

    
    % RT outlier removal
    RTs = LFPIED.RTs;
    
    % ----------- EXCLUDE CONTROL TRIALS (added) -----------
    isNonControl = ~LFPIED.isControl;
    RTs = RTs(isNonControl);
    % ------------------------------------------------------
    
    RTsThreshold   = 10;
    OutlierIndices = RTs >= RTsThreshold;
    RTs            = RTs(~OutlierIndices);
    
    % ONLY 2 epochs used
    IEDtrialsPostOnset   = LFPIED.IEDtrialsPostOnset(:,   isNonControl);
    IEDtrialsPreResponse = LFPIED.IEDtrialsPreResponse(:, isNonControl);
    
    IEDtrialsPostOnset   = IEDtrialsPostOnset(:,   ~OutlierIndices);
    IEDtrialsPreResponse = IEDtrialsPreResponse(:, ~OutlierIndices);
    
    % --------------------- non-IED trials ---------------------
    allTimePoints = IEDtrialsPostOnset + IEDtrialsPreResponse;
    nonIEDIndices = find(all(allTimePoints == 0, 1));
    
    nonIEDtrials_bhvr_measure = RTs(nonIEDIndices);
    nonIEDtrials_bhvr_measure_mean(pt) = mean(nonIEDtrials_bhvr_measure);

    % --------------------- only-in-one-epoch trials ---------------------
    IEDtrialsPostOnset_only   = IEDtrialsPostOnset   .* ~IEDtrialsPreResponse;
    IEDtrialsPreResponse_only = IEDtrialsPreResponse .* ~IEDtrialsPostOnset;

    % per-channel significant mean storage
    IEDTrials_bhvr_measure_MeanPerChan_PostOnset   = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreResponse = nan(nChans,1);

    pVal_PostOnset   = nan(1,nChans);
    pVal_PreResponse = nan(1,nChans);

    NumberofPermutations = 10000;

    for chz = 1:nChans

        % ===================== PostOnset =====================
        IEDTrials_bhvr_measure = IEDtrialsPostOnset_only(chz,:);
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure .* RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);

        if numel(IEDTrials_bhvr_measure) > 0
            pVal_PostOnset(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PostOnset(chz) < alpha
                IEDTrials_bhvr_measure_MeanPerChan_PostOnset(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PostOnset(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure

        % ===================== PreResponse =====================
        IEDTrials_bhvr_measure = IEDtrialsPreResponse_only(chz,:);
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure .* RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);

        if numel(IEDTrials_bhvr_measure) > 0
            pVal_PreResponse(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PreResponse(chz) < alpha
                IEDTrials_bhvr_measure_MeanPerChan_PreResponse(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PreResponse(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure

    end % end channels

    % --------------------- per-patient epoch means ---------------------
    IEDtrials_bhvr_measure_mean(pt,PostOnset)   = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PostOnset);
    IEDtrials_bhvr_measure_mean(pt,PreResponse) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreResponse);

    % --------------------- anatomical counts for % plots ---------------------
    % PostOnset
    ChanIndices = find(~isnan(IEDTrials_bhvr_measure_MeanPerChan_PostOnset));

    % -------------------- ADDED: save significant channel numbers for this patient (PostOnset) --------------------
    sigChanNums_PostOnset{pt} = selectedChans(ChanIndices);
    % -------------------------------------------------------------------------------------------------------------

    timePeriodAnatomicalLoc = anatomicalLocs(selectedChans(ChanIndices));
    for location = 1:length(timePeriodAnatomicalLoc)
        element = timePeriodAnatomicalLoc(location);
        tempIndexInLocs  = ismember(AnatomicalLocsVecPostOnset, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPostOnset(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPostOnset(FoundIndexInLocs,pt) + 1;
        else
            nan_index = find(AnatomicalLocsVecPostOnset == "nan", 1);
            AnatomicalLocsVecPostOnset(nan_index) = element;
            AnatomicalLocsPatientsPostOnset(nan_index,pt) = AnatomicalLocsPatientsPostOnset(nan_index,pt) + 1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear ChanIndices timePeriodAnatomicalLoc

    % PreResponse
    ChanIndices = find(~isnan(IEDTrials_bhvr_measure_MeanPerChan_PreResponse));

    % -------------------- ADDED: save significant channel numbers for this patient (PreResponse) --------------------
    sigChanNums_PreResponse{pt} = selectedChans(ChanIndices);
    % ---------------------------------------------------------------------------------------------------------------

    timePeriodAnatomicalLoc = anatomicalLocs(selectedChans(ChanIndices));
    for location = 1:length(timePeriodAnatomicalLoc)
        element = timePeriodAnatomicalLoc(location);
        tempIndexInLocs  = ismember(AnatomicalLocsVecPreResponse, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPreResponse(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPreResponse(FoundIndexInLocs,pt) + 1;
        else
            nan_index = find(AnatomicalLocsVecPreResponse == "nan", 1);
            AnatomicalLocsVecPreResponse(nan_index) = element;
            AnatomicalLocsPatientsPreResponse(nan_index,pt) = AnatomicalLocsPatientsPreResponse(nan_index,pt) + 1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear ChanIndices timePeriodAnatomicalLoc

    clear IEDtrialsPostOnset IEDtrialsPreResponse
    clear IEDtrialsPostOnset_only IEDtrialsPreResponse_only

end % end patients

%% ===================== build "ALL channels" anatomical reference =====================
for pt = 1:PatientsNum
    % NOTE: anatomicalLocs/selectedChans from last pt are still in workspace;
    % this loop matches your original structure. If you prefer, move this block
    % inside the main loop and store per-pt. Leaving it as-is to minimize changes.

    AnatomicalLocsAll = anatomicalLocs(selectedChans);
    for all = 1:length(AnatomicalLocsAll)
        element = AnatomicalLocsAll(all);
        tempIndexInLocs  = ismember(AnatomicalLocsVecAll, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsAll(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll(FoundIndexInLocs,pt) + 1;
        else
            nan_index = find(AnatomicalLocsVecAll == "nan", 1);
            AnatomicalLocsVecAll(nan_index) = element;
            AnatomicalLocsPatientsAll(nan_index,pt) = AnatomicalLocsPatientsAll(nan_index,pt) + 1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear AnatomicalLocsAll
end

%% ===================== save means =====================
name = 'IEDtrials_bhvr_measure_mean_RTs';
save([outputFolderName name '.mat'],'IEDtrials_bhvr_measure_mean');

% -------------------- ADDED: save significant channel numbers (per patient, per epoch) --------------------
name = 'sigChanNums_RTs_PostOnset_PreResponse';
save([outputFolderName name '.mat'],'sigChanNums_PostOnset','sigChanNums_PreResponse','ptID_list');
% ---------------------------------------------------------------------------------------------------------

%% ===================== preprocessing =====================
IEDtrials_bhvr_measure_mean_PostOnset = IEDtrials_bhvr_measure_mean(:,PostOnset);
IEDtrials_bhvr_measure_mean_PostOnset = IEDtrials_bhvr_measure_mean_PostOnset(~isnan(IEDtrials_bhvr_measure_mean_PostOnset));

IEDtrials_bhvr_measure_mean_PreResponse = IEDtrials_bhvr_measure_mean(:,PreResponse);
IEDtrials_bhvr_measure_mean_PreResponse = IEDtrials_bhvr_measure_mean_PreResponse(~isnan(IEDtrials_bhvr_measure_mean_PreResponse));

%% ===================== visualization: box plot =====================
vec1 = nonIEDtrials_bhvr_measure_mean;
vec2 = IEDtrials_bhvr_measure_mean_PostOnset;
vec3 = IEDtrials_bhvr_measure_mean_PreResponse;

allVecs = [vec1; vec2; vec3];
group   = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1)];

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.12, 0.25], 'Visible', 'on');

boxplot(allVecs, group, ...
    'Labels', {'non-IED', 'PostOnset', 'PreResponse'}, ...
    'Color', 'k', ...
    'Symbol', '');

xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('mean RT', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

set(findobj(gca, 'Type', 'line'), 'LineWidth', 0.5);

% --- overlay jittered points with colors + alpha=0.4 ---
hold on;
jitterAmount = 0.12;
markerSize   = 18;
a = 0.4;

x1 = 1 + (rand(size(vec1)) - 0.5) * jitterAmount;
s1 = scatter(x1, vec1, markerSize, 'filled', ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none');
s1.MarkerFaceAlpha = a;

x2 = 2 + (rand(size(vec2)) - 0.5) * jitterAmount;
s2 = scatter(x2, vec2, markerSize, 'filled', ...
    'MarkerFaceColor', epoch_colors(1,:), 'MarkerEdgeColor', 'none');
s2.MarkerFaceAlpha = a;

x3 = 3 + (rand(size(vec3)) - 0.5) * jitterAmount;
s3 = scatter(x3, vec3, markerSize, 'filled', ...
    'MarkerFaceColor', epoch_colors(2,:), 'MarkerEdgeColor', 'none');
s3.MarkerFaceAlpha = a;

hold off;
set(gca, 'box', 'off', 'tickdir', 'out');

set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
% filename = 'RTs_boxplot';
% saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

%% ===================== p-values across patients =====================
pValuePostOnset_pt   = ranksum(vec1, vec2);
pValuePreResponse_pt = ranksum(vec1, vec3);

outputFileName = fullfile(outputFolderName, 'RTs_pValues.txt');
fileID = fopen(outputFileName, 'w');

fprintf(fileID, 'P-values from Rank Sum Tests:\n');
fprintf(fileID, 'P-value (non-IED vs PostOnset): %.4f\n',   pValuePostOnset_pt);
fprintf(fileID, 'P-value (non-IED vs PreResponse): %.4f\n', pValuePreResponse_pt);

fclose(fileID);

%% ===================== visualization: channel percentage =====================
% cleaning AnatomicalLocsVec
startsWithNaC = startsWith(AnatomicalLocsVecPostOnset, "NaC");
containsLCWM  = contains(AnatomicalLocsVecPostOnset, "Left Cerebral White Matter");
AnatomicalLocsVecPostOnset(startsWithNaC | containsLCWM) = "nan";
clear startsWithNaC containsLCWM

startsWithNaC = startsWith(AnatomicalLocsVecPreResponse, "NaC");
containsLCWM  = contains(AnatomicalLocsVecPreResponse, "Left Cerebral White Matter");
AnatomicalLocsVecPreResponse(startsWithNaC | containsLCWM) = "nan";
clear startsWithNaC containsLCWM

startsWithNaC = startsWith(AnatomicalLocsVecAll, "NaC");
containsLCWM  = contains(AnatomicalLocsVecAll, "Left Cerebral White Matter");
AnatomicalLocsVecAll(startsWithNaC | containsLCWM) = "nan";
clear startsWithNaC containsLCWM

% cleaning AnatomicalLocsPatients based on AnatomicalLocsVec
missingIndices = find(AnatomicalLocsVecPostOnset == "nan");
AnatomicalLocsPatientsPostOnset(missingIndices, :) = [];
AnatomicalLocsVecPostOnset(missingIndices) = [];
clear missingIndices

missingIndices = find(AnatomicalLocsVecPreResponse == "nan");
AnatomicalLocsPatientsPreResponse(missingIndices, :) = [];
AnatomicalLocsVecPreResponse(missingIndices) = [];
clear missingIndices

missingIndicesAll = find(AnatomicalLocsVecAll == "nan");
AnatomicalLocsPatientsAll(missingIndicesAll, :) = [];
AnatomicalLocsVecAll(missingIndicesAll) = [];
clear missingIndicesAll

ChanNumsPostOnset   = sum(AnatomicalLocsPatientsPostOnset, 2);
ChanNumsPreResponse = sum(AnatomicalLocsPatientsPreResponse, 2);
ChanNumsAll         = sum(AnatomicalLocsPatientsAll, 2);

% percentage per region
ChanNumsPercentPostOnset = nan(length(ChanNumsPostOnset),1);
for i = 1:length(ChanNumsPostOnset)
    ChanNumInAll     = find(AnatomicalLocsVecAll == AnatomicalLocsVecPostOnset(i));
    ratioNum         = ChanNumsAll(ChanNumInAll);
    ratioNumPercent  = floor((ChanNumsPostOnset(i)/ratioNum)*100);
    ChanNumsPercentPostOnset(i) = ratioNumPercent;
end
clear ChanNumInAll ratioNum ratioNumPercent

ChanNumsPercentPreResponse = nan(length(ChanNumsPreResponse),1);
for i = 1:length(ChanNumsPreResponse)
    ChanNumInAll     = find(AnatomicalLocsVecAll == AnatomicalLocsVecPreResponse(i));
    ratioNum         = ChanNumsAll(ChanNumInAll);
    ratioNumPercent  = floor((ChanNumsPreResponse(i)/ratioNum)*100);
    ChanNumsPercentPreResponse(i) = ratioNumPercent;
end
clear ChanNumInAll ratioNum ratioNumPercent

%% ===================== vis percentage (ONLY 2 plots, colored, alpha=0.4) =====================
figure('Units', 'normalized', 'Position', [0.1, 0, 0.2, 0.8]);

position1 = [0.22, 0.70, 0.6, 0.18]; % PostOnset
position2 = [0.22, 0.35, 0.6, 0.18]; % PreResponse

threshold = 1;

% --------------------- PostOnset ---------------------
subplot('Position', position1);
visibleIndices = ChanNumsPostOnset > threshold & ChanNumsPercentPostOnset > threshold;

values = ChanNumsPercentPostOnset(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPostOnset(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans  = ChanNumsPostOnset(visibleIndices);
sortedChans  = sortedChans(sortOrder);

b = bar(sortedValues, 0.5, 'FaceColor', epoch_colors(1,:), 'EdgeColor', 'none');
b.FaceAlpha = 0.4;

xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 10, 'FontWeight', 'bold');

text(max(xlim)-0.03*max(xlim), max(ylim), 'PostOnset', ...
    'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

for i = 1:length(sortedValues)
    ChanNumInAll  = find(AnatomicalLocsVecAll == sortedLabels(i));
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%', sortedChans(i));
    text(i, 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

% --------------------- PreResponse ---------------------
subplot('Position', position2);
visibleIndices = ChanNumsPreResponse > threshold & ChanNumsPercentPreResponse > threshold;

values = ChanNumsPercentPreResponse(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPreResponse(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans  = ChanNumsPreResponse(visibleIndices);
sortedChans  = sortedChans(sortOrder);

b = bar(sortedValues, 0.5, 'FaceColor', epoch_colors(2,:), 'EdgeColor', 'none');
b.FaceAlpha = 0.4;

xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 10, 'FontWeight', 'bold');

text(max(xlim)-0.03*max(xlim), max(ylim), 'PreResponse', ...
    'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

for i = 1:length(sortedValues)
    ChanNumInAll  = find(AnatomicalLocsVecAll == sortedLabels(i));
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%',sortedChans(i));
    text(i, 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

% Title + labels + save
annotation('textbox', [0.1, 0.9, 0.9, 0.1], ...
    'String', 'percentage of channels with significantly different IED and non-IED RTs', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'FontSize', 18);

annotation('textbox', [0.12, 0.45, 0.3, 0.06], ...
    'String', 'channels across patients (%)', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'FontSize', 18, 'Rotation', 90);

set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'RTs_percentage';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
