% I have skipped version 4 cause I started with RTs and it was bullshit
% Author: Nill


clear;
clc;
close all;
warning('off','all');

%%

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_acc_2_chunks\';
fileList = dir(fullfile(inputFolderName_IEDdata, '*.LFPIED.mat'));
PatientsNum = length(fileList);

PostResponse = 1;
PreOutcome   = 2;
alpha = 0.05;


%% for 2 time periods (Post-Response, Pre-Outcome) for acc

% these means are over trials
nonIEDtrials_bhvr_measure_mean = nan(PatientsNum,1);
IEDtrials_bhvr_measure_mean = nan(PatientsNum, 2);

AnatoicalLocsNums = 150; % arbitrary number
AnatomicalLocsPatientsPostResponse = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsPreOutcome   = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsAll          = zeros(AnatoicalLocsNums, PatientsNum);

nan_cell_array = repmat({'nan'},AnatoicalLocsNums, 1);
AnatomicalLocsVecPostResponse = string(nan_cell_array);
AnatomicalLocsVecPreOutcome   = string(nan_cell_array);
AnatomicalLocsVecAll          = string(nan_cell_array);

% -------------------- ADDED: store significant channel numbers per patient (per epoch) --------------------
sigChanNums_PostResponse = cell(PatientsNum,1);
sigChanNums_PreOutcome   = cell(PatientsNum,1);
ptID_list                = strings(PatientsNum,1);
% ---------------------------------------------------------------------------------------------------------

for pt = 1:PatientsNum
% for pt = 1:2
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    disp("patient: " + ptID);

    % -------------------- ADDED --------------------
    ptID_list(pt) = string(ptID);
    % -----------------------------------------------

    % IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    nChans = length(LFPIED.selectedChans)-1; % removing the last selected chan

    RTs = LFPIED.RTs;

    anatomicalLocs = LFPIED.anatomicalLocs;

    % ===================== COMBINE LEFT/RIGHT =====================
    % Remove leading "Left " / "Right " so L/R of the same area collapse
    % (works for cell arrays of char OR string arrays)
    anatomicalLocs = regexprep(string(anatomicalLocs), '^(Left|Right)\s+', '');
    % =============================================================

    selectedChans = LFPIED.selectedChans;
    selectedChans = selectedChans(1:end-1);

    % ====================== ONLY include trials that are NOT isControl ======================
    % robustly grab isControl vector if field name differs across files
    if isfield(LFPIED,'isControl')
        isControl = LFPIED.isControl;
    elseif isfield(LFPIED,'isControlTrials')
        isControl = LFPIED.isControlTrials;
    elseif isfield(LFPIED,'controlTrials')
        isControl = LFPIED.controlTrials;
    else
        % if not present, keep all trials (but you asked to filter; this avoids crashing)
        isControl = false(size(RTs));
    end
    isControl = logical(isControl(:)');  % row
    % ======================================================================================

    RTsThreshold = 10;
    OutlierIndices = RTs >= RTsThreshold;

    % keep indices: non-outlier AND NOT control
    keepIdx = ~OutlierIndices & ~isControl;

    RTs = RTs(keepIdx);

    nTrials = length(RTs);
    BankedTrials = LFPIED.BankedTrials;
    BankedTrials = BankedTrials(keepIdx);

    % keep matrices consistent across trials
    IEDtrialsPostResponse = LFPIED.IEDtrialsPostResponse(:, keepIdx);
    IEDtrialsPreOutcome   = LFPIED.IEDtrialsPreOutcome(:,   keepIdx);

    % -------------------------------------------------------------------------------------------------
    % nonIED trials = trials with NO ieds in (PostResponse OR PreOutcome) in ANY chans
    allTimePoints = IEDtrialsPostResponse + IEDtrialsPreOutcome;
    nonIEDIndices = find(all(allTimePoints == 0, 1));

    nonIEDtrials_bhvr_measure = BankedTrials(nonIEDIndices);
    nonIEDtrials_bhvr_measure = double(nonIEDtrials_bhvr_measure);

    nonIEDtrials_bhvr_measure_mean(pt) = mean(nonIEDtrials_bhvr_measure);

    IEDTrials_bhvr_measure_MeanPerChan_PostResponse = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreOutcome   = nan(nChans,1);

    % ONLY-one-epoch masks (PostResponse only, PreOutcome only)
    IEDtrialsPostResponse_only = IEDtrialsPostResponse .* ~IEDtrialsPreOutcome;
    IEDtrialsPreOutcome_only   = IEDtrialsPreOutcome   .* ~IEDtrialsPostResponse;

    pVal_PostResponse = nan(1,nChans);
    pVal_PreOutcome   = nan(1,nChans);

    NumberofPermutations = 10000;

    for chz = 1:nChans

        % PostResponse
        IEDTrials_bhvr_measure = IEDtrialsPostResponse_only(chz,:);
        TrialIndx = find(IEDTrials_bhvr_measure == 1);
        IEDTrials_bhvr_measure(IEDTrials_bhvr_measure == 0) = NaN;
        IEDTrials_bhvr_measure(TrialIndx) = BankedTrials(TrialIndx);
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(~isnan(IEDTrials_bhvr_measure));

        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PostResponse(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PostResponse(chz)< alpha
                IEDTrials_bhvr_measure_MeanPerChan_PostResponse(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PostResponse(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure

        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

        % PreOutcome
        IEDTrials_bhvr_measure = IEDtrialsPreOutcome_only(chz,:);
        TrialIndx = find(IEDTrials_bhvr_measure == 1);
        IEDTrials_bhvr_measure(IEDTrials_bhvr_measure == 0) = NaN;
        IEDTrials_bhvr_measure(TrialIndx) = BankedTrials(TrialIndx);
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(~isnan(IEDTrials_bhvr_measure));

        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PreOutcome(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PreOutcome(chz)< alpha
                IEDTrials_bhvr_measure_MeanPerChan_PreOutcome(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PreOutcome(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure

        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

    end % end for chan

    IEDtrials_bhvr_measure_mean(pt,PostResponse) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PostResponse);
    IEDtrials_bhvr_measure_mean(pt,PreOutcome)   = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreOutcome);

    % ========================= channel percentage by anatomical location =========================

    % PostResponse
    ChanIndices = find(~isnan(IEDTrials_bhvr_measure_MeanPerChan_PostResponse));

    % -------------------- save significant channel numbers for this patient (PostResponse) --------------------
    sigChanNums_PostResponse{pt} = selectedChans(ChanIndices);
    % ---------------------------------------------------------------------------------------------------------------

    timePeriodAnatomicalLoc = anatomicalLocs(selectedChans(ChanIndices));
    for location=1:length(timePeriodAnatomicalLoc)
        element = timePeriodAnatomicalLoc(location);
        tempIndexInLocs = ismember(AnatomicalLocsVecPostResponse, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPostResponse(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPostResponse(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPostResponse == "nan", 1);
            AnatomicalLocsVecPostResponse(nan_index) = element;
            AnatomicalLocsPatientsPostResponse(nan_index,pt) = AnatomicalLocsPatientsPostResponse(nan_index,pt)+1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear ChanIndices timePeriodAnatomicalLoc

    % PreOutcome
    ChanIndices = find(~isnan(IEDTrials_bhvr_measure_MeanPerChan_PreOutcome));

    % -------------------- ADDED: save significant channel numbers for this patient (PreOutcome) --------------------
    sigChanNums_PreOutcome{pt} = selectedChans(ChanIndices);
    % ---------------------------------------------------------------------------------------------------------------

    timePeriodAnatomicalLoc = anatomicalLocs(selectedChans(ChanIndices));
    for location=1:length(timePeriodAnatomicalLoc)
        element = timePeriodAnatomicalLoc(location);
        tempIndexInLocs = ismember(AnatomicalLocsVecPreOutcome, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPreOutcome(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPreOutcome(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPreOutcome == "nan", 1);
            AnatomicalLocsVecPreOutcome(nan_index) = element;
            AnatomicalLocsPatientsPreOutcome(nan_index,pt) = AnatomicalLocsPatientsPreOutcome(nan_index,pt)+1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear ChanIndices timePeriodAnatomicalLoc

    clear IEDtrialsPostResponse IEDtrialsPreOutcome
    clear IEDtrialsPostResponse_only IEDtrialsPreOutcome_only
end

for pt = 1:PatientsNum
    AnatomicalLocsAll = anatomicalLocs(selectedChans);
    for all=1:length(AnatomicalLocsAll)
        element = AnatomicalLocsAll(all);
        tempIndexInLocs = ismember(AnatomicalLocsVecAll, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsAll(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecAll == "nan", 1);
            AnatomicalLocsVecAll(nan_index) = element;
            AnatomicalLocsPatientsAll(nan_index,pt) = AnatomicalLocsPatientsAll(nan_index,pt)+1;
        end
       clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear AnatomicalLocsAll
end

%%

name = 'IEDtrials_bhvr_measure_mean_acc';
save([outputFolderName name '.mat'],'IEDtrials_bhvr_measure_mean');

% -------------------- ADDED: save significant channel numbers (per patient, per epoch) --------------------
name = 'sigChanNums_acc_PostResponse_PreOutcome';
save([outputFolderName name '.mat'],'sigChanNums_PostResponse','sigChanNums_PreOutcome','ptID_list');
% ---------------------------------------------------------------------------------------------------------

%% preprocessing

IEDtrials_bhvr_measure_mean_PostResponse = IEDtrials_bhvr_measure_mean(:,PostResponse);
IEDtrials_bhvr_measure_mean_PostResponse = IEDtrials_bhvr_measure_mean_PostResponse(~isnan(IEDtrials_bhvr_measure_mean_PostResponse));

IEDtrials_bhvr_measure_mean_PreOutcome = IEDtrials_bhvr_measure_mean(:,PreOutcome);
IEDtrials_bhvr_measure_mean_PreOutcome = IEDtrials_bhvr_measure_mean_PreOutcome(~isnan(IEDtrials_bhvr_measure_mean_PreOutcome));

% removing zeros:
IEDtrials_bhvr_measure_mean_PostResponse = IEDtrials_bhvr_measure_mean_PostResponse(IEDtrials_bhvr_measure_mean_PostResponse ~= 0);
IEDtrials_bhvr_measure_mean_PreOutcome   = IEDtrials_bhvr_measure_mean_PreOutcome(IEDtrials_bhvr_measure_mean_PreOutcome ~= 0);

%% visualization

vec1 = nonIEDtrials_bhvr_measure_mean;
vec2 = IEDtrials_bhvr_measure_mean_PostResponse;
vec3 = IEDtrials_bhvr_measure_mean_PreOutcome;

allVecs = [vec1; vec2; vec3];
group = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1)];

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.2, 0.3], 'Visible', 'on');

boxplotHandle = boxplot(allVecs, group, ...
    'Labels', {'non-IED','Post-Response','Pre-Outcome'}, ...
    'Color', 'k', ...
    'Symbol', '');

xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('normalized acc', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

set(findobj(gca, 'Type', 'line'), 'LineWidth', 0.5);

hold on;
jitterAmount = 0.1;
markerSize = 10;
scatter(group + (rand(size(group)) - 0.5) * jitterAmount, allVecs, markerSize, ...
        'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', jitterAmount);

hold off;
set(gca, 'box', 'off', 'tickdir', 'out');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
% filename = 'acc_boxplot';
% saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

%% p_val for bhv measures accross patients

pValuePostResponse_acc_pt = ranksum(vec1,vec2);
pValuePreOutcome_acc_pt   = ranksum(vec1,vec3);

outputFileName = fullfile(outputFolderName, 'acc_pValues.txt');
fileID = fopen(outputFileName, 'w');

fprintf(fileID, 'P-values from Rank Sum Tests:\n');
fprintf(fileID, 'P-value (non-IED vs PostResponse): %.4f\n', pValuePostResponse_acc_pt);
fprintf(fileID, 'P-value (non-IED vs PreOutcome): %.4f\n', pValuePreOutcome_acc_pt);

fclose(fileID);

%% visualization channel percentage

% cleaning AnatomicalLocsVec
% (after left/right removal, "Left Cerebral White Matter" becomes "Cerebral White Matter")
startsWithNaC = startsWith(AnatomicalLocsVecPostResponse, "NaC");
containsCWM   = contains(AnatomicalLocsVecPostResponse, "Cerebral White Matter");
AnatomicalLocsVecPostResponse(startsWithNaC | containsCWM) = "nan";
clear startsWithNaC containsCWM

startsWithNaC = startsWith(AnatomicalLocsVecPreOutcome, "NaC");
containsCWM   = contains(AnatomicalLocsVecPreOutcome, "Cerebral White Matter");
AnatomicalLocsVecPreOutcome(startsWithNaC | containsCWM) = "nan";
clear startsWithNaC containsCWM

startsWithNaC = startsWith(AnatomicalLocsVecAll, "NaC");
containsCWM   = contains(AnatomicalLocsVecAll, "Cerebral White Matter");
AnatomicalLocsVecAll(startsWithNaC | containsCWM) = "nan";
clear startsWithNaC containsCWM

% cleaning AnatomicalLocsPatients based on AnatomicalLocsVec

missingIndices = find(AnatomicalLocsVecPostResponse == "nan");
AnatomicalLocsPatientsPostResponse(missingIndices, :) = [];
AnatomicalLocsVecPostResponse(missingIndices) = [];
clear missingIndices

missingIndices = find(AnatomicalLocsVecPreOutcome == "nan");
AnatomicalLocsPatientsPreOutcome(missingIndices, :) = [];
AnatomicalLocsVecPreOutcome(missingIndices) = [];
clear missingIndices

missingIndicesAll = find(AnatomicalLocsVecAll == "nan");
AnatomicalLocsPatientsAll(missingIndicesAll, :) = [];
AnatomicalLocsVecAll(missingIndicesAll) = [];

ChanNumsPostResponse = sum(AnatomicalLocsPatientsPostResponse, 2);
ChanNumsPreOutcome   = sum(AnatomicalLocsPatientsPreOutcome, 2);
ChanNumsAll          = sum(AnatomicalLocsPatientsAll, 2);

% finding percentage
ChanNumsPercentPostResponse = nan(length(ChanNumsPostResponse),1);
for i = 1:length(ChanNumsPostResponse)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPostResponse(i));
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPostResponse(i)/ratioNum)*100);
    ChanNumsPercentPostResponse(i)=ratioNumPercent;
end
clear ChanNumInAll ratioNum ratioNumPercent

ChanNumsPercentPreOutcome = nan(length(ChanNumsPreOutcome),1);
for i = 1:length(ChanNumsPreOutcome)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPreOutcome(i));
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPreOutcome(i)/ratioNum)*100);
    ChanNumsPercentPreOutcome(i)=ratioNumPercent;
end
clear ChanNumInAll ratioNum ratioNumPercent

%% vis percentage (ONLY Post-Response & Pre-Outcome)
figure('Units', 'normalized', 'Position', [0.1, 0, 0.2, 0.8]);

position1 = [0.22, 0.75, 0.6, 0.18]; % PostResponse
position2 = [0.22, 0.45, 0.6, 0.18]; % PreOutcome

threshold = 5;


epoch_names = {'post-response','pre-outcome'};
epoch_colors = [
    0.85 0.33 0.10
    0.60 0.40 0.80
];

bar_alpha = 0.4;

% -------------------- PostResponse plot --------------------
subplot('Position', position1);
visibleIndices = ChanNumsPostResponse > threshold & ChanNumsPercentPostResponse > threshold;

values = ChanNumsPercentPostResponse(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPostResponse(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPostResponse(visibleIndices);
sortedChans = sortedChans(sortOrder);

h1 = bar(sortedValues, 0.5);
h1.FaceColor = epoch_colors(1,:);
h1.FaceAlpha = bar_alpha;
h1.EdgeColor = 'none';

xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 10, 'FontWeight', 'bold');
text(max(xlim)-0.03*max(xlim), max(ylim), epoch_names{1}, ...
    'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

for i = 1:length(sortedValues)
    % label: number of significant chans (not a percent)
    labelText = sprintf('%d', sortedChans(i));
    text(i, 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

% -------------------- PreOutcome plot --------------------
subplot('Position', position2);
visibleIndices = ChanNumsPreOutcome > threshold & ChanNumsPercentPreOutcome > threshold;

values = ChanNumsPercentPreOutcome(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPreOutcome(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPreOutcome(visibleIndices);
sortedChans = sortedChans(sortOrder);

h2 = bar(sortedValues, 0.5);
h2.FaceColor = epoch_colors(2,:);
h2.FaceAlpha = bar_alpha;
h2.EdgeColor = 'none';

xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 10, 'FontWeight', 'bold');
text(max(xlim)-0.03*max(xlim), max(ylim), epoch_names{2}, ...
    'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

for i = 1:length(sortedValues)
    labelText = sprintf('%d', sortedChans(i));
    text(i, 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

suptitle(" ");
annotation('textbox', [0.1, 0.9, 0.9, 0.1], ...
    'String', 'percentage of channels with different accuracy in IED vs non-IED', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);

annotation('textbox', [0.12, 0.45, 0.3, 0.06], ...
    'String', 'channels across patients (%)', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18, 'Rotation', 90);

set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'acc_percentage';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
