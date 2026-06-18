% ccn 2026
% Author: Nill

clear;
clc;
close all;
warning('off','all');

%%

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_RTs_2_chunks\';
fileList = dir(fullfile(inputFolderName_IEDdata, '*.LFPIED.mat'));
PatientsNum = length(fileList);

%% ONLY these epochs + requested colors

PostOnset    = 1;
PreResponse  = 2;

epoch_names  = {'post-onset','pre-response'};

epoch_colors = [
    0.20 0.45 0.70
    0.30 0.70 0.40
];

alpha = 0.4;         % requested
NumberofPermutations = 10000;

%%

% for 2 time periods for RTs
% these means are over trials
nonIEDtrials_bhvr_measure_mean = nan(PatientsNum,1);
IEDtrials_bhvr_measure_mean    = nan(PatientsNum, 2);

for pt = 1:PatientsNum
% for pt = 1:1
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    disp("patient: " + ptID);

    % IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    nChans = length(LFPIED.selectedChans)-1; % removing the last selected chan

    RTs = LFPIED.RTs;

    isNonControl = ~LFPIED.isControl;
    RTs = RTs(isNonControl);
    % ------------------------------------------------------
    
    RTsThreshold = 10;
    OutlierIndices = RTs >= RTsThreshold;
    RTs = RTs(~OutlierIndices);
    nTrials = length(RTs);
    
    % ONLY these epochs
    IEDtrialsPostOnset   = LFPIED.IEDtrialsPostOnset(:, isNonControl);
    IEDtrialsPreResponse = LFPIED.IEDtrialsPreResponse(:, isNonControl);
    
    IEDtrialsPostOnset   = IEDtrialsPostOnset(:, ~OutlierIndices);
    IEDtrialsPreResponse = IEDtrialsPreResponse(:, ~OutlierIndices);


    % first I need to find the nonIED trials, i mean the trials that did
    % not have ANY IEDs in these epochs!!! in ANY chans!
    allTimePoints = IEDtrialsPostOnset + IEDtrialsPreResponse;
    nonIEDIndices = find(all(allTimePoints == 0, 1));

    % now that I've found nonIED trials I need to have their behavioral
    % mearue vector:
    % for the sake of each patient
    nonIEDtrials_bhvr_measure = RTs(nonIEDIndices);

    % Step 2: Use the (non-zscored here) RTs to compute the mean for non-IED trials
    nonIEDtrials_bhvr_measure_mean(pt) = mean(nonIEDtrials_bhvr_measure);

    IEDTrials_bhvr_measure_MeanPerChan_PostOnset   = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreResponse = nan(nChans,1);

    % now I need to take those IEDtrials that ONLY happened in 1 of the 2
    % time periods and not in the other 1:
    IEDtrialsPostOnset_only   = IEDtrialsPostOnset   .* ~IEDtrialsPreResponse;
    IEDtrialsPreResponse_only = IEDtrialsPreResponse .* ~IEDtrialsPostOnset;

    pVal_PostOnset   = nan(1,nChans);
    pVal_PreResponse = nan(1,nChans);

    for chz = 1:nChans

        % Post-Onset
        IEDTrials_bhvr_measure = IEDtrialsPostOnset_only(chz,:);
        % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure .* RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);

        if ~isempty(IEDTrials_bhvr_measure)
            pVal_PostOnset(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PostOnset(chz) < alpha
                IEDTrials_bhvr_measure_MeanPerChan_PostOnset(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PostOnset(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure

        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

        % Pre-Response
        IEDTrials_bhvr_measure = IEDtrialsPreResponse_only(chz,:);
        % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure .* RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);

        if ~isempty(IEDTrials_bhvr_measure)
            pVal_PreResponse(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PreResponse(chz) < alpha
                IEDTrials_bhvr_measure_MeanPerChan_PreResponse(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PreResponse(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure

    end % end for chan

    IEDtrials_bhvr_measure_mean(pt,PostOnset)   = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PostOnset);
    IEDtrials_bhvr_measure_mean(pt,PreResponse) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreResponse);

end

%% save

name = 'IEDtrials_bhvr_measure_mean_RTs';
save([outputFolderName name '.mat'],'IEDtrials_bhvr_measure_mean');

%% preprocessing

IEDtrials_bhvr_measure_mean_PostOnset = IEDtrials_bhvr_measure_mean(:,PostOnset);
IEDtrials_bhvr_measure_mean_PostOnset = IEDtrials_bhvr_measure_mean_PostOnset(~isnan(IEDtrials_bhvr_measure_mean_PostOnset));

IEDtrials_bhvr_measure_mean_PreResponse = IEDtrials_bhvr_measure_mean(:,PreResponse);
IEDtrials_bhvr_measure_mean_PreResponse = IEDtrials_bhvr_measure_mean_PreResponse(~isnan(IEDtrials_bhvr_measure_mean_PreResponse));

%% visualization

% for RTs only these 2 time periods are important:

vec1_full = nonIEDtrials_bhvr_measure_mean(:);
vec2_full = IEDtrials_bhvr_measure_mean_PostOnset(:);
vec3_full = IEDtrials_bhvr_measure_mean_PreResponse(:);

% =========================
% Filter vectors for boxplot
% =========================
valid1 = ~isnan(vec1_full) & vec1_full ~= 0;
valid2 = ~isnan(vec2_full) & vec2_full ~= 0;
valid3 = ~isnan(vec3_full) & vec3_full ~= 0;

vec1 = vec1_full(valid1);
vec2 = vec2_full(valid2);
vec3 = vec3_full(valid3);

% Concatenate all vectors into a single column vector and create group identifiers
allVecs = [vec1; vec2; vec3];
group = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1)];

% square figure, figsize 4 (inches)
figure('Units','inches','Position',[1 1 4 4], 'Visible','on', 'Color','w');

% Plot the boxplot
boxplot(allVecs, group, ...
    'Labels', {'non-IED', epoch_names{1}, epoch_names{2}}, ...
    'Color', 'k', ...
    'Symbol', '');

xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('mean RT', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Make box plot lines slightly thicker
set(findobj(gca, 'Type', 'line'), 'LineWidth', 0.5);

hold on;

% =========================
% Jitter settings
% =========================
jitterAmount = 0.12;
markerSize = 16;
alpha = 0.4;

% Create jittered x positions for FULL vectors
x1_full = nan(size(vec1_full));
x2_full = nan(size(vec2_full));
x3_full = nan(size(vec3_full));

x1_full(valid1) = 1 + (rand(sum(valid1),1) - 0.5) * jitterAmount;
x2_full(valid2) = 2 + (rand(sum(valid2),1) - 0.5) * jitterAmount;
x3_full(valid3) = 3 + (rand(sum(valid3),1) - 0.5) * jitterAmount;

% =========================
% Very tiny gray lines between same participants
% =========================
nParticipants = min([length(vec1_full), length(vec2_full), length(vec3_full)]);

for i = 1:nParticipants
    xline = [x1_full(i), x2_full(i), x3_full(i)];
    yline = [vec1_full(i), vec2_full(i), vec3_full(i)];

    validLine = ~isnan(xline) & ~isnan(yline) & yline ~= 0;

    if sum(validLine) >= 2
        plot(xline(validLine), yline(validLine), ...
            '-', ...
            'Color', [0.82 0.82 0.82], ...
            'LineWidth', 0.25);
    end
end

% =========================
% Overlay jittered data points
% =========================

% non-IED (gray)
scatter(x1_full(valid1), vec1_full(valid1), markerSize, ...
    'filled', ...
    'MarkerFaceColor', [0.55 0.55 0.55], ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', alpha);

% Post-Onset
scatter(x2_full(valid2), vec2_full(valid2), markerSize, ...
    'filled', ...
    'MarkerFaceColor', epoch_colors(1,:), ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', alpha);

% Pre-Response
scatter(x3_full(valid3), vec3_full(valid3), markerSize, ...
    'filled', ...
    'MarkerFaceColor', epoch_colors(2,:), ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', alpha);

hold off;

set(gca, 'box', 'off', 'tickdir', 'out');

% Final adjustments and saving the plot (square 4x4 inches)
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperPosition', [0 0 4 4], 'PaperSize', [4 4]);

filename = 'RTs_boxplot';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

%% p_val for bhv measures accross patients

pValuePostOnset_rt_pt   = ranksum(vec1,vec2);
pValuePreResponse_rt_pt = ranksum(vec1,vec3);

% Specify the output file name
outputFileName = fullfile(outputFolderName, 'RTs_pValues.txt');

% Open the file for writing
fileID = fopen(outputFileName, 'w');

% Write p-values to the file
fprintf(fileID, 'P-values from Rank Sum Tests:\n');
fprintf(fileID, 'P-value (non-IED vs %s): %.4f\n', epoch_names{1}, pValuePostOnset_rt_pt);
fprintf(fileID, 'P-value (non-IED vs %s): %.4f\n', epoch_names{2}, pValuePreResponse_rt_pt);

% Close the file
fclose(fileID);
