 % In this version I am going to find the number of IEDs in each channel in
 % each trial

% Author: Nill


clear;
clc;
close all;
warning('off','all');


%%

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_19_IEDnonIED_bhvr_analysis\v3\';
outputFolderName_vecs = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_19_IEDnonIED_bhvr_analysis\v3\ITs\';

fileList = dir(fullfile(inputFolderName_IEDdata, '*.LFPIED.mat'));
PatientsNum = length(fileList);



PreOnset1 = 1;
PreOnset2 = 2;
PostOnset = 3;
PreResponse = 4;
PostResponse = 5;
PreOutcome = 6;
alpha = 0.05;
%%

nonIEDtrials_bhvr_measure_mean = nan(PatientsNum,6); 
IEDTrials_bhvr_measure_mean_PreOnset1 = nan(PatientsNum,1);
IEDTrials_bhvr_measure_mean_PreOnset2 = nan(PatientsNum,1);
IEDTrials_bhvr_measure_mean_PostOnset = nan(PatientsNum,1);
IEDTrials_bhvr_measure_mean_PreResponse = nan(PatientsNum,1);
IEDTrials_bhvr_measure_mean_PostResponse = nan(PatientsNum,1);
IEDTrials_bhvr_measure_mean_PreOutcome = nan(PatientsNum,1);

for pt = 1:PatientsNum
% for pt = 1:1
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    %IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    nChans = length(LFPIED.selectedChans)-1; % removing the last selected chan

    RTs = LFPIED.RTs;
    ITs = LFPIED.ITs;
    RTsThreshold = 10;
    OutlierIndices = RTs >= RTsThreshold;
    RTs = RTs(~OutlierIndices);
    ITs = ITs(~OutlierIndices);
    nTrials = length(RTs);

    IEDtrialsPreOnset1 = LFPIED.IEDtrialsPreOnset1(:, ~OutlierIndices);
    IEDtrialsPreOnset2 = LFPIED.IEDtrialsPreOnset2(:, ~OutlierIndices);
    IEDtrialsPostOnset = LFPIED.IEDtrialsPostOnset(:, ~OutlierIndices);
    IEDtrialsPreResponse = LFPIED.IEDtrialsPreResponse(:, ~OutlierIndices);
    IEDtrialsPostResponse = LFPIED.IEDtrialsPostResponse(:, ~OutlierIndices);
    IEDtrialsPreOutcome = LFPIED.IEDtrialsPreOutcome(:, ~OutlierIndices);




    NumberofPermutations = 10000;

    IEDChanPercent_PreOnset1 = nan(nTrials,1);
    IEDChanPercent_PreOnset2 = nan(nTrials,1);
    IEDChanPercent_PostOnset = nan(nTrials,1);
    IEDChanPercent_PreResponse = nan(nTrials,1);
    IEDChanPercent_PostResponse = nan(nTrials,1);
    IEDChanPercent_PreOutcome = nan(nTrials,1);


    for trial = 1:nTrials
    
        
        IEDChanPercent_PreOnset1(trial) = sum(IEDtrialsPreOnset1(:,trial))/nChans;  % percent of channels with IED occuring in them
        IEDChanPercent_PreOnset2(trial) = sum(IEDtrialsPreOnset2(:,trial))/nChans;
        IEDChanPercent_PostOnset(trial) = sum(IEDtrialsPostOnset(:,trial))/nChans;
        IEDChanPercent_PreResponse(trial) = sum(IEDtrialsPreResponse(:,trial))/nChans;
        IEDChanPercent_PostResponse(trial) = sum(IEDtrialsPostResponse(:,trial))/nChans;
        IEDChanPercent_PreOutcome(trial) = sum(IEDtrialsPreOutcome(:,trial))/nChans;
                      
    end % end for trial

    % we are ONLY taking into account IED trials that IEDs occur in more
    % than 5% of channels:
    IEDpercentThreshold = 0;

    IEDChanPercent_PreOnset1(IEDChanPercent_PreOnset1 <= IEDpercentThreshold) = 0;
    IEDChanPercent_PreOnset2(IEDChanPercent_PreOnset2 <= IEDpercentThreshold) = 0;
    IEDChanPercent_PostOnset(IEDChanPercent_PostOnset <= IEDpercentThreshold) = 0;
    IEDChanPercent_PreResponse(IEDChanPercent_PreResponse <= IEDpercentThreshold) = 0;
    IEDChanPercent_PostResponse(IEDChanPercent_PostResponse <= IEDpercentThreshold) = 0;
    IEDChanPercent_PreOutcome(IEDChanPercent_PreOutcome <= IEDpercentThreshold) = 0;

    pValPostOnset = []; pValPostResponse = []; pValPreOnset1 = []; pValPreOnset2 = []; pValPreOutcome = []; pValPreResponse = [];


    %PreOnset1

    if sum(IEDChanPercent_PreOnset1~=0) % if any trials with IEDs are left
        NonIEDTrials_bhvr_measure = ITs(IEDChanPercent_PreOnset1 == 0);
        IEDTrials_bhvr_measure = ITs(IEDChanPercent_PreOnset1 ~= 0);
        pValPreOnset1 = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
        if (pValPreOnset1 <alpha)
            IEDTrials_bhvr_measure_mean_PreOnset1(pt) = mean(IEDTrials_bhvr_measure);
            nonIEDtrials_bhvr_measure_mean(pt,PreOnset1) = mean(NonIEDTrials_bhvr_measure);
        end
        clear NonIEDTrials_bhvr_measure IEDTrials_bhvr_measure
    end

    %PreOnset2
    if sum(IEDChanPercent_PreOnset2~=0) % if any trials with IEDs are left
        NonIEDTrials_bhvr_measure = ITs(IEDChanPercent_PreOnset2 == 0);
        IEDTrials_bhvr_measure = ITs(IEDChanPercent_PreOnset2 ~= 0);
        pValPreOnset2 = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
        if (pValPreOnset2 <alpha)
            IEDTrials_bhvr_measure_mean_PreOnset2(pt) = mean(IEDTrials_bhvr_measure);
            nonIEDtrials_bhvr_measure_mean(pt,PreOnset2) = mean(NonIEDTrials_bhvr_measure);
        end
        clear NonIEDTrials_bhvr_measure IEDTrials_bhvr_measure
    end

        %PostOnset
    if sum(IEDChanPercent_PostOnset~=0) % if any trials with IEDs are left
        NonIEDTrials_bhvr_measure = ITs(IEDChanPercent_PostOnset == 0);
        IEDTrials_bhvr_measure = ITs(IEDChanPercent_PostOnset ~= 0);
        pValPostOnset = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
        if (pValPostOnset <alpha)
            IEDTrials_bhvr_measure_mean_PostOnset(pt) = mean(IEDTrials_bhvr_measure);
            nonIEDtrials_bhvr_measure_mean(pt,PostOnset) = mean(NonIEDTrials_bhvr_measure);
        end
        clear NonIEDTrials_bhvr_measure IEDTrials_bhvr_measure
    end

        %PreResponse

   if sum(IEDChanPercent_PreResponse~=0) % if any trials with IEDs are left 
        NonIEDTrials_bhvr_measure = ITs(IEDChanPercent_PreResponse == 0);
        IEDTrials_bhvr_measure = ITs(IEDChanPercent_PreResponse ~= 0);
        pValPreResponse = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
        if (pValPreResponse <alpha)
            IEDTrials_bhvr_measure_mean_PreResponse(pt) = mean(IEDTrials_bhvr_measure);
            nonIEDtrials_bhvr_measure_mean(pt,PreResponse) = mean(NonIEDTrials_bhvr_measure);
        end
        clear NonIEDTrials_bhvr_measure IEDTrials_bhvr_measure
   end 
        %PostResponse
    
   if sum(IEDChanPercent_PostResponse~=0) % if any trials with IEDs are left         
        NonIEDTrials_bhvr_measure = ITs(IEDChanPercent_PostResponse == 0);
        IEDTrials_bhvr_measure = ITs(IEDChanPercent_PostResponse ~= 0);
        pValPostResponse = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
        if (pValPostResponse <alpha)
            IEDTrials_bhvr_measure_mean_PostResponse(pt) = mean(IEDTrials_bhvr_measure);
            nonIEDtrials_bhvr_measure_mean(pt,PostResponse) = mean(NonIEDTrials_bhvr_measure);
        end
        clear NonIEDTrials_bhvr_measure IEDTrials_bhvr_measure
   end 
        %PreOutcome
 
   if sum(IEDChanPercent_PreOutcome~=0) % if any trials with IEDs are left
        NonIEDTrials_bhvr_measure = ITs(IEDChanPercent_PreOutcome == 0);
        IEDTrials_bhvr_measure = ITs(IEDChanPercent_PreOutcome ~= 0);
        pValPreOutcome = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
        if (pValPreOutcome <alpha)
            IEDTrials_bhvr_measure_mean_PreOutcome(pt) = mean(IEDTrials_bhvr_measure);
            nonIEDtrials_bhvr_measure_mean(pt,PreOutcome) = mean(NonIEDTrials_bhvr_measure);
        end
        clear NonIEDTrials_bhvr_measure IEDTrials_bhvr_measure
   end






end % patient for clause


%% preprocessing

% we need to aggregate the nonIED trials for 6 time periods
nonIEDtrials_bhvr_measure_mean_aggregate = nanmean(nonIEDtrials_bhvr_measure_mean,2);



% Removing NaN elements in each vector:
% Remove NaNs and zeros from nonIEDtrials_bhvr_measure_mean_aggregate
Clean_nonIEDtrials_bhvr_measure_mean_aggregate = nonIEDtrials_bhvr_measure_mean_aggregate(~isnan(nonIEDtrials_bhvr_measure_mean_aggregate) & nonIEDtrials_bhvr_measure_mean_aggregate ~= 0);

% Remove NaNs and zeros from IEDTrials_bhvr_measure_mean_PreOnset1
Clean_IEDTrials_bhvr_measure_mean_PreOnset1 = IEDTrials_bhvr_measure_mean_PreOnset1(~isnan(IEDTrials_bhvr_measure_mean_PreOnset1) & IEDTrials_bhvr_measure_mean_PreOnset1 ~= 0);

% Remove NaNs and zeros from IEDTrials_bhvr_measure_mean_PreOnset2
Clean_IEDTrials_bhvr_measure_mean_PreOnset2 = IEDTrials_bhvr_measure_mean_PreOnset2(~isnan(IEDTrials_bhvr_measure_mean_PreOnset2) & IEDTrials_bhvr_measure_mean_PreOnset2 ~= 0);

% Remove NaNs and zeros from IEDTrials_bhvr_measure_mean_PostOnset
Clean_IEDTrials_bhvr_measure_mean_PostOnset = IEDTrials_bhvr_measure_mean_PostOnset(~isnan(IEDTrials_bhvr_measure_mean_PostOnset) & IEDTrials_bhvr_measure_mean_PostOnset ~= 0);

% Remove NaNs and zeros from IEDTrials_bhvr_measure_mean_PreResponse
Clean_IEDTrials_bhvr_measure_mean_PreResponse = IEDTrials_bhvr_measure_mean_PreResponse(~isnan(IEDTrials_bhvr_measure_mean_PreResponse) & IEDTrials_bhvr_measure_mean_PreResponse ~= 0);

% Remove NaNs and zeros from IEDTrials_bhvr_measure_mean_PostResponse
Clean_IEDTrials_bhvr_measure_mean_PostResponse = IEDTrials_bhvr_measure_mean_PostResponse(~isnan(IEDTrials_bhvr_measure_mean_PostResponse) & IEDTrials_bhvr_measure_mean_PostResponse ~= 0);

% Remove NaNs and zeros from IEDTrials_bhvr_measure_mean_PreOutcome
Clean_IEDTrials_bhvr_measure_mean_PreOutcome = IEDTrials_bhvr_measure_mean_PreOutcome(~isnan(IEDTrials_bhvr_measure_mean_PreOutcome) & IEDTrials_bhvr_measure_mean_PreOutcome ~= 0);




%% saving clean vectors
% I am saving the vectors above, because sometimes I run the whole matlab
% file!!!!! it takes a lot of time to run!!!

name = 'nonIED';
save([outputFolderName_vecs name '.mat'],'Clean_nonIEDtrials_bhvr_measure_mean_aggregate');

name = 'PreOnset1';
save([outputFolderName_vecs name '.mat'],'Clean_IEDTrials_bhvr_measure_mean_PreOnset1');

name = 'PreOnset2';
save([outputFolderName_vecs name '.mat'],'Clean_IEDTrials_bhvr_measure_mean_PreOnset2');

name = 'PostOnset';
save([outputFolderName_vecs name '.mat'],'Clean_IEDTrials_bhvr_measure_mean_PostOnset');

name = 'PreResponse';
save([outputFolderName_vecs name '.mat'],'Clean_IEDTrials_bhvr_measure_mean_PreResponse');

name = 'PostResponse';
save([outputFolderName_vecs name '.mat'],'Clean_IEDTrials_bhvr_measure_mean_PostResponse');

name = 'PreOutcome';
save([outputFolderName_vecs name '.mat'],'Clean_IEDTrials_bhvr_measure_mean_PreOutcome');


%% Visualization
vec1 = Clean_nonIEDtrials_bhvr_measure_mean_aggregate;
vec2 = Clean_IEDTrials_bhvr_measure_mean_PreOnset1;
% vec3 = Clean_IEDTrials_bhvr_measure_mean_PreOnset2;
vec3 = Clean_IEDTrials_bhvr_measure_mean_PostOnset;
vec4 = Clean_IEDTrials_bhvr_measure_mean_PreResponse;
vec5 = Clean_IEDTrials_bhvr_measure_mean_PostResponse;
vec6 = Clean_IEDTrials_bhvr_measure_mean_PreOutcome;


% Concatenate all vectors into a single column vector and create group identifiers
allVecs = [vec1; vec2; vec3; vec4; vec5; vec6];
group = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1); 4*ones(length(vec4),1); 5*ones(length(vec5),1); 6*ones(length(vec6),1) ];

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.2, 0.3], 'Visible', 'off'); 

% Plot the boxplot with black color
boxplotHandle = boxplot(allVecs, group, 'Labels', {'non-IED', 'PreOnset1', 'PreOnset2', 'PreResponse', 'PostResponse', 'PreOutcome'}, 'Color', 'k', 'Symbol', '');
xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
% ylim([0 5]);
ylabel('mean IT', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Make box plot lines thicker
set(findobj(gca, 'Type', 'line'), 'LineWidth', 1.5);

% Overlay jittered data points in gray
hold on;
jitterAmount = 0.1;
scatter(group + (rand(size(group)) - 0.5) * jitterAmount, allVecs, 'jitter', 'on', 'jitterAmount', jitterAmount, 'MarkerEdgeColor', [0.5 0.5 0.5]);


ylimValues = ylim; % Get current y-axis limits
line([1.5 1.5], ylimValues, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Draw dashed line

% Final adjustments and saving the plot
hold off;
set(gca, 'box', 'off', 'tickdir', 'out');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'ITs_boxplot';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

%%


% Group all vectors into one matrix
all_data = [vec1(:); vec2(:); vec3(:); vec4(:); vec5(:); vec6(:); vec7(:)];

% Define a grouping variable to indicate which vector each data point belongs to
group = [ones(size(vec1(:))); 2*ones(size(vec2(:))); 3*ones(size(vec3(:))); ...
         4*ones(size(vec4(:))); 5*ones(size(vec5(:))); 6*ones(size(vec6(:))); 7*ones(size(vec7(:)))];

% Perform one-way ANOVA
[p,tbl,stats] = anova1(all_data, group);

% Perform post-hoc multiple comparison test to see which pairs are different
results = multcompare(stats);

% Display the results of the pairwise comparisons
% results has the format [group1, group2, lowerCI, diff, upperCI, p-value]

comparison_table = array2table(results, ...
    'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDifference', 'UpperCI', 'PValue'});
disp(comparison_table);


%%


% Capture the formatted output
str = formattedDisplayText(comparison_table);

% Remove HTML tags like <strong> using regular expressions
clean_str = regexprep(str, '<.*?>', '');  % Remove any HTML-like tags

% Specify the full path for the file
filePath = fullfile(outputFolderName, 'ITs_stat.txt');

% Write the cleaned text to the file
fileID = fopen(filePath, 'w');
fprintf(fileID, '%s', clean_str);
fclose(fileID);

%% p_val for bhv measures

pValuePreOnset1 = ranksum(Clean_nonIEDtrials_bhvr_measure_mean_aggregate, ...
                          Clean_IEDTrials_bhvr_measure_mean_PreOnset1);

pValuePreOnset2 = ranksum(Clean_nonIEDtrials_bhvr_measure_mean_aggregate, ...
                          Clean_IEDTrials_bhvr_measure_mean_PreOnset2);

pValuePostOnset = ranksum(Clean_nonIEDtrials_bhvr_measure_mean_aggregate, ...
                          Clean_IEDTrials_bhvr_measure_mean_PostOnset);

pValuePreResponse = ranksum(Clean_nonIEDtrials_bhvr_measure_mean_aggregate, ...
                          Clean_IEDTrials_bhvr_measure_mean_PreResponse);

pValuePostResponse = ranksum(Clean_nonIEDtrials_bhvr_measure_mean_aggregate, ...
                          Clean_IEDTrials_bhvr_measure_mean_PostResponse);

pValuePreOutcome = ranksum(Clean_nonIEDtrials_bhvr_measure_mean_aggregate, ...
                          Clean_IEDTrials_bhvr_measure_mean_PreOutcome);

