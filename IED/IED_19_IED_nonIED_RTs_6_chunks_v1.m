 % Across patient I am going to do a perm test on RTs between IED and
 % nonIED trials.
 % I want to check when IEDs happen is different time periods (out of those 6 time periods)
 % in terms of RT.

% Author: Nill


clear;
clc;
close all;
warning('off','all');


%%

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_19_IEDnonIED_bhvr_analysis\v1\';
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


nonIEDtrials_bhvr_measure_mean = nan(PatientsNum,6); % for 6 time periods
% at the end we ONLY need to take into account the trials as IED trials
% that all the 6 rows are nonNan
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
    RTsThreshold = 10;
    OutlierIndices = RTs >= RTsThreshold;

    IEDtrialsPreOnset1 = LFPIED.IEDtrialsPreOnset1;
    IEDtrialsPreOnset2 = LFPIED.IEDtrialsPreOnset2;
    IEDtrialsPostOnset = LFPIED.IEDtrialsPostOnset;
    IEDtrialsPreResponse = LFPIED.IEDtrialsPreResponse;
    IEDtrialsPostResponse = LFPIED.IEDtrialsPostResponse;
    IEDtrialsPreOutcome = LFPIED.IEDtrialsPreOutcome;



    
    pVal_PreOnset1 = nan(1,nChans);
    pVal_PreOnset2 = nan(1,nChans);
    pVal_PostOnset = nan(1,nChans);
    pVal_PreResponse = nan(1,nChans);
    pVal_PostResponse = nan(1,nChans);
    pVal_PreOutcome = nan(1,nChans);



    NumberofPermutations = 10000;

    NonIEDTrials_bhvr_measure_MeanPerChan_PreOnset1 = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreOnset1 = nan(nChans,1);

    
    NonIEDTrials_bhvr_measure_MeanPerChan_PreOnset2 = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreOnset2 = nan(nChans,1);

    
    NonIEDTrials_bhvr_measure_MeanPerChan_PostOnset = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PostOnset = nan(nChans,1);

    
    NonIEDTrials_bhvr_measure_MeanPerChan_PreResponse = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreResponse = nan(nChans,1);

    
    NonIEDTrials_bhvr_measure_MeanPerChan_PostResponse = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PostResponse = nan(nChans,1);

    
    NonIEDTrials_bhvr_measure_MeanPerChan_PreOutcome = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreOutcome = nan(nChans,1);




    for chz = 1:nChans
    
        %PreOnset1
        NonIEDTrials_bhvr_measure = ~IEDtrialsPreOnset1(chz,:); % here RT
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*RTs;
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure(NonIEDTrials_bhvr_measure ~= 0);

        IEDTrials_bhvr_measure = IEDtrialsPreOnset1(chz,:); % here RT
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);
   
        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PreOnset1(chz) = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
            if pVal_PreOnset1(chz)<alpha
                NonIEDTrials_bhvr_measure_MeanPerChan_PreOnset1(chz) = nanmean(NonIEDTrials_bhvr_measure);
                IEDTrials_bhvr_measure_MeanPerChan_PreOnset1(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PreOnset1(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure NonIEDTrials_bhvr_measure


         %PreOnset2
        NonIEDTrials_bhvr_measure = ~IEDtrialsPreOnset2(chz,:); % here RT
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*RTs;
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure(NonIEDTrials_bhvr_measure ~= 0);

        IEDTrials_bhvr_measure = IEDtrialsPreOnset2(chz,:); % here RT
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);
   
        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PreOnset2(chz) = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
            if pVal_PreOnset2(chz)<alpha
                NonIEDTrials_bhvr_measure_MeanPerChan_PreOnset2(chz) = nanmean(NonIEDTrials_bhvr_measure);
                IEDTrials_bhvr_measure_MeanPerChan_PreOnset2(chz) = nanmean(IEDTrials_bhvr_measure);
            end
            else
            pVal_PreOnset2(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure NonIEDTrials_bhvr_measure


        %PostOnset
        NonIEDTrials_bhvr_measure = ~IEDtrialsPostOnset(chz,:); % here RT
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*RTs;
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure(NonIEDTrials_bhvr_measure ~= 0);

        IEDTrials_bhvr_measure = IEDtrialsPostOnset(chz,:); % here RT
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);
   
        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PostOnset(chz) = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
            if pVal_PostOnset(chz)< alpha
                NonIEDTrials_bhvr_measure_MeanPerChan_PostOnset(chz) = nanmean(NonIEDTrials_bhvr_measure);
                IEDTrials_bhvr_measure_MeanPerChan_PostOnset(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PostOnset(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure NonIEDTrials_bhvr_measure

        %PreResponse
        NonIEDTrials_bhvr_measure = ~IEDtrialsPreResponse(chz,:); % here RT
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*RTs;
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure(NonIEDTrials_bhvr_measure ~= 0);

        IEDTrials_bhvr_measure = IEDtrialsPreResponse(chz,:); % here RT
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);
   
        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PreResponse(chz) = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
            if pVal_PreResponse(chz)< alpha
                NonIEDTrials_bhvr_measure_MeanPerChan_PreResponse(chz) = nanmean(NonIEDTrials_bhvr_measure);
                IEDTrials_bhvr_measure_MeanPerChan_PreResponse(chz) = nanmean(IEDTrials_bhvr_measure);
            end
            else
            pVal_PreResponse(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure NonIEDTrials_bhvr_measure


         %PostResponse
        NonIEDTrials_bhvr_measure = ~IEDtrialsPostResponse(chz,:); % here RT
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*RTs;
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure(NonIEDTrials_bhvr_measure ~= 0);

        IEDTrials_bhvr_measure = IEDtrialsPostResponse(chz,:); % here RT
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);
   
        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PostResponse(chz) = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
            if pVal_PostResponse(chz)< alpha
                NonIEDTrials_bhvr_measure_MeanPerChan_PostResponse(chz) = nanmean(NonIEDTrials_bhvr_measure);
                IEDTrials_bhvr_measure_MeanPerChan_PostResponse(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PostResponse(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure NonIEDTrials_bhvr_measure


         %PreOutcome
        NonIEDTrials_bhvr_measure = ~IEDtrialsPreOutcome(chz,:); % here RT
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure.*RTs;
        NonIEDTrials_bhvr_measure = NonIEDTrials_bhvr_measure(NonIEDTrials_bhvr_measure ~= 0);

        IEDTrials_bhvr_measure = IEDtrialsPreOutcome(chz,:); % here RT
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*(~OutlierIndices); % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);
   
        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PreOutcome(chz) = permutationTest(IEDTrials_bhvr_measure,NonIEDTrials_bhvr_measure, NumberofPermutations);
            if pVal_PreOutcome(chz)< alpha
                NonIEDTrials_bhvr_measure_MeanPerChan_PreOutcome(chz) = nanmean(NonIEDTrials_bhvr_measure);
                IEDTrials_bhvr_measure_MeanPerChan_PreOutcome(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PreOutcome(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure NonIEDTrials_bhvr_measure
                    
                      
    
    end % end for chan

    nonIEDtrials_bhvr_measure_mean(pt,PreOnset1) = nanmean(NonIEDTrials_bhvr_measure_MeanPerChan_PreOnset1);
    IEDTrials_bhvr_measure_mean_PreOnset1(pt) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreOnset1);

    nonIEDtrials_bhvr_measure_mean(pt,PreOnset2) = nanmean(NonIEDTrials_bhvr_measure_MeanPerChan_PreOnset2);
    IEDTrials_bhvr_measure_mean_PreOnset2(pt) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreOnset2);

    nonIEDtrials_bhvr_measure_mean(pt,PostOnset) = nanmean(NonIEDTrials_bhvr_measure_MeanPerChan_PostOnset);
    IEDTrials_bhvr_measure_mean_PostOnset(pt) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PostOnset);

    nonIEDtrials_bhvr_measure_mean(pt,PreResponse) = nanmean(NonIEDTrials_bhvr_measure_MeanPerChan_PreResponse);
    IEDTrials_bhvr_measure_mean_PreResponse(pt) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreResponse);

    nonIEDtrials_bhvr_measure_mean(pt,PostResponse) = nanmean(NonIEDTrials_bhvr_measure_MeanPerChan_PostResponse);
    IEDTrials_bhvr_measure_mean_PostResponse(pt) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PostResponse);

    nonIEDtrials_bhvr_measure_mean(pt,PreOutcome) = nanmean(NonIEDTrials_bhvr_measure_MeanPerChan_PreOutcome);
    IEDTrials_bhvr_measure_mean_PreOutcome(pt) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreOutcome);



end
%% preprocessing

% we need to aggregate the nonIED trials for 6 time periods
nonIEDtrials_bhvr_measure_mean_aggregate = nanmean(nonIEDtrials_bhvr_measure_mean,2);



% Removing NaN elements in each vector:

Clean_nonIEDtrials_bhvr_measure_mean_aggregate = nonIEDtrials_bhvr_measure_mean_aggregate(~isnan(nonIEDtrials_bhvr_measure_mean_aggregate));

Clean_IEDTrials_bhvr_measure_mean_PreOnset1 = IEDTrials_bhvr_measure_mean_PreOnset1(~isnan(IEDTrials_bhvr_measure_mean_PreOnset1));

Clean_IEDTrials_bhvr_measure_mean_PreOnset2 = IEDTrials_bhvr_measure_mean_PreOnset2(~isnan(IEDTrials_bhvr_measure_mean_PreOnset2));

Clean_IEDTrials_bhvr_measure_mean_PostOnset = IEDTrials_bhvr_measure_mean_PostOnset(~isnan(IEDTrials_bhvr_measure_mean_PostOnset));

Clean_IEDTrials_bhvr_measure_mean_PreResponse = IEDTrials_bhvr_measure_mean_PreResponse(~isnan(IEDTrials_bhvr_measure_mean_PreResponse));

Clean_IEDTrials_bhvr_measure_mean_PostResponse = IEDTrials_bhvr_measure_mean_PostResponse(~isnan(IEDTrials_bhvr_measure_mean_PostResponse));

Clean_IEDTrials_bhvr_measure_mean_PreOutcome = IEDTrials_bhvr_measure_mean_PreOutcome(~isnan(IEDTrials_bhvr_measure_mean_PreOutcome));
%% saving clean vectors
% I am saving the vectors above, because sometimes I run the whole matlab
% file!!!!! it takes a lot of time to run!!!
outputFolderName_vecs = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_19_IEDnonIED_bhvr_analysis\v1\RTs\';

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


%% Visualization
vec1 = Clean_nonIEDtrials_bhvr_measure_mean_aggregate;
vec2 = Clean_IEDTrials_bhvr_measure_mean_PreOnset1;
vec3 = Clean_IEDTrials_bhvr_measure_mean_PreOnset2;
vec4 = Clean_IEDTrials_bhvr_measure_mean_PostOnset;
vec5 = Clean_IEDTrials_bhvr_measure_mean_PreResponse;
vec6 = Clean_IEDTrials_bhvr_measure_mean_PostResponse;
vec7 = Clean_IEDTrials_bhvr_measure_mean_PreOutcome;


% Concatenate all vectors into a single column vector and create group identifiers
allVecs = [vec1; vec2; vec3; vec4; vec5; vec6; vec7];
group = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1); 4*ones(length(vec4),1); 5*ones(length(vec5),1); 6*ones(length(vec6),1); 7*ones(length(vec7),1) ];

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.2, 0.3], 'Visible', 'off'); 

% Plot the boxplot with black color
boxplotHandle = boxplot(allVecs, group, 'Labels', {'non-IED', 'PreOnset1', 'PreOnset2', 'PostOnset', 'PreResponse', 'PostResponse', 'PreOutcome'}, 'Color', 'k', 'Symbol', '');
xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
% ylim([0 5]);
ylabel('mean RT', 'FontSize', 14, 'FontWeight', 'bold');
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
filename = 'RTs_boxplot';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

%%

% Assuming the vectors are already defined as vec1, vec2, ..., vec7

% Group all vectors into one matrix
all_data = [vec1(:); vec2(:); vec3(:); vec4(:); vec5(:); vec6(:); vec7(:)];

% Define a grouping variable to indicate which vector each data point belongs to
group = [ones(size(vec1(:))); 2*ones(size(vec2(:))); 3*ones(size(vec3(:))); ...
         4*ones(size(vec4(:))); 5*ones(size(vec5(:))); 6*ones(size(vec6(:))); ...
         7*ones(size(vec7(:)))];

% Perform one-way ANOVA
[p,tbl,stats] = anova1(all_data, group);

% Perform post-hoc multiple comparison test to see which pairs are different
results = multcompare(stats);

% Display the results of the pairwise comparisons
% results has the format [group1, group2, lowerCI, diff, upperCI, p-value]
disp('Pairwise comparisons and p-values:');
comparison_table = array2table(results, ...
    'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDifference', 'UpperCI', 'PValue'});
disp(comparison_table);

