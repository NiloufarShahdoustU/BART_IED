% In this code I am taking a look at the percentage of IEDs in channels. 


% AUTHOR: Nill

clear;
clc;
close all;
warning('off','all');
%% loading neural and event data
 
% getting the numbers of different patients


inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\IEDdata';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\SingleNeuron\output\5_IEDs_chan_percentage_across_pt';
inputFolderName_IEDtrials_info = '\\155.100.91.44\d\Data\Nill\BART\IEDTrials';
fileList = dir(fullfile(inputFolderName_IEDtrials_info, '*.IEDTrials.mat'));
PatientsNum = length(fileList);
SignalLength = 3001;
binWidth = 50;

xTicks = strings(1, PatientsNum);



%% Second visulization: looking at number of channels per each anatomical location in each patient


AnatoicalLocsNums = 150; % arbitrary number

AnatomicalLocsPatientsAll = zeros(AnatoicalLocsNums, PatientsNum);

nan_cell_array = repmat({'nan'},AnatoicalLocsNums, 1);
AnatomicalLocsVecAll = string(nan_cell_array);


% here we need to fill AnatomicalLocsPatients which is the final matrix
% that we're gonna use for visualization and the other one is AnatomicalLocsVec
% which is a vector of strings that has all the anatomical locs consistent
% with the AnatomicalLocsPatients

for pt = 1:PatientsNum
% for pt = 1:1
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    xTicks(pt) = ptID;
    disp("patient: " + ptID);
    IEDdata = [inputFolderName_IEDtrials '\' ptID '.IEDdata.mat'];
    load(IEDdata);
    IED_timepoints = IEDdata.IED_timepoints;
    nChans = size(IED_timepoints,1);
    nTrials = size(IED_timepoints,3);
    anatomicalLocs = IEDdata.anatomicalLocs;
    anatomicalLocs = string(anatomicalLocs);
    selectedChans = IEDdata.selectedChans;
    selectedChans = selectedChans(1:end-1);


    IEDtrials = [inputFolderName_IEDtrials_info '\' ptID '.IEDtrials.mat'];
    load(IEDtrials);
    IEDtrials = IEDTrialsInfo.IEDtrials;
    IEDChansAcrossTrials = any(IEDtrials, 2);


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




clear AnatomicalLocsAll_IED AnatomicalLocsAll IEDChansAcrossTrials
end


%% cleaning AnatomicalLocsVec

% Find elements that start with 'NaC'


startsWithNaC4 = startsWith(AnatomicalLocsVecAll, "NaC");
AnatomicalLocsVecAll(startsWithNaC4) = "nan";

%% cleaning AnatomicalLocsPatients based on AnatomicalLocsVec


missingIndicesAll = find(AnatomicalLocsVecAll == "nan");
AnatomicalLocsPatientsAll(missingIndicesAll, :) = [];
AnatomicalLocsVecAll(missingIndicesAll) = [];



%% visualization 1
% visualization of all channels
% so far we've had the visualization of the percentage of channels that are
% included in the study with their percentage. 




figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.8], 'Visible', 'off');
threshold = 1.5; 


ChanNumsAll = sum(AnatomicalLocsPatientsAll, 2);
ChanNumsAllPercent = (ChanNumsAll/sum(ChanNumsAll))*100;
visibleIndices = ChanNumsAllPercent > threshold;
values = ChanNumsAllPercent(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend'); % Sort values in descending order
sortedLabels = AnatomicalLocsVecAll(visibleIndices);
sortedLabels = sortedLabels(sortOrder); % Sort labels according to the sorted values
bar(sortedValues, 'FaceColor', [0.329, 0.329, 0.306]); 
xticks(1:length(sortedValues)); 
xticklabels(sortedLabels);
ylabel('channels(%)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('brain area', 'FontSize', 14, 'FontWeight', 'bold');

title('percentage of each channel across all patients');
% 
% Add percentage values
for i = 1:length(sortedValues)
    text(i, sortedValues(i), sprintf('%.2f', sortedValues(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

ax = gca;
ax.XAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';



% Remove top and right borders
set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

clear values sortedValues sortOrder sortedLabels


% save
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'IEDChans_percentage_almost_correct';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');





%% visualization 2
% now I am finding the results in another way. first I will calculate the
% percentage of the number of channels in each area in each patient and
% then I will find the mean of percentages, in this way I have both mean
% and std and then I need to have error bar for each bar in the barplot. 


figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.8], 'Visible', 'on');
threshold = 1.5;

columnSums = sum(AnatomicalLocsPatientsAll, 1);
normalizedMatrixByChanNum = (AnatomicalLocsPatientsAll ./ columnSums) * 100;

ChanNumsAllPercent_v2 = mean(normalizedMatrixByChanNum, 2);

numSamples = size(normalizedMatrixByChanNum, 2);
rowStd = std(normalizedMatrixByChanNum, 0, 2);
ChanNumsAllPercent_v2_sem = rowStd / sqrt(numSamples);
visibleIndices = ChanNumsAllPercent_v2 > threshold;
values = ChanNumsAllPercent_v2(visibleIndices);
errors = ChanNumsAllPercent_v2_sem(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecAll(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedErrors = errors(sortOrder);

barHandle = bar(sortedValues, 'FaceColor', [0.329, 0.329, 0.306]);

hold on;

x = 1:length(sortedValues);
% Set the lower error bar values to NaN
er = errorbar(x, sortedValues, NaN(size(sortedErrors)), sortedErrors, 'k', 'LineStyle', 'none');
er.CapSize = 10;  % Customize cap size if needed
er.LineWidth = 1.5;  % Customize line width if needed

xticks(x);
xticklabels(sortedLabels);
ylabel('channels(%)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('brain area', 'FontSize', 14, 'FontWeight', 'bold');
title('Percentage of each channel across all patients');

ax = gca;
ax.XAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';


% Add percentage values on top of bars at the same height
yPosition = 0.5; % Adjust this value as needed to set the position
for i = 1:length(sortedValues)
    text(i, yPosition, sprintf('%.2f', sortedValues(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'white', FontSize=12, FontWeight='bold');
end


set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'IEDChans_percentage_correct';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
