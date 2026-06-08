% In this code I am taking a look at the number of IEDs in each channel
% AUTHOR: Nill

clear;
clc;
close all;
warning('off','all');
%% loading neural and event data

inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\IEDdata';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\SingleNeuron\output\6_IED_sum_chans';
inputFolderName_IEDtrials_info = '\\155.100.91.44\d\Data\Nill\BART\IEDTrials';
fileList = dir(fullfile(inputFolderName_IEDtrials_info, '*.IEDTrials.mat'));
PatientsNum = length(fileList);
SignalLength = 3001;
xTicks = strings(1, PatientsNum);

%% Second visulization: looking at number of channels per each anatomical location in each patient


AnatoicalLocsNums = 150; % arbitrary number

AnatomicalLocsPatientsAll = zeros(AnatoicalLocsNums, PatientsNum); % this will contain the number of IEDs in a specific channel 
AnatomicalLocsPatientsAll_chans = zeros(AnatoicalLocsNums, PatientsNum);
% normalized to the number of trials for each patient.

nan_cell_array = repmat({'nan'},AnatoicalLocsNums, 1);
AnatomicalLocsVecAll = string(nan_cell_array);


for pt = 1:PatientsNum
% for pt = 2:2
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
    IED_timepoints = IED_timepoints(1:end-1,:,:);


    IEDtrials = [inputFolderName_IEDtrials_info '\' ptID '.IEDtrials.mat'];
    load(IEDtrials);
    IEDtrials = IEDTrialsInfo.IEDtrials;


    AnatomicalLocsAll = anatomicalLocs(selectedChans);
    for all=1:length(AnatomicalLocsAll)
        element = AnatomicalLocsAll(all);
        tempIndexInLocs = ismember(AnatomicalLocsVecAll, element);
        FoundIndexInLocs = find(tempIndexInLocs);


        slice = squeeze(IED_timepoints(all, :, :)); % Extract the 2D slice 
        non_nan_values = slice(~isnan(slice)); % Find non-NaN values
        sum_value = 0;
        if ~isempty(non_nan_values)
            sum_value = (sum(non_nan_values)/(nTrials*SignalLength))*10000; % Calculate the sum of IEDs for this channel and find normalized to the number of trials
        end

        if ~isempty(FoundIndexInLocs)
            
            AnatomicalLocsPatientsAll(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll(FoundIndexInLocs,pt) + sum_value;
            AnatomicalLocsPatientsAll_chans(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll_chans(FoundIndexInLocs,pt) + 1;
        else
            nan_index = find(AnatomicalLocsVecAll == "nan", 1);
            AnatomicalLocsVecAll(nan_index) = element;
            AnatomicalLocsPatientsAll(nan_index,pt) = AnatomicalLocsPatientsAll(nan_index,pt) + sum_value;
            AnatomicalLocsPatientsAll_chans(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll_chans(FoundIndexInLocs,pt) + 1;
            
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element slice non_nan_values sum_value
    end


clear AnatomicalLocsAll_IED AnatomicalLocsAll 
end


% cleaning AnatomicalLocsVec

% Find elements that start with 'NaC'


startsWithNaC4 = startsWith(AnatomicalLocsVecAll, "NaC");
AnatomicalLocsVecAll(startsWithNaC4) = "nan";

% cleaning AnatomicalLocsPatients based on AnatomicalLocsVec


missingIndicesAll = find(AnatomicalLocsVecAll == "nan");
AnatomicalLocsPatientsAll(missingIndicesAll, :) = [];
AnatomicalLocsPatientsAll_chans(missingIndicesAll, :) = [];
AnatomicalLocsVecAll(missingIndicesAll) = [];


%% visualization 1/3

% Original code to compute necessary values
AllChansSum = sum(AnatomicalLocsPatientsAll_chans, 2);
ZeroIndices = find(AllChansSum == 0);
AllChansSum(ZeroIndices) = [];
AnatomicalLocsPatientsAll(ZeroIndices,:) = [];
AnatomicalLocsPatientsAll_chans(ZeroIndices,:) = [];

threshold = 0.3;
NormalizedMatrix = (AnatomicalLocsPatientsAll ./ AllChansSum) * 100;
AnatomicalLocsNumberOfIEDs = mean(NormalizedMatrix, 2);

visibleIndices_IED = AnatomicalLocsNumberOfIEDs > threshold;
values_IED = AnatomicalLocsNumberOfIEDs(visibleIndices_IED);
[sortedValues_IED, sortOrder_IED] = sort(values_IED, 'descend');
sortedLabels = AnatomicalLocsVecAll(visibleIndices_IED);
sortedLabels = sortedLabels(sortOrder_IED);

% Calculate standard error for error bars
numSamples = size(NormalizedMatrix, 2);
rowStd = std(NormalizedMatrix, 0, 2);
ChanNumsAllPercent_v2_sem = rowStd / sqrt(numSamples);
errors_IED = ChanNumsAllPercent_v2_sem(visibleIndices_IED);
sortedErrors = errors_IED(sortOrder_IED);

% Create the bar plot
figure('Units', 'normalized', 'Position', [0.3, 0.2, 0.4, 0.5], 'Visible', 'on');
barHandle = bar(sortedValues_IED, 'FaceColor', [0.329, 0.329, 0.306]);
hold on;

% Add error bars
x = 1:length(sortedValues_IED);
er = errorbar(x, sortedValues_IED, NaN(size(sortedErrors)), sortedErrors, 'k', 'LineStyle', 'none');
er.CapSize = 10;
er.LineWidth = 1.5;

% Customize axes and labels
xticks(x);
xticklabels(sortedLabels);
ylabel('channels(%)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('brain area', 'FontSize', 14, 'FontWeight', 'bold');
title({'Number of IEDs in each channel across patients'; 'normalized by the # of trials and # of channels'}, 'FontSize', 16, 'FontWeight', 'bold');

% Modify x-tick and y-tick labels font size and weight
ax = gca;
ax.XAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';

% Add percentage values on top of bars
yPosition = 0.05;
for i = 1:length(sortedValues_IED)
    text(i, yPosition, sprintf('%.2f', sortedValues_IED(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
end

% Remove top and right borders
set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

% Save the figure
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'IEDChans_percentage_normalize_trial_channels';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');




%% Another version of Visualization 2/3:
% now in this section I am going to create another figure of the data that
% is normalized by the number of Trials ONLY:


% Compute necessary values
AllChansSum = sum(AnatomicalLocsPatientsAll_chans, 2);
ZeroIndices = find(AllChansSum == 0);
AllChansSum(ZeroIndices) = [];
AnatomicalLocsPatientsAll(ZeroIndices,:) = [];
AnatomicalLocsPatientsAll_chans(ZeroIndices,:) = [];

threshold = 0.1;
NormalizedMatrix = AnatomicalLocsPatientsAll; % not normalizing based on the number of channels

AnatomicalLocsNumberOfIEDs = mean(NormalizedMatrix, 2);

visibleIndices_IED = AnatomicalLocsNumberOfIEDs > threshold;
values_IED = AnatomicalLocsNumberOfIEDs(visibleIndices_IED);
[sortedValues_IED, sortOrder_IED] = sort(values_IED, 'descend');
sortedLabels = AnatomicalLocsVecAll(visibleIndices_IED);
sortedLabels = sortedLabels(sortOrder_IED);

% Calculate standard error for error bars
numSamples = size(NormalizedMatrix, 2);
rowStd = std(NormalizedMatrix, 0, 2);
ChanNumsAllPercent_v2_sem = rowStd / sqrt(numSamples);
errors_IED = ChanNumsAllPercent_v2_sem(visibleIndices_IED);
sortedErrors = errors_IED(sortOrder_IED);

% Create the bar plot with error bars
figure('Units', 'normalized', 'Position', [0.3, 0.2, 0.4, 0.5], 'Visible', 'on');
barHandle = bar(sortedValues_IED, 'FaceColor', [0.329, 0.329, 0.306]);
hold on;

% Add error bars
x = 1:length(sortedValues_IED);
er = errorbar(x, sortedValues_IED, NaN(size(sortedErrors)), sortedErrors, 'k', 'LineStyle', 'none');
er.CapSize = 10;
er.LineWidth = 1.5;

% Customize axes and labels
xticks(x);
xticklabels(sortedLabels);
ylabel('channels(%)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('brain area', 'FontSize', 14, 'FontWeight', 'bold');
title({'Number of IEDs in each channel across patients'; 'normalized by the # of trials'}, 'FontSize', 16, 'FontWeight', 'bold');

% Modify x-tick and y-tick labels font size and weight
ax = gca;
ax.XAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';

% Add percentage values on top of bars
yPosition = 0.02;
for i = 1:length(sortedValues_IED)
    text(i, yPosition, sprintf('%.2f', sortedValues_IED(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
end

% Remove top and right borders
set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

% Save the figure
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'IEDChans_percentage_normalize_trials';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');




%% Another version of Visualization 3/3:
% now in this section I am going to create another figure of the data that
% is normalized by the number of Channels ONLY:

clear;
clc;
close all;
warning('off','all');
% loading neural and event data

inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\IEDdata';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\SingleNeuron\output\6_IED_sum_chans';
inputFolderName_IEDtrials_info = '\\155.100.91.44\d\Data\Nill\BART\IEDTrials';
fileList = dir(fullfile(inputFolderName_IEDtrials_info, '*.IEDTrials.mat'));
PatientsNum = length(fileList);
SignalLength = 3001;
xTicks = strings(1, PatientsNum);

% Second visulization: looking at number of channels per each anatomical location in each patient


AnatoicalLocsNums = 150; % arbitrary number

AnatomicalLocsPatientsAll = zeros(AnatoicalLocsNums, PatientsNum); % this will contain the number of IEDs in a specific channel 
AnatomicalLocsPatientsAll_chans = zeros(AnatoicalLocsNums, PatientsNum);
% normalized to the number of trials for each patient.

nan_cell_array = repmat({'nan'},AnatoicalLocsNums, 1);
AnatomicalLocsVecAll = string(nan_cell_array);

for pt = 1:PatientsNum
% for pt = 2:2
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
    IED_timepoints = IED_timepoints(1:end-1,:,:);


    IEDtrials = [inputFolderName_IEDtrials_info '\' ptID '.IEDtrials.mat'];
    load(IEDtrials);
    IEDtrials = IEDTrialsInfo.IEDtrials;


    AnatomicalLocsAll = anatomicalLocs(selectedChans);
    for all=1:length(AnatomicalLocsAll)
        element = AnatomicalLocsAll(all);
        tempIndexInLocs = ismember(AnatomicalLocsVecAll, element);
        FoundIndexInLocs = find(tempIndexInLocs);


        slice = squeeze(IED_timepoints(all, :, :)); % Extract the 2D slice 
        non_nan_values = slice(~isnan(slice)); % Find non-NaN values
        sum_value = 0;
        if ~isempty(non_nan_values)
            sum_value = (sum(non_nan_values)/SignalLength)*100; % Calculate the sum of IEDs for this channel and find normalized to the number of trials
        end

        if ~isempty(FoundIndexInLocs)
            
            AnatomicalLocsPatientsAll(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll(FoundIndexInLocs,pt) + sum_value;
            AnatomicalLocsPatientsAll_chans(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll_chans(FoundIndexInLocs,pt) + 1;
        else
            nan_index = find(AnatomicalLocsVecAll == "nan", 1);
            AnatomicalLocsVecAll(nan_index) = element;
            AnatomicalLocsPatientsAll(nan_index,pt) = AnatomicalLocsPatientsAll(nan_index,pt) + sum_value;
            AnatomicalLocsPatientsAll_chans(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll_chans(FoundIndexInLocs,pt) + 1;
            
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element slice non_nan_values sum_value
    end


clear AnatomicalLocsAll_IED AnatomicalLocsAll 
end


% cleaning AnatomicalLocsVec

% Find elements that start with 'NaC'


startsWithNaC4 = startsWith(AnatomicalLocsVecAll, "NaC");
AnatomicalLocsVecAll(startsWithNaC4) = "nan";

% cleaning AnatomicalLocsPatients based on AnatomicalLocsVec


missingIndicesAll = find(AnatomicalLocsVecAll == "nan");
AnatomicalLocsPatientsAll(missingIndicesAll, :) = [];
AnatomicalLocsPatientsAll_chans(missingIndicesAll, :) = [];
AnatomicalLocsVecAll(missingIndicesAll) = [];


%% visualization 3/3

% Summarize channels and remove zero data points
AllChansSum = sum(AnatomicalLocsPatientsAll_chans, 2);
ZeroIndices = find(AllChansSum == 0);
AllChansSum(ZeroIndices) = [];
AnatomicalLocsPatientsAll(ZeroIndices,:) = [];
AnatomicalLocsPatientsAll_chans(ZeroIndices,:) = [];

threshold = 0.6;

% Normalize and calculate mean values
NormalizedMatrix = (AnatomicalLocsPatientsAll ./ AllChansSum) * 100;
AnatomicalLocsNumberOfIEDs = mean(NormalizedMatrix, 2);

visibleIndices_IED = AnatomicalLocsNumberOfIEDs > threshold;
values_IED = AnatomicalLocsNumberOfIEDs(visibleIndices_IED);
[sortedValues_IED, sortOrder_IED] = sort(values_IED, 'descend');
sortedLabels = AnatomicalLocsVecAll(visibleIndices_IED);
sortedLabels = sortedLabels(sortOrder_IED);

% Calculate standard error for error bars
numSamples = size(NormalizedMatrix, 2);
rowStd = std(NormalizedMatrix, 0, 2);
ChanNumsAllPercent_v2_sem = rowStd / sqrt(numSamples);
errors_IED = ChanNumsAllPercent_v2_sem(visibleIndices_IED);
sortedErrors = errors_IED(sortOrder_IED);

% Create the bar plot
figure('Units', 'normalized', 'Position', [0.3, 0.2, 0.4, 0.5], 'Visible', 'on');
barHandle = bar(sortedValues_IED, 'FaceColor', [0.329, 0.329, 0.306]);
hold on;

% Add error bars
x = 1:length(sortedValues_IED);
er = errorbar(x, sortedValues_IED, NaN(size(sortedErrors)), sortedErrors, 'k', 'LineStyle', 'none');
er.CapSize = 10;
er.LineWidth = 1.5;

% Customize axes and labels
xticks(x);
xticklabels(sortedLabels);
ylabel('channels(%)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('brain area', 'FontSize', 14, 'FontWeight', 'bold');
title({'Number of IEDs in each channel across patients'; 'normalized by the # of channels'}, 'FontSize', 16, 'FontWeight', 'bold');

% Modify x-tick and y-tick labels font size and weight
ax = gca;
ax.XAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';

% Remove top and right borders
set(gca, 'box', 'off');
set(gca, 'TickDir', 'out');

% Save the figure
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'IEDChans_percentage_normalize_channels';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');



