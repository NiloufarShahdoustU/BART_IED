% Across patient pVal visualization for different patients and different
% locations, now I need to save them in order to visualize them
% Author: Nill

clear;
clc;
close all;
warning('off','all');



%% loading neural and event data
 
% getting the numbers of different patients

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\IEDdata_ACs';
fileList = dir(fullfile(inputFolderName_IEDdata, '*.IEDdata.mat'));
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\AES2024\IED6_8_10_output';
Fnew = 500;
PatientsNum = length(fileList);

%% first visulization: looking at different channels with different pvals
figure;
for pt = 1:PatientsNum
% for pt = 1:7
    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1}; 
    disp(' ');
    disp(['Processing patient ID: ' patientID]);
    ptID = patientID;
    matFile_IEDdata = [inputFolderName_IEDdata '\' ptID '.IEDdata.mat'];
    load(matFile_IEDdata);



    RTSampleSize_NonIED = IEDdata.RTSampleSize_NonIED;
    RTSampleSize_preStim = IEDdata.RTSampleSize_preStim;
    RTSampleSize_periStim = IEDdata.RTSampleSize_periStim;
    RTSampleSize_postStim = IEDdata.RTSampleSize_postStim;
    anatomicalLocs = IEDdata.anatomicalLocs;
    selectedChans = IEDdata.selectedChans;
    selectedChans = selectedChans(1:end-1);
    pVal_preStimRT = IEDdata.pVal_preStimRT;
    pVal_periStimRT = IEDdata.pVal_periStimRT;
    pVal_postStimRT = IEDdata.pVal_postStimRT;
    pVal_preStimRT_filtered = IEDdata.pVal_preStimRT_filtered;
    pVal_periStimRT_filtered = IEDdata.pVal_periStimRT_filtered;
    pVal_postStimRT_filtered = IEDdata.pVal_postStimRT_filtered;

    nChans = length(pVal_preStimRT);

    



    % Create subplot for Pre Stim RT
    subplot(1,3,1);
    scatter(1:length(pVal_preStimRT_filtered), pVal_preStimRT_filtered, 50, 'o', 'filled', 'DisplayName', 'Pre Stim RT');
    xlabel('Channel number', 'FontSize', 18);
    xlim([0 120]);
    ylabel('P-values', 'FontSize', 18);
    title('Pre Stim RT (1:balloon onset)', 'FontSize', 20);
    hold on;
    % legend('show');


    % for i = 1:length(pVal_preStimRT)
    %     if isnan(pVal_preStimRT_filtered(i))
    %         continue;
    %     end
    %     text(i, pVal_preStimRT_filtered(i), ['channel = ', anatomicalLocs{selectedChans(i)}, ' NonIED, IED number of samples = ', num2str(RTSampleSize_NonIED(i)), ', ', num2str(RTSampleSize_preStim(i))], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    % 
    % end


    % Create subplot for Peri Stim RT
    subplot(1,3,2);
    scatter(1:length(pVal_periStimRT_filtered), pVal_periStimRT_filtered, 50, 'o', 'filled', 'DisplayName', 'Peri Stim RT');
    xlabel('Channel number', 'FontSize', 18);
    xlim([0 120]);
    ylabel('P-values', 'FontSize', 18);
    title('Peri Stim RT (balloon onset-100ms:balloon onset+100ms)', 'FontSize', 20);
    hold on;
    % legend('show');
    % 
    % for i = 1:length(pVal_periStimRT)
    %     if isnan(pVal_periStimRT_filtered(i))
    %         continue;
    %     end
    %     text(i, pVal_periStimRT_filtered(i), ['channel = ', anatomicalLocs{selectedChans(i)}, ' NonIED, IED number of samples = ', num2str(RTSampleSize_NonIED(i)), ', ', num2str(RTSampleSize_periStim(i))], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    % 
    % end


    % Create subplot for Post Stim RT
    subplot(1,3,3);
    scatter(1:length(pVal_postStimRT_filtered), pVal_postStimRT_filtered, 50, 'o', 'filled', 'DisplayName', 'Post Stim RT');
    xlabel('Channel number', 'FontSize', 18);
    xlim([0 120]);
    ylabel('P-values', 'FontSize', 18);
    title('Post Stim RT (balloon onset: inflation start)', 'FontSize', 20);
    hold on;
    % legend('show');
    % for i = 1:length(pVal_postStimRT)
    %     if isnan(pVal_postStimRT_filtered(i))
    %         continue;
    %     end
    %        text(i, pVal_postStimRT_filtered(i), ['channel = ', anatomicalLocs{selectedChans(i)}, ' NonIED, IED number of samples = ', num2str(RTSampleSize_NonIED(i)), ', ', num2str(RTSampleSize_postStim(i))], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    % 
    % 
    % end
    % 


end
close all;

%% Second visulization: looking at number of channels per each anatomical location in each patient


AnatoicalLocsNums = 150; % arbitrary number
AnatomicalLocsPatientsPre = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsPeri = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsPost = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsAll = zeros(AnatoicalLocsNums, PatientsNum);

nan_cell_array = repmat({'nan'},AnatoicalLocsNums, 1);
AnatomicalLocsVecPre = string(nan_cell_array);
AnatomicalLocsVecPeri = string(nan_cell_array);
AnatomicalLocsVecPost = string(nan_cell_array);
AnatomicalLocsVecAll = string(nan_cell_array);


% here we need to fill AnatomicalLocsPatients which is the final matrix
% that we're gonna use for visualization and the other one is AnatomicalLocsVec
% which is a vector of strings that has all the anatomical locs consistent
% with the AnatomicalLocsPatients

for pt = 1:PatientsNum

    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1}; 
    disp(' ');
    disp(['Processing patient ID: ' patientID]);
    ptID = patientID;
    matFile_IEDdata = [inputFolderName_IEDdata '\' ptID '.IEDdata.mat'];
    load(matFile_IEDdata);



    RTSampleSize_NonIED = IEDdata.RTSampleSize_NonIED;
    RTSampleSize_preStim = IEDdata.RTSampleSize_preStim;
    RTSampleSize_periStim = IEDdata.RTSampleSize_periStim;
    RTSampleSize_postStim = IEDdata.RTSampleSize_postStim;
    anatomicalLocs = IEDdata.anatomicalLocs;
    anatomicalLocs = string(anatomicalLocs);
    selectedChans = IEDdata.selectedChans;
    selectedChans = selectedChans(1:end-1);
    pVal_preStimRT = IEDdata.pVal_preStimRT;
    pVal_periStimRT = IEDdata.pVal_periStimRT;
    pVal_postStimRT = IEDdata.pVal_postStimRT;
    pVal_preStimRT_filtered = IEDdata.pVal_preStimRT_filtered;
    pVal_periStimRT_filtered = IEDdata.pVal_periStimRT_filtered;
    pVal_postStimRT_filtered = IEDdata.pVal_postStimRT_filtered;
    nChans = length(pVal_preStimRT);


    ChanIndicesPre = find(~isnan(pVal_preStimRT_filtered));
    AnatomicalLocsPre = anatomicalLocs(selectedChans(ChanIndicesPre));
    for pre=1:length(AnatomicalLocsPre)
        element = AnatomicalLocsPre(pre);
        tempIndexInLocs = ismember(AnatomicalLocsVecPre, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPre(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPre(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPre == "nan", 1);
            AnatomicalLocsVecPre(nan_index) = element;
            AnatomicalLocsPatientsPre(nan_index,pt) = AnatomicalLocsPatientsPre(nan_index,pt)+1;
        end
        clear ChanIndicesPre FoundIndexInLocs nan_index
    end

    ChanIndicesPeri = find(~isnan(pVal_periStimRT_filtered));
    AnatomicalLocsPeri = anatomicalLocs(selectedChans(ChanIndicesPeri));
    for peri=1:length(AnatomicalLocsPeri)
        element = AnatomicalLocsPeri(peri);
        tempIndexInLocs = ismember(AnatomicalLocsVecPeri, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPeri(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPeri(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPeri == "nan", 1);
            AnatomicalLocsVecPeri(nan_index) = element;
            AnatomicalLocsPatientsPeri(nan_index,pt) = AnatomicalLocsPatientsPeri(nan_index,pt)+1;
        end
        clear ChanIndicesPeri FoundIndexInLocs nan_index
        
    end



    ChanIndicesPost = find(~isnan(pVal_postStimRT_filtered));
    AnatomicalLocsPost = anatomicalLocs(selectedChans(ChanIndicesPost));
    for post=1:length(AnatomicalLocsPost)
        element = AnatomicalLocsPost(post);
        tempIndexInLocs = ismember(AnatomicalLocsVecPost, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPost(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPost(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPost == "nan", 1);
            AnatomicalLocsVecPost(nan_index) = element;
            AnatomicalLocsPatientsPost(nan_index,pt) = AnatomicalLocsPatientsPost(nan_index,pt)+1;
        end
        clear ChanIndicesPost FoundIndexInLocs nan_index
    end


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
        clear FoundIndexInLocs nan_index tempIndexInLocs
    end

end


%% cleaning AnatomicalLocsVec

% Find elements that start with 'NaC'
startsWithNaC1 = startsWith(AnatomicalLocsVecPre, "NaC");
% Replace elements starting with 'NaC' with missing values
AnatomicalLocsVecPre(startsWithNaC1) = "nan";

startsWithNaC2 = startsWith(AnatomicalLocsVecPeri, "NaC");
AnatomicalLocsVecPeri(startsWithNaC2) = "nan";

startsWithNaC3 = startsWith(AnatomicalLocsVecPost, "NaC");
AnatomicalLocsVecPost(startsWithNaC3) = "nan";

startsWithNaC4 = startsWith(AnatomicalLocsVecAll, "NaC");
AnatomicalLocsVecAll(startsWithNaC4) = "nan";



%% cleaning AnatomicalLocsPatients based on AnatomicalLocsVec

missingIndicesPre = find(AnatomicalLocsVecPre == "nan");
AnatomicalLocsPatientsPre(missingIndicesPre, :) = [];
AnatomicalLocsVecPre(missingIndicesPre) = [];


missingIndicesPeri = find(AnatomicalLocsVecPeri == "nan");
AnatomicalLocsPatientsPeri(missingIndicesPeri, :) = [];
AnatomicalLocsVecPeri(missingIndicesPeri) = [];


missingIndicesPost = find(AnatomicalLocsVecPost == "nan");
AnatomicalLocsPatientsPost(missingIndicesPost, :) = [];
AnatomicalLocsVecPost(missingIndicesPost) = [];

missingIndicesAll = find(AnatomicalLocsVecAll == "nan");
AnatomicalLocsPatientsAll(missingIndicesAll, :) = [];
AnatomicalLocsVecAll(missingIndicesAll) = [];

%%

ChanNumsPre = sum(AnatomicalLocsPatientsPre, 2);
ChanNumsPeri = sum(AnatomicalLocsPatientsPeri, 2);
ChanNumsPost = sum(AnatomicalLocsPatientsPost, 2);
ChanNumsAll = sum(AnatomicalLocsPatientsAll, 2);


%% vis
figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);

% this part is for spacing between subplots
position1 = [0.13, 0.68, 0.775, 0.18];
position2 = [0.13, 0.39, 0.775, 0.18];
position3 = [0.13, 0.1, 0.775, 0.18]; 

% Define a threshold for visible bars
threshold = 4; % Set this to the minimum value that should be visible

% Pre plot
subplot('Position',position1);
% Filter values and labels
visibleIndices = ChanNumsPre > threshold;
values = ChanNumsPre(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend'); % Sort values in descending order
sortedLabels = AnatomicalLocsVecPre(visibleIndices);
sortedLabels = sortedLabels(sortOrder); % Sort labels according to the sorted values
b = bar(sortedValues);
xticks(1:length(sortedValues)); % Adjust tick positions to match the sorted order
xticklabels(sortedLabels);
ylabel('number of channels');
ylim([1 max(sortedValues)+5]);
title('pre');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((sortedValues(i)/ratioNum)*100);
    labelText = sprintf('%d / %d \n (%d%%)', sortedValues(i), ratioNum, ratioNumPercent);
    text(i, sortedValues(i)+ 0.1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

% Peri plot
subplot('Position',position2);
% Filter values and labels
visibleIndices = ChanNumsPeri > threshold;
values = ChanNumsPeri(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPeri(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
b = bar(sortedValues);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
ylabel('number of channels');
ylim([1 max(sortedValues)+5]);
title('peri');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((sortedValues(i)/ratioNum)*100);
    labelText = sprintf('%d / %d \n (%d%%)', sortedValues(i), ratioNum, ratioNumPercent);
    text(i, sortedValues(i)+ 0.1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

% Post plot
subplot('Position',position3);
% Filter values and labels
visibleIndices = ChanNumsPost > threshold;
values = ChanNumsPost(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPost(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
b = bar(sortedValues);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
xlabel('brain area');
ylabel('number of channels');
ylim([1 max(sortedValues)+5]);
title('post');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((sortedValues(i)/ratioNum)*100);
    labelText = sprintf('%d / %d \n (%d%%)', sortedValues(i), ratioNum, ratioNumPercent);
    text(i, sortedValues(i)+ 0.1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

box off;

suptitle("Channels by brain area with significantly different Accuracies from Non-IED trials across patients");


%saving screensize pdf
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
% saveas(gcf,'ACs.pdf')


%% finding percentage
% in this part I need to find out the percentage of the number of channels
% in an specific area in the analysis, among all channels in that brain
% area. I want to use it in the bar plot, I want to sort the bar plot based
% on the percentage and not the number of channels

ChanNumsPrePercent = nan( length(ChanNumsPre),1);
for i = 1:length(ChanNumsPre)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPre(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPre(i)/ratioNum)*100);
    ChanNumsPrePercent(i)=ratioNumPercent;
end


ChanNumsPeriPercent = nan(length(ChanNumsPeri),1);
for i = 1:length(ChanNumsPeri)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPeri(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPeri(i)/ratioNum)*100);
    ChanNumsPeriPercent(i)=ratioNumPercent;
end


ChanNumsPostPercent = nan( length(ChanNumsPost),1);
for i = 1:length(ChanNumsPost)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPost(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPost(i)/ratioNum)*100);
    ChanNumsPostPercent(i)=ratioNumPercent;
end


%% this is the same barplot, but sorted based on the percentage of channels

close all;


figure('Units', 'normalized', 'Position', [0.1, 0, 0.28, 0.9]);

position1 = [0.2, 0.88, 0.775, 0.18];
position2 = [0.2, 0.55, 0.775, 0.18];
position3 = [0.2, 0.2, 0.775, 0.18]; 

threshold = 2; % Set to the minimum value that should be visible

% Pre plot
subplot('Position', position1);
visibleIndices = ChanNumsPre > threshold  & ChanNumsPrePercent > 4;

values = ChanNumsPrePercent(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPre(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPre(visibleIndices);
sortedChans = sortedChans(sortOrder);
bar(sortedValues, 0.5, 'FaceColor', [0.329, 0.329, 0.306]);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% set(gca, 'YTickLabel', get(gca, 'YTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
ylim([0 15]);
text(max(xlim)-0.03*max(xlim), max(ylim), 'pre stimulus onset', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%% \n %d/%d', sortedValues(i), sortedChans(i), NumOfAllChans);
    text(i, sortedValues(i) + 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

% Peri plot
subplot('Position', position2);
visibleIndices = ChanNumsPeri > threshold  & ChanNumsPeriPercent > 4;
values = ChanNumsPeriPercent(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPeri(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPeri(visibleIndices);
sortedChans = sortedChans(sortOrder);
bar(sortedValues, 0.5, 'FaceColor', [0.329, 0.329, 0.306]);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% set(gca, 'YTickLabel', get(gca, 'YTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
ylim([0 max(sortedValues)+5]);
% ylabel('channels across patients (%)',  'FontSize', 16, 'FontWeight', 'bold', 'Position', [-1.4, 0, 0]); 
text(max(xlim)-0.03*max(xlim), max(ylim), 'peri stimulus onset', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%% \n %d/%d', sortedValues(i), sortedChans(i), NumOfAllChans);
    text(i, sortedValues(i) +  1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

% Post plot
subplot('Position', position3);
visibleIndices = ChanNumsPost > threshold & ChanNumsPostPercent > 4;
values = ChanNumsPostPercent(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPost(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPost(visibleIndices);
sortedChans = sortedChans(sortOrder);
bar(sortedValues, 0.5, 'FaceColor', [0.329, 0.329, 0.306]);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% set(gca, 'YTickLabel', get(gca, 'YTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
ylim([0 15]);
xlabel('brain area', 'FontSize', 16, 'FontWeight', 'bold');
text(max(xlim)-0.03*max(xlim), max(ylim), 'post stimulus onset', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%% \n %d/%d', sortedValues(i), sortedChans(i), NumOfAllChans);
    text(i, sortedValues(i) + 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;

% Title and saving
suptitle(" ");
annotation('textbox', [0.1, 0.85, 0.9, 0.1], 'String', 'percentage of channels with significantly different IED and non-IED accuraciess', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
annotation('textbox', [0.11, 0.3, 0.8, 0.06], 'String', 'channels across patients (%)', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18, 'Rotation', 90);

set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);

filename = 'ACs_percentage';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');


%%
close all;

nonIEDtrialsRTs_mean = nan(PatientsNum,1);
IEDtrials_preStimRT_mean = nan(PatientsNum,1);
IEDtrials_periStimRT_mean = nan(PatientsNum,1);
IEDtrials_postStimRT_mean = nan(PatientsNum,1);


for pt = 1:PatientsNum

    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1}; 
    ptID = patientID;
    matFile_IEDdata = [inputFolderName_IEDdata '\' ptID '.IEDdata.mat'];
    load(matFile_IEDdata);

    selectedChans = IEDdata.selectedChans;
    selectedChans = selectedChans(1:end-1);
    NonIEDTrialsRTs =  IEDdata.NonIEDTrialsRTs;
    IEDtrials_preStimRT = IEDdata.IEDtrials_preStimRT;
    IEDtrials_periStimRT = IEDdata.IEDtrials_periStimRT;
    IEDtrials_postStimRT = IEDdata.IEDtrials_postStimRT;   
    pVal_preStimRT_filtered = IEDdata.pVal_preStimRT_filtered;
    pVal_periStimRT_filtered = IEDdata.pVal_periStimRT_filtered;
    pVal_postStimRT_filtered = IEDdata.pVal_postStimRT_filtered;

    nChans = length(selectedChans);
    
    nanIndicesPre = isnan(pVal_preStimRT_filtered);
    IEDtrials_preStimRT_filtered = IEDtrials_preStimRT;
    IEDtrials_preStimRT_filtered(nanIndicesPre, :) = NaN;

    nanIndicesPeri = isnan(pVal_periStimRT_filtered);
    IEDtrials_periStimRT_filtered = IEDtrials_periStimRT;
    IEDtrials_periStimRT_filtered(nanIndicesPeri, :) = NaN;

    nanIndicesPost = isnan(pVal_postStimRT_filtered);
    IEDtrials_postStimRT_filtered = IEDtrials_postStimRT;
    IEDtrials_postStimRT_filtered(nanIndicesPost, :) = NaN;    

    nonIEDtrialsRTs_mean(pt) = nanmean(nanmean(NonIEDTrialsRTs));
    IEDtrials_preStimRT_mean(pt) = nanmean(nanmean(IEDtrials_preStimRT_filtered));
    IEDtrials_periStimRT_mean(pt) = nanmean(nanmean(IEDtrials_periStimRT_filtered));
    IEDtrials_postStimRT_mean(pt) = nanmean(nanmean(IEDtrials_postStimRT_filtered));

end


% Removing NaN elements in each vector:
notNanIndexNonIED = ~isnan(nonIEDtrialsRTs_mean);
nonIEDtrialsRTs_mean = nonIEDtrialsRTs_mean(notNanIndexNonIED);

notNanIndexPre = ~isnan(IEDtrials_preStimRT_mean);
IEDtrials_preStimRT_mean = IEDtrials_preStimRT_mean(notNanIndexPre);

notNanIndexPeri = ~isnan(IEDtrials_periStimRT_mean);
IEDtrials_periStimRT_mean = IEDtrials_periStimRT_mean(notNanIndexPeri);

notNanIndexPost = ~isnan(IEDtrials_postStimRT_mean);
IEDtrials_postStimRT_mean = IEDtrials_postStimRT_mean(notNanIndexPost);

% Perform rank sum tests
[pValuePre, ~] = ranksum(nonIEDtrialsRTs_mean, IEDtrials_preStimRT_mean);
[pValuePeri, ~] = ranksum(nonIEDtrialsRTs_mean, IEDtrials_periStimRT_mean);
[pValuePost, ~] = ranksum(nonIEDtrialsRTs_mean, IEDtrials_postStimRT_mean);

% Visualization:
figure;

vec1 = nonIEDtrialsRTs_mean;
vec2 = IEDtrials_preStimRT_mean;
vec3 = IEDtrials_periStimRT_mean;
vec4 = IEDtrials_postStimRT_mean;

% Concatenate all vectors into a single column vector and create group identifiers
allVecs = [vec1; vec2; vec3; vec4];
group = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1); 4*ones(length(vec4),1)];

% Plot the boxplot with black color
boxplotHandle = boxplot(allVecs, group, 'Labels', {'non-IED', 'pre', 'peri', 'post'}, 'Color', 'k', 'Symbol', '');
xlabel('non-IED trials vs IED trials in three time periods', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0 max(allVecs)+0.2]);
ylabel('mean acccuracy across patients (s)', 'FontSize', 14, 'FontWeight', 'bold');
title('mean accuracy in three time periods of stimulus onset', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Make box plot lines thicker
set(findobj(gca, 'Type', 'line'), 'LineWidth', 1.5);

% Overlay jittered data points in gray
hold on;
jitterAmount = 0.1;
scatter(group + (rand(size(group)) - 0.5) * jitterAmount, allVecs, 'jitter', 'on', 'jitterAmount', jitterAmount, 'MarkerEdgeColor', [0.5 0.5 0.5]);

% Add significance markers if applicable
significanceThreshold = 0.05; % Adjust if needed
if pValuePre < significanceThreshold
    text(2, max(allVecs) * 1.1, '*', 'FontSize', 24, 'HorizontalAlignment', 'center');
end
if pValuePeri < significanceThreshold
    text(3, max(allVecs) * 1.1, '*', 'FontSize', 24, 'HorizontalAlignment', 'center');
end
if pValuePost < significanceThreshold
    text(4, max(allVecs) * 1.1, '*', 'FontSize', 24, 'HorizontalAlignment', 'center');
end

ylimValues = ylim; % Get current y-axis limits
line([1.5 1.5], ylimValues, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Draw dashed line

% Final adjustments and saving the plot
hold off;
set(gca, 'box', 'off', 'tickdir', 'out');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'ACs_boxplot';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');



