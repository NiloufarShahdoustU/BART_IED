% In this code I am taking a look at the percentage of IED trials for all
% the patients in order to have a better understanding of the data. 
% data.
% AUTHOR: Nill

clear;
clc;
close all;
warning('off','all');
%% loading neural and event data
 
% getting the numbers of different patients


inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\IEDTrials';
fileList = dir(fullfile(inputFolderName_IEDtrials, '*.IEDTrials.mat'));
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\SingleNeuron\output\3_IED_percentages';
PatientsNum = length(fileList);

IED_percentage_vector = nan(1, PatientsNum);
xTicks = strings(1, PatientsNum);

for pt = 1:PatientsNum
    
% for pt = 3:3
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    xTicks(pt) = ptID;
    disp("patient: " + ptID);
    IEDtrials = [inputFolderName_IEDtrials '\' ptID '.IEDtrials.mat'];
    load(IEDtrials);
    IEDtrials = IEDTrialsInfo.IEDtrials;
    IEDtrialsAcrossLFPChans = any(IEDtrials, 1); % we're gonna take the trials that an IED occured on any LFP channel
    IED_percentage_vector(pt) = (sum(IEDtrialsAcrossLFPChans)/length(IEDtrialsAcrossLFPChans) )*100;

    clear IEDtrials IEDtrialsAcrossLFPChans

end % pt for


%% plotting bart chart of IEDs percentage accross patients
% Assuming IED_percentage_vector and xTicks are defined previously

figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.99], 'Visible','off');
bar(IED_percentage_vector, 'FaceColor', [1 0.75 0.79]);
xticks(1:length(xTicks)); 
xticklabels(xTicks);
xlabel('patients', 'FontSize',14);
ylabel('IED trials (%)','FontSize',14);

% Calculate the mean of the IED percentage vector
mean_value = mean(IED_percentage_vector);
std_value = std(IED_percentage_vector);

hold on;
yline(mean_value, 'r--', 'LineWidth', 2); 
hold off;
text(length(IED_percentage_vector) - 1, mean_value+0.01, ['Mean: ', num2str(mean_value)], 'VerticalAlignment', 'bottom');
title('IED percentage over trials for each patient', FontSize=30);


% saving

set(gca, 'box', 'off', 'tickdir', 'out');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'IEDTrials_percentage';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
%% debug part


