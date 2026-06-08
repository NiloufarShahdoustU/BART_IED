% In this code I am taking a look at the percentage of IEDs in time
% AUTHOR: Nill

clear;
clc;
close all;
warning('off','all');
%% loading neural and event data
 
% getting the numbers of different patients


inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\IEDdata';
fileList = dir(fullfile(inputFolderName_IEDtrials, '*.IEDdata.mat'));
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\SingleNeuron\output\4_IEDs_in_time_percentage';
PatientsNum = length(fileList);
SignalLength = 3001;
binWidth = 50;

IED_percentage_vector = nan(1, PatientsNum);
xTicks = strings(1, PatientsNum);
AllPtAllChAllTr = nan(PatientsNum, SignalLength); % each row is one patient data and is the 
                                                 % mean of number of IEDs on all trials and
                                                 % all channels
%%
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
    MeanofIEds = nan(1, SignalLength);
    % Calculate the mean and normalize for each element in the second dimension
    for j = 1:SignalLength
        slice = squeeze(IED_timepoints(:, j, :)); % Extract the 2D slice for the j-th element of the second dimension
        non_nan_values = slice(~isnan(slice)); % Find non-NaN values
        
        if ~isempty(non_nan_values)
            sum_value = sum(non_nan_values); % Calculate the mean of non-NaN values
            normalized_mean = sum_value / (nChans * nTrials); % Normalize by the product of the first and last dimensions
            MeanofIEds(j) = normalized_mean; % Store the result in the final vector
        end
    end
    AllPtAllChAllTr(pt,:) = MeanofIEds*100*100; %percentage



    clear IEDtrials IEDtrialsAcrossLFPChans

end % pt for
%% remove nans
nanRemovedAllPtAllChAllTr = AllPtAllChAllTr;
nanRemovedAllPtAllChAllTr(isnan(nanRemovedAllPtAllChAllTr)) = 0;


%% save

% numBins = round(size(nanRemovedAllPtAllChAllTr, 2) / binWidth);
MyVec = mean(nanRemovedAllPtAllChAllTr, 1);
timeVector = linspace(binWidth, size(nanRemovedAllPtAllChAllTr, 2), SignalLength);
sem_vec = std(nanRemovedAllPtAllChAllTr, 0, 1) / sqrt(size(nanRemovedAllPtAllChAllTr, 1));

figure('Units', 'normalized',  'Visible','on');
plotMeanAndStd(MyVec, sem_vec, timeVector, 'red', 100);

set(gca, 'box', 'off', 'tickdir', 'out');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'IEDTrials_percentage_inTime';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');


