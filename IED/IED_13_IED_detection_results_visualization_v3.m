% In this version I am using the result of the newest version of IED
% detection that are saved in IEDdata_bad_chans_removed_v2 
% hopefully it's good


% AUTHOR: Nill
%%
clear;
clc;
close all;
warning('off','all');

%%

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\IEDdata_bad_chans_removed_v2';
inputFolderName_LFPmat = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed\LFPmat_bad_chans_removed';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_13_sanity_check_v3\';
fileList = dir(fullfile(inputFolderName_IEDdata, '*.IEDdata.mat'));
PatientsNum = length(fileList);
% I want to have 25 figures, i.e. 25 samples of LFP visualizations and I
% want them to be in 5x5 order
nFigures = 25;
nRows = 5;
nCols = 5;

%%


for pt = 1:PatientsNum
% for pt = 56:60
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    %IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.IEDdata.mat'];
    load(IEDdata);
    IED_timepoints = IEDdata.IED_timepoints;
    IED_trials = IEDdata.IEDtrials;
    
    % LFP read
    LFPmat = [inputFolderName_LFPmat '\' ptID '.LFPmat.mat'];
    load(LFPmat);
    LFPmat = LFPmatStruct_new.LFPmat;

    % Find the indices of 1 elements in order to find the chan and trial
    % that an IED happened in them.
    [rowIndices, colIndices] = find(IED_trials == 1);
    % Combine the indices into a nx2 matrix. I put chan index in first col
    % and trial index in the second call
    IEDChanTrial_indices = [rowIndices, colIndices];

    % I want to choose 25 random number between 1 and the size of
    % IEDChanTrial_indices. 

    RandomChanTrial = randperm(size(IEDChanTrial_indices,1), nFigures);
    
    figure('Units', 'normalized', 'Position', [0, 0, 0.95, 0.99], 'Visible', 'off'); 
    for indice=1:nFigures
        Chan = IEDChanTrial_indices(RandomChanTrial(indice), 1);
        Trial = IEDChanTrial_indices(RandomChanTrial(indice),2);

        subplot(nRows,nCols, indice);
        LFPData = squeeze(LFPmat(Chan,:,Trial));
        % plot(smooth(LFPData,10)); hold on;
        plot(LFPData);hold on;
        IEDTimes = find(squeeze(IED_timepoints(Chan,:,Trial))== 1);
        scatter(IEDTimes, LFPData(IEDTimes), 300, 'r', 'x', 'LineWidth', 2);
        xlim([1,3001]);
        title(['chan = ' num2str(Chan) ', trial = ' num2str(Trial)]);
        set(gca, 'box', 'off');
        set(gca, 'TickDir', 'out');
        hold off;

        clear LFPData IEDTimes
    end
   sgtitle(['patient = ' char(string(ptID))]);

   clear IEDdata IED_timepoints IED_trials LFPmat rowIndices colIndices
   clear IEDChanTrial_indices RandomChanTrial


   % save
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = string(ptID) + '_IED_detection_sanity_check';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
end


%% debug

