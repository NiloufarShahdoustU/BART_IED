% In this version I am using the result of the newest version of IED
% detection on the chuck of data that I have!


% AUTHOR: Nill
%%
clear;
clc;
close all;
warning('off','all');

%%

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';

fileList = dir(fullfile(inputFolderName_IEDdata, '*.LFPIED.mat'));
PatientsNum = length(fileList);
% I want to have 25 figures, i.e. 25 samples of LFP visualizations and I
% want them to be in 5x5 order
nFigures = 25;
nRows = 5;
nCols = 5;

%%

for pt = 1:PatientsNum
% for pt = 1:1
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    %IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    %% PreOnset1
    IED_timepointsPreOnset1 = LFPIED.IED_timepointsPreOnset1;
    IEDtrialsPreOnset1 = LFPIED.IEDtrialsPreOnset1;
    LFPmatPreOnset1 = LFPIED.LFPmatPreOnset1;


    % Find the indices of 1 elements in order to find the chan and trial
    % that an IED happened in them.
    [rowIndices, colIndices] = find(IEDtrialsPreOnset1 == 1);
    % Combine the indices into a nx2 matrix. I put chan index in first col
    % and trial index in the second call
    IEDChanTrial_indices = [rowIndices, colIndices];

    % I want to choose 25 random number between 1 and the size of
    % IEDChanTrial_indices. 
    if(~isempty(IEDChanTrial_indices))
        if (size(IEDChanTrial_indices,1) >= nFigures)
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), nFigures);
        else
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), size(IEDChanTrial_indices,1));
            nFigures = size(IEDChanTrial_indices,1);
            % here we change the number of figures to a smaller number than
            % 25 and then at the end we need to change it back to 25;
        end
        
        figure('Units', 'normalized', 'Position', [0, 0, 0.95, 0.99], 'Visible', 'off'); 
        for indice=1:nFigures
            Chan = IEDChanTrial_indices(RandomChanTrial(indice), 1);
            Trial = IEDChanTrial_indices(RandomChanTrial(indice),2);
    
            subplot(nRows,nCols, indice);
            LFPData = squeeze(LFPmatPreOnset1(Chan,:,Trial));
            % plot(smooth(LFPData,10)); hold on;
            plot(LFPData);hold on;
            IEDTimes = find(squeeze(IED_timepointsPreOnset1(Chan,:,Trial))== 1);
            scatter(IEDTimes, LFPData(IEDTimes), 300, 'r', 'x', 'LineWidth', 2);
            xlim([1,500]);
            title(['chan = ' num2str(Chan) ', trial = ' num2str(Trial)]);
            set(gca, 'box', 'off');
            set(gca, 'TickDir', 'out');
            hold off;
    
            clear LFPData IEDTimes
        end
       sgtitle(['patient = ' char(string(ptID)) ', period = PreOnset1']);
    
    
       % save
        outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_13_sanity_check_v5\preOnset1\';
        set(gcf, 'Units', 'inches');
        screenposition = get(gcf, 'Position');
        set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
        filename = string(ptID) + '_IED_detection_sanity_check';
        saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
        nFigures = 25;
    end 
    clear rowIndices colIndices RandomChanTrial IEDChanTrial_indices outputFolderName

    %% PreOnset2
    IED_timepointsPreOnset2 = LFPIED.IED_timepointsPreOnset2;
    IEDtrialsPreOnset2 = LFPIED.IEDtrialsPreOnset2;
    LFPmatPreOnset2 = LFPIED.LFPmatPreOnset2;


    % Find the indices of 1 elements in order to find the chan and trial
    % that an IED happened in them.
    [rowIndices, colIndices] = find(IEDtrialsPreOnset2 == 1);
    % Combine the indices into a nx2 matrix. I put chan index in first col
    % and trial index in the second call
    IEDChanTrial_indices = [rowIndices, colIndices];

    % I want to choose 25 random number between 1 and the size of
    % IEDChanTrial_indices. 
    if(~isempty(IEDChanTrial_indices))
        if (size(IEDChanTrial_indices,1) >= nFigures)
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), nFigures);
        else
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), size(IEDChanTrial_indices,1));
            nFigures = size(IEDChanTrial_indices,1);
        end

        figure('Units', 'normalized', 'Position', [0, 0, 0.95, 0.99], 'Visible', 'off'); 
        for indice=1:nFigures
            Chan = IEDChanTrial_indices(RandomChanTrial(indice), 1);
            Trial = IEDChanTrial_indices(RandomChanTrial(indice),2);

            subplot(nRows,nCols, indice);
            LFPData = squeeze(LFPmatPreOnset2(Chan,:,Trial));
            % plot(smooth(LFPData,10)); hold on;
            plot(LFPData);hold on;
            IEDTimes = find(squeeze(IED_timepointsPreOnset2(Chan,:,Trial))== 1);
            scatter(IEDTimes, LFPData(IEDTimes), 300, 'r', 'x', 'LineWidth', 2);
            xlim([1,500]);
            title(['chan = ' num2str(Chan) ', trial = ' num2str(Trial)]);
            set(gca, 'box', 'off');
            set(gca, 'TickDir', 'out');
            hold off;

            clear LFPData IEDTimes
        end
       sgtitle(['patient = ' char(string(ptID)) ', period = PreOnset2']);


       % save
        outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_13_sanity_check_v5\PreOnset2\';
        set(gcf, 'Units', 'inches');
        screenposition = get(gcf, 'Position');
        set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
        filename = string(ptID) + '_IED_detection_sanity_check';
        saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
        nFigures = 25;
    end 

    clear rowIndices colIndices RandomChanTrial IEDChanTrial_indices outputFolderName

    %% PostOnset
    IED_timepointsPostOnset = LFPIED.IED_timepointsPostOnset;
    IEDtrialsPostOnset = LFPIED.IEDtrialsPostOnset;
    LFPmatPostOnset = LFPIED.LFPmatPostOnset;


    % Find the indices of 1 elements in order to find the chan and trial
    % that an IED happened in them.
    [rowIndices, colIndices] = find(IEDtrialsPostOnset == 1);
    % Combine the indices into a nx2 matrix. I put chan index in first col
    % and trial index in the second call
    IEDChanTrial_indices = [rowIndices, colIndices];

    % I want to choose 25 random number between 1 and the size of
    % IEDChanTrial_indices. 
    if(~isempty(IEDChanTrial_indices))
        if (size(IEDChanTrial_indices,1) >= nFigures)
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), nFigures);
        else
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), size(IEDChanTrial_indices,1));
            nFigures = size(IEDChanTrial_indices,1);
        end

        figure('Units', 'normalized', 'Position', [0, 0, 0.95, 0.99], 'Visible', 'off'); 
        for indice=1:nFigures
            Chan = IEDChanTrial_indices(RandomChanTrial(indice), 1);
            Trial = IEDChanTrial_indices(RandomChanTrial(indice),2);

            subplot(nRows,nCols, indice);
            LFPData = squeeze(LFPmatPostOnset(Chan,:,Trial));
            % plot(smooth(LFPData,10)); hold on;
            plot(LFPData);hold on;
            IEDTimes = find(squeeze(IED_timepointsPostOnset(Chan,:,Trial))== 1);
            scatter(IEDTimes, LFPData(IEDTimes), 300, 'r', 'x', 'LineWidth', 2);
            xlim([1,500]);
            title(['chan = ' num2str(Chan) ', trial = ' num2str(Trial)]);
            set(gca, 'box', 'off');
            set(gca, 'TickDir', 'out');
            hold off;

            clear LFPData IEDTimes
        end
       sgtitle(['patient = ' char(string(ptID)) ', period = PostOnset']);


       % save
        outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_13_sanity_check_v5\PostOnset\';
        set(gcf, 'Units', 'inches');
        screenposition = get(gcf, 'Position');
        set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
        filename = string(ptID) + '_IED_detection_sanity_check';
        saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
        nFigures = 25;
    end

    clear rowIndices colIndices RandomChanTrial IEDChanTrial_indices outputFolderName


     %% PreResponse
    IED_timepointsPreResponse = LFPIED.IED_timepointsPreResponse;
    IEDtrialsPreResponse = LFPIED.IEDtrialsPreResponse;
    LFPmatPreResponse = LFPIED.LFPmatPreResponse;


    % Find the indices of 1 elements in order to find the chan and trial
    % that an IED happened in them.
    [rowIndices, colIndices] = find(IEDtrialsPreResponse == 1);
    % Combine the indices into a nx2 matrix. I put chan index in first col
    % and trial index in the second call
    IEDChanTrial_indices = [rowIndices, colIndices];

    % I want to choose 25 random number between 1 and the size of
    % IEDChanTrial_indices. 
    if(~isempty(IEDChanTrial_indices))
        if (size(IEDChanTrial_indices,1) >= nFigures)
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), nFigures);
        else
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), size(IEDChanTrial_indices,1));
            nFigures = size(IEDChanTrial_indices,1);
        end

        figure('Units', 'normalized', 'Position', [0, 0, 0.95, 0.99], 'Visible', 'off'); 
        for indice=1:nFigures
            Chan = IEDChanTrial_indices(RandomChanTrial(indice), 1);
            Trial = IEDChanTrial_indices(RandomChanTrial(indice),2);

            subplot(nRows,nCols, indice);
            LFPData = squeeze(LFPmatPreResponse(Chan,:,Trial));
            % plot(smooth(LFPData,10)); hold on;
            plot(LFPData);hold on;
            IEDTimes = find(squeeze(IED_timepointsPreResponse(Chan,:,Trial))== 1);
            scatter(IEDTimes, LFPData(IEDTimes), 300, 'r', 'x', 'LineWidth', 2);
            xlim([1,500]);
            title(['chan = ' num2str(Chan) ', trial = ' num2str(Trial)]);
            set(gca, 'box', 'off');
            set(gca, 'TickDir', 'out');
            hold off;

            clear LFPData IEDTimes
        end
       sgtitle(['patient = ' char(string(ptID)) ', period = PreResponse']);


       % save
        outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_13_sanity_check_v5\PreResponse\';
        set(gcf, 'Units', 'inches');
        screenposition = get(gcf, 'Position');
        set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
        filename = string(ptID) + '_IED_detection_sanity_check';
        saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
        nFigures = 25;
    end 

    clear rowIndices colIndices RandomChanTrial IEDChanTrial_indices outputFolderName


 %% PostResponse
    IED_timepointsPostResponse = LFPIED.IED_timepointsPostResponse;
    IEDtrialsPostResponse = LFPIED.IEDtrialsPostResponse;
    LFPmatPostResponse = LFPIED.LFPmatPostResponse;


    % Find the indices of 1 elements in order to find the chan and trial
    % that an IED happened in them.
    [rowIndices, colIndices] = find(IEDtrialsPostResponse == 1);
    % Combine the indices into a nx2 matrix. I put chan index in first col
    % and trial index in the second call
    IEDChanTrial_indices = [rowIndices, colIndices]; % this might be 0 sometimes!!

    % I want to choose 25 random number between 1 and the size of
    % IEDChanTrial_indices. 
    if(~isempty(IEDChanTrial_indices))
        if (size(IEDChanTrial_indices,1) >= nFigures)
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), nFigures);
        else
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), size(IEDChanTrial_indices,1));
            nFigures = size(IEDChanTrial_indices,1);
        end

        figure('Units', 'normalized', 'Position', [0, 0, 0.95, 0.99], 'Visible', 'off'); 
        for indice=1:nFigures
            Chan = IEDChanTrial_indices(RandomChanTrial(indice), 1);
            Trial = IEDChanTrial_indices(RandomChanTrial(indice),2);

            subplot(nRows,nCols, indice);
            LFPData = squeeze(LFPmatPostResponse(Chan,:,Trial));
            % plot(smooth(LFPData,10)); hold on;
            plot(LFPData);hold on;
            IEDTimes = find(squeeze(IED_timepointsPostResponse(Chan,:,Trial))== 1);
            scatter(IEDTimes, LFPData(IEDTimes), 300, 'r', 'x', 'LineWidth', 2);
            xlim([1,500]);
            title(['chan = ' num2str(Chan) ', trial = ' num2str(Trial)]);
            set(gca, 'box', 'off');
            set(gca, 'TickDir', 'out');
            hold off;

            clear LFPData IEDTimes
        end
       sgtitle(['patient = ' char(string(ptID)) ', period = PostResponse']);


       % save
        outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_13_sanity_check_v5\PostResponse\';
        set(gcf, 'Units', 'inches');
        screenposition = get(gcf, 'Position');
        set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
        filename = string(ptID) + '_IED_detection_sanity_check';
        saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
        nFigures = 25;
    end 
    clear rowIndices colIndices RandomChanTrial IEDChanTrial_indices outputFolderName


     %% PreOutcome
    IED_timepointsPreOutcome = LFPIED.IED_timepointsPreOutcome;
    IEDtrialsPreOutcome = LFPIED.IEDtrialsPreOutcome;
    LFPmatPreOutcome = LFPIED.LFPmatPreOutcome;


    % Find the indices of 1 elements in order to find the chan and trial
    % that an IED happened in them.
    [rowIndices, colIndices] = find(IEDtrialsPreOutcome == 1);
    % Combine the indices into a nx2 matrix. I put chan index in first col
    % and trial index in the second call
    IEDChanTrial_indices = [rowIndices, colIndices];

    % I want to choose 25 random number between 1 and the size of
    % IEDChanTrial_indices. 
    if(~isempty(IEDChanTrial_indices))
        if (size(IEDChanTrial_indices,1) >= nFigures)
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), nFigures);
        else
            RandomChanTrial = randperm(size(IEDChanTrial_indices,1), size(IEDChanTrial_indices,1));
            nFigures = size(IEDChanTrial_indices,1);
        end

        figure('Units', 'normalized', 'Position', [0, 0, 0.95, 0.99], 'Visible', 'off'); 
        for indice=1:nFigures
            Chan = IEDChanTrial_indices(RandomChanTrial(indice), 1);
            Trial = IEDChanTrial_indices(RandomChanTrial(indice),2);

            subplot(nRows,nCols, indice);
            LFPData = squeeze(LFPmatPreOutcome(Chan,:,Trial));
            % plot(smooth(LFPData,10)); hold on;
            plot(LFPData);hold on;
            IEDTimes = find(squeeze(IED_timepointsPreOutcome(Chan,:,Trial))== 1);
            scatter(IEDTimes, LFPData(IEDTimes), 300, 'r', 'x', 'LineWidth', 2);
            xlim([1,500]);
            title(['chan = ' num2str(Chan) ', trial = ' num2str(Trial)]);
            set(gca, 'box', 'off');
            set(gca, 'TickDir', 'out');
            hold off;

            clear LFPData IEDTimes 
        end
       sgtitle(['patient = ' char(string(ptID)) ', period = PreOutcome']);


       % save
        outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_13_sanity_check_v5\PreOutcome\';
        set(gcf, 'Units', 'inches');
        screenposition = get(gcf, 'Position');
        set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
        filename = string(ptID) + '_IED_detection_sanity_check';
        saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
        nFigures = 25;
    end

    clear rowIndices colIndices RandomChanTrial IEDChanTrial_indices outputFolderName

    %% all done



   clear IEDdata IED_timepoints IED_trials LFPmat rowIndices colIndices
   clear IEDChanTrial_indices RandomChanTrial
end


%% debug

% aaa = sum(IEDtrialsPreOnset1(:));