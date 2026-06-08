% In this code I am taking a look at the percentage of IED trials for all
% the patients for all 6 time chunks in order to have a better understanding of the data. 
% data.
% AUTHOR: Nill

clear;
clc;
close all;
warning('off','all');
%% loading neural and event data
 
% getting the numbers of different patients


inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
fileList = dir(fullfile(inputFolderName_IEDtrials, '*.LFPIED.mat'));
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_0_outputs\v2\';
PatientsNum = length(fileList);

IED_percentage_vector = nan(6, PatientsNum); % 6 time chunks
xTicks = strings(1, PatientsNum);
PreOnset1 = 1;
PreOnset2 = 2;
PostOnset = 3;
PreResponse = 4;
PostResponse = 5;
PreOutcome = 6;
%%

for pt = 1:PatientsNum
% for pt = 1:1
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    xTicks(pt) = ptID;
    disp("patient: " + ptID);
    LFPIED = [inputFolderName_IEDtrials '\' ptID '.LFPIED.mat'];
    load(LFPIED);

    % preOnet1
    IEDtrials = LFPIED.IEDtrialsPreOnset1;
    IEDtrialsAcrossLFPChans = any(IEDtrials, 1); % we're gonna take the trials that an IED occured on any LFP channel
    if (~isempty(IEDtrialsAcrossLFPChans))
        IED_percentage_vector(PreOnset1, pt) = (nansum(IEDtrialsAcrossLFPChans)/length(IEDtrialsAcrossLFPChans) )*100;
    end
    clear IEDtrials IEDtrialsAcrossLFPChans

    % preOnet2
    IEDtrials = LFPIED.IEDtrialsPreOnset2;
    IEDtrialsAcrossLFPChans = any(IEDtrials, 1); % we're gonna take the trials that an IED occured on any LFP channel
    if (~isempty(IEDtrialsAcrossLFPChans))
        IED_percentage_vector(PreOnset2, pt) = (nansum(IEDtrialsAcrossLFPChans)/length(IEDtrialsAcrossLFPChans) )*100;
    end
    clear IEDtrials IEDtrialsAcrossLFPChans


    % PostOnset
    IEDtrials = LFPIED.IEDtrialsPostOnset;
    IEDtrialsAcrossLFPChans = any(IEDtrials, 1); % we're gonna take the trials that an IED occured on any LFP channel
    if (~isempty(IEDtrialsAcrossLFPChans))
        IED_percentage_vector(PostOnset, pt) = (nansum(IEDtrialsAcrossLFPChans)/length(IEDtrialsAcrossLFPChans) )*100;
    end
    clear IEDtrials IEDtrialsAcrossLFPChans

    %PreResponse
    IEDtrials = LFPIED.IEDtrialsPreResponse;
    IEDtrialsAcrossLFPChans = any(IEDtrials, 1); % we're gonna take the trials that an IED occured on any LFP channel
    if (~isempty(IEDtrialsAcrossLFPChans))
        IED_percentage_vector(PreResponse, pt) = (nansum(IEDtrialsAcrossLFPChans)/length(IEDtrialsAcrossLFPChans) )*100;
    end
    clear IEDtrials IEDtrialsAcrossLFPChans

    %PostResponse
    IEDtrials = LFPIED.IEDtrialsPostResponse;
    IEDtrialsAcrossLFPChans = any(IEDtrials, 1); % we're gonna take the trials that an IED occured on any LFP channel
    if (~isempty(IEDtrialsAcrossLFPChans))
        IED_percentage_vector(PostResponse, pt) = (nansum(IEDtrialsAcrossLFPChans)/length(IEDtrialsAcrossLFPChans) )*100;
    end
    clear IEDtrials IEDtrialsAcrossLFPChans

    %PreOutcome
    IEDtrials = LFPIED.IEDtrialsPreOutcome;
    IEDtrialsAcrossLFPChans = any(IEDtrials, 1); % we're gonna take the trials that an IED occured on any LFP channel
    if (~isempty(IEDtrialsAcrossLFPChans))
        IED_percentage_vector(PreOutcome, pt) = (nansum(IEDtrialsAcrossLFPChans)/length(IEDtrialsAcrossLFPChans) )*100;
    end
    clear IEDtrials IEDtrialsAcrossLFPChans

end % pt for


%% plotting bart chart of IEDs percentage accross patients
% Assuming IED_percentage_vector and xTicks are defined previously

% PreOnset1
    temp_vector = squeeze(IED_percentage_vector(PreOnset1, :));
    figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.99], 'Visible','off');
    bar(temp_vector, 'FaceColor', [1 0.75 0.79]);
    xticks(1:length(xTicks)); 
    xticklabels(xTicks);
    xlabel('patients', 'FontSize',14);
    ylabel('IED trials (%)','FontSize',14);
    
    % Calculate the mean of the IED percentage vector
    mean_value = mean(temp_vector);
     
    
    hold on;
    yline(mean_value, 'r--', 'LineWidth', 2); 
    hold off;
    text(length(temp_vector) - 1, mean_value+0.01, ['Mean: ', num2str(mean_value)], 'VerticalAlignment', 'bottom');
    title('IED percentage over trials for each patient - PreOnset1', FontSize=30);
    
    
    % saving
    
    set(gca, 'box', 'off', 'tickdir', 'out');
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = '0_IEDTrials_percentage_PreOnset1';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
    clear temp_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PreOnset2
    temp_vector = squeeze(IED_percentage_vector(PreOnset2, :));
    figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.99], 'Visible','off');
    bar(temp_vector, 'FaceColor', [1 0.75 0.79]);
    xticks(1:length(xTicks)); 
    xticklabels(xTicks);
    xlabel('patients', 'FontSize',14);
    ylabel('IED trials (%)','FontSize',14);
    
    % Calculate the mean of the IED percentage vector
    mean_value = mean(temp_vector);
     
    
    hold on;
    yline(mean_value, 'r--', 'LineWidth', 2); 
    hold off;
    text(length(temp_vector) - 1, mean_value+0.01, ['Mean: ', num2str(mean_value)], 'VerticalAlignment', 'bottom');
    title('IED percentage over trials for each patient - PreOnset2', FontSize=30);
    
    
    % saving
    
    set(gca, 'box', 'off', 'tickdir', 'out');
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = '0_IEDTrials_percentage_PreOnset2';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
    clear temp_vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PostOnset
    temp_vector = squeeze(IED_percentage_vector(PostOnset, :));
    figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.99], 'Visible','off');
    bar(temp_vector, 'FaceColor', [1 0.75 0.79]);
    xticks(1:length(xTicks)); 
    xticklabels(xTicks);
    xlabel('patients', 'FontSize',14);
    ylabel('IED trials (%)','FontSize',14);
    
    % Calculate the mean of the IED percentage vector
    mean_value = mean(temp_vector);
     
    
    hold on;
    yline(mean_value, 'r--', 'LineWidth', 2); 
    hold off;
    text(length(temp_vector) - 1, mean_value+0.01, ['Mean: ', num2str(mean_value)], 'VerticalAlignment', 'bottom');
    title('IED percentage over trials for each patient - PostOnset', FontSize=30);
    
    
    % saving
    
    set(gca, 'box', 'off', 'tickdir', 'out');
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = '0_IEDTrials_percentage_PostOnset';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
    clear temp_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PreResponse
    temp_vector = squeeze(IED_percentage_vector(PreResponse, :));
    figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.99], 'Visible','off');
    bar(temp_vector, 'FaceColor', [1 0.75 0.79]);
    xticks(1:length(xTicks)); 
    xticklabels(xTicks);
    xlabel('patients', 'FontSize',14);
    ylabel('IED trials (%)','FontSize',14);
    
    % Calculate the mean of the IED percentage vector
    mean_value = mean(temp_vector);
     
    
    hold on;
    yline(mean_value, 'r--', 'LineWidth', 2); 
    hold off;
    text(length(temp_vector) - 1, mean_value+0.01, ['Mean: ', num2str(mean_value)], 'VerticalAlignment', 'bottom');
    title('IED percentage over trials for each patient - PreResponse', FontSize=30);
    
    
    % saving
    
    set(gca, 'box', 'off', 'tickdir', 'out');
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = '0_IEDTrials_percentage_PreResponse';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
    clear temp_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PostResponse
    temp_vector = squeeze(IED_percentage_vector(PostResponse, :));
    figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.99], 'Visible','off');
    bar(temp_vector, 'FaceColor', [1 0.75 0.79]);
    xticks(1:length(xTicks)); 
    xticklabels(xTicks);
    xlabel('patients', 'FontSize',14);
    ylabel('IED trials (%)','FontSize',14);
    
    % Calculate the mean of the IED percentage vector
    mean_value = mean(temp_vector);
     
    
    hold on;
    yline(mean_value, 'r--', 'LineWidth', 2); 
    hold off;
    text(length(temp_vector) - 1, mean_value+0.01, ['Mean: ', num2str(mean_value)], 'VerticalAlignment', 'bottom');
    title('IED percentage over trials for each patient - PostResponse', FontSize=30);
    
    
    % saving
    
    set(gca, 'box', 'off', 'tickdir', 'out');
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = '0_IEDTrials_percentage_PostResponse';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
    clear temp_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PreOutcome
    temp_vector = squeeze(IED_percentage_vector(PreOutcome, :));
    figure('Units', 'normalized', 'Position', [0.1, 0, 0.7, 0.99], 'Visible','off');
    bar(temp_vector, 'FaceColor', [1 0.75 0.79]);
    xticks(1:length(xTicks)); 
    xticklabels(xTicks);
    xlabel('patients', 'FontSize',14);
    ylabel('IED trials (%)','FontSize',14);
    
    % Calculate the mean of the IED percentage vector
    mean_value = mean(temp_vector);
     
    
    hold on;
    yline(mean_value, 'r--', 'LineWidth', 2); 
    hold off;
    text(length(temp_vector) - 1, mean_value+0.01, ['Mean: ', num2str(mean_value)], 'VerticalAlignment', 'bottom');
    title('IED percentage over trials for each patient - PreOutcome', FontSize=30);
    
    
    % saving
    
    set(gca, 'box', 'off', 'tickdir', 'out');
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = '0_IEDTrials_percentage_PreOutcome';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
    clear temp_vector


%% debug:
% aaaaa  = squeeze(LFPIED.LFPmatPreOnset1(35,:,100:150));
% plot(aaaaa)


% aa = sum(LFPIED.IEDtrialsPreOutcome(:));
