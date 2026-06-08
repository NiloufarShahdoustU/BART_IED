% this is writtern for each individual patient.


clear;
clc;
close all;
warning('off','all');
%% loading neural and event data

% getting the numbers of different patients
inputFolderName = '\\155.100.91.44\d\Data\Nill\BART\IEDtrialsAcrossLFPchans_ITs_RTs_regression';
fileList = dir(fullfile(inputFolderName, '*.RegData.mat'));
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_15_output_regression\';
PatientsNum = length(fileList);


pre = 1;
peri = 2;
post = 3;

for pt = 1:PatientsNum
% for pt = 1:1    
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    IEDtrialsAcrossLFPChansRegData = [inputFolderName '\' ptID '.RegData.mat'];
    load(IEDtrialsAcrossLFPChansRegData);
    IEDtrialsAcrossLFPChans = IEDtrialsAcrossLFPChansRegData.IEDtrialsAcrossLFPChans_PrePeriPost;
    RTs = IEDtrialsAcrossLFPChansRegData.RTs;
    ITs = IEDtrialsAcrossLFPChansRegData.ITs;

    % Remove reaction times that are more than 10.
    ReactTimeThreshold = 10;
    OutlierIndices = RTs >= ReactTimeThreshold;
    RTs = RTs(~OutlierIndices);
    ITs = ITs(~OutlierIndices);
    IEDtrialsAcrossLFPChans = IEDtrialsAcrossLFPChans(:, ~OutlierIndices);
    IEDtrialsAcrossLFPChans_pre = IEDtrialsAcrossLFPChans(pre,:);
    IEDtrialsAcrossLFPChans_peri = IEDtrialsAcrossLFPChans(peri,:);
    IEDtrialsAcrossLFPChans_post = IEDtrialsAcrossLFPChans(post,:);
    
    % pre
    IEDtrialsAcrossLFPChans_preNonZero_idx = find(IEDtrialsAcrossLFPChans_pre~=0);
    IEDtrialsAcrossLFPChans_preNonZero = IEDtrialsAcrossLFPChans_pre(IEDtrialsAcrossLFPChans_preNonZero_idx);
    RTs_pre = RTs(IEDtrialsAcrossLFPChans_preNonZero_idx);
    ITs_pre = ITs(IEDtrialsAcrossLFPChans_preNonZero_idx);

    % peri
    IEDtrialsAcrossLFPChans_periNonZero_idx = find(IEDtrialsAcrossLFPChans_peri~=0);
    IEDtrialsAcrossLFPChans_periNonZero = IEDtrialsAcrossLFPChans_peri(IEDtrialsAcrossLFPChans_periNonZero_idx);
    RTs_peri = RTs(IEDtrialsAcrossLFPChans_periNonZero_idx);
    ITs_peri = ITs(IEDtrialsAcrossLFPChans_periNonZero_idx);

    % post
    IEDtrialsAcrossLFPChans_postNonZero_idx = find(IEDtrialsAcrossLFPChans_post~=0);
    IEDtrialsAcrossLFPChans_postNonZero = IEDtrialsAcrossLFPChans_post(IEDtrialsAcrossLFPChans_postNonZero_idx);
    RTs_post = RTs(IEDtrialsAcrossLFPChans_postNonZero_idx);
    ITs_post = ITs(IEDtrialsAcrossLFPChans_postNonZero_idx);

    % Create figure and subplots
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.7], 'Visible', 'off'); 
    
    [ha, pos] = tight_subplot(2, 3, [0.05 0.05], [0.15 0.05], [0.05 0.05]);
    
    % First row of subplots
    axes(ha(1));
    regressionPlot_function(RTs_pre, IEDtrialsAcrossLFPChans_preNonZero);
    title('pre', 'FontSize', 14, 'FontWeight', 'bold');
    text(-0.1, 0.5, 'RT(s)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 18, 'FontWeight', 'bold');
    
    axes(ha(2));
    regressionPlot_function(RTs_peri, IEDtrialsAcrossLFPChans_periNonZero);
    title('peri', 'FontSize', 14, 'FontWeight', 'bold');
    
    axes(ha(3));
    regressionPlot_function(RTs_post, IEDtrialsAcrossLFPChans_postNonZero);
    title('post', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Second row of subplots
    axes(ha(4));
    regressionPlot_function(ITs_pre, IEDtrialsAcrossLFPChans_preNonZero);
    text(-0.1, 0.5, 'IT(s)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
    
    axes(ha(5));
    regressionPlot_function(ITs_peri, IEDtrialsAcrossLFPChans_periNonZero);
    text(0.5, -0.15, 'number of channels', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'Rotation', 0, 'FontSize', 14, 'FontWeight', 'bold');
    
    axes(ha(6));
    regressionPlot_function(ITs_post, IEDtrialsAcrossLFPChans_postNonZero);
    
    % Add the overall title
    sgtitle(['patient = ' char(string(ptID))], 'FontSize', 18, 'FontWeight', 'bold');
    
    % Adjust positions of subplots
    for ih = 1:length(ha)
        pos = get(ha(ih), 'Position');
        pos(2) = pos(2) - 0.05; % shift down by 0.03 units
        set(ha(ih), 'Position', pos);
    end

    % Save figure
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf, 'Position');
    set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
    filename = string(ptID) + '_regression';
    saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

    clear IEDtrialsAcrossLFPChans RTs ITs IEDtrialsAcrossLFPChans_pre IEDtrialsAcrossLFPChans_peri IEDtrialsAcrossLFPChans_post
    close all;
end % pt for


