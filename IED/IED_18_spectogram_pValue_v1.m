% In this section I am going to read all the pValues (chans, freqs, time)
% and save the pValue figures. 


%TODO: you need to run this file after when you are done with IED_17. 
%%
clear;
clc;
close all;
%%

inputFolderName = '\\155.100.91.44\d\Data\Nill\BART\IED_nonIED_freq_time_pValues';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_18_pValue_v1\';
fileList = dir(fullfile(inputFolderName, '*.pValue.mat'));
PatientsNum = length(fileList);
%%



% for pt = 1:PatientsNum
for pt = 50:50
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);


    pValueMatrix = [inputFolderName '\' ptID '.pValue.mat'];
    load(pValueMatrix);
    nChans = size(pValueMatrix,1);
    nSignal = size(pValueMatrix,3);
    t = 1:nSignal;

    % filter the pValueMatrix and put 0 for pVal>=0.05
    mask_ge_003 = pValueMatrix >= 0.05;
    mask_lt_003 = pValueMatrix < 0.05;

    pValueMatrix(mask_ge_003) = 0;
    pValueMatrix(mask_lt_003) = 1;

    % bear in mind that some values in pValueMatrix are Nans! because there
    % are some channels that do not have any IEDs!

    % iterate over channels and create figures:
    for chan=1:nChans
    % for chan=64:64
    % disp("chan: " + chan);
        %visualize the pValMatrix for that channel:
        pValForChan = squeeze(pValueMatrix(chan,:,:));
        %check if this channel has IED (if it's all nan it does not have
        %any IEDs):
        
        if(~all(isnan(pValForChan), 'all'))
            % disp('not empty');
            figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.5], 'Visible', 'off'); 
            aaa = pcolor(t, sFreqs, pValForChan);
             % aaa = surf(t, sFreqs, pValForChan);
             set(gca, 'YScale', 'log');
             view(2);
            shading flat;
            axis square;       
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title(['patient = ' char(string(ptID)) ', chan = ' char(string(chan))], 'FontSize', 14, 'FontWeight', 'bold');
            
            % Save figure
            filename = string(ptID) + '_chan' + char(string(chan));
            saveas(gcf, fullfile(outputFolderName, filename), 'pdf');
        else
            % disp('empty');
        end
        clear pValForChan 
    end

    clear mask_ge_003 mask_lt_003;

end

%%
% pValForChan = squeeze(pValueMatrix(64,:,:));
% aaa = sum(pValForChan(:));