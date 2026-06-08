% In this code I am creating the data for regression of number of channels
% and ITs and RTs:

% data.
% AUTHOR: Nill

clear;
clc;
close all;
warning('off','all');
%% loading neural and event data
 
% getting the numbers of different patients

inputFolderName_IEDtrials = '\\155.100.91.44\d\Data\Nill\BART\IEDTrials';
inputFolderName_timepoints = '\\155.100.91.44\d\Data\Nill\BART\IEDdata';
fileList = dir(fullfile(inputFolderName_timepoints, '*.IEDdata.mat'));
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\IEDtrialsAcrossLFPchans_ITs_RTs_regression\';
PatientsNum = length(fileList);
balloonOnsetTime = 2000;
PlotRows = 2;
StimOnsetTime = 1001; % this is actually ~2000 in ms like the variable balloonOnsetTime. 
binWidth = 100;

preStim = StimOnsetTime-375:StimOnsetTime-125; %- 750 to -250 ms to stim onset
periStim = StimOnsetTime-125:StimOnsetTime+125; %-250 to 250 ms to stim onset
postStim = StimOnsetTime+125:StimOnsetTime+375; % 250 to 750 to stim onset

pre = 1;
peri = 2;
post = 3;
for pt = 1:PatientsNum
% for pt = 1:1
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    IEDtrials = [inputFolderName_IEDtrials '\' ptID '.IEDtrials.mat'];
    load(IEDtrials);
    IEDtrials = IEDTrialsInfo.IEDtrials;
    IEDtrialsAcrossLFPChans = any(IEDtrials, 1); % we're gonna take the trials that an IED occured on any LFP channel
    nChans = size(IEDtrials,1);
    nTrials = size(IEDtrials,2);




    IEDdata = [inputFolderName_timepoints '\' ptID '.IEDdata.mat'];
    load(IEDdata);
    
    IED_timepoints = IEDdata.IED_timepoints;
    RTs = IEDdata.RTs;
    ITs = IEDdata.ITs;


    IEDtrialsAcrossLFPChans_PrePeriPost = zeros(3,nTrials); %first row is pre, second is peri, and third is post. 

   
    numberOfIEDs = 1;

    for myTrial = 1:nTrials
        slicePre = IED_timepoints(:, preStim, myTrial);
            slicePre(isnan(slicePre)) = 0; % Replace NaN with 0
            IEDsInPreIndice = find(any(slicePre, 1));
            if sum(slicePre(:)) > numberOfIEDs
                IEDtrialsAcrossLFPChans_PrePeriPost(pre, myTrial) = sum(slicePre(:));
            end

        slicePeri = IED_timepoints(:, periStim, myTrial);
            slicePeri(isnan(slicePeri)) = 0; % Replace NaN with 0
            IEDsInPeriIndice = find(any(slicePeri, 1));
            if sum(slicePeri(:)) > numberOfIEDs
                IEDtrialsAcrossLFPChans_PrePeriPost(peri, myTrial) = sum(slicePeri(:));
            end

        slicePost = IED_timepoints(:, postStim, myTrial);
            slicePost(isnan(slicePost)) = 0; % Replace NaN with 0
            IEDsInPostIndice = find(any(slicePost, 1));
            if sum(slicePost(:)) > numberOfIEDs
                IEDtrialsAcrossLFPChans_PrePeriPost(post, myTrial) = sum(slicePost(:));
            end

        
    end

    % here I am going to create a matrix with the same size as IEDtrialsAcrossLFPChans_PrePeriPost
    % and I am going to put 1 in the places where the number of IEDs are
    % max

        IEDtrialsAcrossLFPChans_PrePeriPost_accepted = zeros(size(IEDtrialsAcrossLFPChans_PrePeriPost));

        % Loop through each column to find the max value and its indices
        for col = 1:size(IEDtrialsAcrossLFPChans_PrePeriPost, 2)
            maxNum = max(IEDtrialsAcrossLFPChans_PrePeriPost(:, col));
            if maxNum > 0
                maxIdx = find(IEDtrialsAcrossLFPChans_PrePeriPost(:, col) == maxNum);
                IEDtrialsAcrossLFPChans_PrePeriPost_accepted(maxIdx, col) = 1;
            end
        end



    IEDtrialsAcrossLFPChansRegData.IEDtrialsAcrossLFPChans_PrePeriPost = IEDtrialsAcrossLFPChans_PrePeriPost;
    IEDtrialsAcrossLFPChansRegData.ITs = ITs;
    IEDtrialsAcrossLFPChansRegData.RTs = RTs;

    save([outputFolderName ptID '.RegData.mat'],'IEDtrialsAcrossLFPChansRegData');



    clear Spikes microLabels inclChans Timepoints IEDdata IEDtrials

end % pt for



%% debug part




% aaa = sum(slicePre(:));