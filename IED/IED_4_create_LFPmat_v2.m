% I'm gonna create epoched LFPmat data for each subject
% author: Nill


clear;
clc;
close all;



%% loading neural and event data
 
% getting the numbers of different patients

inputFolderName = '\\155.100.91.44\d\Data\Nill\BART\bhvStruct_Nill_made';
fileList = dir(fullfile(inputFolderName, '*.bhvStruct.mat'));
outputFolderName = '\\155.100.91.44\d\Data\Nill\BART\LFPmat\'; 

for pt = 1:length(fileList)
    data = load(fullfile(inputFolderName, fileList(pt).name));
   
    fileNameParts = strsplit(fileList(pt).name, '.');
    patientID = fileNameParts{1}; 
    disp(' ');
    disp(['Processing patient ID: ' patientID]);

    ptID = patientID;

    nevList = dir(['\\155.100.91.44\d\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);
    
    if length(nevList)>1
        error('many nev files available for this patient. Please specify...')
    elseif length(nevList)<1
        error('no nev files found...')
    else
        nevFile = fullfile(nevList.folder,nevList.name);
    end
    [trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesBART_2(ptID);
    
    % load and define triggers from nevFle
    NEV = openNEV(nevFile,'overwrite');
    trigs = NEV.Data.SerialDigitalIO.UnparsedData;
    trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
    
    % load neural data
    [nevPath,nevName,nevExt] = fileparts(nevFile);
    NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));
    
    % data parameters
    selectedChans = find(isECoG);
    nChans = length(selectedChans);
    nSamps = size(NSX.Data,2);
    Fs = NSX.MetaTags.SamplingFreq;
    
    % resampling LFP at Fnew sampling frequency
    
    
    notchFilter = true;
    Fnew = 500;
    for ch = 1:nChans
        if notchFilter
            [b1,a1] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
            tmp(ch,:) = filtfilt(b1,a1,resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs));
            
            [b2,a2] = iirnotch(120/(Fnew/2),(120/(Fnew/2))/50);
            data2K(ch,:) = filtfilt(b2,a2,tmp(ch,:));
        else
            data2K(ch,:) = resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs);
        end
        clear tmp 
         
    end
    clear NSX NEV a1 a2 b1 b2
    Fs = Fnew;
    
    % common average re-referencing.
    %   data2K = remove1stPC(data2K);
    
    % timing parameters.
    pre = 2;
    post = 4;
    tSec = linspace(-pre,post,Fs*(pre+post)+1);
    
    % TF parameters
    params.fpass = [1 250]; % pick a value that's 50 higher for wavelets
    params.Fs = Fs;
    params.dBconversion = false;
    params.normalized = false;
    params.theoreticalNorm = false;
    params.baseline = true;
    
    % picking pre-bird baseline period
    if params.baseline
        baselineType = 'preCue';
        switch baselineType
            case {'preCue'}
                bP = [-2.5 -1.5];
            case {'preTask'}
                % normalizing based on the mean spectrum from a pre-task baseline epoch.
                secsPreTask = 50;
                if (ppData.Event.trigTimes(1)/3e4)<secsPreTask
                    fprintf('only %d seconds before first trigger.',ppData.Event.trigTimes(1)./3e4);
                    % TODO:: then do spectrum
                else
                    % TODO:: do spectrum.
                end
        end
    end
    
    % task parameters in chronological order..
    % There aren't any trigs that == 4
    balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
    inflateTimes = trigTimes(trigs==23);
    balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
    if length(balloonTimes)>length(inflateTimes)
        balloonTimes(end) = [];
    end
    
    
    % balloon color colormap - [yellow, orange, read, gray]
    cMap(1,:) = [1 0.9 0];
    cMap(2,:) = [1 0.5 0];
    cMap(3,:) = [1 0 0];
    cMap(4,:) = [0.5 0.5 0.5];
    
    % task parameters in chronological order..
    respTimes = trigTimes(trigs==23);
    outcomeTimes = trigTimes(trigs==25 | trigs==26);
    outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
    nTrials = length(outcomeType);
    
    
    if length(balloonType)>nTrials; balloonType=balloonType(1:end-1); end
    
    
    % TODO:: don't analyz 'NaC' trodes...
     
    % epoching data
    LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
    % plotting bank/pop responses
    for ch2 = 1:nChans
        % epoch the spectral data for each channel.
        for tt = 1:nTrials
            updateUser('finished spectral calculations',tt,50,nTrials);
            % epoch the data here [channels X samples X trials]
            LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
        end
    end

    LFPmatStruct.LFPmat = LFPmat;
    LFPmatStruct.anatomicalLocs = anatomicalLocs;
    LFPmatStruct.balloonType = balloonType;
    LFPmatStruct.balloonTimes = balloonTimes;
    LFPmatStruct.nTrials = nTrials;
    LFPmatStruct.selectedChans = selectedChans;
    LFPmatStruct.trodeLabels = trodeLabels;



    

    save([outputFolderName ptID '.LFPmat.mat'],'LFPmatStruct')

    clear ch2 tt LFPmat LFPmatStruct ch data2K balloonTimes balloonType anatomicalLocs nTrials trodeLabels 

end

%%  debug
