function [] = IEDwaves_sEEG_minima_Lo(ptID,minPkProm,visualBadChans,CARflag)
% IEDWAVES_SEEG looks for IED traveling waves on sEEG
%
% takes patient Identifier and minimum peak prominence as inputs.
%
% now taking bad channels as well.

% EHS20210713
tic

% first finding data
topDir = ['D:\Data\' ptID '\Cadwell\'];
dataFiles = subdir([topDir 'fullDurationExport.edf']);

% loading data.
load(fullfile('D:\','Data',ptID,'Imaging','Registered','ChannelMap2.mat'))

segmentSize = 3600; % in s
Fo = 2e3; % in samples TODO:: double check that this is correct.

% looping over data files.
for fl = 1:length(dataFiles)
    % just looping over 24 segments [a priori file duration]
    for sg = 1:24 % one day's worth of recording??
        clear ptResults

        try
            if sg==1
                [hdr,D] = fastEDFRead('File',dataFiles(fl).name,'Range',[1 Fo*segmentSize]);
            else
                [hdr,D] = fastEDFRead('File',dataFiles(fl).name,'Range',[(sg-1)*Fo*segmentSize (sg)*Fo*segmentSize]);
            end
            fileFlag = true;
        catch
            fprintf('\n unable to load file segment %d, for file %d, for patient: %s using fastEDFread trying slow readedf...',sg,fl,ptID)
            %             try
            %                 if sg==1
            %                     [hdr,D] = fastEDFRead('File',dataFiles(fl).name,'Range',[1 Fo*segmentSize]);
            %                 else
            %                     [hdr,D] = fastEDFRead('File',dataFiles(fl).name,'Range',[(sg-1)*Fo*segmentSize (sg)*Fo*segmentSize]);
            %                 end
            %                 fileFlag = true;
            %             catch
            fileFlag = false;
            %             fprintf('\n unable to load file segment %d, for file %d, for patient: %s using any of the edf readers...',sg,fl,ptID)
            %             end

        end

        if fileFlag
            % getting the start time of the EDF
            fileStart = datetime(str2double(['20' hdr.startdate(7:8)]),...
                str2double(hdr.startdate(4:5)),str2double(hdr.startdate(1:2)),...
                str2double(hdr.starttime(1:2)),str2double(hdr.starttime(4:5)),...
                str2double(hdr.starttime(7:8)));

            % resampling ECoG
            Fs = 500;

            %% reordering electrode locations and label maps to fit ChannelMap2
            % taking the union of the channelMaps.
            [CMunion,CMcommonvalues,CMrepeatedValues] = union(ChannelMap1,ChannelMap2,'stable');

            % just building an interpreter matrix in which the first column is the
            % electrode identitites in ChannelMap1 and the second columd is the
            % electrode identities in ChannelMap2.
            CMtranslationMatrix = [ChannelMap1(CMcommonvalues) ChannelMap2(CMcommonvalues)];
            CMtranslationMask = ~isnan(CMtranslationMatrix(:,2));
            CMtranslationMatrix = CMtranslationMatrix(CMtranslationMask,:);

            [~,CM1order] = sort(CMtranslationMatrix(:,1));
            [~,CM2order] = sort(CMtranslationMatrix(:,2));

            tmp = LabelMap(CMcommonvalues);
            tmp = tmp(CMtranslationMask);
            macroLabels = tmp(CM1order);
            clear tmp

            % sorting ElecXYZProj to be in CM2 order.
            ElecXYZProj = ElecXYZProj(CM2order,:);
            ElecXYZProj(isnan(ElecXYZProj(:,1)),:) = [];

            % and now sorting the data
            tmpdats = resample(double(D),Fs,Fo)'; %D(:,goodTrodes)

            nChans = length(macroLabels);
            % [channels X samples] @ 500 Hz

            % need to remove noisy channels first. The major issue appears
            % to be 60 Hz, so will start there.
            [b,a] = butter(4,[20 40]./(Fs/2));
            [b2,a2] = butter(4,[1 50]./(Fs/2));
            for ch = 1:nChans
                updateUser('filtering channels',ch,20,nChans)

                LoSignal(ch,:) = filtfilt(b2,a2,tmpdats(ch,:));

                IEDsignal(ch,:) = filtfilt(b,a,tmpdats(ch,:));

            end % looping over channels.

            % common average rereferencing
            if CARflag
                ECoG = remove1stPC(LoSignal);
            else
                ECoG = LoSignal;
            end
            clear tmpdats D tmpECoG                % not the most efficient, but hopefully no memory errors

            %% IED detection: detecting peaks in beta power.
            % [20191001]:: This min peak prominence measure works great for
            % detecting IEDs from beta power.
            % for doing the full spectrogram.
            %             IEDsignal = smooth(squeeze(mean(dF.Sft(closestTrodeToUMA,dF.fHz>20 & dF.fHz<40,:)))...
            %                 ./max(squeeze(mean(dF.Sft(closestTrodeToUMA,dF.fHz>20 & dF.fHz<40,:))))*100,50);

            % [20210907] OK, in order to try and detect fewer, more
            % porminent IEDs, I'm going to detect peaks across all
            % channels, and retain those IEDs that recruit several
            % channels.

            % just filtering in the beta band for now. cons: may lead to
            % spurious detections.

            % detecting IEDs from peaks in median beta power.
            %             minPkProm = 15;
            [IEDmaximaVals,IEDmaxima,IEDwidths,~] = findpeaks(mean(smoothdata(IEDsignal','gaussian',20)'),'MinPeakProminence',minPkProm);

            % TODO:: reject large channels based on maxima here?
            maxOutliers = outliers(IEDmaximaVals,[0 75],2);
            IEDmaxima(maxOutliers) = [];
            IEDmaximaVals(maxOutliers) = [];
            IEDwidths(maxOutliers) = [];

            % excluding IEDs within a half second on either side of each
            % data segment.
            IEDmaxima((IEDmaxima-(Fs/2))<=0) = [];
            IEDmaxima((IEDmaxima+(Fs/2))>length(ECoG)) = [];

            % TODO:: realign data to its maximum decrease or voltage
            % minima, as I'm getting two peaks in minima... [20210813]
            % [20210907] actually, this only matters for comparing among
            % groups of IEDs, so not essential right now.

            %% looping over IED detections.
            for ied = 1:length(IEDmaxima)

                % looping over channels to filter data
                % [20210727] as of today. I'm going to apply different
                % amounts of smoothing to each type of data. 10 samples for
                % the full BW data and 25 samples for the beta peaks.
                for ch2 = 1:nChans
                    % full BW data
                    IEDdata(ied).data(ch2,:) = smoothdata(ECoG(ch2,IEDmaxima(ied)-(Fs/2):IEDmaxima(ied)+(Fs/2))...
                        -mean(ECoG(ch2,IEDmaxima(ied)-(Fs/2):IEDmaxima(ied)+(Fs/2)),2)','gaussian',10)';

                    % skipping beta measures for now.
                    %                     % beta
                    %                     IEDdata(ied).beta(ch2,:) = abs(hilbert(IEDsignal(ch2,IEDmaxima(ied)-(Fs/2):IEDmaxima(ied)+(Fs/2))));
                    %                     % beta phase
                    %                     IEDdata(ied).phiBeta(ch2,:) = angle(hilbert(IEDsignal(ch2,IEDmaxima(ied)-(Fs/2):IEDmaxima(ied)+(Fs/2))));

                    % finding timings of IED minima and beta power peaks
                    try
                        [IEDmins(ch2),IEDmindices(ch2)] = min(IEDdata(ied).data(ch2,(Fs/2)-round(IEDwidths(ied)):(Fs/2)+round(IEDwidths(ied))));
                        % [betaMaxes(ch2),betaMaxdices(ch2)] = max(IEDdata(ied).beta(ch2,(Fs/2)-round(IEDwidths(ied)):(Fs/2)+round(IEDwidths(ied))));
                    catch
                        [IEDmins(ch2),IEDmindices(ch2)] = min(IEDdata(ied).data(ch2,150:351));
                        % [betaMaxes(ch2),betaMaxdices(ch2)] = max(IEDdata(ied).beta(ch2,200:301));
                    end

                    % Note the grid code using vector calculus methods is in ECoG_IEDs
                end

                % grouping data
                IEDdata(ied).file = dataFiles(fl).name;
                IEDdata(ied).segment = sg;
                IEDdata(ied).segmentSize = segmentSize;
                IEDdata(ied).minPkProminence = minPkProm;
                IEDdata(ied).mindices = IEDmindices;
                % IEDdata(ied).betaMaxdices = betaMaxdices;


                %% distances
                % [20210727] OK, I'm going to look at distance from largest
                % beta peak.
                % [maxVal,sourceChan] = max(betaMaxes);
                % [20210907] switching to largest minimum
                [~,sourceChan] = min(IEDmins);
                [~,channelMindices] = sort(IEDmins);

                % Distances from stim to response electrodes
                LR = '5';
                % euclidean distances
                D_euclidean = sqrt(sum((ElecXYZProj.^2) - repmat(ElecXYZProj(sourceChan,:),nChans,1),2));
                %[~,ascendingDidcs_euclidean] = sort(D_euclidean);

                % number of connections
                numConnsFile = fullfile('D:\Data\',ptID,'Imaging','Probtrackx',sprintf('FA_%s_LR_%s',ptID,LR),'fdt_network_matrix');
                waytotalFile = fullfile('D:\Data\',ptID,'Imaging','Probtrackx',sprintf('FA_%s_LR_%s',ptID,LR),'waytotal');
                load(numConnsFile)
                load(waytotalFile)

                % need to adjust for number of good channels.
                if strcmp(ptID,'CS202007')
                    D_probCon = log10(fdt_network_matrix(sourceChan,1:end/2)./waytotal(1:end/2)');
                else
                    D_probCon = log10(fdt_network_matrix(sourceChan,:)./waytotal');
                end
                D_probCon(isinf(D_probCon)) = 0;
                %[~,ascendingDidcs_numCons] = sort(D_numCons);

                % path length.
                pathLenFile = fullfile('D:\Data\',ptID,'Imaging','Probtrackx',sprintf('FA_%s_LR_%s',ptID,LR),'fdt_network_matrix_lengths');
                load(pathLenFile)

                % need to adjust for number of good channels.
                if strcmp(ptID,'CS202007')
                    D_pathLen = fdt_network_matrix_lengths(sourceChan,1:end/2);
                else
                    D_pathLen = fdt_network_matrix_lengths(sourceChan,:);
                end
                %[~,ascendingDidcs_pathLen] = sort(D_pathLen);

                % now redefining channel numbers and labels.
                inclChans = ChannelMap2(~isnan(ChannelMap2));
                nInclChans = length(inclChans);
                chanLabels = LabelMap(~isnan(ChannelMap2));

                %%
                % updating user every 50 IEDs.
                updateUser('processed IED',ied,50,length(IEDmaxima))

                % starting a structure for data storage
                ptResults.pt = ptID;
                ptResults.chanLabel{ied} = chanLabels(sourceChan);
                ptResults.chanNumber(ied) = sourceChan;
                ptResults.fileNumber = fl;
                ptResults.segmentNumber = sg;
                ptResults.edfHeader = hdr;


                %% runniing models.
                % linear linear models for distance versus IED mins.
                % times 2 to convert to seconds.
                %                 betaMaxLM = fitlm(D_euclidean,betaMaxdices*2);      %,'Exclude',maxDescentLocs==min(maxDescentLocs),'Weights',ampWeights);
                minsLM = fitlm(D_euclidean,IEDmindices*2);                      %,'Exclude',maxesLocs==min(maxesLocs,'Weights',ampWeights));

                %                 betaMaxLM_numCons = fitlm(D_probCon,betaMaxdices*2);      %,'Exclude',maxDescentLocs==min(maxDescentLocs),'Weights',ampWeights);
                minsLM_numCons = fitlm(D_probCon,IEDmindices*2);                      %,'Exclude',maxesLocs==min(maxesLocs,'Weights',ampWeights));

                %                 betaMaxLM_pathLen = fitlm(D_pathLen,betaMaxdices*2);      %,'Exclude',maxDescentLocs==min(maxDescentLocs),'Weights',ampWeights);
                minsLM_pathLen = fitlm(D_pathLen,IEDmindices*2);                      %,'Exclude',maxesLocs==min(maxesLocs,'Weights',ampWeights));

                % hierarchical models
                %                 distTbl = table(2*betaMaxdices',2*IEDmindices',D_euclidean,D_probCon',D_pathLen',...
                %                     'VariableNames',{'beta','minima','euclidean','probCon','pathLen'});
                %                 betaLMM = fitglme(distTbl,'beta ~ 1 + probCon*pathLen');
                %                 minsLMM = fitglme(distTbl,'minima ~ 1 + probCon*pathLen');


                %% saving models.
                % LMs
                ptResults.FileStartFromEDF = fileStart;
                ptResults.detectionTimesSec_ReFileStart = (segmentSize*(sg-1))+(IEDmaxima./Fs);
                %                 ptResults.betaMaxLM(ied).LM = betaMaxLM;
                ptResults.minsLM(ied).LM = minsLM;
                %                 ptResults.betaMaxLM_probCon(ied).LM = betaMaxLM_numCons;
                ptResults.minsLM_probCon(ied).LM = minsLM_numCons;
                %                 ptResults.betaMaxLM_pathLen(ied).LM = betaMaxLM_pathLen;
                ptResults.minsLM_pathLen(ied).LM = minsLM_pathLen;
                % LMMs
                %                 ptResults.betaMaxLMM(ied).LM = betaLMM;
                %                 ptResults.minsLMM(ied).LM = minsLMM;
                % IED data.
                ptResults.IEDdata(ied) = IEDdata(ied);

                % color map
                delayMap = cool(nChans);

                %% visualize data for each IED

                % plotting
                if (minsLM.Coefficients{2,3}>=2.5 || (minsLM_numCons.Coefficients{2,3}<=-2.5 && minsLM_pathLen.Coefficients{2,3}>=2.5))
                    plotIED = true;
                else
                    plotIED = false;
                end
                if plotIED
                    figure(ied)
                    ax1 = subplot(4,2,1);
                    hold on
                    for ch3 = 1:nChans
                        plot(linspace(0,1,length(IEDdata(ied).data)),IEDdata(ied).data(channelMindices(ch3),:),'color',delayMap(channelMindices(ch3),:),'linewidth',0.5)
                    end
                    hold off
                    axis tight
                    title(sprintf('patient: %s, file: %d, segment: %d, IED number: %d',ptID,fl,sg,ied))
                    ylabel('voltage (uV)')
                    xlabel('time (s)')

                    %                     ax2 = subplot(4,2,2);
                    %                     hold on
                    %                     for ch4 = 1:nChans
                    %                         plot(linspace(0,1,length(IEDdata(ied).beta)),IEDdata(ied).beta(ch3,:),'color',delayMap(betaMaxdices(ch3),:),'linewidth',0.5)
                    %                     end
                    %                     hold off
                    %                     axis tight
                    %                     ylabel('voltage (uV)')
                    %                     xlabel('time (s)')

                    ax3 = subplot(4,2,3);
                    mP = plotAdded(minsLM);
                    mP(1).Marker = '.';
                    mP(1).MarkerEdgeColor = 'k';
                    legend off
                    axis tight
                    ylabel('timing of IED minima')
                    xlabel('euclidean distance')
                    if minsLM.Coefficients{2,4}<0.05 && ~isinf(minsLM.Coefficients{2,3}) && sign(minsLM.Coefficients{2,3})>0 %only looking at positive effects.
                        title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',minsLM.DFE,minsLM.Coefficients{2,3},minsLM.Coefficients{2,4}))
                        ptResults.sigMinModel(ied) = true;
                        ptResults.signMinModel(ied) = sign(minsLM.Coefficients{2,3});
                    else
                        title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',minsLM.DFE,minsLM.Coefficients{2,3},minsLM.Coefficients{2,4}))
                        ptResults.sigMinModel(ied) = false;
                        ptResults.signMinModel(ied) = sign(minsLM.Coefficients{2,3});
                    end

                    %                     ax4 = subplot(4,2,4);
                    %                     bP = plotAdded(betaMaxLM);
                    %                     bP(1).Marker = '.';
                    %                     bP(1).MarkerEdgeColor = 'k';
                    %                     legend off
                    %                     axis tight
                    %                     ylabel('timing of beta maxima')
                    %                     xlabel('euclidean distance')
                    %                     if betaMaxLM.Coefficients{2,4}<0.05 && ~isinf(betaMaxLM.Coefficients{2,3}) && sign(betaMaxLM.Coefficients{2,3})>0
                    %                         title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',betaMaxLM.DFE,betaMaxLM.Coefficients{2,3},betaMaxLM.Coefficients{2,4}))
                    %                         ptResults.sigBetaModel(ied) = true;
                    %                         ptResults.signBetaModel(ied) = sign(betaMaxLM.Coefficients{2,3});
                    %                     else
                    %                         title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',betaMaxLM.DFE,betaMaxLM.Coefficients{2,3},betaMaxLM.Coefficients{2,4}))
                    %                         ptResults.sigBetaModel(ied) = false;
                    %                         ptResults.signBetaModel(ied) = sign(betaMaxLM.Coefficients{2,3});
                    %                     end

                    ax5 = subplot(4,2,5);
                    mP = plotAdded(minsLM_numCons);
                    mP(1).Marker = '.';
                    mP(1).MarkerEdgeColor = 'k';
                    legend off
                    axis tight
                    ylabel('timing of IED minima')
                    xlabel('number of connections')
                    if minsLM_numCons.Coefficients{2,4}<0.05 && ~isinf(minsLM_numCons.Coefficients{2,3}) && sign(minsLM_numCons.Coefficients{2,3})>0
                        title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',minsLM_numCons.DFE,minsLM_numCons.Coefficients{2,3},minsLM_numCons.Coefficients{2,4}))
                        ptResults.sigMinModel_numCons(ied) = true;
                        ptResults.signMinModel_numCons(ied) = sign(minsLM_numCons.Coefficients{2,3});
                    else
                        title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',minsLM_numCons.DFE,minsLM_numCons.Coefficients{2,3},minsLM_numCons.Coefficients{2,4}))
                        ptResults.sigMinModel_numCons(ied) = false;
                        ptResults.signMinModel_numCons(ied) = sign(minsLM_numCons.Coefficients{2,3});
                    end

                    %                     ax6 = subplot(4,2,6);
                    %                     bP = plotAdded(betaMaxLM_numCons);
                    %                     bP(1).Marker = '.';
                    %                     bP(1).MarkerEdgeColor = 'k';
                    %                     legend off
                    %                     axis tight
                    %                     ylabel('timing of beta maxima')
                    %                     xlabel('number of connections')
                    %                     if betaMaxLM_numCons.Coefficients{2,4}<0.05 && ~isinf(betaMaxLM_numCons.Coefficients{2,3}) && sign(betaMaxLM_numCons.Coefficients{2,3})>0
                    %                         title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',betaMaxLM_numCons.DFE,betaMaxLM_numCons.Coefficients{2,3},betaMaxLM_numCons.Coefficients{2,4}))
                    %                         ptResults.sigBetaModel_numCons(ied) = true;
                    %                         ptResults.signBetaModel_numCons(ied) = sign(betaMaxLM_numCons.Coefficients{2,3});
                    %                     else
                    %                         title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',betaMaxLM_numCons.DFE,betaMaxLM_numCons.Coefficients{2,3},betaMaxLM_numCons.Coefficients{2,4}))
                    %                         ptResults.sigBetaModel_numCons(ied) = false;
                    %                         ptResults.signBetaModel_numCons(ied) = sign(betaMaxLM_numCons.Coefficients{2,3});
                    %                     end

                    ax7 = subplot(4,2,7);
                    mP = plotAdded(minsLM_pathLen);
                    mP(1).Marker = '.';
                    mP(1).MarkerEdgeColor = 'k';
                    legend off
                    axis tight
                    ylabel('timing of IED minima')
                    xlabel('path length')
                    if minsLM_pathLen.Coefficients{2,4}<0.05 && ~isinf(minsLM_pathLen.Coefficients{2,3}) && sign(minsLM_pathLen.Coefficients{2,3})>0
                        title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',minsLM_pathLen.DFE,minsLM_pathLen.Coefficients{2,3},minsLM_pathLen.Coefficients{2,4}))
                        ptResults.sigMinModel_pathLen(ied) = true;
                        ptResults.signMinModel_pathLen(ied) = sign(minsLM_pathLen.Coefficients{2,3});
                    else
                        title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',minsLM_pathLen.DFE,minsLM_pathLen.Coefficients{2,3},minsLM_pathLen.Coefficients{2,4}))
                        ptResults.sigMinModel_pathLen(ied) = false;
                        ptResults.signMinModel_pathLen(ied) = sign(minsLM_pathLen.Coefficients{2,3});
                    end

                    %                     ax8 = subplot(4,2,8);
                    %                     bP = plotAdded(betaMaxLM_pathLen);
                    %                     bP(1).Marker = '.';
                    %                     bP(1).MarkerEdgeColor = 'k';
                    %                     legend off
                    %                     axis tight
                    %                     ylabel('timing of beta maxima')
                    %                     xlabel('path length')
                    %                     if betaMaxLM_pathLen.Coefficients{2,4}<0.05 && ~isinf(betaMaxLM_pathLen.Coefficients{2,3}) && sign(betaMaxLM_pathLen.Coefficients{2,3})>0
                    %                         title(sprintf('significant: (t(%d) = %.2f, p = %.2f)',betaMaxLM_pathLen.DFE,betaMaxLM_pathLen.Coefficients{2,3},betaMaxLM_pathLen.Coefficients{2,4}))
                    %                         ptResults.sigBetaModel_pathLen(ied) = true;
                    %                         ptResults.signBetaModel_pathLen(ied) = sign(betaMaxLM_pathLen.Coefficients{2,3});
                    %                     else
                    %                         title(sprintf('NOT significant: (t(%d) = %.2f, p = %.2f)',betaMaxLM_pathLen.DFE,betaMaxLM_pathLen.Coefficients{2,3},betaMaxLM_pathLen.Coefficients{2,4}))
                    %                         ptResults.sigBetaModel_pathLen(ied) = false;
                    %                         ptResults.signBetaModel_pathLen(ied) = sign(betaMaxLM_pathLen.Coefficients{2,3});
                    %                     end

                    % saving figure
                    halfMaximize(ied,'left')
                    if CARflag
                        saveas(ied,fullfile('D:\','Figs','Elliot','IEDs','sEEG','car',sprintf('%s_file%d_segment%d_IED%d_distanceFromChannel%d.pdf',ptID,fl,sg,ied,sourceChan)))
                    else
                        saveas(ied,fullfile('D:\','Figs','Elliot','IEDs','sEEG','noreref',sprintf('%s_file%d_segment%d_IED%d_distanceFromChannel%d.pdf',ptID,fl,sg,ied,sourceChan)))
                    end
                    close(ied)
                end

            end % loop over ieds
            % TODO :: save results
            try
                save(['D:\Data\Elliot\sEEG_IED_stats\' ptID '\' sprintf('%s_file%d_segment%d_linearAndMixedEffectsModels.mat',ptID,fl,sg)],'ptResults', '-v7.3')
            catch
                fprintf('\nno results for patient %s, file %d, segment %d...',ptID,fl,sg)
            end
        end
        %         made this try/catch into the above fileFlag statement.
        %         catch
        %             fprintf('\n could not read data for segment %d in file %d\n',sg,fl)
        %             pause(10)
        %         end
    end % looping over segments
end % looping over files.


A = toc;
fprintf('\nCode took %.2f hours for patient %s\n',A/3600,ptID)

%
% if isfield('ptResults','sigMinModel')
%     % IED minima results.
%     fprintf('\n%d of %d IEDs significant using timing of IED minima models (Euclidean Distance)\n',sum(ptResults.sigMinModel),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of IED minima models (Number of Connections)\n',sum(ptResults.sigMinModel_numCons),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of IED minima models (Path Length)\n',sum(ptResults.sigMinModel_pathLen),length(IEDdata))
%
%     fprintf('\n%d of %d IEDs significant using timing of IED minima models (Euclidean Distance & Number of Connections)\n',sum(ptResults.sigMinModel),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of IED minima models (Euclidean Distance & Path Length)\n',sum(ptResults.sigMinModel & ptResults.sigMinModel_numCons),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of IED minima models (Number of Connections & Path Length)\n',sum(ptResults.sigMinModel_numCons & ptResults.sigMinModel_pathLen),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of IED minima models (All three dsitance metrics)\n',sum(ptResults.sigMinModel & ptResults.sigMinModel_numCons & ptResults.sigMinModel_pathLen),length(IEDdata))
%
%     fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
%
%     % beta maxima results.
%     fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Euclidean Distance)\n',sum(ptResults.sigBetaModel),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Number of Connections)\n',sum(ptResults.sigBetaModel_numCons),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Path Length)\n',sum(ptResults.sigBetaModel_pathLen),length(IEDdata))
%
%     fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Euclidean Distance & Number of Connections)\n',sum(ptResults.sigBetaModel),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Euclidean Distance & Path Length)\n',sum(ptResults.sigBetaModel & ptResults.sigBetaModel_numCons),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (Number of Connections & Path Length)\n',sum(ptResults.sigBetaModel_numCons & ptResults.sigBetaModel_pathLen),length(IEDdata))
%     fprintf('\n%d of %d IEDs significant using timing of beta power maxima models (All three dsitance metrics)\n',sum(ptResults.sigBetaModel & ptResults.sigBetaModel_numCons & ptResults.sigBetaModel_pathLen),length(IEDdata))
%
%     fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
%     %     fprintf('\n     These electrodes had highest beta power:\n\n %s',cell2mat(ptResults.chanLabel))
%
% else
%     fprintf('\nno results for patient %s',ptID)
% end









%% code for looking at different types of distance.
%                 LR = '5';
%                 switch whichDist
%                     case {'euclidean'}
%                         D = sqrt(sum(repmat(ElecXYZProj(closestTrodeToUMA,:),nChans,1)-ElecXYZProj.^2,2));
%                         [~,ascendingDidcs] = sort(D);
%                     case {'numconnections'}
%                         % TODO:: add structural connetivtiy distances...
%                         numConnsFile = fullfile('D:\Data\',ptID,'Imaging','Probtrackx',sprintf('FA_%s_LR_%s',ptID,LR),'fdt_network_matrix');
%                         load(numConnsFile)
%
%                         % need to adjust for number of good channels.
%                         D = fdt_network_matrix_lengths(closestTrodeToUMA,:);
%                         [~,ascendingDidcs] = sort(D);
%
%                     case {'pathlength'}
%                         pathLenFile = fullfile('D:\Data\',ptID,'Imaging','Probtrackx',sprintf('FA_%s_LR_%s',ptID,LR),'fdt_network_matrix_lengths');
%                         load(pathLenFile)
%
%                         % need to adjust for number of good channels.
%                         D = fdt_network_matrix_lengths(:,:);
%                         [~,ascendingDidcs] = sort(D);
%
%                         % now redefining channel numbers and labels.
%                         inclChans = ChannelMap2(~isnan(ChannelMap2));
%                         nInclChans = length(inclChans);
%                         chanLabels = LabelMap(~isnan(ChannelMap2));
%
%                         % TODO:: May have to remove bad channels still...
%                 end


