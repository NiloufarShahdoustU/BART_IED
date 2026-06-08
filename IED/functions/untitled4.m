clc;
clear;
close all;

sigBHFEEGcorr = struct();
unitIDs_Pts = {};
clusterTbl = [];
moreImpulsive = false(0,1);
MNItrodeLocs_allPTs = [];
trodeLabels_allPts = {''};
anatomicalLocs_allPts = {''};
ptIDs_all = {''};

noWMpts = {};

[ptArray,bhvStruct,hazEEG] = BARTnumbers;
% looping over patients
nPts = length(ptArray);


%% plotting the flat brain once.
% visualization flag: could interpolate between pial, inflated, and flat
close all;
whichVis = '2D';
[f100,vFlat,vSurf,fFlat,fSurf]  = visSurfMNI(whichVis);
hold on;

% [f300,vFlat,vSurf,fFlat,fSurf]  = visSurfMNI(whichVis);


% Only those patients who have EEG or micros.
% pts = [1:nPts];
pts = [1:2];


% looping over patients.
for pt = fliplr(pts)



    % clearing some variables so they aren't accidentally combined across patients.
    % clear  MNItrodeLocs sEEGverts els els2 ChanMap ElecXYZMNIProj ElecAtlasProjRaw labelsFromChanMap

    % which patient?
    ptID = ptArray{pt};

    % which directory
    BARTdir = (['\\155.100.91.44\d\Data\preProcessed\BART_preprocessed\' ptID '\Data']);
    [trodeLabelsNeuro,isECoG,~] = ptTrodesBART(ptID);
    load(['\\155.100.91.44\d\Data\preProcessed\BART_preprocessed\' ptID '\Imaging\Registered\ChannelMap.mat']);

    labelsFromChanMap = trodeLabelsNeuro(isECoG);


        for els = length(labelsFromChanMap):-1:1
            
            if exist('ElecXYZMNIProj','var') % intermediate GUI variable names...
                CHANNEL = ChannelMap1(strcmp(LabelMap,labelsFromChanMap(els)));
                if ~isempty(CHANNEL) && ~isnan(CHANNEL)
                    MNItrodeLocs(els,:) = ElecXYZMNIProj(CHANNEL,:);
                    trodeLabels{els} = ElecAtlasProj{CHANNEL,1};
                else
                    MNItrodeLocs(els,:) = nan(1,3);
                    trodeLabels{els} = blanks(1);
                end
            elseif exist('ChanMap','var') % old GUI
                CHANNEL = ChanMap.ChannelMap1(strcmp(ChanMap.LabelMap,labelsFromChanMap(els)));

                if ~isempty(CHANNEL) && ~isnan(CHANNEL)
                    MNItrodeLocs(els,:) = ChanMap.ElecXYZMNIProj(CHANNEL,:);
                    trodeLabels{els} = ChanMap.ElecNMMProj{CHANNEL};
                else
                    MNItrodeLocs(els,:) = nan(1,3);
                    trodeLabels{els} = blanks(1);
                end

            else % current GUI variables.
                CHANNEL = ChannelMap1(strcmp(LabelMap,labelsFromChanMap(els)));
                if ~isempty(CHANNEL) && ~isnan(CHANNEL)
                    MNItrodeLocs(els,:) = ElecAtlasProjRaw(CHANNEL,:);
                    trodeLabels{els} = ElecNMMProj{CHANNEL};
                else
                    MNItrodeLocs(els,:) = nan(1,3);
                    trodeLabels{els} = blanks(1);
                end
                if (size(trodeLabels{els},1)==1 && size(trodeLabels{els},2)==1)
                    keyboard
                end
            end

                sEEGverts(els,:) = vFlat(els,:);
        end
        for i=1:length(labelsFromChanMap)
            scatter(sEEGverts(i,1),sEEGverts(i,2),20,rgb('red'),'filled')
        end
end
                    