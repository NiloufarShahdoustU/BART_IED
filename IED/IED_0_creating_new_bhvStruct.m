%% in this script I will create new bhv struct:

clear;
clc;
close all;
%%

inputFolderName = '\\155.100.91.44\d\Data\Nill\BART\bhvStruct';
fileList = dir(fullfile(inputFolderName, '*.bhvStruct.mat'));

for i = 1:length(fileList)
    data = load(fullfile(inputFolderName, fileList(i).name));
   
    fileNameParts = strsplit(fileList(i).name, '.');
    patientID = fileNameParts{1};
    disp(['Processing patient ID: ' patientID]);
    BARTbehavior(patientID)
    
end