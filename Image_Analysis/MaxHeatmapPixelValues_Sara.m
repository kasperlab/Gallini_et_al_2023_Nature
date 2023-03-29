%%
%{
Version 20201116 - Sara
MatLab R2018a
This program will allow the user to select a .tif file within the same
folder and isolate the average MFI values based on X [woundMFI], which
approximates the values distance from the wound, and organizes the mean,
standard deviation and number into an excel file.
%}

clear;  % Clear workspace
clc;    % Clear command window

%% USER INPUT TO SELECT IMAGE FILE TO IMPORT

%{
    Creates a directory [files] of .tif files in current folder and
    isolates the indivdual image file names into [fileList]
%}
files       = dir('*.tif');
fileList    = struct2cell(files);
fileList    = fileList(1,:).';

%{
    Creates a dialog box listing [image] for user selection and creates
    the integer [image] based on the selection.
%}
[image tf] = listdlg('PromptString','Select a file to analyze.',...
    'SelectionMode','single',...
    'ListString',fileList,...
    'ListSize',[300,200],...
    'Name','File Selection');

%{
    If user fails to select a [image] then [tf]=0 and an error dialog is
    generated to stop the rest of the script. If [tf]=1 then a [image] was
    selected and the script proceeds.
%}
if tf == 0
    kill = errordlg('No file selected.','File Error');
end

%% IMPORT SELECTED IMAGE FILE & ANALYZE

% Import selected image file name into 2D array [imageRaw]
fileName                = files(image).name;
info                    = imfinfo(fileName);
imageRaw                = double(imread(fileName));

% Remove 0 values from [imageRaw]
imageRaw(imageRaw==0)   = NaN;

%{
Creates a 4 column 2D array [woundMFI]:
    Column 1: Relative distance from wound (X)
    Column 2: Mean of column
    Column 3: Standard deviation of column
    Column 4: Number of non-NaN elements in column
%}
woundMFI                = [(1:1:length(imageRaw))'... 
                           (mean(imageRaw,1,'omitnan'))'...
                           (std(imageRaw,1,'omitnan'))'...
                           (sum(~isnan(imageRaw),1))'];

% Saves [woundMFI] to excel
xlswrite([fileName,'.xlsx'],woundMFI(:,:));
