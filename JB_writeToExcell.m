function [ output_args ] = JB_writeToExcell(savingDir,filename,data,count)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

startingDir = cd;
baseFileName = strcat(savingDir,'\',filename); %save old figure
 
cd(savingDir);
xlswrite(filename,data,count)
e = actxserver('Excel.Application'); % # open Activex server
ewb = e.Workbooks.Open(baseFileName); % # open file (enter full path!)
ewb.Worksheets.Item(1).Name = DATAavg.name; % # rename 1st sheet
ewb.Save % # save to the same file
ewb.Close(false)
e.Quit

end

