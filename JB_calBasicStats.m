function [meanData,stdData,countData,semData] = JB_calBasicStats(data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


meanData = nanmean(data);
    stdData = nanstd(data);
    countData = sum(~isnan(data),1);
    semData = stdData./(sqrt(countData));


end

