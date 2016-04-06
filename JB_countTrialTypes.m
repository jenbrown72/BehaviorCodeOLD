function [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(tempdata)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
Hit = sum(tempdata==1);
FA = sum(tempdata==2);
Miss = sum(tempdata==3);
CR = sum(tempdata==4);
performance =  (Hit+CR)/(length(tempdata));
licking =  (Hit+FA)/(Miss+CR+Hit+FA);

    if any(Hit+Miss==0) || any(FA+CR==0) %is if both go and no go data is being analysed, if yes calculate dprime
dprime = NaN;
    else
dprime = JB_dPrime(Hit,Miss,CR, FA);
    end

end

