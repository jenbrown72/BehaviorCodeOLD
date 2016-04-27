function [out] = JB_dPrime(Hit,Miss,CR, FA)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    %        %D' - discrimitability index
    %when H=F, then d' =0
    %Highest Possible d' (greatest sensitivity) is 6.93, the effective limit
    %(using .99 and .01) 4.65, typical values are up to 2 and 69% correct for
    %both different and same trials corresponds to a d' of 1.0.
    
% see if any value is a 0, if so increase each by 1.

if (~any(Hit)|| ~any(Miss) || ~any(CR) || ~any(FA))  
    Hit = Hit+1;
    Miss = Miss+1;
    CR = CR+1;
    FA = FA+1;
end

%X = norminv(P,mu,sigma): mu= mean, sigma = standard deviation 

  out=norminv((Hit/( Hit+Miss)),0,1)-norminv((FA/(FA+CR)),0,1);
  
  % [X,XLO,XUP] = norminv(P,mu,sigma,pcov,alpha)

end

