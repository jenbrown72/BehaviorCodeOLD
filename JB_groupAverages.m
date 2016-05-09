function [AllDATA] = JB_groupAverages(AllDATA,listToAnalyse,condition)

 % JB_groupAverages : Selects already analysed data from independent mouse IDs
%for group analysis. Outputs the matrix AllDATA
 
 %   [] = JB_groupAverages(AllDATA,listToAnalyse,condition)
 %   AllDATA = []: start a new matrix called AllDATA, or input AllDATA and save to an already generated matrix
 %   listToAnalyse = Mouse ID to select from
 %   condition 1 = first data set, or 2 = second data set. e.g. control (1) and
 %   stim (2)

 % Examples:
  %   [AllDATA] = JB_groupAverages([],listToAnalyse,1);
 %   [AllDATA] = JB_groupAverages(AllDATA,listToAnalyse,2);
 %   [AllDATA] = JB_groupAverages([],{'JB022_RL_';'JB022_R_'},1);
 %
 % 

 %   Used Lists:
% eg:
% listToAnalyse = {'JB022_RL_';'SG020_RC_';'JB016_R_';'JB017_RL_';'JB017_RC_';'JB017_NT_'}
% listToAnalyse = {'JB022_R_';'JB022_NT_';'SG020_RL_';'SG020_NT_';'JB016_CE_';'JB017_R_';'JB017_L_'}
% listToAnalyse = {'SG022_RC_';'SG022_RL_'};
% listToAnalyse = {'JB022_RL_';'JB022_R_';'JB022_NT_';'SG020_RL_';'SG020_RC_';'SG020_NT_'};
% listToAnalyse = {'SG022_RC_';'SG022_RL_';'SG023_L_'};
%  listToAnalyse = {'JB017_RL_';'JB017_RC_';'JB017_R_';'JB017_NT_'};
% listToAnalyse = {'JB022_RL_';'JB022_R_';'JB022_NT_';'SG020_RL_';'SG020_RC_';'SG020_NT_';'JB016_R_';'JB016_CE_';'JB016_NT_';'JB017_RL_';'JB017_RC_';'JB017_R_';'JB017_L_';'JB017_NT_'};
% listToAnalyse = {'JB022_RL_';'JB022_R_';'SG020_RL_';'SG020_RC_';'SG020_NT_';'JB016_R_';'JB016_CE_';'JB017_RL_';'JB017_RC_';'JB017_R_';'JB017_L_';'JB017_NT_'};
%listToAnalyse = {'JB016_CE_';'JB016_L_';'JB016_NT_';'JB016_R_';'JB016_RL_'};


%concat if from different directories:  C = [AllDATA{1,1}.data AllDATA2{1,1}.data]


currDirectory = cd;

for h = 1:length(listToAnalyse) %for each mouse ID folder
    
    cd(char(strcat(currDirectory,'\',listToAnalyse{h},'\txtFiles\data')));  %go to the specific mouse folder
    
    %Add here whish functions you want to run
    [basicPropertiesToPlot,possibleAngles] = JB_basicBehaviorProperties(0,0); %run basic properties
     [DATA] = JB_plotSegments(basicPropertiesToPlot,0)
%   [DATAtrim] = JB_plotTrimmingTest(basicPropertiesToPlot,1,0);
%     [DATA] = JB_plotOptogenetics(basicPropertiesToPlot,possibleAngles,1);

    
    %save data into structure
 AllDATA{condition}.data{h} = DATA;
%     AllDATA{condition}.data{h} = basicPropertiesToPlot;
%      AllDATA{condition}.data{h}.name = listToAnalyse{h};
    cd(currDirectory);
end

end

