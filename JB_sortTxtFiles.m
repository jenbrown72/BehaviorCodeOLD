

function [] = JB_sortTxtFiles()
% JB_sortTxtFiles 
%   Script to sort incoming txt files into their appropriate mouse ID
%   folder. Make sure the below directories are appropriate for the
%   computer and files being used

%Where are the txt files being loaded and sorted to
outputdir = 'C:\Users\adesniklab\Documents\BehaviorRawData\MouseIDs\';
cd(outputdir);
idDir = dir; %get a list of all MouseID folders

%Where are the unsorted txt files stored
txtdir = 'C:\Users\adesniklab\Documents\BehaviorRawData\TxtFiles';
cd(txtdir);
txtFiles = dir('*.txt'); %get a list of all data txt files

display(' ')
disp(['Sorting: ' num2str(length(txtFiles)) ' file'])
display(' ')

for j = 1:length(txtFiles); % go through and extract data

    filed=0; %Keep track of whether we have filed away the txtFile or not
        
        if (txtFiles(j).bytes == 0) % ignor files with 0 bytes.
            filed=1;
            continue
        elseif  strncmp(txtFiles(j).name,'_',1) %ignor files with no name
            filed=1;
            continue
        else
            
    for k = 1:length(idDir); 
        
        if(any(strcmp(idDir(k).name,{'.','..'}))) % ignor . .. dir
            continue
        else
            %match mouse ID with appropriate mouse ID folder
            match = strcmp(idDir(k).name,txtFiles(j).name(1:length(idDir(k).name)));
            if match==1; %see if file already exsists in the mouse ID folder
                tempIDdir = char(strcat(outputdir,idDir(k).name)); %temp directory to mouse ID folder
                cd(tempIDdir); % go to mouse ID folder
                tempIDtxtFiles = dir('*.txt'); %get a list of current txtfiles in that folder
                cd(txtdir);
                if(any(strcmp(txtFiles(j).name, {tempIDtxtFiles.name}))); %see if the file already exsists
                    filed = 1;
                    continue; % if yes, dont copy and move on
                else
                    cd(txtdir);
                    disp(['copying: ' txtFiles(j).name , '  to  ', idDir(k).name])
                    copyfile(txtFiles(j).name,tempIDdir) %copy txt file to appropriate Mouse ID folder
                    filed=1;
                    
                end
                
            else
                
                continue
                
            end
        end
        
    end
    
    if  (filed==1)
        
        continue
        
    else % if the file has not been sufficiently dealt with give an error message
        disp(['file not dealt with: ' txtFiles(j).name ])
        prompt = {'Create new mouse ID folder:','Modify txt file name'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'',txtFiles(j).name };
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);  
        newMouseIDfolder = answer(1);
        newMouseIDname = answer(2);
        
        if ~(strcmp(newMouseIDname{1},txtFiles(j).name)) % if you updated the mouse name, rename this in txtFiles
            movefile(txtFiles(j).name, newMouseIDname{1});   
        end
        
        if isempty(newMouseIDfolder{1}) % if you didnt have to create a new folder, go on
            continue;     
        else
            
            newfolderName = char(newMouseIDfolder); %if you created a new Mouse ID folder
            
            if (any(strcmp(newMouseIDfolder, {idDir.name}))); % Check to make sure whether or not the folder already exsists
                % if folder does exsist - save this txt file in it
                currFolder = char(strcat(outputdir, newfolderName));
                cd(txtdir);
                copyfile(newMouseIDname{1},currFolder);
                disp(['copying: ' newMouseIDname{1} ' to ', newfolderName])
                filed=1;
            else
                cd(outputdir);
                str = strcat(newfolderName,'\txtFiles\data');
                mkdir(str);
                disp(['New folder made : ' (newfolderName) ]);
                addpath(newfolderName);
                tempName = char(strcat(newfolderName, '*'));
                idDirTemp = dir(tempName); % find this new folder and add it to the list of folders
                idDir = [idDir; idDirTemp];
                currFolder = char(strcat(outputdir, newfolderName));

                disp(['copying: ' newMouseIDname{1} ' to ', newfolderName])
                cd(txtdir);
                
                copyfile(newMouseIDname{1},currFolder)
                filed=1;
            end
        end
    end
    
        end
    pause(0.05);

end

cd(outputdir);
display(' ')
disp('Finished Loading Data ')
display(' ')

