function [] = JB_sortTxtFilesAddBox()

txtdir = cd('C:\Users\adesniklab\Documents\BehaviorRawData\TxtFiles'); % directory where the sorted txt files will be copied to.
txtdirSaved = dir('*.txt'); % get a list of current saved txt files in output directory

currDir = cd('C:\Users\adesniklab\Desktop\TabletDownloads'); % define the directory for matlab to search in
tabletFolderDir = dir('Tablet*');

for i = 1:length(tabletFolderDir); % go through and extract data from each tablet
    %    i=1;
    display(' ')
    disp(['Sorting Data from: ' tabletFolderDir(i).name ''])
    display(' ')
    
    s = char(strcat(currDir,'/',tabletFolderDir(i).name,'/DataBuffer/')); % this will take you to the correct direction containing the data
    cd(s); % go into this directory
    txtFiles = dir('*.txt'); % create list of all txt files in current directory
    
    for j = 1:length(txtFiles); % for all txt files - sort
        
        if (txtFiles(j).bytes == 0) % delete files with 0 bytes.

            delete(txtFiles(j).name);
            continue
            
        elseif  strncmpi(txtFiles(j).name,'_',1) %if there is no file name

            delete(txtFiles(j).name);
            continue
            
            % now see if the txt files say which tablet they were recorded
            % on - if not, add this to the end of the file name.
            
        elseif ~(strcmp(txtFiles(j).name(end-7:end-5),'Box'))
            newname = char(strcat(txtFiles(j).name(1:end-4),'_Box',tabletFolderDir(i).name(end),'.txt'));
            movefile(txtFiles(j).name,newname);
            
        end
        
        
        
    end
    
    clear txtFiles;
    txtFiles = dir('*.txt'); % Update list of all txt files in current directory
    
    for j = 1:length(txtFiles); %copy any new files to txtFiles directory
        
        if ~(any(strcmp(txtFiles(j).name,{txtdirSaved.name})));
            
            copyfile(txtFiles(j).name,txtdir);
            disp(['copying: ' txtFiles(j).name]);
            
        else
            continue
        end
        
    end
    
end

display(' ')
disp('-----Finished Sorting and Saving Data-----');
display(' ')
cd(currDir);

end
