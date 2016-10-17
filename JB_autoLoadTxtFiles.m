function [] = JB_autoLoadTxtFiles()

mouseIDfiles = dir('*_');
currcd = cd;

display(' ')
disp(['loading files from: ' num2str(length(mouseIDfiles)) ' MouseID files'])
display(' ')

for j = 1:length(mouseIDfiles); % go through and extract data

    filed=0; %Keep track of whether we have filed away the txtFile or not
       
    tempMouseIDsFolder = char(strcat(currcd, '\', mouseIDfiles(j).name));
    
    cd(tempMouseIDsFolder);
    JB_loadTxtFile;
    
    display(' ')
disp(['loaded Txt Files from: ' (mouseIDfiles(j).name) ''])
display(' ')
    cd(currcd)
    
end