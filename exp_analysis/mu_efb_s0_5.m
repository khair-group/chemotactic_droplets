%%% Loops through all sub-directories in the parent directory
%%% and stores the number of BrOct neighbors in contact with each
%%% EFB droplet as a time-series.

dist_scale_fac=3.26; %microns per pixel
broct_dia=65; %microns
efb_dia=110; %microns
dt=2; %seconds per frame
trim_len=1; % which time-frame to start reading data from. Setting it equal
            % to 1 means data will be read from the beginning, with no
            % frames discarded.

contact_dist=((broct_dia+efb_dia)/2)/dist_scale_fac;

%%%% specification based on BrOct or EFB droplets
broct_baseFileName = 'new_broct_frames_cell_196_392.mat';
efb_baseFileName = 'new_efb_frames_cell_196_392.mat';
ndrop_baseFileName = 'n_drops_info.dat';
spec_prefix='efb_based';


myFolder='broct_65micron_dia_efb_110micron_dia';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end

d=dir(myFolder);
dfolders = d([d(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));

chk={dfolders.name};
raw_folder_list=string(chk);
TF=contains(raw_folder_list,"N_");
folder_list=raw_folder_list(TF);

cstats=[];

for k = 1 : length(folder_list) %loop through folders and perform the desired operation
%for k = 5 : 5 %loop through folders and perform the desired operation

    neighb_list={};
    fullFolderName = fullfile(myFolder,folder_list(k));
    fprintf(1, 'Now reading %s\n', fullFolderName);
    broct_fullFileName = fullfile(fullFolderName, broct_baseFileName);
    efb_fullFileName = fullfile(fullFolderName, efb_baseFileName);
    ndrop_fullFileName = fullfile(fullFolderName, ndrop_baseFileName); 

    broct_temp=load(broct_fullFileName);
    efb_temp=load(efb_fullFileName);
    ndrop_temp=load(ndrop_fullFileName)

    broct_data=broct_temp.broct_frames_cell;
    efb_data=efb_temp.efb_frames_cell;
    n_broct = ndrop_temp(1);
    n_efb = ndrop_temp(2); 

    num_frames=length(efb_temp.efb_frames_cell);
    num_frame_list=1:num_frames;

    [smallest_efb_drops] = find_smallest_num_drops(efb_data)
    [smallest_broct_drops] = find_smallest_num_drops(broct_data)

    for nd=1:smallest_efb_drops 
        tag_id=nd
        [num_neighbor,nc_nn] = ret_efb_contact_dynamics(broct_data,efb_data,num_frames,tag_id,trim_len,contact_dist);
        neighb_list{end+1}=num_neighbor;
    end

    tmp_op_name=sprintf('neighb_info.mat');
    tmp_op_name=strcat(spec_prefix,'_',tmp_op_name);
    combined_op_name=fullfile(fullFolderName,tmp_op_name);
    save(combined_op_name,'neighb_list'); 
    
end




