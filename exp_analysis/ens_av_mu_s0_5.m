%%% Performs the ensemble-average over the neighbor-timeseries of
%%% EFB droplets and calculates the ensemble average. 

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
nlist_baseFileName='efb_based_neighb_info.mat';

time_cutoff = 20; % minimum length of time for which max neighbor needs to persist, so that it gets counted during ensemble averaging

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

    max_neighb_list={};
    fullFolderName = fullfile(myFolder,folder_list(k));
    fprintf(1, 'Now reading %s\n', fullFolderName);
    ndrop_fullFileName = fullfile(fullFolderName, ndrop_baseFileName); 
    nlist_fullFileName = fullfile(fullFolderName, nlist_baseFileName);

    ndrop_temp=load(ndrop_fullFileName);

    n_broct = ndrop_temp(1);
    n_efb = ndrop_temp(2); 

    samp = importdata(nlist_fullFileName);
    size_samp = size(samp);
    smallest_efb_drops = size_samp(2);

    for nd=1:smallest_efb_drops 
        tag_id=nd
        num_neighbor = samp{tag_id};
        max_neighb=max(num_neighbor);
        c=nnz(num_neighbor==max_neighb);
        while 1
            if(c<time_cutoff)
                max_neighb=max_neighb-1;
                c=nnz(num_neighbor==max_neighb);
            else
                break;
            end
        end
        max_neighb_list{end+1}=max_neighb;
    end

    avg_max_nn=mean(cell2mat(max_neighb_list));
    std_max_nn=std(cell2mat(max_neighb_list));
    std_err_nn=std_max_nn/sqrt(length(cell2mat(max_neighb_list)));


    cstats=[cstats;n_broct./n_efb n_broct n_efb avg_max_nn std_err_nn];   
    
end

tmp_op_name=sprintf('twindow_%d_mu_info.mat',time_cutoff);
tmp_op_name=strcat(spec_prefix,'_',tmp_op_name);
op_name=fullfile(myFolder,tmp_op_name);

save(op_name,'cstats');




