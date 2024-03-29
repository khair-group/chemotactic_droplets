

d_efb=0.4;
d_bo=0.4;
contact_dist=1.05*((d_efb+d_bo)/2);


skip_fac=1;
dt=1e-1;
L=35.;
trim_len=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myFolder='eq_dia';
% myFolder = '/Users/kailashamr/Desktop/kailasham/work_updates/two_dimensional_autophoretic_disks/video_gen/pe5_dat_files';
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


% Get a list of all files in the folder with the desired file name pattern.
% filePattern = fullfile(myFolder, 'L10_chemo_0_hs_1_N400_pos_ubar_0.05_D_r_0.001_omeg_10_tsteps_100000_dt_0.1000.mat'); % Change to whatever pattern you need.
filePattern = fullfile(myFolder, 'L35*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
chk={theFiles.name};
str_inp=string(chk);
sorted_fname_list=(str_inp);
%sorted_fname_list=natsortfiles(str_inp);

op_folder='eq_dia';

%%%%% process the instantaneous velocity information
%%%%% to be passed on to the image generation routine
%%%%% so that the arrow is drawn appropriately.

fullFileName_cstats_op = fullfile(theFiles(1).folder, 'no_trim_contact_stats_ens_avg_efb_drops.mat');
trace_op_name = fullfile(theFiles(1).folder, 'maxBrOct_neighb_trace_op.mat');

cstats=[];

temp_op_name='maxBrOct_neighb_trace_';

for k = 1 : length(theFiles)
% for k = 1 : 1
    k
    baseFileName = sorted_fname_list(k);
    efb_index_FileName = strcat('EFB_indices_',baseFileName);
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    efb_index_fullFileName = fullfile(theFiles(k).folder, efb_index_FileName);
    
    efb_indices=importdata(efb_index_fullFileName);
    tag_id=efb_indices;
    
    bo_index_FileName = strcat('BO_indices_',baseFileName);
    bo_index_fullFileName = fullfile(theFiles(k).folder, bo_index_FileName);
    bo_indices=importdata(bo_index_fullFileName);
    
    fprintf(1, 'Now reading %s\n', fullFileName);
    a=importdata(fullFileName);
    [n_bo,n_efb,avg_nc,std_nc,avg_max_nn,std_max_nn,nc_nn,num_neighbor_coll] = calc_save_contact_dyn_stats_ens_avg_efb_based(a,efb_indices,bo_indices,contact_dist,baseFileName,trim_len,L,dt);
    
    cstats=[cstats; n_bo n_efb avg_nc,std_nc,avg_max_nn,std_max_nn];

    temp_op_name=strcat(temp_op_name,baseFileName);
    trace_op_name=fullfile(theFiles(1).folder,temp_op_name);
    save(trace_op_name,'num_neighbor_coll');

end



