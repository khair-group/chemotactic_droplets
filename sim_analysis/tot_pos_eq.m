

d_efb=0.4;
d_bo=0.4;


skip_fac=1;
dt=1e-1;
L=35.;
trim_len=1;
dist_scale_fac=1;

%%%%% specify here wheteher the code has to process BrOct based or
%%%%% EFB based positions

d_cut=0.02;
%d_cut=0.05*1*(d_bo+d_efb);
spec_prefix='tot';

%%%%%%%%%
%%%%%%%%%



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

cstats=[];

for k = 1 : length(theFiles)
% for k = 1 : 1
    k
    baseFileName = sorted_fname_list(k);
    efb_index_FileName = strcat('EFB_indices_',baseFileName);
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    efb_index_fullFileName = fullfile(theFiles(k).folder, efb_index_FileName);
    efb_indices=importdata(efb_index_fullFileName);
    
    bo_index_FileName = strcat('BO_indices_',baseFileName);
    bo_index_fullFileName = fullfile(theFiles(k).folder, bo_index_FileName);
    bo_indices=importdata(bo_index_fullFileName);
    
    efb_expr='EFB(.+?)_R';
    efb_matchStr = regexp(baseFileName,efb_expr,'match');
    n_efb=sscanf(efb_matchStr,'EFB%d_R');
    
    bo_expr='BO(.+?)_N';
    bo_matchStr = regexp(baseFileName,bo_expr,'match');
    n_bo=sscanf(bo_matchStr,'BO%d_N');
    
    fprintf(1, 'Now reading %s\n', fullFileName);
    a=importdata(fullFileName);
    
    efb_pos=a(:,efb_indices,:);
    bo_pos=a(:,bo_indices,:);
    
    %%%%%% specifying name of output file %%%%%
    tmp_op_name=sprintf('qt_dcut_%4.3f.dat',d_cut);
    tmp_op_name=strcat(spec_prefix,'N_BO',num2str(n_bo),'_N_EFB',num2str(n_efb),'_',tmp_op_name);
    op_name=fullfile(op_folder,tmp_op_name);
    
    [hv_series_op] = multi_traj_hvsd(a,trim_len,dist_scale_fac,d_cut);
    [findat] = calc_qt_mean_std_err(hv_series_op,dt);
    writematrix(findat,op_name,'delimiter','tab');
end

