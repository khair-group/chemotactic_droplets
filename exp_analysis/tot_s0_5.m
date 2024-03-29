%%% Loops through all sub-directories in the parent directory
%%% and performs Q(t) calculation

dist_scale_fac=3.26; %microns per pixel
broct_dia=65; %microns
efb_dia=110; %microns
dt=2; %seconds per frame
trim_len=1; % which time-frame to start reading data from. Setting it equal
            % to 1 means data will be read from the beginning, with no
            % frames discarded.

%%%% specification based on BrOct or EFB droplets
d_cut=0.05*(broct_dia+efb_dia);
broct_baseFileName = 'broct_inp.xlsx';
efb_baseFileName='efb_inp.xlsx';
spec_prefix='tot';


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

for k = 1 : length(folder_list) %loop through folders and perform the desired operation

    fullFolderName = fullfile(myFolder,folder_list(k));
%    fprintf(1, 'Now reading %s\n', fullFolderName);
    broct_fullFileName = fullfile(fullFolderName, broct_baseFileName);
    efb_fullFileName = fullfile(fullFolderName, efb_baseFileName);

    fprintf(1, 'Now reading %s\n', folder_list(k));
    
    tmp_op_name=sprintf('qt_dcut_%4.2f_microns.dat',d_cut);
    tmp_op_name=strcat(spec_prefix,'_',tmp_op_name);
    tmp_op_name=strcat(folder_list(k),'_',tmp_op_name);
    op_name=fullfile(fullFolderName,tmp_op_name);

%    fprintf(1, 'Output file name: %s\n', op_name);
 
    
    tmp_op_name=sprintf('multip_hvsd_dcut_%4.2f_microns.mat',d_cut);
    tmp_op_name=strcat(spec_prefix,'_',tmp_op_name);
    combined_op_name=fullfile(fullFolderName,tmp_op_name);
    
    [broct_cell_mat,n_bo] = parse_drop_based_inp(broct_fullFileName);
    [efb_cell_mat,n_efb] = parse_drop_based_inp(efb_fullFileName);
    
    n_drops=n_efb+n_bo;
    tot_cell_mat=[broct_cell_mat,efb_cell_mat];
    
    all_process_qt_inp(tot_cell_mat,n_drops,dist_scale_fac,dt,trim_len,d_cut,combined_op_name,op_name); 
    
end

