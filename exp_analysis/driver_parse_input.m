%%% This file processes Excel sheets (or csv files) that contain the
%%% positions of the BrOct or EFB droplets at various time-instances.
%%% It produces as output .mat files which contain 392 cells. 
%%% The number of cells equals the number of frames in the recorded
%%% video. Each cell contains the positions of the positions of droplets
%%% in that frame.

%%% The format of the Excel sheet is as follows: 
%%% x-pixels, y-pixels, frame number, drop number


dt=1/0.5;
pos_fac=3.26;

offset=0; % this is due to a feature in the video processing tool. Indicates
          % how much the starting frame is offset from 1.

broct_inp='broct_65micron_dia_efb_110micron_dia/N_BrOct_by_N_EFB_20/broct_inp.xlsx';
[broct_frames_cell,broct_num_frames] = parse_frame_inp(broct_inp,offset);

efb_inp='broct_65micron_dia_efb_110micron_dia/N_BrOct_by_N_EFB_20/efb_inp.xlsx';
[efb_frames_cell,efb_num_frames] = parse_frame_inp(efb_inp,offset);

save('broct_65micron_dia_efb_110micron_dia/N_BrOct_by_N_EFB_20/broct_frames_cell_196_392','broct_frames_cell');
save('broct_65micron_dia_efb_110micron_dia/N_BrOct_by_N_EFB_20/efb_frames_cell_196_392','efb_frames_cell');