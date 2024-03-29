function [write_val] = single_traj_hvsd(pos_i,dist_scale_fac,d_cut)
% Calculates H[d_cut-r(t)] for a single drop trajectory
% where H(...) is the Heaviside function, d_cut is the cut-off distance
% and r(t) is the displacement undergone by the droplet in that timestep.

% Output is a two-column matrix. The first column is the frame number, 
% and the second column is the output of the Heaviside function at that
% frame.

r_drp_i=vecnorm(pos_i,2,2);
disp_drp_i=diff(r_drp_i);
disp_drp_i=[0;disp_drp_i];
%%% converting pixels to microns
%%% while working with simulation data, dist_scale_fac=1

disp_drp_i=disp_drp_i.*dist_scale_fac;
disp_drp_i=abs(disp_drp_i);
hvsd_arg=d_cut-disp_drp_i;
chk=heaviside(hvsd_arg);
num_frame=height(pos_i);
dum_fnum=1:num_frame;
frame_num=dum_fnum';
write_val=[frame_num chk];


end

