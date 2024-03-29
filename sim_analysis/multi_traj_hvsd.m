function [hv_series_op] = multi_traj_hvsd(pos_mat,trim_len,dist_scale_fac,d_cut)
% Function calls "single_traj_hvsd.m" for each of the "n_drops"
% trajectories, and returns a cell-list as the output.
% Each entry in the cell-list a two-column matrix returned
% by "single_traj_hvsd.m"


hv_series_op={};
chk=size(pos_mat);
n_drops=chk(2);

for nval=1:n_drops
    pos_nval=squeeze(pos_mat(trim_len:end,nval,:));
    [write_val] = single_traj_hvsd(pos_nval,dist_scale_fac,d_cut);
    hv_series_op{end+1}=(write_val);
end

end

