function [hv_series_op] = multi_traj_hvsd(cell_mat,n_drops,trim_len,dist_scale_fac,d_cut)
% Function calls "single_traj_hvsd.m" for each of the "n_drops"
% trajectories, and returns a cell-list as the output.
% Each entry in the cell-list a two-column matrix returned
% by "single_traj_hvsd.m"


hv_series_op={};

for i=1:n_drops
    drp_i=cell_mat{i};
    pos_i=drp_i(trim_len:end,1:2);
    [write_val] = single_traj_hvsd(pos_i,dist_scale_fac,d_cut);
    hv_series_op{end+1}=(write_val);
end

end

