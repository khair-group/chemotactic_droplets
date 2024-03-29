function [] = all_process_qt_inp(cell_mat,n_drops,dist_scale_fac,dt,trim_len,d_cut,combined_op_name,op_name)

[hv_series_op] = multi_traj_hvsd(cell_mat,n_drops,trim_len,dist_scale_fac,d_cut);
[findat] = calc_qt_mean_std_err(hv_series_op,dt);

writematrix(findat,op_name,'delimiter','tab');
% save(combined_op_name,'hv_series_op');
end



