function [n_bo,n_efb,avg_nc,std_nc,avg_max_nn,std_max_nn,nc_nn,num_neighbor_coll] = calc_save_contact_dyn_stats_ens_avg_efb_based(a,efb_indices,bo_indices,contact_dist,baseFileName,trim_len,L,dt)

%%% Processes the simulation trajectories stored in .mat files, and
%%% calculates and saves the contact dynamics statistics

%%% In this case, tag_id is a list of EFB drop indices. This routine
%%% performs an ensemble average calculation.

efb_expr='EFB(.+?)_R';
efb_matchStr = regexp(baseFileName,efb_expr,'match');
n_efb=sscanf(efb_matchStr,'EFB%d_R');
bo_expr='BO(.+?)_N';
bo_matchStr = regexp(baseFileName,bo_expr,'match');
n_bo=sscanf(bo_matchStr,'BO%d_N');

pbc_flag=1;

length_tag=length(efb_indices);

nc_nn_list={};
max_neighb_list={};

num_neighbor_coll = [];

%%%% extract the positions of the EFB and BrOct droplets
%%%% separately. This would be useful for performing contact 
%%%% dynamics on a particular species.

efb_pos=a(:,efb_indices,:);
bo_pos=a(:,bo_indices,:);


for i=1:ceil(length_tag)
    [num_neighbor,nc_nn] = ret_contact_dynamics_efb_based(efb_pos,bo_pos,i,L,trim_len,dt,pbc_flag,contact_dist);
    nc_nn_list{end+1}=nc_nn;
    max_neighb_list{end+1}=max(num_neighbor);
    num_neighbor_coll = [num_neighbor_coll; num_neighbor];
end

avg_nc=mean(cell2mat(nc_nn_list));
std_nc=std(cell2mat(nc_nn_list));

avg_max_nn=mean(cell2mat(max_neighb_list));
std_max_nn=std(cell2mat(max_neighb_list));



end
