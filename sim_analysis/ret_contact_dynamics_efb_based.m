function [num_neighbor,nc_nn] = ret_contact_dynamics_efb_based(efb_pos,bo_pos,tag_id,L,trim_len,dt,pbc_flag,contact_dist)

%Return a timeseries pertaining to the contact dynamics of tagged EFBparticle

% setting opt=2 ensures that EFB-BrOct distances are calculated and stored 
% in a matrix D. Number of rows in the matrix is equal to the number of 
% EFB droplets.
opt=2;


%Return a timeseries pertaining to the contact dynamics of tagged particle
tsteps=length(efb_pos);
tser=dt*(trim_len:tsteps);
num_neighbor=zeros(1,length(tser));

ct=1;

tag_id

for i=trim_len:tsteps
    t_inst=i;
    pos_br=squeeze(bo_pos(t_inst,:,:));
    pos_efb=squeeze(efb_pos(t_inst,:,:));
    [D] = dist_calc_bidisp_mix_sim_data(pos_br,pos_efb,L,pbc_flag,opt);
    D_tag=D(tag_id,:);
    idx=find(D_tag<=contact_dist & D_tag>0);
    num_neighbor(ct)=length(idx);
    ct=ct+1;
end

chk_flag=(diff(num_neighbor)~=0);
nc_nn=sum(chk_flag); % number of changes in the number of neighbors


end

