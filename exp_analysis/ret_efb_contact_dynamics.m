function [num_neighbor,nc_nn] = ret_efb_contact_dynamics(broct_data,efb_data,num_frames,tag_id,trim_len,contact_dist)
%Return a timeseries pertaining to the contact dynamics of tagged EFBparticle

% setting opt=2 ensures that EFB-BrOct distances are calculated and stored 
% in a matrix D. Number of rows in the matrix is equal to the number of 
% EFB droplets.
opt=2;

num_neighbor=zeros(1,num_frames);

ct=1;

for i=trim_len:num_frames
    pos_br=broct_data{i};
    pos_efb=efb_data{i};
    [D] = dist_calc_bidisp_mix(pos_br,pos_efb,opt);
    D_tag=D(tag_id,:);
    idx=find(D_tag<=contact_dist & D_tag>0);
    % idx_ser{end+1}=idx;
    num_neighbor(ct)=length(idx);
    ct=ct+1;
end

chk_flag=(diff(num_neighbor)~=0);
nc_nn=sum(chk_flag); % number of changes in the number of neighbors
end

