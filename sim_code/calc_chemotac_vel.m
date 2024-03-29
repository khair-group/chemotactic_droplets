function [v_chem] = calc_chemotac_vel(xval,yval,N,n,A,M,L,pm_ind,inert_ind)
% Takes in the inter-particle distances at a given time-frame.
% Returns a list of particles that are within a distance "contact_dist"
% of each other, i.e., those that are overlapping

% d_matrix = squareform(D); %Matrix representation for the distance between particles

v_chem=zeros(N,2);

for j = 1:length(inert_ind)
    for i=1:length(pm_ind)
        pairwise_vel_ij=ret_chemo_vel_ij(xval(pm_ind(i)),yval(pm_ind(i)),xval(inert_ind(j)),yval(inert_ind(j)),n,A(pm_ind(i)),M(inert_ind(j)),L);
        v_chem(inert_ind(j),:)=v_chem(inert_ind(j),:)+pairwise_vel_ij;
    end
end  


for i=1:length(pm_ind)
    for j = 1:length(inert_ind)
        pairwise_vel_ji=ret_chemo_vel_ij(xval(inert_ind(j)),yval(inert_ind(j)),xval(pm_ind(i)),yval(pm_ind(i)),n,A(inert_ind(j)),M(pm_ind(i)),L);
        v_chem(pm_ind(i),:)=v_chem(pm_ind(i),:)+pairwise_vel_ji;
    end
end  


for i=1:length(pm_ind)
    for j = 1:length(pm_ind)
        if(i~=j)
            pairwise_vel_ji=ret_chemo_vel_ij(xval(pm_ind(j)),yval(pm_ind(j)),xval(pm_ind(i)),yval(pm_ind(i)),n,A(pm_ind(j)),M(pm_ind(i)),L);
            v_chem(pm_ind(i),:)=v_chem(pm_ind(i),:)+pairwise_vel_ji;
        end
    end
end  



end

