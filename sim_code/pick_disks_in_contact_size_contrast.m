function [v_hs,contact_list] = pick_disks_in_contact_size_contrast(xval,yval,N,D,...
                                                     hs_flag,dt,L,EFB_ind,BO_ind,...
                                                     rad_EFB,rad_BO)
% Takes in the inter-particle distances at a given time-frame.
% Returns a list of particles that are within a distance "contact_dist"
% of each other, i.e., those that are overlapping

pot_strength=(0.5);
contact_list=[];
v_hs=zeros(N,2);

d_EE=2*rad_EFB;
d_BB=2*rad_BO;
d_BE=rad_EFB+rad_BO;

contact_dist=max([d_EE,d_BB,d_BE]); 


% contact_dist=min([d_EE,d_BB,d_BE]); 

% This is just a first-approximation
% depending on the identities of the particles that intrude upon the space
% of their neighbors, the value for contact_dist has to be suitably
% modified.


%%%%%% first step: resolve BO-BO contacts %%%%%%%%
contact_dist=d_BB;
sftness=0.01*contact_dist;
M_near = D; %Matrix representation for the distance between particles
[l1,l2]=find((M_near+sftness<contact_dist) & (M_near>0));

if(hs_flag)    
    for ind = 1:(length(BO_ind))
        i=BO_ind(ind);
        list = l1(l2==i);
        if ~isempty(list)
            n_list=length(list);
            for k=1:n_list
                contact_list=[contact_list,list'];
                if(ismember(list(k),BO_ind))
                    [v_hs] = modular_vhs_gen(i,list,k,xval,yval,M_near,...
                        v_hs,contact_dist,L,dt,pot_strength);
                end
            end
        end
    end
end


%%%%% second step: resolve EFB-EFB contacts %%%%%%%%

contact_dist=d_EE;
sftness=0.01*contact_dist;
M_near = D; %Matrix representation for the distance between particles
[l1,l2]=find((M_near+sftness<contact_dist) & (M_near>0));

if(hs_flag)    
    for ind = 1:(length(EFB_ind))
        i=EFB_ind(ind);
        list = l1(l2==i);
        if ~isempty(list)
            n_list=length(list);
            for k=1:n_list
                contact_list=[contact_list,list'];
                if(ismember(list(k),EFB_ind))
                    [v_hs] = modular_vhs_gen(i,list,k,xval,yval,M_near,...
                        v_hs,contact_dist,L,dt,pot_strength);
                end
            end
        end
    end
end





%%%%% third step: resolve EFB-BO contacts %%%%%%%%%


contact_dist=d_BE;
sftness=0.01*contact_dist;
M_near = D; %Matrix representation for the distance between particles
[l1,l2]=find((M_near+sftness<contact_dist) & (M_near>0));

if(hs_flag)    
    for ind = 1:(length(EFB_ind))
        i=EFB_ind(ind);
        list = l1(l2==i);
        if ~isempty(list)
            n_list=length(list);
            for k=1:n_list
                contact_list=[contact_list,list'];
                if(ismember(list(k),BO_ind))
                    [v_hs] = modular_vhs_gen(i,list,k,xval,yval,M_near,...
                        v_hs,contact_dist,L,dt,pot_strength);
                end
            end
        end
    end
end

if(hs_flag)    
    for ind = 1:(length(BO_ind))
        i=BO_ind(ind);
        list = l1(l2==i);
        if ~isempty(list)
            n_list=length(list);
            for k=1:n_list
                contact_list=[contact_list,list'];
                if(ismember(list(k),EFB_ind))
                    [v_hs] = modular_vhs_gen(i,list,k,xval,yval,M_near,...
                        v_hs,contact_dist,L,dt,pot_strength);
                end
            end
        end
    end
end
        

end

