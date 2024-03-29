%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Backbone for Vicsek model simulation downloaded on 4th Oct 2022 from
%%% https://www.mathworks.com/matlabcentral/fileexchange/64208-vicsek-model-simulation
%%%
%%% RK's modifications: 
%%%
%%% [] Implementing the Heyes-Melrose algorithm so that particles behave
%%% as hard spheres.
%%%
%%% [] Introducing chemotactic interactions between particles.
%%%


L=35; % we consider a periodic box of size L x L


N_EFB=[90,150,300,445,600,700,800,850,900];

N=990*ones(1,length(N_EFB));

N_BO=N-N_EFB;

rad_EFB=0.2;
contrast_fac=2.; %size contrast
rad_BO=contrast_fac*rad_EFB; 

rho=N/(L*L);


tsteps=5e1;
dt=1e-1;

pbc_flag=1;
hs_flag=1;
chemo_coup_flag=1;



dummy=size(N);
len_N=dummy(2);

mu=0.; % magnitude of speed; same for all particles
D_r=0.; %rotational diffusivity
bet=0.; %rate constant for Poisson process, speed fluctuation
del=0.; %half-width of uniform distribution from which the speeds are picked randomly

a1_val=sqrt(0.1);
a2_val=-sqrt(0.1);
m1_val=sqrt(0.1);
m2_val=0.1*sqrt(0.1);
                               
for i=1:len_N    
    
    n=1.; % concentration varies as c~r^{-n};
            
    pos=zeros(tsteps+1,N(i),2); %stores x and y position of particles
    theta_all=zeros(tsteps+1,N(i),1);
    A=zeros(1,N(i)); %surface activity of each disk
    M=zeros(1,N(i)); %mobility coefficient for each disk
    
    
    BO_ind = 1:N_BO(i);
    A(BO_ind)=a1_val;
    M(BO_ind)=m1_val;
    EFB_ind=N_BO(i)+1:N(i);
    A(EFB_ind)=a2_val;
    M(EFB_ind)=m2_val;
    
    
    [pos,theta_all] = init_pos_ori(pos,theta_all,L,N(i));
                                      
       
    vmag0=mu;
    vel_all=zeros(tsteps+1,N(i),2); %stores components of the velocity
    vel_all(1,BO_ind,1)=vmag0*cos(theta_all(1,BO_ind,1));
    vel_all(1,BO_ind,2)=vmag0*sin(theta_all(1,BO_ind,1));
    theta_all(EFB_ind)=0.;
           
    
    [pos,theta_all,vel_all] = function_for_sim(L,N(i),mu,D_r,...
        pos,theta_all,vel_all,...
        dt,tsteps,pbc_flag,hs_flag,...
        bet,del,rad_EFB,rad_BO,chemo_coup_flag,A,M,n,BO_ind,EFB_ind);
    
   
    op_name=sprintf('sample_op_folder/L%d_chemo_%d_hs_%d_N_BO%d_N_EFB%d_R_BO_eq_%.2f_R_EFB_pos_vmag_%.2f_del_%.4f_D_r_%g_beta_%.2f_tsteps_%d_dt_%4.4f.mat',L,chemo_coup_flag,hs_flag,N_BO(i),N_EFB(i),contrast_fac,mu,del,D_r,bet,tsteps,dt);
    save(op_name,'pos','-v7.3');
    op_name=sprintf('sample_op_folder/EFB_indices_L%d_chemo_%d_hs_%d_N_BO%d_N_EFB%d_R_BO_eq_%.2f_R_EFB_pos_vmag_%.2f_del_%.4f_D_r_%g_beta_%.2f_tsteps_%d_dt_%4.4f.mat',L,chemo_coup_flag,hs_flag,N_BO(i),N_EFB(i),contrast_fac,mu,del,D_r,bet,tsteps,dt);
    save(op_name,'EFB_ind','-v7.3');
    op_name=sprintf('sample_op_folder/BO_indices_L%d_chemo_%d_hs_%d_N_BO%d_N_EFB%d_R_BO_eq_%.2f_R_EFB_pos_vmag_%.2f_del_%.4f_D_r_%g_beta_%.2f_tsteps_%d_dt_%4.4f.mat',L,chemo_coup_flag,hs_flag,N_BO(i),N_EFB(i),contrast_fac,mu,del,D_r,bet,tsteps,dt);
    save(op_name,'BO_ind','-v7.3');

end


