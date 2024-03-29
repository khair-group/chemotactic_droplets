function [vc_ij] = ret_chemo_vel_ij(x_i,y_i,x_j,y_j,n,A_i,M_j,L)
%Returns the chemotactic velocity of the j-th particle due to the i-th
%particle
vc_ij=zeros(1,2);
dx=x_j-x_i;
dy=y_j-y_i;

if (dx > 0.5*L)
    dx = dx - L;
end
if (dx < (-0.5*L))
    dx = dx + L;
end
% similarly for the y-coordinate
if (dy > 0.5*L)
    dy = dy - L;
end
if (dy < (-0.5*L))
    dy = dy + L;
end
phi=atan2(dy,dx);
rmag=sqrt(dx.^2+dy.^2);
vc_ij(1)=n*A_i*M_j*(rmag*cos(phi))./(rmag.^(n+2));
vc_ij(2)=n*A_i*M_j*(rmag*sin(phi))./(rmag.^(n+2));
end

