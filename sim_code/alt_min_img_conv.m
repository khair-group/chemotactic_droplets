function [D] = alt_min_img_conv(x,y,L,pbc_flag)
% This is an implementation of the minimum-image convention.
% see: https://en.wikipedia.org/wiki/Periodic_boundary_conditions

nlim=length(x);
D=zeros(nlim,nlim);
% D_x=zeros(nlim,nlim);
% D_y=zeros(nlim,nlim);

for i=1:nlim
    for j=1:nlim
        dx=x(j)-x(i);
        dy=y(j)-y(i);
        if (pbc_flag==1)
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
        end
        D(i,j) = sqrt((dx*dx)+(dy*dy));
%         D_x(i,j)=dx;
%         D_y(i,j)=dy;
    end
end

    
    
    
end





