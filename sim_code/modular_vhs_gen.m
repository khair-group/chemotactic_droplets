function [v_hs] = modular_vhs_gen(i,list,k,xval,yval,M_near,v_hs,...
    contact_dist,L,dt,pot_strength)
%calculation of the steric repulsion velocity
    dy=yval(list(k))-yval(i);
    dx=xval(list(k))-xval(i);
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
    dist_temp=M_near(i,list(k));
    v_hs(i,1)=v_hs(i,1)+((pot_strength/dt)*(dist_temp-contact_dist)*cos(phi));
    v_hs(i,2)=v_hs(i,2)+((pot_strength/dt)*(dist_temp-contact_dist)*sin(phi));
end

