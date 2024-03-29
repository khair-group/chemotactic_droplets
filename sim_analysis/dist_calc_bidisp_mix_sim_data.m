function [D] = dist_calc_bidisp_mix_sim_data(pos_br,pos_efb,L,pbc_flag,opt)
% Returns the distance matrix between BrOct-BrOct (opt=1), 
% EFB-BrOct (opt=2), and EFB-EFB (opt=3)
% droplets

switch opt
    case 1 % BrOct-BrOct    
        [D]=alt_min_img_conv(pos_br(:,1),pos_br(:,2),L,pbc_flag);
        
    case 2 % EFB-BrOct
        for i=1:length(pos_efb)
            for j=1:length(pos_br)
                dx=pos_br(j,1)-pos_efb(i,1);
                dy=pos_br(j,2)-pos_efb(i,2);
                if(pbc_flag==1)
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
            end
        end
                  
        
    case 3 % EFB-EFB
        [D]=alt_min_img_conv(pos_efb(:,1),pos_efb(:,2),L,pbc_flag);
end

end

