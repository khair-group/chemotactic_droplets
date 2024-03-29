function [D] = dist_calc_bidisp_mix(pos_br,pos_efb,opt)
% Returns the distance matrix between BrOct-BrOct (opt=1), 
% EFB-BrOct (opt=2), and EFB-EFB (opt=3)
% droplets

switch opt
    case 1 % BrOct-BrOct
        [list_of_D] = pdist([pos_br(:,1) pos_br(:,2)],'euclidean');
        D=squareform(list_of_D);        
        
    case 2 % EFB-BrOct
        for i=1:length(pos_efb)
            for j=1:length(pos_br)
                dx=pos_br(j,1)-pos_efb(i,1);
                dy=pos_br(j,2)-pos_efb(i,2);
                D(i,j) = sqrt((dx*dx)+(dy*dy));
            end
        end          
        
    case 3 % EFB-EFB
        [list_of_D] = pdist([pos_efb(:,1) pos_efb(:,2)],'euclidean');
        D=squareform(list_of_D); 
end

end

