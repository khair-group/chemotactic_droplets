function [cell_mat,n_drops] = parse_drop_based_inp(fname)
%%%% reads in an Excel sheet where different numbers of 
%%%% frames have been recorded for the different droplets, 
%%%% and creates a cell-array. Each element of the cell-array
%%%% contains the {x,y,t} information of one droplet.
%%%% These elements need not be equal in length.

a=readtable(fname);
raw_inp_mat=table2array(a);
inp_mat=rmmissing(raw_inp_mat);

len=size(inp_mat);
num_lines=len(1);

sw_flag=0;
num_switch=0;

st_ptr=1; %points to the start of data-series for a droplet
cell_mat=cell(1,0);
lin_chk=0;

for i=1:(num_lines-1)
    cur_drop_ind=inp_mat(i,4);
    next_drop_ind=inp_mat(i+1,4);
    end_ptr=i;
    if (cur_drop_ind ~= next_drop_ind)
        sw_flag=1;
        num_switch=num_switch+1; %keeps track of number of droplets
        tmp_arr=inp_mat(st_ptr:end_ptr,1:3);
        st_ptr=i+1; % moves the pointer to start of data for next droplet
        cell_mat{end+1}=tmp_arr;
        lin_chk=lin_chk+length(cell_mat{num_switch});
    else
        sw_flag=0;
    end
end

tmp_arr=inp_mat(st_ptr:end,1:3);
cell_mat{end+1}=tmp_arr;

n_drops=num_switch+1;

end

