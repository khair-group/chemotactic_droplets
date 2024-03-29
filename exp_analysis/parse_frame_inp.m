function [frames_cell,num_frames] = parse_frame_inp(fname,offset)
%%%% reads in an Excel sheet where different numbers of 
%%%% frames have been recorded for the different droplets, 
%%%% and creates a cell-array. The number of elements in
%%%% the cell array equals the maximum number of frames. 
%%%% Each element of the cell-array
%%%% contains the {x,y} information of one frame.

a=readtable(fname);
raw_inp_mat=table2array(a);
inp_mat=rmmissing(raw_inp_mat);


len=size(inp_mat);
num_lines=len(1);



num_frames=max(inp_mat(:,3))-offset;
% num_frames=(max(inp_mat(:,3))-min(inp_mat(:,3)))+1;


frames_cell=cell(1,num_frames);

for i=1:num_lines
    i
    chk=inp_mat(i,3)-offset;
    frames_cell{chk}=[frames_cell{chk};inp_mat(i,1:2)];
end

end

