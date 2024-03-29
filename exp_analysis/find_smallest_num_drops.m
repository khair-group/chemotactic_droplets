function [smallest_num_drops] = find_smallest_num_drops(pos_inp)
% Find the smallest number of drops amongst the various frames
% recorded experimentally. This is the number over which 
% statistics will be averaged.

smallest_num_drops=inf;

len_val=length(pos_inp);

for i=1:len_val
    num_drops=length(pos_inp{i});
    if(num_drops<smallest_num_drops)
        smallest_num_drops=num_drops;
    end
end


end

