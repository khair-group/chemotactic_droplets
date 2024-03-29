function [x_op,theta_op] = init_pos_ori(x_inp,theta_inp,L,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
x_inp(1,:,1)=L*rand(1,N);
x_inp(1,:,2)=L*rand(1,N);
theta_inp(1,:,1)=2*pi*(rand(1,N)-0.5); % randomly distributed
x_op=x_inp;
theta_op=theta_inp;
end

