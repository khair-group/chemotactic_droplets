function [theta_op] = ret_theta_gauss_dist(N,dt)
theta_op=normrnd(0,sqrt(dt),[1,N]);
end

