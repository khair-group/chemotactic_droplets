function [findat] = calc_qt_mean_std_err(hv_inp,dt)
% Goes over various Heaviside values for multiple drops
% and returns the mean and standard error in the self-overlap function 
% Q(t). See Fig. 2D and associated discussion in 
% https://doi.org/10.1073/pnas.2216497120 for explanation

n_drops=length(hv_inp);
file_count=0;
minlen=inf;

for nval=1:n_drops
    nval
    file_count=file_count+1;
    templen=length(hv_inp{nval});
    if templen<minlen
        minlen=templen;
    end
end

xvals=ones(n_drops,minlen);
yvals=ones(n_drops,minlen);

for i=1:n_drops
    afile=hv_inp{i};
    xvals(i,:)=afile(1:minlen,1);
    yvals(i,:)=afile(1:minlen,2);
end


findat=ones(minlen,3);
findat(:,1)=dt*xvals(1,:)';
findat(:,2)=mean(yvals,1);
findat(:,3)=std(yvals,1)./sqrt(n_drops);

end

