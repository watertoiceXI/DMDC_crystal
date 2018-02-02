function [dim,ddim] = estimateDimension(d,epsilon)
% Estimates dimension as derivative of log double sum with respect to log
% epsilon
    ds1 = sum(sum(exp(-d/(2*epsilon*(1+1e-4)))));
    ds2 = sum(sum(exp(-d/(2*epsilon*(1-1e-4)))));
    %%% dlog(kernelsum)/dlog(epsilon) ~ (log(ds1)-log(ds2))/(log(e1)-log(e2))
    dim = 2*log(ds1/ds2)/log((1+1e-4)/(1-1e-4));
    %dim = log(ds1/ds2)/1e-4;
    if (nargout == 2)
        ds = sum(sum(exp(-d/(2*epsilon))));
        %%% (log(ds1)+log(ds2)-2log(ds))/(log(e1)+log(e2)-2log(e))
        %%% d^2log(kernelsum)/d^2log(epsilon) ~ epsilon *d^2log(kernelsum)/depsilon + dim
        ddim = 2*log(ds1*ds2/ds^2)/1e-8 + dim;
       
    end
    
end

