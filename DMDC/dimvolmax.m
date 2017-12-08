function [d,vol,volB,curv,holes,t] = dimvolmax(spectrum,sigma)

    s=-log(spectrum)/sigma;

    if (size(s,1)==1) s = s'; end
    
    ss = @(t) exp(1).^(-repmat(s,1,length(t)).*repmat(t,length(s),1));    
    dim = @(t) -2*t.*sum(repmat(s,1,length(t)).*ss(t),1)./sum(ss(t),1);
    
   
    t = fminsearch(dim,1);
    
%        keyboard
    %size(t)
    x=sqrt(t);
    d=-dim(t);
    ss = ss(t);
    
    c2 = d*(d-1)*x.^(d-2) - (4*d+2)*s*x.^d + 4*(s.^2).*(x.^(d+2));
    c2 = c2.*ss;
    c2 = ((4*pi)^(d/2))*sum(c2)/2;
    
    c1 = ((4*pi)^(d/2))*sum((d*x.^(d-1) - 2*s.*x.^(d+1)).*ss) - 2*c2*x;
    c0 = ((4*pi)^(d/2))*sum(x.^d.*ss) - c1*x-c2*x.^2;

    vol=c0;
    volB=-2*c1/sqrt(pi);
    curv=3*c2/(2*pi);
    holes=1-(curv/d);
    
end

