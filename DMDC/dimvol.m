function [dim,vol,volB,curv,holes] = dimvol(s,t)
    
    ss=exp(-s.*t);
    x=sqrt(t);
    dim = 2*x^2*sum(s.*ss)/sum(ss);
    d=dim;
    
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

