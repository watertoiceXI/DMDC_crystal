function [y,q,b,d,sigma,dim,vol,volB,curv,holes] = DiffusionMapping(x,k,nvars,sigma,alpha)
fprintf('GPU Diffusion Maps...                 ');
tic;
N = size(x,1); %number of points

%[ds,inds] = knn(x,x,k,0);
[ds,inds] = knnCPU(x,x,k);

if (sigma == 0)
    %     sigma = min(min(ds(:,3:32).^2))
    sigma = mean(mean(ds(:,2:2).^2))
    toc
    % sigma = mean(mean(ds(:,end).^2))/2  %%% Larger sigma if previous is small
end
zr=10;
%    sigma=ds;
if zr>0
    if (sigma == -1)
        normalize =sqrt(ds(:,end)/2);
        normalize(normalize == 0) = 1;
        ds = ds./repmat(normalize,1,k);
        ds = ds./normalize(inds);
    end
    
    ds = exp(-ds.^2/(4*abs(sigma)));      %%% Standard kernel
    %ds = max(1-ds.^2/(10*abs(sigma)),0); %%% Quadratic kernel
    
    if (sigma == -1)
        sigma = mean(mean(repmat(normalize,1,k).*normalize(inds)))^2
    end
    
    d = sparse(reshape(double(inds'),N*k,1),repmat(1:N,k,1),reshape(double(ds'),N*k,1),N,N,N*k);
    
    d = (d+d')/2;
    
    if (alpha > 0)
        Dinv = spdiags((sum(d)').^(-alpha),0,N,N);
        d = Dinv*d*Dinv;
    end
    
    Dinv = spdiags((sum(d)').^(-1/2),0,N,N);
    Drt = spdiags((sum(d)').^(1/2),0,N,N);
    d = Dinv*d*Dinv;
    
      opts.maxiter = 200;
    [y,b] = eigs(d,nvars+1,'lm',opts);
    
    b = abs(diag(b));
    [~,perm] = sort(b,'descend');
    b = b(perm);
    if (nargout>5)
        [dim,vol,volB,curv,holes] = dimvolmax(b(1:end),sigma);
    end
    b=diag(b);
    y = y(:,perm);
    q = Dinv*y;
    
    y = Drt*y;
    toc;
else
    y=0;
    q=0;
    b=0;
    d=0;
    dim=0;
    vol=0;
    volB=0;
    curv=0;
    holes=0;
     sigma=ds;
end
end


