function [subims,subinds,Z] = GetAllSubims(A,n,m,t,timeindex)

    %%% A collection of images
    %%% Returns numSubsToGet subimages that are n-by-m-by-t
    
    [N,M,T]=size(A);        %%% full image size

    %%%% Collect the indices of the subimages for fast extraction
    
    [xs,ys]=meshgrid(0:N-n,0:M-m);
    xs=xs(:);
    ys=ys(:);
    numsubs=length(xs);
    
    [X,Y,Z]=meshgrid(1:n,1:m,1:t);
    X=X(:);Y=Y(:);Z=Z(:);
    %keyboard;
    X = repmat(X,1,numsubs) + repmat(xs',length(X),1);
    Y = repmat(Y,1,numsubs) + repmat(ys',length(Y),1);
    Z = repmat(Z,1,numsubs) + timeindex-1;

    subinds = sub2ind([N M T],X(:),Y(:),Z(:)); 

    subims  = reshape(A(subinds)',n*m*t,numsubs);

end

