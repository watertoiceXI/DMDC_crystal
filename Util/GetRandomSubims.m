function [subims,subinds,Z] = GetRandomSubims(A,n,m,t,numSubsToGet)

    %%% A collection of images
    %%% Returns numSubsToGet subimages that are n-by-m-by-t
    
    [N,M,T]=size(A);        %%% full image size

    %%%% Collect the indices of the subimages for fast extraction
    
    [X,Y,Z]=meshgrid(1:n,1:m,1:t);
    X=X(:);Y=Y(:);Z=Z(:);

    X = repmat(X,1,numSubsToGet);
    Y = repmat(Y,1,numSubsToGet);
    Z = repmat(Z,1,numSubsToGet);
    
    for q = 1:numSubsToGet
        X(:,q) = X(:,q) + randi([0 N-n],1);
        Y(:,q) = Y(:,q) + randi([0 M-m],1);
        Z(:,q) = Z(:,q) + randi([0 T-t],1);
    end
    
    subinds = sub2ind([N M T],X(:),Y(:),Z(:)); 

    subims  = reshape(A(subinds)',n*m*t,numSubsToGet);

end

