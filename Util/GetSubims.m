function [subims,subinds] = GetSubims(A,subimSize)

    if (nargin<2)
        subimSize=4;
    end
    n=subimSize;        %%% Subimage size
    N=size(A,1);        %%% full image size

    %%%% Collect the indices of the subimages for fast extraction

    nsubs = N-n+1;
    subinds = zeros(nsubs^2,n^2);
    for i = 1:nsubs
        for j = 1:nsubs
            [X,Y]=meshgrid((1:n)+i-1,(1:n)+j-1);
            subinds((N-n+1)*(i-1) + j,:) = sub2ind([N N],X(:),Y(:));        
        end
    end

    subims = A(subinds)';

end

