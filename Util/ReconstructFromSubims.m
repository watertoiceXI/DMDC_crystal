function Arecon = ReconstructFromSubims(subims,A,subimSize)

    if (nargin<2)
        subimSize=4;
    end
    n=subimSize;        %%% Subimage size
    N=size(A,1);        %%% full image size
    Arecon = zeros(N);
    pixelCounts = zeros(N);

    %%%% Collect the indices of the subimages for fast extraction

    nsubs = N-n+1;
    subinds = zeros(nsubs^2,n^2);
    for i = 1:nsubs
        for j = 1:nsubs
            Arecon((1:n)+i-1,(1:n)+j-1) = Arecon((1:n)+i-1,(1:n)+j-1) + reshape(subims(:,(j-1)*nsubs+i),subimSize,subimSize);
            pixelCounts((1:n)+i-1,(1:n)+j-1) = pixelCounts((1:n)+i-1,(1:n)+j-1) + 1;
        end
    end
    
    Arecon = Arecon./pixelCounts;
    Arecon = Arecon';

end

