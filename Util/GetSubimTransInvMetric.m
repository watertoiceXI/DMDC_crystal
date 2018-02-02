function d = GetSubimTransInvMetric(A,subimSize)

    if (nargin<2)
        subimSize=4;
    end
    n=subimSize;        %%% Subimage size
    N=size(A,1);        %%% full image size

    %%%% Collect the indices of the subimages for fast extraction
    
    pixelInds = reshape(1:subimSize^2,subimSize,subimSize);
    allInds = zeros(subimSize^2,subimSize^2);
    for i=1:subimSize       
        for j=1:subimSize
            allInds(:,(i-1)*subimSize + j) = pixelInds(:);
            pixelInds = circshift(pixelInds,1,2);
        end
        pixelInds = circshift(pixelInds,1,1);
    end
    
    
    nsubs = N-n+1;
    subinds = zeros(nsubs^2,n^2);
    for i = 1:nsubs
        for j = 1:nsubs
            [X,Y]=meshgrid((1:n)+i-1,(1:n)+j-1);
            subinds(nsubs*(i-1) + j,:) = sub2ind([N N],X(:),Y(:));        
        end
    end

    subims = A(subinds);

    nsubs = N-n+1;
    d = zeros(nsubs^2,nsubs^2);
    for i = 1:nsubs^2
        s1 = repmat(subims(i,:),[nsubs^2 1 subimSize^2]);
        s2 = reshape(subims(:,allInds),nsubs^2,subimSize^2,subimSize^2);
        d(i,:) = sqrt(squeeze(min(sum((s1-s2).^2,2),[],3)));       
    end
    

%     nsubs = N-n+1;
%     d = zeros(nsubs^2,nsubs^2);
%     for i = 1:nsubs^2
%         for j = 1:nsubs
%             s1 = A((1:n)+i-1,(1:n)+j-1);
%             s2 = A((1:n)+i-1,(1:n)+j-1);
%             d(i,j) = sqrt(min(sum((repmat(s1(:),1,subimSize^2)-s2(allInds)).^2)));       
%         end
%     end

end

