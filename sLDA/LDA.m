function [U,L] = LDA(data,labels,alpha,nvars)
% LDA - Linear Discriminant Analysis projection of a labeled data set
% data      n-by-N matrix representing a set of N points in R^n
% labels    1-by-N vector of categorical labels containing C unique values  

        [n,~]=size(data);
        
        if (nargin<3)
            alpha = 0;
        end
        if (nargin<4)
            nvars = 10;
        end
        nvars = min(nvars,n-1);
 
        [C,~,labels] = unique(labels);
        C = length(C);
        
        categoryMeans = zeros(n,C);
        for i = 1:C
            %%% Find all the indices with label i
            iinds = find(labels==i); 
            %%% Find the mean of the data points with the i-th label
            categoryMeans(:,i) = mean(data(:,iinds),2);
            %%% Subtract the i-th mean from each data point with label i
            data(:,iinds) = data(:,iinds) - repmat(categoryMeans(:,i),1,length(iinds));
        end

        sigmab = cov(categoryMeans');   %%% Cov matrix of the means
        sigma = cov(data');
        sigma = (1-alpha)*diag(diag(sigma))+alpha*mean(diag(sigma))*eye(size(sigma));             %%% Cov matrix of the centered data
        %sigma=eye(size(sigma));
        [U,L]=eigs(sigmab,sigma,nvars); %%% Solve the generalized eigenvalue problem
                                        %%% A*u=lambda*B*u, u=eigenvector, lambda =eigenvalue
                                        %%% A = sigmab, B=sigma
    
        [L,sinds]=sort(abs(diag(L)),'descend'); %%% Sort the eigenvalues
        U = U(:,sinds);

end



