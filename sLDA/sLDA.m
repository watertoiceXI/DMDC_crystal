function [U,L] = sLDA(data,target,slices,alpha,nvar)
% LDA - Linear Discriminant Analysis projection of a labeled data set
% data      n-by-N matrix representing a set of N points in R^n
% target    1-by-N vector containing the target function values
% slices    number of slices to cut the target into

        if (nargin<3) slices = 30;  end
        if (nargin<4) alpha = 0;    end
        if (nargin<5) nvar = 10;    end

        N=length(target);
       
        labels = zeros(1,N);    %%% create a list of sorted labels, evenly slicing the number of data points
        l = ceil(N/slices);
        for i = 1:slices-1
           labels((i-1)*l+1:i*l) = i;
        end
        labels((slices-1)*l:end) = slices;

        [~,sinds]=sort(target);         %%% find the permutation of indices (sinds) that sorts the target
        sinv(sinds)=1:length(sinds);    %%% find the inverse of that permutation
        labels = labels(sinv);          %%% apply the inverse permutation to the sorted labels
 
        [U,L]=LDA(data,labels,alpha,nvar);
        
        for i = 1:size(U,2)
            U(:,i) = U(:,i)/norm(U(:,i));
        end

end

