function [u,l,qest,epsilon,dim] = DiffusionMap(x,nvars,k,tuningMethod)
%%% Conformally Invariant Diffusion Map (CIDM)
%%% Inputs
    %%% x       - n-by-N data set with N data points in R^n
    %%% nvars   - number of eigenfunctions/eigenvalues to compute
    %%% k       - number of nearest neighbors to use
    %%% tuningMethod - methods of autotuning epsilon
    
%%% Outputs
    %%% u       - Eigenfunctions of the generator/Laplacian
    %%% l       - Eigenvalues
    %%% qest    - Sampling measure
    %%% epsilon - scale, derived from the tuning kernel sum
    %%% dim     - intrinsic dimension, derived from the tuning kernel sum

    N = size(x,1); %number of points
    
    if (nargin<4)   tuningMethod = 0;   end
    if (nargin<3)   k=ceil(log(N)^2);   end
    if (nargin<2)   nvars = 2*k;        end

    [d,inds] = pdist2(x,x,'euclidean','smallest',k);

    %%% CkNN/CIDM normalized distances
    d = d.^2;

    %%% Tune epsilon and estimate dimension
    if (tuningMethod == 0)  
        epsilon = mean(mean(d(2:k,:)))
    else
        epsilon = tuneEpsilon(d,tuningMethod)
    end
    dim = estimateDimension(d,epsilon);

    %%% Exponential kernel
    d = exp(-d./(2*epsilon));
    d = sparse(reshape(double(inds),N*k,1),repmat(1:N,k,1),reshape(double(d),N*k,1),N,N,N*k)';
    d = (d+d')/2;

    %%% CIDM Normalization
    qest = full(sum(d,2));
    Dinv = spdiags(qest.^(-1),0,N,N);
    d = Dinv*d*Dinv;
    Dinv = spdiags((full(sum(d,2))).^(-1/2),0,N,N);
    d = Dinv*d*Dinv;
    d = (d+d')/2;

    if (nvars > 0)
    
        [u,l] = eigs(d,nvars);

        l = (diag(l));
        [~,perm] = sort(abs(l),'descend');
        l = abs(log(l(perm)))/epsilon^2;
        l=diag(l);
        u = Dinv*u(:,perm);
        
    end

    qest = qest/N/((2*pi*epsilon)^(dim/2));
    
end




