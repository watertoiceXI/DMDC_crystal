function [X,s] = metricMDS(D,p)
%%% Algorithm 1.2 Metric MDS
%%% Inputs: D - N-by-N matrix if distances D_{ij} = d(y_i,y_j)
%           p - percentage of variance to maintain (or dimension if p>1)
%%% Output: X - m-by-N matrix with columns x_i in R^m for i=1,...,N

D = D.^2;                            % Square the distances
D = -(1/2)*DoubleCenter(D);          % Double center to convert D into a Gram matrix
D = (D+D')/2;
[U,S]=eigs(D,p);
[S,sinds] = sort(diag(S),'descend');
X = diag(sqrt(S(1:p)))*U(:,sinds(1:p))';