function [X,S,U,mu] = PCA(Y,p)
%%% Algorithm 1.1 Principal Component Analysis (PCA)
%%% Inputs: Y - n-by-N matrix with columns y_i in R^n for i=1,...,N
%           p - projection dimension
%%% Output: X - m-by-N matrix with columns x_i in R^m for i=1,...,N

mu = mean(Y,2);
Y = Y - repmat(mu,1,size(Y,2));         % Center the data
[U,S,~] = svds(Y,p);                    % Compute the SVD
S=diag(S).^2;                           % Extract variance of principal components
X = U(:,1:p)'*Y;                        % Project onto first m principal components
