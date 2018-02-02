function A = DoubleCenter(A)
% Double Center the Matrix A
% subtract the corresponding row and column averages from each entry and
% add back the total average

A = A - repmat(mean(A),size(A,1),1) - repmat(mean(A,2),1,size(A,2)) + mean(mean(A));

end