function clusterIndicators = BerryCluster2(q)
%%%%%%

%q=data.q(:,1:8);

n = size(q,2);
A = zeros(n,n);

% while (cond(A)>1000)
%     A(:,ceil(rand*n)) = q(ceil(rand*size(q,1)),:)';
% end
A(:,1) = q(ceil(rand*size(q,1)),:)';

for i = 1:n-1
    [~,ind] = min(sum((q*A(:,1:i)).^2,2));
    A(:,i+1) = q(ind,:)';
end
   
FinalA = A;
lastC = cond(A);

for i = 1:5

    A(:,1) = q(ceil(rand*size(q,1)),:)';

    for i = 1:n-1
        [~,ind] = min(sum((q*A(:,1:i)).^2,2));
        A(:,i+1) = q(ind,:)';
    end

    if (cond(A)<lastC)
        lastC = cond(A);
        FinalA = A;
    end

end

clusterIndicators = q*pinv(FinalA)';

plot(clusterIndicators(:,1:n))
