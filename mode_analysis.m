function clusterIndicators = mode_analysis(q)
%%%%%%
n = size(q,2);
A = zeros(n,n);

A(:,1) = q(ceil(rand*size(q,1)),:)';

for i = 1:n-1
    [~,ind] = min(sum((q*A(:,1:i)).^2,2));
    A(:,i+1) = q(ind,:)';
end
   
FinalA = A;
lastC = cond(A);

insufficient = true;
nModes = 5;
corr_thr = 0.9;

while insufficient
    for i = 1:nModes
    
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
    insufficient = false;
    for i = 1:nModes
        for j = i:nModes
            if abs(corr(clusterIndicators(i,:),clusterIndicators(j,:)))>corr_thr
                insufficient = true;
            end
        end
    end
    nModes = nModes-1;
end

stds = std(clusterIndicators,[],2);

nonlin = abs(clusterIndicators) > 5*stds;

end
