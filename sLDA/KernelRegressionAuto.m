function [Yrecon,insampleRMSE,outsampleRMSE,insampleCORR,outsampleCORR, c] = KernelRegressionAuto(X,Y,trainingInds,reg,k,k2)
%   X - n-by-N matrix representing N data points in R^n to be mapped to Y
%   Y - 1-by-N matrix
%   trainingInds - indices of the data points to be used for training,
%                   should be a vector of numbers between 1 and N
%   reg          - regularization parameter for the regression, input -1 to autotune (increases runtime by a factor of 10)
%   k            - number of nearest neighbors in the kernel matrix, should be as high as possible but smaller k will save memory
%   k2           - number of nearest neighbors for determining the bandwidth
%   Yrecon       - reconstruction of Y from the regression
%   c            - third letter of the alphabet
    
    [n,N]=size(X);

    if (nargin < 3) trainingInds=1:floor(N/2);  end
    if (nargin < 4) reg = -1;                   end
    if (nargin < 5) k=min(min(400,N),length(trainingInds));               end
    if (nargin < 6) k2=min(30,k);               end
    
    if (length(trainingInds) == 1)
        trainingInds = 1:trainingInds;
    end

    [d,inds]=pdist2(X(:,trainingInds)',X','euclidean','smallest',k);
    

    %%% CkNN Normalization
    dk=mean(d(2:k2,:));
    d = (d.^2)./(repmat(dk,k,1).*dk(inds));
    dk = mean(mean(d(2:k2,:)));
    
    %%% RBF Kernel
    d = exp(-d/dk);
    d = sparse(reshape(double(inds),N*k,1),repmat(1:N,k,1),reshape(double(d),N*k,1),length(trainingInds),N,N*k)';

    testingInds =setdiff(1:N,trainingInds);
    
    D = spdiags(1./sum(d(:,trainingInds),2),0,N,N);
    d = D*d;
    
    K  = d(trainingInds,trainingInds);
    Kz = d(testingInds,trainingInds);
    y  = Y(:,trainingInds)';
    fz = Y(:,testingInds)';
    I = speye(length(trainingInds));
    
    if (reg == -1)

        range = logspace(-14,0,8);

        bestdel=0;
        bestout=inf;

        for del = range

            Kreg = (K'*K+del*I);
            c = Kreg\(K'*y);
            Yrecon=(d*c)';
            out = sqrt(mean(mean((Y(:,testingInds)-Yrecon(:,testingInds)).^2)));

            if (out<bestout)
                bestout=out;
                bestdel=del;
            end
        end

        deldel=inf;
        del = bestdel;
        prevout = bestout;
        LM = 1e-6; %%% levenberg marquardt parameter

        iter = 0;

        while (abs(deldel/del)>1e-2)&&(iter<10)

            iter = iter+1;
            Kreg = (K'*K+del*I);
            c = Kreg\(K'*y);
            v = Kreg\c;
            vv = Kreg\v;

            dL = fz'*Kz*v - v'*Kreg*Kz'*Kz*v;
            ddL = v'*Kz'*Kz*v + del*v'*v - 2*fz'*Kz*vv + 2*v'*Kreg*Kz'*Kz*vv;

            deldel = dL/(ddL+LM);
            deltemp = del - deldel;

            Kreg = (K'*K+deltemp*I);
            c = Kreg\(K'*y);
            Yrecon=(d*c)';
            out = sqrt(mean(mean((Y(:,testingInds)-Yrecon(:,testingInds)).^2)));%+deltemp*sum(c.^2);

            while ((out>prevout)||(deltemp<=0))&&(iter<10)
                iter = iter+1;
                LM=1000*LM;

                deldel = dL/(ddL+LM);
                deltemp = del - deldel;

                Kreg = (K'*K+deltemp*I);
                c = Kreg\(K'*y);
                Yrecon=(d*c)';
                out = sqrt(mean(mean((Y(:,testingInds)-Yrecon(:,testingInds)).^2)));%+deltemp*sum(c.^2);

            end

            if (out<bestout)
                bestout=out;
                bestdel=deltemp;
            end

            %[del out LM]
            prevout = out;
            del=deltemp;
            LM = LM/1000;

        end
    
        sum(d(:,trainingInds),2);
        
    else
        
        del = eigs(d,1)*reg;
        
    end
    
    Kreg = (K'*K+del*I);
    c = Kreg\(K'*y);
    Yrecon=(d*c)';
    
    
    insampleRMSE = sqrt(mean(mean((Y(:,trainingInds)-Yrecon(:,trainingInds)).^2)));
    insampleCORR = corr(Y(:,trainingInds)',Yrecon(:,trainingInds)');
    
    outsampleRMSE = sqrt(mean(mean((Y(:,testingInds)-Yrecon(:,testingInds)).^2)));
    outsampleCORR = corr(Y(:,testingInds)',Yrecon(:,testingInds)');
    
    

end

