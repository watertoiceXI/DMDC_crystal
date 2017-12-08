
cpuTimes = 0;
gpuTimes = 0;

mErrVec = 0;
maxErrVec = 0;



        
for n=8
        d=1;
        k = d*256;
        N = 9879;    %number of points
        M = (2048);    %dimension
        [N M]

        X = randn(N,M);

        tic;
        ds = zeros(N,k);
        ii = zeros(N,k);
        parfor i = 1:256
            %dtemp = sum(bsxfun(@minus,X,X(i,:)).^2,2);
            dtemp = sum(bsxfun(@mydist,X,X(i,:)),2);
            [sortedd,inds] = sort(dtemp);
            ds(i,:) = dtemp(inds(1:k)).^(.5);
            ii(i,:) = inds(1:k);
        end
        cpuTimes(n,d) = toc;
        
        
        tic;
        %[dss,i] = knnGPU(single(X),single(X),k,0);
        [dss,iind] = knn(X,X,k,0);
        numNaNs(n,d) = sum(sum(isnan(dss)));
        gpuTimes(n,d) = toc;
        
        mErr = 0;
        maxErr = 0;
        for j = 1:256
            dtemp = sqrt(sum(bsxfun(@minus,X(iind(j,:),:),X(j,:)).^2,2));
            mErr = mErr + (mean(abs(dtemp'-dss(j,:))))/N;
            if (max(abs(dtemp'-dss(j,:))) > maxErr)
                maxErr = max(abs(dtemp'-dss(j,:)));
            end
        end
        
        mErrVec(n,d) = mErr;
        maxErrVec(n,d) = maxErr;

        errors(n,d) = sum(sum(abs(iind-int32(ii)) > 0));
        inds = find(abs(iind-int32(ii)) == 0);
        meanErr(n,d) = mean(mean(abs(ds(inds)-dss(inds))));
     
end
        


        



