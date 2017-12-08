function [ds,inds] = knn(R,Q,k,gpuID)
% Find the k nearest neighbors for each element of the query set Q among
% the points in the reference set R

% Q is Nxd where N is the number of query points and d is the dimension
% R is Mxd where M is the number of reference points and d is the dimension
% gpuID is the ID of the gpuDevice to use or -1 to use the CPU.
% ds is Nxk and contains the distances to the k nearest neighbors.
% inds is Nxk and contains the indices of the k nearest neighbors.
% the dimension of the points must be less than 2^13

% %     
%     gpuID=-1
    if (gpuID>=0)
        
        maxN = 2^10*256/(2^ceil(log(k)/log(2)));
        
        NQ = size(Q,1);
        NR = size(R,1);
        
        nqi = ceil(NQ/maxN);
        nqr = ceil(NR/maxN);
        
        if (nqi*nqr == 1)
            
            [ds,inds] = knnGPU(single(R),single(Q),k,gpuID);
                
        else
        
            ds = [];
            inds = [];

            for i = 1:nqi

                dsi = [];
                indsi = [];

                if (i<nqi)
                    Qinds = (1:maxN)+maxN*(i-1);
                else
                    Qinds = (maxN*(i-1)+1):NQ;
                end

                for r = 1:nqr

                    if (r<nqr)
                        Rinds = (1:maxN)+maxN*(r-1);
                    else
                        Rinds = (maxN*(r-1)+1):NR;
                    end
%                     keyboard
% size(single(R(Rinds,:)))
% size(single(Q(Qinds,:)))
% % size(min(k,length(Rinds)))     
[dst,indst] = knnGPU(single(R(Rinds,:)),single(Q(Qinds,:)),min(k,length(Rinds)),gpuID);

                    dsi = [dsi dst];
                    indsi = [indsi indst+(r-1)*maxN];

                end

                [dst,indst] = sort(dsi');
                dst = dst';
                indst = indst';

                ds = [ds; dst(:,1:k)];
                inds = [inds; indsi(sub2ind(size(indsi),repmat((1:size(indsi,1))',1,k),indst(:,1:k)))];

            end
            
        end       
    
    else
    
        N = size(Q,1);
        
        ds = zeros(N,k);
        inds = zeros(N,k);

        for i = 1:N            %find initial distances

            dtemp = sum(bsxfun(@mydist,R,Q(i,:)),2);
            [dst,indst] = sort(dtemp);

            ds(i,:) = sqrt(dst(1:k));
            inds(i,:) = indst(1:k);

        end
        
        
    end


end

