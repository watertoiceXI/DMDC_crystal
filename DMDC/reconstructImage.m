function [ reconIm ] = reconstructImage(info,timeIndex,componentIndices)
   info.q(timeIndex,ninds)
    ninds = componentIndices;
    reconIm = info.q(timeIndex,ninds)*info.b(ninds,ninds).^(1/info.sigma)*info.components(:,ninds)';
    reconIm = reshape(reconIm,[size(info.ims,2),size(info.ims,3)]);


end

