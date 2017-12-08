function [x,meanState,preCompression,originalImageIndices,meanDistPerFrame] = DelayCoordinates(macroVars,inds,segments,delays,kappa,compressionDim,setBlockSize)

%INPUTS : macroVars are TxN dimensional where T represents time and N
%represents the macroscopic state.

fprintf('Delay embedding and precompression... ');
tic;

if (isempty(segments))
    segments = [0 length(inds)];
else
    segments = [0 segments];
end

if (nargin<7)
    setBlockSize = 2^2;   %process this many states per block
end

numStates = length(inds)-(length(segments)-1)*max(delays);
stateSize = size(macroVars,2);
embeddDim = stateSize*length(delays);
preCompression = randn(embeddDim,compressionDim);
[preCompression,~] = qr(preCompression,0);

if (size(preCompression,2) < compressionDim)
    
    compressionDim = size(preCompression,2);
end
x = zeros(numStates,embeddDim);

meanState = [];
originalImageIndices = zeros(numStates,1);

for jj = 1:length(segments)-1
    segStart = segments(jj);
    numStatesSeg = segments(jj+1)-max(delays)-segStart;
    % Split into blocks to conserve memory and maintain speed
    blockSize = setBlockSize;   %process this many states pe r block
    
    numBlocks = floor(numStatesSeg/blockSize);
    
    for i = 1:numBlocks
        blockStart = max(delays) + segStart + (i-1)*blockSize;
        fullState = zeros(blockSize,embeddDim);
        %
        for j = 1:length(delays)
            fullState(:,stateSize*(j-1) + (1:stateSize)) = exp(-delays(j)*kappa)*macroVars(blockStart + (1:blockSize) - delays(j),:);
           
        end
        if (isempty(meanState))
            meanState = mean(fullState);
        end
        
        x(blockStart + (1:blockSize) - jj*max(delays),:) = (double(fullState)-repmat(meanState,size(fullState,1),1));
        originalImageIndices(blockStart + (1:blockSize) - jj*max(delays)) = blockStart + (1:blockSize) - delays(j);
    end
   
    blockStart = max(delays) + segStart + numBlocks*blockSize;
    if (numStatesSeg > blockStart - segStart - max(delays))
        
        blockSize = numStatesSeg-numBlocks*blockSize;
        fullState = zeros(blockSize,embeddDim);
        for j = 1:length(delays)
            fullState(:,stateSize*(j-1) + (1:stateSize)) = exp(-delays(j)*kappa)*macroVars(blockStart + (1:blockSize) - delays(j),:);
        end
        if (isempty(meanState))
            meanState = mean(fullState);
        end
        x(blockStart + (1:blockSize) - jj*max(delays),:) = (double(fullState)-repmat(meanState,size(fullState,1),1));
        originalImageIndices(blockStart + (1:blockSize) - jj*max(delays)) = blockStart + (1:blockSize) - delays(j);
    end
    %     distPerFrame(jj) = mean(sqrt(sum((x(2:end,:)-x(1:end-1,:)).^2,2)));
    
    xsegstart = segStart + max(delays)+1-jj*max(delays);
   
%     xsegend = blockStart+blockSize - jj*max(delays) % crashes
   
    xsegend = blockStart - jj*max(delays);  %fix
    
    distPerFrame(jj) = mean(sqrt(sum((x(xsegstart+1:xsegend,:)-x(xsegstart:xsegend-1,:)).^2,2)));
%     if size(find(isnan(distPerFrame)),2)>0
%         keyboard
%     end
end

meanDistPerFrame = (mean(distPerFrame));
% keyboard
x = x/meanDistPerFrame;

toc;


end

