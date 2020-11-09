function [DimVol,q,b,VBAutoDimension] = DMDC_simple(sig,delays,numdims_out,segments)

if (nargin<3)
    numdims_out = 16;
end
if (nargin<4)
    segments=[];
end

kappa = 1/max(delays);
k = 256;
sigma = 0;
alpha = 1/2;
videoTitle = 'tempVid.mpeg';
pr = 0;
k2 = k/4; 
inds = 1:size(sig,1); 
NoVariab = 256;

% EXAMPLE CALLS
% info = DISO6(ims,inds,segments,1,0,1:8,8,1,256,0,1/2);
% info = DISO6(ims,inds,segments,1,0,1:8);


%%%%% EXAMPLE IMAGE LOADING ROUTINES            %%%%%%%%%%%%%%%%%%%%

% Load images with a prefix followed by digits indicating the time
%ims = loadImages(fileprefix,fileextension,numdigits,inds,xoffset,yoffset,width,height)
%segments = [];

% Single segment of data, all in a single folder
%[ims,inds,segments,numLoaded] = loadFolder('TestData',512,'tif');

% Multiple data segments, each in a subfolder of single "super folder"
%[ims,inds,segments,numLoaded] = loadSuperFolder('TestData',512,'tif');
% global gpuID
% gpuID = -1


micromodel = 1; %Random Projections


%%%%% INITIALIZE DISO PARAMETERS                %%%%%%%%%%%%%%%%%%%%






%%%%%% TIME DELAY EMBEDDING & PRECOMPRESSION    %%%%%%%%%%%%%%%%%%%%

compressionDim = 2048;
% compressionDim = 512;
micromodelNumber = micromodel;
%size(macroVars)

 [x,meanState,preCompression,originalImageIndices,meanDistPerFrame] = DelayCoordinates(sig,inds,segments,delays,kappa,compressionDim);


%%%%%% GPU ADAPTED DIFFUSION MAPS               %%%%%%%%%%%%%%%%%%%%




[q,b,sigma,quest,VBAutoDimension] = VBAuto(x,k,k2,numdims_out,1);

% ims=0;


%%%%%% DMDC COMPONENTS                          %%%%%%%%%%%%%%%%%%%%

hh=diag(b);
%     [DLGb,VLgb] =DimT(hh,infon);
DimVol = zeros(2,1);
NoVariab = min(NoVariab,length(hh));
for ii=5:NoVariab
    [DLGb(ii), VLgb(ii)] = dimvolmax(hh(1:ii),sigma);
end

ok=find(abs(diff(DLGb(:)))<0.01);
ok=4+find(abs(diff(DLGb(5:end)))<0.00001);
if size(ok)~=0
    mn=ok;
    mm=mn(1);
    DimVol(1)=DLGb(mm);
    DimVol(2)=VLgb(mm);
end


% YouGotMail
end



