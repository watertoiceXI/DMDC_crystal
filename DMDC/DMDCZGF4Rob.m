function info = DMDCZGF4Rob(ims,inds,segments,micromodel,microCoords,delays)
numdims = 16;
kappa = 1/max(delays);
k = 256;
sigma = 0;
alpha = 1/2;
videoTitle = 'tempVid.mpeg';
pr = 0;
k2 = k/4;
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



%%%%% INITIALIZE DISO PARAMETERS                %%%%%%%%%%%%%%%%%%%%


info.ims = ims;
info.height = size(ims,2);
info.width = size(ims,3);
if (inds == 0) 
    inds = 1:size(ims,1); 
end
info.inds = inds;
info.segments = segments;
info.delays = delays;
info.numdims = numdims;
info.kappa = kappa;
info.k = k;
info.sigma = sigma;
info.alpha = alpha;
info.videoTitle = videoTitle;
info.pr=pr;


%%%%%% MICROSCOPIC MODEL                        %%%%%%%%%%%%%%%%%%%%

info.micromodelNumber = micromodel;
%     [info.macroVars,info.microModel] = buildMicroModel(info,micromodel,microCoords);

[info.macroVars,info.microModel] = buildMicroModel(info,micromodel,microCoords);

% keyboard

%%%%%% TIME DELAY EMBEDDING & PRECOMPRESSION    %%%%%%%%%%%%%%%%%%%%

info.compressionDim = 2048;
% info.compressionDim = 512;
info.micromodelNumber = micromodel;
%size(info.macroVars)

 [info.x,info.meanState,info.preCompression,info.originalImageIndices,info.meanDistPerFrame] = DelayCoordinates(info.macroVars,info.inds,info.segments,info.delays,info.kappa,info.compressionDim);


%%%%%% GPU ADAPTED DIFFUSION MAPS               %%%%%%%%%%%%%%%%%%%%




[info.q,info.b,info.sigma,info.quest,info.VBAutoDimension] = VBAuto(info.x,k,k2,numdims,1);

% info.ims=0;


%%%%%% DMDC COMPONENTS                          %%%%%%%%%%%%%%%%%%%%

[ info.components ] = DMDC(info,0,0,videoTitle);

% YouGotMail
end



