
clc;
close all;

N = 5000;
c = 15;
delays = 1:8;
kappa = 1/max(delays);
compressionDim = 2048; %ignore, mostly
k = 256;
sigma = 0;
k2 = k/4; 
numdims_out = 5;

sig = randn(c,N);

% PCA (just for visualization)
[X,s,U,mu]=PCA(sig,floor(c/4));

% Diffusion Map
[u,l,qest]=CIDM(sig');

% Diffusion map "modes"
A = sig*u;



 
Y = randn(1,N);

% Try to regress the voltages from the PCA coordiantes of the macroscopic
% variables.  Need to compare this to trying to regress directly from the
% full images.
figure;plot(Y);hold on;%plot(X');
AA=X'\Y;
plot(X'*AA);


[xD,meanState,preCompression,originalImageIndices,meanDistPerFrame] = DelayCoordinates(sig,1:N,[],delays,kappa,compressionDim);
[info.q,info.b,info.sigma,info.quest,info.VBAutoDimension] = VBAuto(xD,k,k2,numdims_out,1);
