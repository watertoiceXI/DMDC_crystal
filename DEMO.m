
clc;
close all;
%addpath('Util');

if (~exist('videoData'))
    load('LC_18.10--08-Dec-2017.mat')
    % for i = 1:1000
    %     imagesc(squeeze(videoData(400:3:1600,1:3:2000,1,i)));drawnow;
    % end
    videoData = squeeze(videoData(400:3:1600,1:3:2000,1,:));
end

subsize = 20;
t=1;
N=5000;
subims = GetRandomSubims(videoData,subsize,subsize,t,N);

% Convert from integers to double precision
subims = double(subims);

% PCA (just for visualization)
[coords,s,U,mu]=PCA(subims,16);

% Diffusion Map
[u,l,qest]=CIDM(subims');

% Diffusion map "modes"
A = subims*u;

% Represent each full image as a function and compute its Fourier
% coefficients (via convolution with diffusion modes) giving macroscopic
% coordinates
numims = size(videoData,3);
nummodes = size(A,2);
macroCoords = zeros(numims,nummodes);
tic;
for i = 1:numims
    for j = 1:nummodes
        
        macroCoords(i,j) = exp(-l(j,j)*1e4)*sum(sum(conv2(squeeze(videoData(:,:,i)),reshape(A(:,j),subsize,subsize),'valid')));
        
    end
end
toc;

% PCA of full images
X=PCA(macroCoords',9);
 
% Normalize voltages and PCA variables
V=voltage(1:1000:end);
V=V-mean(V);
V=V/std(V);
X=X./repmat(std(X')',1,1000);

% Try to regress the voltages from the PCA coordiantes of the macroscopic
% variables.  Need to compare this to trying to regress directly from the
% full images.
figure;plot(V);hold on;%plot(X');
AA=X'\V;
plot(X'*AA);


%%%%%%%% Plots

%%% Diffusion modes
figure;
for i = 1:16
    subplot(4,4,i);imagesc(reshape(A(:,i),20,20));
end

%%% Convolution with diffusion modes
figure;
imagesc(squeeze(videoData(:,:,1)));
figure;
for i = 1:16
    subplot(4,4,i);
    imagesc(conv2(squeeze(videoData(:,:,1)),reshape(A(:,i),20,20),'valid'));
end

%%% PCA modes
figure;
subplot(4,4,1);
imagesc(reshape(mu(1:subsize^2),subsize,subsize));
for i = 2:16
    subplot(4,4,i);
    imagesc(reshape(U(1:subsize^2,i-1),subsize,subsize));
end

%%% Convolution with PCA modes
figure;
imagesc(squeeze(videoData(:,:,1)));
figure;
for i = 1:16
    subplot(4,4,i);
    imagesc(conv2(squeeze(videoData(:,:,1)),reshape(U(:,i),20,20),'valid'));
end


 L=30; % number of subimages in each direction
% 
 figure;plot_flattened_dataset(coords(1:2,:),subims(1:subsize^2,:),L);
% figure;plot_flattened_dataset(coords(3:4,:),subims(1:subsize^2,:),L);
% figure;plot_flattened_dataset(coords(5:6,:),subims(1:subsize^2,:),L);
% figure;plot_flattened_dataset(coords(7:8,:),subims(1:subsize^2,:),L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


