function [A,dataToyIm]=DEMOrossler()

aparams = .1:.05:.21;

P = length(aparams);
T = 5*2048;
modes = 96;
delays = 500;
kappa = 1/125;
k = 128*2;
sigma = 0.00;
bb=50;
D=500;
pr=1;
NoVariab = 256;
W=50;
A = zeros(T,W,W);
vDefect = 1/2;
xcoords(1:W,1) = 1:W;
xcoords = repmat(xcoords,1,W);
ycoords = xcoords';

numdims_out = 16;
kappa = 1/max(delays);
k = 256;
sigma = 0;
alpha = 1/2;
pr = 0;
k2 = k/4;

for p = 1:P
    
    p
    aparam=aparams(p)
    xx = generateRossler(T,aparam,0,0);
    
    A = xx(:,1:2);
    %xx=func(p).xx;
    %     figure
%     for i=1:T
%         
%         tt=xx(i,1);
%         x2 = W/2;
%         y2 = W/2 + (W/2)*.1*((-1)*vDefect*tt);
%         A(i,:,:) = squeeze(A(i,:,:)) + exp(1).^((-(xcoords-x2).^2-(ycoords-y2).^2)/10);
%         
%         %          imagesc(squeeze(A(i,:,:)),[0 1]);
%         % %     %
%         %             drawnow;
%     end
    
    [DimVol,q,b,VBAutoDimension]  = DMDC_simple(A,delays,numdims_out);
    %
    %
    figure(78);
    hold all
    plot(DimVol(1),DimVol(2),'*')
    
end
%%
kappa = 1/125;
k = 128*2;
sigma = 0.00;
bb=50;
D=500;
pr=1;
NoVariab = 256;
W=50;
A = zeros(T,W,W);
vDefect = 1/2;
xcoords(1:W,1) = 1:W;
xcoords = repmat(xcoords,1,W);
ycoords = xcoords';


figure
for D=50:50:5000
    for i=1:T
        
        tt=xx(i,1);
        x2 = W/2;
        y2 = W/2 + (W/2)*.1*((-1)*vDefect*tt);
        A(i,:,:) = squeeze(A(i,:,:)) + exp(1).^((-(xcoords-x2).^2-(ycoords-y2).^2)/10);
        
        %          imagesc(squeeze(A(i,:,:)),[0 1]);
        % %     %
        %             drawnow;
    end
    
    infon  = DMDCZGF(A,1:size(xx,1),size(xx,1),1,[],0:bb:D,NoVariab,4/D,256,0.1,1/2,0,pr);
    %
    %
    hh=diag(infon.b);
    %     [DLGb,VLgb] =DimT(hh,infon);
    [DimVol] = dimmod(infon, NoVariab);
    figure(78);
    hold all
    plot(DimVol(1),DimVol(2),'*')
    
    figure(79);
    hold all
    plot( DimVol(1),DimVol(2),'.')
    dataToyDelay(D).dim=DimVol(1);
    dataToyDelay(D).vol=DimVol(2);
    dataToyDelay(D).parameter = aparam;
    dataToyDelay(D).hh = hh;
    dataToyDelay(D).D  = D;
    dataToyDelay(D).sigma  = infon.sigs;
    
end



