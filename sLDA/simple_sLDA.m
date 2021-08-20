function [insampleRMSE,insampleCORR,outsampleRMSE,outsampleCORR] = simple_sLDA(X,Y,nvars,delays,show,trainingInds)
% X [time, channel]
% Y [time, 1]
% slices  10:10:80 Just a sweep to see what works best. 
% alpha 0 "could simply look at the variance, or weigh
% things different. 
% sLDA_dim given
[N,n]=size(X);
if (nargin < 4) 
    delays = 0;
end
if (nargin < 6) 
    trainingInds=1:floor(N/2);  
end
if (nargin < 5) 
    show=false;  
end
if (nargin < 3) 
    nvars=10;  
end
slices = 10;
alpha = 0;

if delays ~= 0
    for d=delays
        X=[X(:,1:end-(d-1));X(:,d:end)];
    end
    Y=Y(max(delays):end);
end
[U,L] = sLDA(X(:,trainingInds),Y(trainingInds),slices,alpha,nvars);
X = (X'*U)';
[Yrecon,insampleRMSE,outsampleRMSE,insampleCORR,outsampleCORR, c] = KernelRegressionAuto(X,Y,trainingInds);
if show
    insampleRMSE
    insampleCORR
    outsampleRMSE
    outsampleCORR
    figure();
    plot(Y,'b')
    hold on;
    plot((Yrecon-mean(Yrecon))*std(Y)/std(Yrecon)+mean(Yrecon),'r')
end
end