function [insampleRMSE,insampleCORR,outsampleRMSE,outsampleCORR] = auto_Regr(LFP,predict_step,fhi,flo,fpred,nvars,delays)
fs = 1000;
%building the cwt
[cwt_X,f] = cwt(LFP,'morse',fs);
[~,arghi] = min(abs(f-fhi));
[~,arglo] = min(abs(f-flo));
[~,argpred] = min(abs(f-fpred));
fprintf(['actual hi: ', num2str(f(arghi)),' actual lo: ',num2str(f(arglo)), 'actual predict: ',num2str(f(argpred))]);
X = abs(cwt_X(arghi:arglo,1:(end-predict_step)));
Y = abs(cwt_X(argpred,predict_step:end));
%now, recording variables for the loop. If the RSME is not less than the
%std, then we're worse than just guessing the mean.... 
thr_RSME = std(Y);
fprintf(['Threshold RSME: ', num2str(thr_RSME)]);
bestn = 0;
best_RSME=1000000000;
%now, looping over the parameters
for n=nvars
    [insampleRMSE,insampleCORR,outsampleRMSE,outsampleCORR] = simple_sLDA(X,Y,n,delays);
    if insampleRMSE<thr_RSME
        if outsampleRMSE<best_RSME
            bestn = n;
            best_RSME = outsampleRMSE;
        end
    end
end

%[insampleRMSE,insampleCORR,outsampleRMSE,outsampleCORR] = ReservoirRegression(X,Y);


if bestn == 0
    fprintf('Failed!!!!')
    return
end
fprintf(['Best Nvar: ',num2str(bestn)])
[insampleRMSE,insampleCORR,outsampleRMSE,outsampleCORR] = simple_sLDA(X,Y,bestn,delays,true);
