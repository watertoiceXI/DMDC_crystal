function [epsilon] = tuneEpsilon(d,method)
% Different methods of estimating the dimension
        if (method == 1)
            
            %%% minimize: epsilons.^4.*ddims.^2./dims.^4
            %%% find the smallest epsilon with derivative of dimension = 0
            epsilon = min(d(d>0))/10;
            f = @(x) leftMostFlatCriterion(d,x);
            epsilon = SimpleNewton(f,epsilon);
            
        elseif (method == 2)
            %%% look for smallest epsilon before loss of numerical
            %%% precision in the kernel sums
            dim0=0;dim=0;
            epsilon = max(max(d));
            while (abs(dim-dim0)<1e-1)
                epsilon=epsilon/2;
                dim0 = 2*(log(sum(sum(exp(-d(d~=0)/(2*epsilon*(1+1e-4))))))-log(sum(sum(exp(-d(d~=0)/(2*epsilon))))))/(log(epsilon*(1+1e-4))-log(epsilon));
                dim = 2*(log(1+sum(sum(exp(-d(d~=0)/(2*epsilon*(1+1e-4))))))-log(1+sum(sum(exp(-d(d~=0)/(2*epsilon))))))/(log(epsilon*(1+1e-4))-log(epsilon));
            end
        
        elseif (method == 3)
            
            %%% maximize the dimension
            epsilon = min(d(d>0))
            f = @(x) estimateDimension(d,x);
            epsilon = SimpleNewton(f,epsilon);
    
        end
        
        
        plot=11;
        if (plot)
            epsilons = 2.^(-30:.1:10);
            epsilons = epsilon*2.^(-5:.01:5);
            dims=zeros(1,length(epsilons));
            ddims=zeros(1,length(epsilons));
            for i=1:length(epsilons)
                [dims(i),ddims(i)]=estimateDimension(d,epsilons(i));
            end 
            figure(plot);hold off;
            semilogx(epsilons,dims);hold on;
            semilogx(epsilons,ddims);
            semilogx(epsilons,epsilons.^4.*abs(ddims).^2./dims.^4);
            dim = estimateDimension(d,epsilon);
            semilogx(epsilon,dim,'o');
            ylim([-1 5]);
        end
end

function val = leftMostFlatCriterion(d,epsilon)
    [dim,ddim]=estimateDimension(d,epsilon);
    val = epsilon^4*ddim^2/dim^4;
end

