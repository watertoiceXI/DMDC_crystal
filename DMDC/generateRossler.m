function [ xx ] = generateRossler(T,aparam,dnoise,onoise)
%%% aparam = .1 for periodic, increase up to .4 for chaos

    global paramset; 
    paramset = [aparam dnoise];     
    

    xx = zeros(T,3);

    state = [rand+1;0;0];
    
    %%%% throw out the initial transient
    for i = 1:10000
        substeps = 1;    
        dt = .05/substeps;
        for j = 1:substeps
            k1=dt*RosslerVectorField(state);
            k2=dt*RosslerVectorField(state+k1/2);
            k3=dt*RosslerVectorField(state+k2/2);
            k4=dt*RosslerVectorField(state+k3);
            state=state+k1/6+k2/3+k3/3+k4/6 + sqrt(dnoise)*randn(size(state));
        end
    end

    for i = 1:T
        substeps = 1;    
        dt = .05/substeps;
        for j = 1:substeps
            k1=dt*RosslerVectorField(state);
            k2=dt*RosslerVectorField(state+k1/2);
            k3=dt*RosslerVectorField(state+k2/2);
            k4=dt*RosslerVectorField(state+k3);
            state=state+k1/6+k2/3+k3/3+k4/6 + sqrt(dnoise)*randn(size(state));
        end       
        xx(i,:) = state + sqrt(onoise)*randn(size(state));
    end

end

function dx = RosslerVectorField(state)
    x=state;
    global paramset;
    a = paramset(1);
    b = 0.2;
    c = 5.7;

    dx = zeros(size(x));
    dx(1,:) = -x(2,:) - x(3,:);
    dx(2,:) = x(1,:) + a*x(2,:);
    dx(3,:) = b + x(3,:)*(x(1,:)-c);
end