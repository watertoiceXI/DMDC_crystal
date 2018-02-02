function x = SimpleNewton(f,x)
%SimpleNewton, minimize f wrt x

    y=f(x);
    miny=abs(y);
    minx=x;
    for i = 1:10
       x = x - y*1e-4/(f(x+1e-4)-y);
       y=f(x);
       if (abs(y)<miny)
           miny=abs(y);
           minx=x;
       end  
    end
    if (abs(y)>miny)
        x=minx;
    end
    
end

