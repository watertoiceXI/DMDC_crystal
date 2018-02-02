

%%% flat torus, real 2d data

t0=2*pi*(1/10:1/10:1);
[t,p]=meshgrid(t0,t0);t=t(:);p=p(:);

x=[cos(t');sin(t');cos(p');sin(p')];
[u,l,qest]=CIDM(x');

%%% no extrinsic curvature circle

t=2*pi*(1/100:1/100:1);
x=[cos(t);sin(t)];
[u,l,qest]=CIDM(x');

%%% bad embedding circle, fake 2d data

x=[cos(t);sin(t);cos(4*t);sin(4*t)];
[u,l,qest]=CIDM(x');