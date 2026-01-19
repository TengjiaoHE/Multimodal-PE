function [ClenCurtis_H,wcc] = ClenshawCurtis (n,H)
a = H(1);
b = H(end);
xk = cos((0:1:n)*pi/n);
%         xk = cos((0.5:1:n-0.5)*pi/n);
ClenCurtis_H = (xk+1)*(b-a)/2+a;
[wf1,wf2,wcc] = fejer(n);
wcc(n+1) = wcc(1);