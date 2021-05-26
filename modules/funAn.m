% Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, 2018-11-29 16:38
% Adaptive Simpson quadrature formula to calculate the v_perp integral
% for ring beam dispersion relation.
% To do: Speed up using Gauss integral
% 18-12-05 08:50 We have test this version with loss-cone version at
% drift bi-Maxwellian limit, and find the runtime not change too much, 
% e.g., 97s -> 106s. Thus, the computation speed and accuracy of this
% function is not a big issue at this moment.
function An=funAn(n,a,b,c)
na=length(a); An=0.*a;
ym=5.0; 
% ym=10.0; 
% reltol=1e-15; abstol=1e-10;
ymin=max(0,b-ym);
% ymin=0;
ymax=ym+b;
for ja=1:na
    f=@(y) besselj(n,a(ja)*y).^2.*exp(-(y-b).^2).*(y-c);  
    An(ja)=integral(f,ymin,ymax);
%     An(ja)=integral(f,ymin,ymax,'RelTol',reltol,'AbsTol',abstol);
%     An(ja)=integral(f,ymin,ymax,'RelTol',reltol);
end