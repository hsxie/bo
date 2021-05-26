% 18-12-23 23:05 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Normalization data of omega_n and k_n for BO

% if(imus(1)==0 || iem==0)
%   kn=1/lambdaDs(1); % normalized k by kn, modify here to use other kn
%   wn=wps1;
% else
%   kn=1/max(abs(rhocs));
%   wn=min(abs(wcs));
% end

% % 18-12-19 00:21 for Xie2014 Darwin model
% wn=wcs1;
% % cwp=sqrt(c2)/wcs1; % c/omega_{p1}
% kn=1/cwp; % normalized k by kn, modify here to use other kn

% kn=1/lambdaD; % 18-12-16 17:52 to compare with Umeda12
% kn=1/lambdaDs(1); % 18-12-17 19:37 to compare with TaoX's ES3D loss cone
% wn=wps1; % 18-12-17 20:26 to compare with unmagnetized ES3D version

lambdaD=sqrt(1./sum(1./lambdaDs.^2));
wn=abs(wcs(1)); % normalized omega by wn, modify here to use other wn
% kn=cwp; % normalized k by kn, modify here to use other kn

kn=1/lambdaD;
