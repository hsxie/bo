% 18-12-16 08:58, Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China.
% This file set the BO (PDRK) kernel matrix elements, of electromagnetic
% 3D case.
%
% This new version support: theta=pi/2, kz=0; theta=0, kx=0; and also kz<0.
%
% 18-12-18 21:21 We try two methods to calculate b11, b22, b33 terms, it
% seems that imethod=1 works better, which convergent easily for small N.
% To understand why.
%
% 18-12-28 21:16 A bug is fixed for EM3D-U p13 & p31 terms
% 19-01-29 12:33 Another bug is fixed for EM3D-U p13 terms
% 20-10-25 21:40 Another bug is fixed for EM3D-M p32 vdsy term
% 20-11-14 08:32 Update to can handle the calculation of density and
% velocity perturbation/polarization of each species.

% % -- main program begin to set the matrix elements --

if(kx==0) % 18-12-16 14:59, to remove the singularity when kx=0. To update.
  kx=1e-30;
end
% if(kz==0) % not needed
%   kz=1e-30;
% end

k=sqrt(kz^2+kx^2);
asab=kx*rhocsab*sqrt(2); % as=kx*vtps./wcs, note a sqrt(2), 18-11-29 16:34

kvtsab2(1,:)=((kx.*vtpsab(1,:)).^2+(kz.*vtzs).^2);
kvtsab2(2,:)=((kx.*vtpsab(2,:)).^2+(kz.*vtzs).^2);
kvtsab=sqrt(kvtsab2); % define new (k*vts)^2 for unmagnetized species

M=sparse(NN,NN);
snj=0;

if(kz<0) % 18-12-22 00:54
  czjj=-czj;
else
  czjj=czj;
end

% initialize
% b11=0; b12=0; b13=0; b21=0; b22=0; b23=0; b31=0; b32=0; b33=0;
% 20-11-14 10:35 update b_ij to b_{ij,s} for calcualte dJ_s
b11s=zeros(1,S); b12s=0.*b11s; b13s=0.*b11s; b21s=0.*b11s; b22s=0.*b11s; 
b23s=0.*b11s; b31s=0.*b11s; b32s=0.*b11s; b33s=0.*b11s;
% csnj=zeros(1,3*SNJ);
csnj=zeros(1,SNJ); % 18-12-17 22:42
b11snj=csnj.*0; b12snj=csnj.*0; b13snj=csnj.*0;
b21snj=csnj.*0; b22snj=csnj.*0; b23snj=csnj.*0;
b31snj=csnj.*0; b32snj=csnj.*0; b33snj=csnj.*0;

Anabb=zeros(1,2); Anab0=0.*Anabb; Bnabb=0.*Anabb; Bnab0=0.*Anabb;
Cnabb=0.*Anabb; Cnab0=0.*Anabb;

% 18-12-18 21:17 !!!%%%
% imethod=1; % seems =0 not works well or need large N to convergent, % 19-11-30 12:56 especially for ion Bernstein wave
imethod=0; 

for s=1:S % species
 if(ismember(s,jds))
  Nsss=1:Ns(s); % unmagnetized species, for different sigma
 else
  Nsss=-Nss(s):Nss(s); % magnetized species, for different N
 end
 for n=Nsss % sum_N or sum_sigma
  
  for j=1:J % poles of Z(zeta)
   snj=snj+1;
   if(ismember(s,jds)) % unmagnetized species
       
    iab=n; % for unmagnetized species, n is sigma
    if(iab>2)
      disp('Some thing goes wrong, please check Ns, Nss and Nsss.');
    end
    % Note: kvtsab always > 0, thus we do not need special treat kx, kz <0
    csnj(snj)=czj(j)*kvtsab(iab,s)+kz*vdsz(s)+kx*vdsx(s)-1i*nus(s);
    cnj=csnj(snj);
      
    tmp=rsab(iab,s)*wps2(s)*2*bzj(j)/cnj;
    
    p11snj=-1i*nus(s)*((kx*vtpsab(iab,s))^2*czj(j)^2+kx*vdsx(s)*...
      kvtsab(iab,s)*czj(j)+0.5*(kz*vtzs(s))^2)/kvtsab2(iab,s) + ...
      (k*kx*vtpsab(iab,s)^2)^2/kvtsab(iab,s)^3*czj(j)^3 + ...
      2*k^2*kx*vtpsab(iab,s)^2/kvtsab2(iab,s)*vdsx(s)*czj(j)^2 + ...
      (0.5*vtpsab(iab,s)^2*k^2*kz^2*vtzs(s)^2/kvtsab(iab,s)^3+...
      k^2*vdsx(s)^2/kvtsab(iab,s)+(lmdTab(iab,s)-1)*kx^2*kz^2*...
      vtpsab(iab,s)^4/kvtsab(iab,s)^3 )*czj(j) + (lmdTab(iab,s)-1)*...
      kx*kz^2/kvtsab2(iab,s)*vtpsab(iab,s)^2*vdsx(s);
  
    p12snj=vdsy(s)*(k^2*kx*vtpsab(iab,s)^2*czj(j)^2 + ...
      k^2*vdsx(s)*kvtsab(iab,s)*czj(j) + 0.5*(lmdTab(iab,s)-1)*...
      kx*kz^2*vtpsab(iab,s)^2)/kvtsab2(iab,s);
  
    p21snj=vdsy(s)*(-1i*nus(s)*kx*kvtsab(iab,s)*czj(j) + k^2*kx*...
      vtpsab(iab,s)^2*czj(j)^2 + k^2*vdsx(s)*kvtsab(iab,s)*czj(j) +...
      0.5*(lmdTab(iab,s)-1)*kx*kz^2*vtpsab(iab,s)^2)/kvtsab2(iab,s);
  
    p22snj=-0.5*1i*nus(s)+(0.5*vtpsab(iab,s)^2+vdsy(s)^2)*...
      k^2/kvtsab(iab,s)*czj(j);

%     p13snj=-1i*nus(s)*kz*(kx*vtpsab(iab,s)^2*czj(j)^2+vdsx(s)*...
%       kvtsab(iab,s)*czj(j)-0.5*kx*vtpsab(iab,s)^2)/kvtsab2(iab,s) + ...
%       k^2*kx*kz*vtzs(s)^2*vtpsab(iab,s)^2/kvtsab(iab,s)^3*czj(j)^3 + ...
%       k^2*(kx*vtpsab(iab,s)^2*vdsz(s)+kz*vtzs(s)^2*vdsz(s)...
%       )/kvtsab2(iab,s)*czj(j)^2 + (k^2*vdsx(s)*vdsz(s)*kvtsab2(iab,s)-...
%       0.5*k^2*kx*kz*vtzs(s)^2*vtpsab(iab,s)^2-0.5*(lmdTab(iab,s)-1)*...
%       (kx^2*vtpsab(iab,s)^2-kz^2*vtzs(s)^2)*kx*kz*vtpsab(iab,s)^2 ...
%       )/kvtsab(iab,s)^3*czj(j) + 0.5*(lmdTab(iab,s)-1)*kx*kz*...
%       vtpsab(iab,s)^2*(kx*vdsz(s)-kx*vdsx(s))/kvtsab2(iab,s);
  
    % fixed a bug, 18-12-28 20:01
    p13snj=-1i*nus(s)*kz*(kx*vtpsab(iab,s)^2*czj(j)^2+vdsx(s)*...
      kvtsab(iab,s)*czj(j)-0.5*kx*vtpsab(iab,s)^2)/kvtsab2(iab,s) + ...
      k^2*kx*kz*vtzs(s)^2*vtpsab(iab,s)^2/kvtsab(iab,s)^3*czj(j)^3 + ...
      k^2*(kx*vtpsab(iab,s)^2*vdsz(s)+kz*vtzs(s)^2*vdsx(s)... % a bug fixed vdsz -> vdsx
      )/kvtsab2(iab,s)*czj(j)^2 + (k^2*vdsx(s)*vdsz(s)*kvtsab2(iab,s)-...
      0.5*k^2*kx*kz*vtzs(s)^2*vtpsab(iab,s)^2-0.5*(lmdTab(iab,s)-1)*...
      (kx^2*vtpsab(iab,s)^2-kz^2*vtzs(s)^2)*kx*kz*vtpsab(iab,s)^2 ...
      )/kvtsab(iab,s)^3*czj(j) + 0.5*(lmdTab(iab,s)-1)*kx*kz*...
      vtpsab(iab,s)^2*(kz*vdsz(s)-kx*vdsx(s))/kvtsab2(iab,s); % 19-01-29 11:42 kx*vdsz -> kz*vdsz
  
%     p31snj=-1i*nus(s)*kx*(kz*vtzs(s)^2*czj(j)^2+vdsz(s)*kvtsab(iab,s)* ...
%       czj(j)-0.5*kz*vtzs(s)^2)/kvtsab2(iab,s) + ...
%       k^2*kx*kz*vtzs(s)^2*vtpsab(iab,s)^2/kvtsab(iab,s)^3*czj(j)^3 + ...
%       k^2*(kx*vtpsab(iab,s)^2*vdsz(s)+kz*vtzs(s)^2*vdsz(s)...
%       )/kvtsab2(iab,s)*czj(j)^2 + (k^2*vdsx(s)*vdsz(s)*kvtsab2(iab,s)-...
%       0.5*k^2*kx*kz*vtzs(s)^2*vtpsab(iab,s)^2-0.5*(lmdTab(iab,s)-1)*...
%       (kx^2*vtpsab(iab,s)^2-kz^2*vtzs(s)^2)*kx*kz*vtpsab(iab,s)^2 ...
%       )/kvtsab(iab,s)^3*czj(j) + 0.5*(lmdTab(iab,s)-1)*kx*kz*...
%       vtpsab(iab,s)^2*(kx*vdsz(s)-kx*vdsx(s))/kvtsab2(iab,s);
  
    p31snj=-1i*nus(s)*kx*(kz*vtzs(s)^2*czj(j)^2+vdsz(s)*kvtsab(iab,s)* ...
      czj(j)-0.5*kz*vtzs(s)^2)/kvtsab2(iab,s) + ...
      k^2*kx*kz*vtzs(s)^2*vtpsab(iab,s)^2/kvtsab(iab,s)^3*czj(j)^3 + ...
      k^2*(kx*vtpsab(iab,s)^2*vdsz(s)+kz*vtzs(s)^2*vdsx(s)... % a bug fixed vdsz -> vdsx
      )/kvtsab2(iab,s)*czj(j)^2 + (k^2*vdsx(s)*vdsz(s)*kvtsab2(iab,s)-...
      0.5*k^2*kx*kz*vtzs(s)^2*vtpsab(iab,s)^2-0.5*(lmdTab(iab,s)-1)*...
      (kx^2*vtpsab(iab,s)^2-kz^2*vtzs(s)^2)*kx*kz*vtpsab(iab,s)^2 ...
      )/kvtsab(iab,s)^3*czj(j) + 0.5*(lmdTab(iab,s)-1)*kx*kz*...
      vtpsab(iab,s)^2*(kz*vdsz(s)-kx*vdsx(s))/kvtsab2(iab,s); % a bug fixed kx -> kz

    p23snj=vdsy(s)*(-1i*nus(s)*kz/kvtsab(iab,s)*czj(j)+...
      (k^2*kz*vtzs(s)^2*czj(j)^2 + k^2*vdsz(s)*kvtsab(iab,s)*czj(j) - ...
      0.5*(lmdTab(iab,s)-1)*kx^2*kz*vtpsab(iab,s)^2)/kvtsab2(iab,s));

    p32snj=vdsy(s)*(...
      (k^2*kz*vtzs(s)^2*czj(j)^2 + k^2*vdsz(s)*kvtsab(iab,s)*czj(j) - ...
      0.5*(lmdTab(iab,s)-1)*kx^2*kz*vtpsab(iab,s)^2)/kvtsab2(iab,s));
  
    p33snj=-1i*nus(s)*(kz^2*vtzs(s)^2*czj(j)^2+kz*vdsz(s)*kvtsab(iab,s)* ...
      czj(j)+0.5*kx^2*vtpsab(iab,s)^2)/kvtsab2(iab,s) + ...
      k^2*kz^2*vtzs(s)^4/kvtsab(iab,s)^3*czj(j)^3 + ...
      2*k^2*kz*vtzs(s)^2*vdsz(s)/kvtsab2(iab,s)*czj(j)^2 + ...
      (k^2*vdsz(s)^2*kvtsab2(iab,s)+...
      0.5*k^2*kx^2*vtzs(s)^2*vtpsab(iab,s)^2-(lmdTab(iab,s)-1)*...
      kx^2*kz^2*vtzs(s)^2*vtpsab(iab,s)^2 ...
      )/kvtsab(iab,s)^3*czj(j) - (lmdTab(iab,s)-1)*kx^2*kz*...
      vtpsab(iab,s)^2*vdsz(s)/kvtsab2(iab,s);
  
    b11snj(snj)=tmp*p11snj;
%     b11=b11-tmp*p11snj; % the extra wps2(s) term will be add in final step
    b11s(s)=b11s(s)-tmp*p11snj; % 20-11-14 10:41 update
    
    b12snj(snj)=tmp*p12snj;
    b12s(s)=b12s(s)-tmp*p12snj; %
    
    b21snj(snj)=tmp*p21snj;
    b21s(s)=b21s(s)-tmp*p21snj; %
    
    b22snj(snj)=tmp*p22snj;
    b22s(s)=b22s(s)-tmp*p22snj; % the extra wps2(s) term will be add in final step
    
    b13snj(snj)=tmp*p13snj;
    b13s(s)=b13s(s)-tmp*p13snj; %
    
    b31snj(snj)=tmp*p31snj;
    b31s(s)=b31s(s)-tmp*p31snj; %
    
    b23snj(snj)=tmp*p23snj;
    b23s(s)=b23s(s)-tmp*p23snj; %
    
    b32snj(snj)=tmp*p32snj;
    b32s(s)=b32s(s)-tmp*p32snj; %
    
    b33snj(snj)=tmp*p33snj;
    b33s(s)=b33s(s)-tmp*p33snj; % the extra wps2(s) term will be add in final step
      
   else % magnetized species
    
    % 
    % 18-12-20 17:08 Note the czjj(j) term in csnj, which is to support 
    % also kz<0 in Z function, the previous one only correct for kz>0.
    % csnj(snj)=czj(j)*kz*vtzs(s)+kz*vdsz(s)+kx*vdsx(s)-1i*nus(s)+n*wcs(s);
    csnj(snj)=czjj(j)*kz*vtzs(s)+kz*vdsz(s)+kx*vdsx(s)-1i*nus(s)+n*wcs(s);
    cnj=csnj(snj);
    
    for iab=1:nrsab(s) % sum loss cone
        
     if(j==1) % 18-12-17 00:08 only need calculate once
         
      if(vdsr(s)==0) % if no ring beam, still use Besseli
       as2=0.5*asab(iab,s)^2;
       Gamn=besseli(n,as2,1);
       Gamnp=(besseli(n+1,as2,1)+besseli(n-1,as2,1)-2*besseli(n,as2,1))/2;
       Anabb(iab)=4.0/Asab(iab,s)*Gamn/2.0;
       Anab0(iab)=Anabb(iab);
       Bnabb(iab)=4.0/Asab(iab,s)*Gamnp*asab(iab,s)/4.0;
       Bnab0(iab)=Bnabb(iab);
       Cnabb(iab)=4.0/Asab(iab,s)*(Gamn*n^2/as2/4-Gamnp*as2/2.0);
       Cnab0(iab)=Cnabb(iab);
      else
       % 18-12-02 10:00, calculate An, Bn, Cn using numerical integral
       Anabb(iab)=4.0/Asab(iab,s)*funAn(n,asab(iab,s),bsab(iab,s),bsab(iab,s));
       Anab0(iab)=4.0/Asab(iab,s)*funAn(n,asab(iab,s),bsab(iab,s),0);
       Bnabb(iab)=4.0/Asab(iab,s)*funBn(n,asab(iab,s),bsab(iab,s),bsab(iab,s));
       Bnab0(iab)=4.0/Asab(iab,s)*funBn(n,asab(iab,s),bsab(iab,s),0);
       Cnabb(iab)=4.0/Asab(iab,s)*funCn(n,asab(iab,s),bsab(iab,s),bsab(iab,s));
       Cnab0(iab)=4.0/Asab(iab,s)*funCn(n,asab(iab,s),bsab(iab,s),0);
      end
      
      % 18-12-18 17:05
      if(imethod==1)
      p11snj_b=(n*wcs(s)/kx+vdsx(s))*n*wcs(s)/kx*Anabb(iab)/vtpsab(iab,s)^2;
      p22snj_b=(Cnabb(iab)+1i*vdsy(s)*Bnabb(iab)/vtpsab(iab,s));
%       p33snj_b=Anab0(iab);
      p33snj_b=Anab0(iab)/2; % 18-12-18 23:55
%       b11=b11-rsab(iab,s)*wps2(s)*p11snj_b;
%       b22=b22-rsab(iab,s)*wps2(s)*p22snj_b;
%       b33=b33-rsab(iab,s)*wps2(s)*p33snj_b;
      b11s(s)=b11s(s)-rsab(iab,s)*wps2(s)*p11snj_b; % 20-11-14 10:44 update
      b22s(s)=b22s(s)-rsab(iab,s)*wps2(s)*p22snj_b;
      b33s(s)=b33s(s)-rsab(iab,s)*wps2(s)*p33snj_b;
      end

     end
    
    tmp=rsab(iab,s)*wps2(s)*bzj(j)/cnj;
    
    p11snj=(n*wcs(s)/kx+vdsx(s))*(n*wcs(s)/kx*(n*wcs(s)+kx*vdsx(s)-...
      1i*nus(s))*Anabb(iab)/vtpsab(iab,s)^2 + Anab0(iab)*(n*wcs(s)/kx+vdsx(s)...
      )*kz*czjj(j)/vtzs(s));
  
    p12snj=(n*wcs(s)/kx+vdsx(s))*(vdsy(s)*(n*wcs(s)*Anabb(iab)/vtpsab(iab,s)^2 ...
      +Anab0(iab)*kz*czjj(j)/vtzs(s)) + 1i*((n*wcs(s)-1i*nus(s))*Bnabb(iab)/...
      vtpsab(iab,s)+Bnab0(iab)*vtpsab(iab,s)/vtzs(s)*kz*czjj(j)));
  
    p21snj=n*wcs(s)/kx*(n*wcs(s)+kx*vdsx(s)-1i*nus(s))*(-1i* ...
      vtpsab(iab,s)*Bnabb(iab)+vdsy(s)*Anabb(iab))/vtpsab(iab,s)^2 + ...
      (n*wcs(s)/kx+vdsx(s))/vtzs(s)*(-1i*vtpsab(iab,s)*Bnab0(iab)+...
      vdsy(s)*Anab0(iab))*kz*czjj(j);
  
    p22snj=vdsy(s)*n*wcs(s)*(-1i*vtpsab(iab,s)*Bnabb(iab)+vdsy(s)*Anabb(iab) ...
      )/vtpsab(iab,s)^2 +vdsy(s)*(-1i*vtpsab(iab,s)*Bnab0(iab)+...
      vdsy(s)*Anab0(iab))/vtzs(s)*kz*czjj(j) + (n*wcs(s)-1i*nus(s))*(Cnabb(iab)+...
      1i*vdsy(s)*Bnabb(iab)/vtpsab(iab,s)) + vtpsab(iab,s)/vtzs(s)*(...
      vtpsab(iab,s)*Cnab0(iab)+1i*vdsy(s)*Bnab0(iab))*kz*czjj(j);
  
    p13snj=(n*wcs(s)/kx+vdsx(s))*(n*wcs(s)*(czjj(j)*vtzs(s)+...
      vdsz(s))*Anabb(iab)/vtpsab(iab,s)^2 + Anab0(iab)*(kz*czjj(j)^2+(...
      kz*vdsz(s)-1i*nus(s))/vtzs(s)*czjj(j)));
    
    p31snj=n*wcs(s)/kx*(n*wcs(s)+kx*vdsx(s)-1i*nus(s))*(czjj(j)*vtzs(s)+...
      vdsz(s))*Anabb(iab)/vtpsab(iab,s)^2 +Anab0(iab)*(n*wcs(s)/kx+vdsx(s))*(czjj(j)+ ...
      vdsz(s)/vtzs(s))*kz*czjj(j);
    
    p23snj=(vdsy(s)*Anabb(iab)-1i*vtpsab(iab,s)*Bnabb(iab))*(czjj(j)*vtzs(s)+vdsz(s) ...
      )*n*wcs(s)/vtpsab(iab,s)^2 + (vdsy(s)*Anab0(iab)-1i*vtpsab(iab,s)*Bnab0(iab) ...
      )*(kz*czjj(j)+(kz*vdsz(s)-1i*nus(s))/vtzs(s))*czjj(j);
  
    % fixed a bug, 20-10-25 21:15
    p32snj=vdsy(s)*(n*wcs(s)*Anabb(iab)*czjj(j)*vtzs(s)/vtpsab(iab,s)^2+...
      Anab0(iab)*kz*czjj(j)^2) + vdsy(s)*vdsz(s)*(n*wcs(s)*Anabb(iab)/ ... % kz^2 -> kz
      vtpsab(iab,s)^2+Anab0(iab)*kz*czjj(j)/vtzs(s)) + 1i*((n*wcs(s)-1i*nus(s))* ...
      Bnabb(iab)*vtzs(s)/vtpsab(iab,s)+Bnab0(iab)*kz*vtpsab(iab,s)*czjj(j))*(czjj(j)+...
      vdsz(s)/vtzs(s));
%     p32snj=vdsy(s)*(n*wcs(s)*Anabb(iab)*czjj(j)*vtzs(s)/vtpsab(iab,s)^2+...
%       Anab0(iab)*kz^2*czjj(j)^2) + vdsy(s)*vdsz(s)*(n*wcs(s)*Anabb(iab)/ ...
%       vtpsab(iab,s)^2+Anab0(iab)*kz*czjj(j)/vtzs(s)) + 1i*((n*wcs(s)-1i*nus(s))* ...
%       Bnabb(iab)*vtzs(s)/vtpsab(iab,s)+Bnab0(iab)*kz*vtpsab(iab,s)*czjj(j))*(czjj(j)+...
%       vdsz(s)/vtzs(s));  % 20-10-25 21:15 kz^2 should be kz
  
    p33snj=( n*wcs(s)*(vtzs(s)^2*czjj(j)+vtzs(s)*vdsz(s))*Anabb(iab)/ ...
      vtpsab(iab,s)^2 + Anab0(iab)*czjj(j)*(czjj(j)*kz*vtzs(s)+(kz*vdsz(s)-...
      1i*nus(s))) )*(czjj(j)+ vdsz(s)/vtzs(s));
    
    b11snj(snj)=b11snj(snj)+tmp*p11snj;
%     b11=b11-tmp*p11snj; % the extra wps2(s) term will be add in final step
    b11s(s)=b11s(s)-tmp*p11snj; %  20-11-14  10:45 update
    
    b12snj(snj)=b12snj(snj)+tmp*p12snj;
%     b12=b12-b12snj(snj); % wrong if nrsab=2! 18-12-18 09:05
%     b12=b12-tmp*p12snj;
    b12s(s)=b12s(s)-tmp*p12snj; %  20-11-14  10:45 update
    
    b21snj(snj)=b21snj(snj)+tmp*p21snj;
    b21s(s)=b21s(s)-tmp*p21snj; %
    
    b22snj(snj)=b22snj(snj)+tmp*p22snj;
    b22s(s)=b22s(s)-tmp*p22snj; % the extra wps2(s) term will be add in final step
    
    b13snj(snj)=b13snj(snj)+tmp*p13snj;
    b13s(s)=b13s(s)-tmp*p13snj; %
    
    b31snj(snj)=b31snj(snj)+tmp*p31snj;
    b31s(s)=b31s(s)-tmp*p31snj; %
    
    b23snj(snj)=b23snj(snj)+tmp*p23snj;
    b23s(s)=b23s(s)-tmp*p23snj; %
    
    b32snj(snj)=b32snj(snj)+tmp*p32snj;
    b32s(s)=b32s(s)-tmp*p32snj; %
    
    b33snj(snj)=b33snj(snj)+tmp*p33snj;
    b33s(s)=b33s(s)-tmp*p33snj; % the extra wps2(s) term will be add in final step
      
    end
   end
  end
 end
 
 % extra term for diag elements, 20-11-14 10:49 update
 if(imethod==0) % seems need large N to convergent
     b11s(s)=b11s(s)-wps2(s);
     b22s(s)=b22s(s)-wps2(s);
     b33s(s)=b33s(s)-wps2(s);
 else % for unmagnetized species
     b11s(s)=b11s(s)-wps2(intersect(s,jds)); % to check
     b22s(s)=b22s(s)-wps2(intersect(s,jds));
     b33s(s)=b33s(s)-wps2(intersect(s,jds));
 end

end

% % extra term for diag elements, 2018-12-16 10:55
% if(imethod==0) % seems need large N to convergent
%  b11=b11-sum(wps2);
%  b22=b22-sum(wps2);
%  b33=b33-sum(wps2);
% else % for unmagnetized species, 2018-12-18 00:30
%  b11=b11-sum(wps2(jds));
%  b22=b22-sum(wps2(jds));
%  b33=b33-sum(wps2(jds));  
% end

for snj=1:SNJ % set the eigen matrix
  jjx=snj+0*SNJ1;
  jjy=snj+1*SNJ1;
  jjz=snj+2*SNJ1;
  % v_snjx
  M=M+sparse(jjx,jjx,csnj(snj),NN,NN)+...
    sparse(jjx,SNJ3+1,b11snj(snj),NN,NN)+...
    sparse(jjx,SNJ3+2,b12snj(snj),NN,NN)+...
    sparse(jjx,SNJ3+3,b13snj(snj),NN,NN);

  % v_snjy
  M=M+sparse(jjy,jjy,csnj(snj),NN,NN)+...
    sparse(jjy,SNJ3+1,b21snj(snj),NN,NN)+...
    sparse(jjy,SNJ3+2,b22snj(snj),NN,NN)+...
    sparse(jjy,SNJ3+3,b23snj(snj),NN,NN);

  % v_snjz
  M=M+sparse(jjz,jjz,csnj(snj),NN,NN)+...
    sparse(jjz,SNJ3+1,b31snj(snj),NN,NN)+...
    sparse(jjz,SNJ3+2,b32snj(snj),NN,NN)+...
    sparse(jjz,SNJ3+3,b33snj(snj),NN,NN);

end

% E(J), J_{x,y,z}=j_{x,y,z}+sum(v_snj{x,y,z})
tp=-1;
jj=(0*SNJ1+1):(1*SNJ1); ii=0.*jj+SNJ3+1; M=M+sparse(ii,jj,tp,NN,NN);
jj=(1*SNJ1+1):(2*SNJ1); ii=0.*jj+SNJ3+2; M=M+sparse(ii,jj,tp,NN,NN);
jj=(2*SNJ1+1):(3*SNJ1); ii=0.*jj+SNJ3+3; M=M+sparse(ii,jj,tp,NN,NN);

% % jx(E), jy(E), jz(E)
% M=M+sparse(1*SNJ1,SNJ3+1,b11,NN,NN)+...
%   sparse(1*SNJ1,SNJ3+2,b12,NN,NN)+...
%   sparse(1*SNJ1,SNJ3+3,b13,NN,NN)+...
%   sparse(2*SNJ1,SNJ3+1,b21,NN,NN)+...
%   sparse(2*SNJ1,SNJ3+2,b22,NN,NN)+...
%   sparse(2*SNJ1,SNJ3+3,b23,NN,NN)+...
%   sparse(3*SNJ1,SNJ3+1,b31,NN,NN)+...
%   sparse(3*SNJ1,SNJ3+2,b32,NN,NN)+...
%   sparse(3*SNJ1,SNJ3+3,b33,NN,NN);

% jsx(E), jsy(E), jsz(E), 20-11-14 11:02 separate each species
for s=1:S
  M=M+sparse(1*SNJ1-S+s,SNJ3+1,b11s(s),NN,NN)+...
    sparse(1*SNJ1-S+s,SNJ3+2,b12s(s),NN,NN)+...
    sparse(1*SNJ1-S+s,SNJ3+3,b13s(s),NN,NN)+...
    sparse(2*SNJ1-S+s,SNJ3+1,b21s(s),NN,NN)+...
    sparse(2*SNJ1-S+s,SNJ3+2,b22s(s),NN,NN)+...
    sparse(2*SNJ1-S+s,SNJ3+3,b23s(s),NN,NN)+...
    sparse(3*SNJ1-S+s,SNJ3+1,b31s(s),NN,NN)+...
    sparse(3*SNJ1-S+s,SNJ3+2,b32s(s),NN,NN)+...
    sparse(3*SNJ1-S+s,SNJ3+3,b33s(s),NN,NN);
end

% E(B)
M=M+sparse(SNJ3+1,SNJ3+5,c2*kz,NN,NN)+...
  sparse(SNJ3+2,SNJ3+4,-c2*kz,NN,NN)+...
  sparse(SNJ3+2,SNJ3+6,c2*kx,NN,NN)+...
  sparse(SNJ3+3,SNJ3+5,-c2*kx,NN,NN);

% B(E)
M=M+sparse(SNJ3+4,SNJ3+2,-kz,NN,NN)+...
  sparse(SNJ3+5,SNJ3+1,kz,NN,NN)+...
  sparse(SNJ3+5,SNJ3+3,-kx,NN,NN)+...
  sparse(SNJ3+6,SNJ3+2,kx,NN,NN);

% for Darwin model, use lambda*MB*X=M*X, 18-12-17 22:55
if(iem==2) % to check
 MB=speye(NN);
 % modify the omega*E equation
 MB=MB+sparse(SNJ3+1,SNJ3+1,kx^2/k^2,NN,NN)+...
  sparse(SNJ3+3,SNJ3+1,kx*kz/k^2,NN,NN)+...
  sparse(SNJ3+1,SNJ3+3,kx*kx/k^2,NN,NN)+...
  sparse(SNJ3+3,SNJ3+3,kz^2/k^2,NN,NN)-...
  sparse(SNJ3+1,SNJ3+1,1,NN,NN)-...
  sparse(SNJ3+2,SNJ3+2,1,NN,NN)-...
  sparse(SNJ3+3,SNJ3+3,1,NN,NN);
 
end

% % -- main program end of set the matrix elements --
