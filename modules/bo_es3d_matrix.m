% 18-12-16 08:58, Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China.
% This file set the BO (PDRK) kernel matrix elements, of electrostatic
% 3D case. This version support also kz<=0.

% % -- main program begin to set the matrix elements --

k=sqrt(kz^2+kx^2);
asab=kx*rhocsab*sqrt(2); % as=kx*vtps./wcs, note a sqrt(2), 18-11-29 16:34 

kvtsab2(1,:)=((kx.*vtpsab(1,:)).^2+(kz.*vtzs).^2);
kvtsab2(2,:)=((kx.*vtpsab(2,:)).^2+(kz.*vtzs).^2);
kvtsab=sqrt(kvtsab2); % define new (k*vts)^2 for unmagnetized species

M=sparse(NN,NN);
snj=0;

% initialize
csnj=zeros(1,SNJ);
bsnj=csnj.*0;

if(kz<0) % 18-12-22 00:54
  czjj=-czj;
else
  czjj=czj;
end

Anabb=zeros(1,2); Anab0=0.*Anabb;
for s=1:S % species
 if(ismember(s,jds))
  Nsss=1:Ns(s); % unmagnetized species, for different sigma
%   Nsss=1:Nss(s); % i.e., should be 1 not nrsab?
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
     
    bsnj(snj)=2*rsab(iab,s)*wps2(s)*bzj(j)*czj(j)/kvtsab(iab,s);
    
   else % magnetized species
    
    csnj(snj)=czjj(j)*kz*vtzs(s)+kz*vdsz(s)+kx*vdsx(s)-1i*nus(s)+n*wcs(s);
    cnj=csnj(snj);
    
    for iab=1:nrsab(s) % sum loss cone
     if(j==1) % 18-12-17 00:08 only need calculate once
         
      if(vdsr(s)==0) % if no ring beam, still use Besseli
       as2=0.5*asab(iab,s)^2;
       Gamn=besseli(n,as2,1);
       %Gamnp=(besseli(n+1,as2,1)+besseli(n-1,as2,1)-2*besseli(n,as2,1))/2;
       Anabb(iab)=4.0/Asab(iab,s)*Gamn/2.0;
       Anab0(iab)=Anabb(iab);
      else
       Anabb(iab)=4.0/Asab(iab,s)*funAn(n,asab(iab,s),bsab(iab,s),bsab(iab,s));
       Anab0(iab)=4.0/Asab(iab,s)*funAn(n,asab(iab,s),bsab(iab,s),0);
      end
      
     end
     
     bsnj(snj)=bsnj(snj)+rsab(iab,s)*wps2(s)*bzj(j)*( Anab0(iab)*kz/vtzs(s)*czjj(j)+...
         Anabb(iab)*n*wcs(s)/vtpsab(iab,s)^2 )/k^2;
 
    end
    
   end
  end
 end
end

for snj=1:SNJ % set the eigen matrix
    
  % (n_snj, n_snj)
  M=M+sparse(snj,snj,csnj(snj),NN,NN);
  
  % (n_snj, E)
  M=M+sparse(snj,NN,bsnj(snj),NN,NN);

  % (E, n_snj)
  M=M+sparse(NN,snj,-csnj(snj),NN,NN);

end

% (E, E)
M=M+sparse(NN,NN,-sum(bsnj),NN,NN);

% % -- main program end of set the matrix elements --
