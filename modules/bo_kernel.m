% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2014-06-01 17:00
% bo_em3d.m, Plasma Dispersion Relastion solver, kinetic, EM3D,
% bi-Maxwellian equilibrium distribution with parallel drift.
% Transform to matrix eigenvalue problem lambda*X=M*X.
% J-pole approximation for Z(zeta)=sum(b_j/(zeta-c_j))
%
% sparse, eigs
%
% 2014-08-22 01:09, support T_perp/T_para \neq 1,  vds \neq 0
%
% Ref:
%  [Xie2016] H. S. Xie & Y. Xiao, PDRK: A General Kinetic Dispersion
%    Relation Solver for Magnetized Plasma, Plasma Science and Technology,
%    Vol.18, No.2, p97 (2016). (Also 2014arXiv, note: have several typos)
%
% Documents/codes/Erratum: 
%  http://hsxie.me/codes/pdrk/, https://github.com/hsxie/pdrk
%
% 18-09-28 12:42 We find that omega_cs=-qB/m in Miyamoto2004, not the
% standard omega_cs=qB/m. This sign difference affect several terms. This
% version we have fixed the bugs of the sign of omega_cs.
% 18-10-01 09:43 based on new derivations, fixed a bug in b33 term, where a
% wcs^2 is missed.
% 16:01 Benchmark with Gary1993 Fig.7.4 and several other cases ok.
%
% Ref (typos/bugs fixed):
%  H. S. Xie, Detailed Derivations of PDRK-EM3D Equations, 2018-10-03 (
%  10 pages).
% 
% 18-10-03 13:14 This version should be bugs free now.
% 18-10-13 10:37 update to with loss-cone distribution, which thus can
% handle all the same input cases as WHAMP [Ronnmark1982] code
% 18-10-27 The major update version to github, of loss-cone, EM3D, ES3D.
% 18-12-19 10:49 A second major update to a most powerful version, with
% bi-Maxwellian, parallel drift vdz, penpendicular drift across magnetic
% field (vdx, vdy), ring beam vdr, loss cone, and Krook collision. EM3D,
% ES3D, and Darwin. Unmagnetized and magnetized.
% ------------------------------------------------------------------------
% This file is the kernel part of BO (PDRK) code.
% Modify sp and nw if you only need one or several solutions around the
% initial guess value(es) wg0.
% 18-12-19 16:46 add to support sp=2
% 18-12-20 17:02 update the matrix.m, to support kz<0

wwp=zeros(npa,npb,nw);
[ppa,ppb]=ndgrid(pa,pb);
kxx=0.*ppa; % for store kx
kzz=0.*ppa; % for store kz
kk=0.*ppa; % for store k
tt=0.*ppa; % for store theta
%tic;
for jpa = 1:npa % scan 1st parameter a
% you can try 'parfor' to parallization scan:
%  5. If you need parallize the scan of k or theta, please modify the
% kernel.m to use 'parfor', by copy the matrix.m to embed to kernel.m.
  jpa
  
  for jpb = 1:npb % scan 2nd parameter b
    
    % ------- to update 2018-10-18 20:10 --------
    if(iloga==0)
      par(ipa)=pa(jpa);
    else
      par(ipa)=10^(pa(jpa));
    end
    if(ipa~=ipb)
      if(ilogb==0)
        par(ipb)=pb(jpb);
      else
        par(ipb)=10^(pb(jpb));  
      end
    end
    
    k=abs(par(1))*kn; % k only >=0
    theta=par(2); % the angle between k and B0
    if(ipa>2 && ipb>2) % !! to update
      kz=par(3)*kn;
      kx=par(4)*kn;
    else
      kz=cos(theta*pi/180)*k;
      kx=sin(theta*pi/180)*k;
    end
    k=sqrt(kx*kx+kz*kz); % update k
    theta=angle(kz+1i*kx)*180/pi; % update theta
    % !! if you want scan other parameters, modify here, e.g.,
    % betasz=par(5); % or vds(1)=par(5); % etc
    % ---------------
    
    if(iem==1 || iem==2) % electromagnetic run, or Darwin run
      run bo_em3d_matrix; % modify here if you use 'parfor'
    elseif(iem==0) % electrostatic run
      run bo_es3d_matrix;
    elseif(iem==3) % electromagnetic fluid run, 20-12-23 07:27
      run bo_fluid_matrix;
    end
    
    % calculate the solution use either eig() or eigs()
    if(sp==0) % solve all the roots, slow for large N
      % d=eig(full(M));
      % it's found in pdrf (fluid version) that vpa(eig) may more accurate
      if(iem==2) % 18-12-17 22:59 for Darwin model
%        d0=eig(full(M),full(MB));d=double(d0);
       d0=vpa(eig(full(M),full(MB)),16);d=double(d0);
%        d0=vpa(eig(full(MB),full(M)),16);d=double(d0);
      else
%        d0=vpa(eig(full(M)),16);d=double(d0);
       d0=eig(full(M));d=double(d0);
      end
    else % solve only nw solutions, fast for large N, similar to WHAMP
      % d=eigs(M,nw,wg);
      if((iem==1 || iem==2 || iem==3) && iout==2 && icalp==1) % for bo_output.m
        % calcute also the polarizations
        disp(['jpa=',num2str(jpa),', jpb=',num2str(jpb)]);
        eps=1i*1e-8; % 18-10-06 10:29
%         wg=wws(jpa,jpb,jpl)*wn+eps; % add eps to avoid singular warning
        wg=(wws(jpa,jpb,jpl)+eps)*wn; % 20-11-14 08:46 update
        %d0=vpa(eigs(M,nw,wg),16);d=double(d0);
        if(iem==2) % 18-12-17 22:59 for Darwin model
         [V,D]=eigs(M,MB,nw,wg);
        else
         [V,D]=eigs(M,nw,wg);
        end
        
        %[V0,D0]=vpa(eigs(M,nw,wg)); V=double(V0); D=double(D0);
        dEB=V((NN-5):NN,1); % polarization data (dE,dB)
        wd=D(1,1); % wd=omega solved from eigen matrix
        
        % To do: update the normalization of polarizations
        if(abs(dEB(1))<0.1*sqrt(abs(dEB(1))^2+abs(dEB(2))^2+abs(dEB(3))^2))
            tmp=dEB(3); % 21-04-24 20:26 if dEx is small, normalize by dEz, to test
        else
            tmp=dEB(1);
        end
        dV=V(:,1)/tmp;
        dEB=dEB/tmp;
%         dV=V(:,1)/dEB(1); % 20-11-15 10:10
%         dEB=dEB/dEB(1); % 18-10-19 17:04
%         ctmp=sqrt(real(dEB(1))^2+real(dEB(2))^2+real(dEB(3))^2)+...
%             1i*sqrt(imag(dEB(1))^2+imag(dEB(2))^2+imag(dEB(3))^2);
%         ctmp=sqrt(real(dEB(1))^2+real(dEB(2))^2+real(dEB(3))^2)+...
%             imag(dEB(1))^2+imag(dEB(2))^2+imag(dEB(3))^2; % wrong, 20-12-15 17:21
        ctmp=sqrt(real(dEB(1))^2+real(dEB(2))^2+real(dEB(3))^2+...
            imag(dEB(1))^2+imag(dEB(2))^2+imag(dEB(3))^2); % 20-12-15 17:21
        dEB=dEB/ctmp; % normalized (E,B) % 18-10-18 07:14
        dV=dV/ctmp; % normalized V % 20-11-15 10:10
        wws2(jpa,jpb,jpl)=wd/wn; % store new wws, in case
        
        dEx=dEB(1); % Ex
        dEy=dEB(2); % Ey
        dEz=dEB(3); % Ez
        dBx=dEB(4); % Bx
        dBy=dEB(5); % By
        dBz=dEB(6); % Bz
        
        % npolaf=50; npolas=30;
        % Output of polarization quantities Pola(:,:,:,ipola) for ipola=
        % 1 - Ex (complex)
        % 2 - Ey (complex)
        % 3 - Ez (complex)
        % 4 - Bx (complex)
        % 5 - By (complex)
        % 6 - Bz (complex)
        % 7 - UE = electric field energy density (J/m^3)
        % 8 - UB = magnetic field energy density (J/m^3)
        % 9 - UE/(UE+UB) = fraction of field energy in the electric field
        % 10 - |Ex|^2/|E|^2 = fraction of electric field energy in Ex
        % 11 - |Ey|^2/|E|^2 = fraction of electric field energy in Ey
        % 12 - |Ez|^2/|E|^2 = fraction of electric field energy in Ez
        %     = |E_parallel|^2/|E|^2
        % 13 - | k dot E |^2 / ( k^2 * |E|^2 ) = a measure of how electrostatic the
        %     electric field is
        % 14 - | k dot E |^2 / ( | k dot E |^2 + | k cross E |^2 ) = another measure
        %     of how electrostatic the electric field is
        % 15 - |Bx|^2/|B|^2 = fraction of magnetic field energy in Bx
        % 16 - |By|^2/|B|^2 = fraction of magnetic field energy in By
        % 17 - |Bz|^2/|B|^2 = fraction of magnetic field energy in Bz
        %     = |B_parallel|^2/|B|^2
        % 18 - Magnetic ellipticity = -1/0/1 for left-hand polarized, linear
        %     polarized, or right hand polarized waves
        % 19 - thetaB = Angle of major axis of the magnetic ellipse
        % 20 - |By/Bx|
        % 21 - angle(By/Bx)
        % npolaf-3 - kx
        % npolaf-2 - kz
        % npolaf-1 - k
        % npolaf - frequency (complex)
        % npolaf+1 - Jx
        % npolaf+2 - Jy
        % npolaf+3 - Jz
        % npolaf+3+npolas*(s-1)+1:7 - Jsx, Jsy, Jsz, ns, vsx, vsy, vsz
  
        % 20-11-15 16:58 Calculate other polarization quantities.
        % partly copy from R. Denton's version
        dExs=dEx*conj(dEx); % Ex^2, here 's' means square
        dEys=dEy*conj(dEy); % Ey^2
        dEzs=dEz*conj(dEz); % Ez^2
        dEs=dExs+dEys+dEzs; % E^2
        dUE=0.5*dEs*epsilon0*0.5; % electric energy; Note: extra factor 0.5 for time average
        dBxs=dBx*conj(dBx); % Bx^2
        dBys=dBy*conj(dBy); % By^2
        dBzs=dBz*conj(dBz); % Bz^2
        dBs=dBxs+dBys+dBzs; % B^2
        dUB=0.5*dBs/mu0*0.5; % magnetic energy; Note: extra factor 0.5 for time average
        
        Pola(jpa,jpb,jpl,1:6)=dEB; % store dEB
        Pola(jpa,jpb,jpl,7)=dUE;
        Pola(jpa,jpb,jpl,8)=dUB;
        Pola(jpa,jpb,jpl,9)=dUE/(dUE+dUB); % UE/(UE+UB)
        Pola(jpa,jpb,jpl,10)=dExs/dEs; % Ex^2/E^2
        Pola(jpa,jpb,jpl,11)=dEys/dEs; % Ey^2/E^2
        Pola(jpa,jpb,jpl,12)=dEzs/dEs; % Ez^2/E^2=E_parallel^2/E^2
        
        % kdotEs = kx^2*dExs + kz^2*dEzs; % |k dot E|^2 from Denton, wrong?
        kdotdEs=(kx*dEx+kz*dEz)*conj(kx*dEx+kz*dEz); % |k dot E|^2
        % kcrossdEs=abs(sqrt((-kz*dEy)^2+(kz*dEx-kx*dEz)^2+(kx*dEy)^2)); % to check
        kcrossdEs=abs((-kz*dEy)^2)+abs((kz*dEx-kx*dEz)^2)+abs((kx*dEy)^2);
        kdEs=(kx^2+kz^2)*dEs; % k^2 * |E|^2
        Pola(jpa,jpb,jpl,13)=kdotdEs/kdEs; % |k dot E|^2/(k^2 * |E|^2)
        Pola(jpa,jpb,jpl,14)=kdotdEs/(kdotdEs+kcrossdEs); % |k dot E|^2 / ( kdotE + |k cross E|^2 )
        Pola(jpa,jpb,jpl,15)=dBxs/dBs; % Bx^2/B^2
        Pola(jpa,jpb,jpl,16)=dBys/dBs; % By^2/B^2
        Pola(jpa,jpb,jpl,17)=dBzs/dBs; % Bz^2/B^2=B_parallel^2/B^2
        
        % ellipticity and angle of B major axis 20-12-24 19:57
        alphaB=dBy/dBx; % By/Bx
        epsilonB=(abs(1-1i*alphaB)-abs(1+1i*alphaB))/(abs(1-1i*alphaB)+...
            abs(1+1i*alphaB)+1e-100);
        absalphaB=abs(alphaB);
        phiB=angle(alphaB);
        thetaA=0.5*acot(-(cos(2*phiB)+1/absalphaB^2 )/( sin(2*phiB) ) );
        if(cos(2*thetaA)+absalphaB^2*cos(2*(phiB+thetaA))<0)
            thetaA=thetaA+pi/2;
        end
        thetaB=atan((absalphaB*cos(thetaA+phiB))/cos(thetaA));
        
        Pola(jpa,jpb,jpl,18)=epsilonB; % ellipticity
        Pola(jpa,jpb,jpl,19)=thetaB; % angle of B major axis
        Pola(jpa,jpb,jpl,20)=absalphaB; % |By/Bx|
        Pola(jpa,jpb,jpl,21)=phiB; % angle(By/Bx)
        
        Pola(jpa,jpb,jpl,npolaf-3)=kx; % kx without normalized
        Pola(jpa,jpb,jpl,npolaf-2)=kz; % kz without normalized
        Pola(jpa,jpb,jpl,npolaf-1)=sqrt(kx^2+kz^2); % k without normalized
        Pola(jpa,jpb,jpl,npolaf)=wd; % frequency without normalized
        
        dJsx=zeros(1,S); dJsy=0.*dJsx; dJsz=0.*dJsx;
        if(iem==1 || iem==2) % kinetic EM version
          % 20-11-15 00:32 new feature to calculate dJs_xyz of each species
          for s=1:S
            indjs=[ids0(s):1:ids1(s),ids2(s)]; % index [v_snj, js] for species #s
            dJsx(s)=-sum(dV(0*SNJ1+indjs))*(1i*epsilon0);
            dJsy(s)=-sum(dV(1*SNJ1+indjs))*(1i*epsilon0);
            dJsz(s)=-sum(dV(2*SNJ1+indjs))*(1i*epsilon0);
          end
          dJx=sum(dJsx); dJy=sum(dJsy); dJz=sum(dJsz);
          
          % calculate density & velocity polarization of each species from
          % dJs_xyz, where the belove equations are from the comparison with
          % fluid model, i.e., dJs=qs*(dns*vds0+ns0*dvs)
          tmps=(kx*dJsx+kz*dJsz)./wd;
          dns=tmps./qs;
          
          withdn=1; % 2020-12-16 12:27
          dvsx=(dJsx-withdn*tmps.*vdsx)./(qs.*ns0);
          dvsy=(dJsy-withdn*tmps.*vdsy)./(qs.*ns0);
          dvsz=(dJsz-withdn*tmps.*vdsz)./(qs.*ns0);
        elseif(iem==3) % multi-fluid EM version, 20-12-24 20:29
          for s=1:S
            dns(s)=dV(4*(s-1)+1);
            dvsx(s)=dV(4*(s-1)+2);
            dvsy(s)=dV(4*(s-1)+3);
            dvsz(s)=dV(4*(s-1)+4);
          end
          dJsx=qs.*(dns.*vdsx+ns0.*dvsx);
          dJsy=qs.*(dns.*vdsy+ns0.*dvsy);
          dJsz=qs.*(dns.*vdsz+ns0.*dvsz);
          dJx=sum(dJsx); dJy=sum(dJsy); dJz=sum(dJsz);
        end
        
        % store dJ_xyz (3), dJ_sxyz (3*S), dn_s (S), dv_sxyz (3*S)
        Pola(jpa,jpb,jpl,npolaf+1)=dJx;
        Pola(jpa,jpb,jpl,npolaf+2)=dJy;
        Pola(jpa,jpb,jpl,npolaf+3)=dJz;
        for s=1:S
          idps=npolaf+3+npolas*(s-1);
          Pola(jpa,jpb,jpl,idps+1)=dJsx(s);
          Pola(jpa,jpb,jpl,idps+2)=dJsy(s);
          Pola(jpa,jpb,jpl,idps+3)=dJsz(s);
          Pola(jpa,jpb,jpl,idps+4)=dns(s);
          Pola(jpa,jpb,jpl,idps+5)=dvsx(s);
          Pola(jpa,jpb,jpl,idps+6)=dvsy(s);
          Pola(jpa,jpb,jpl,idps+7)=dvsz(s);
        end

        if(iem==1 || iem==2) % kinetic EM version
          % 20-11-15 04:23 calculate 3-by-3 tensors sigmas, Qs, K, D(w,k) etc
          % This part can also be used to solve the dispersion realtion by
          % conventional root finding method |D(w,k)|=0 using such as fsolve().
          jtensor=1; % jtensor=1, calculate these tensors, otherwise not
          if(jtensor==1)
            sigmawks=zeros(3,3,S);
            for s=1:S
              inds=ids0(s):1:ids1(s);
              sigmawks(1,1,s)=b11s(s)/wd+sum(b11snj(inds)./(wd-csnj(inds)));
              sigmawks(1,2,s)=b12s(s)/wd+sum(b12snj(inds)./(wd-csnj(inds)));
              sigmawks(1,3,s)=b13s(s)/wd+sum(b13snj(inds)./(wd-csnj(inds)));
              sigmawks(2,1,s)=b21s(s)/wd+sum(b21snj(inds)./(wd-csnj(inds)));
              sigmawks(2,2,s)=b22s(s)/wd+sum(b22snj(inds)./(wd-csnj(inds)));
              sigmawks(2,3,s)=b23s(s)/wd+sum(b23snj(inds)./(wd-csnj(inds)));
              sigmawks(3,1,s)=b31s(s)/wd+sum(b31snj(inds)./(wd-csnj(inds)));
              sigmawks(3,2,s)=b32s(s)/wd+sum(b32snj(inds)./(wd-csnj(inds)));
              sigmawks(3,3,s)=b33s(s)/wd+sum(b33snj(inds)./(wd-csnj(inds)));
            end
            sigmawks=-1i*epsilon0*sigmawks; % conductivity tensor of each species s
            
            Qwks=-sigmawks/(1i*wd*epsilon0);
            sigmawk=sum(sigmawks,3); % conductivity tensor
            Qwk=sum(Qwks,3);
            Kwk=eye(3)+Qwk; % dielectric tensor
            Dwk=Kwk+(kron([kx,0,kz],[kx;0;kz])/(kx^2+kz^2)-eye(3))*((kx^2+...
                kz^2)*c2/wd^2+(iem==2)); % dispersion tensor, D(w,k){\cdot}E=0
            
            % 20-11-15 05:19 check each relations
            dBx0=-kz*dEy/wd; % w*E=kxB
            dBy0=(kz*dEx-kx*dEz)/wd;
            dBz0=kx*dEy/wd;
            if(iem==1)
              dJx0=(wd*dEx-c2*kz*dBy0)*(1i*epsilon0);
              dJy0=(wd*dEy-c2*(kx*dBz0-kz*dBx0))*(1i*epsilon0);
              dJz0=(wd*dEz+c2*kx*dBy0)*(1i*epsilon0);
            elseif(iem==2) % Darwin model
              dJx0=(wd*kx*(kx*dEx+kz*dEz)/(kx^2+kz^2)-c2*kz*dBy0)*(1i*epsilon0);
              dJy0=(wd*0*dEy-c2*(kx*dBz0-kz*dBx0))*(1i*epsilon0);
              dJz0=(wd*kz*(kx*dEx+kz*dEz)/(kx^2+kz^2)+c2*kx*dBy0)*(1i*epsilon0);
            end
            
            dJs1=zeros(3,S); dJsx1=zeros(1,S); dJsy1=dJsx1; dJsz1=dJsx1;
            for s=1:S
              % dJs1(:,s)=squeeze(sigmawks(s))*[dEx;dEy;dEz];  % wrong!
              dJs1(:,s)=squeeze(sigmawks(:,:,s))*[dEx;dEy;dEz]; % 20-11-15 12:17
              dJsx1(s)=dJs1(1,s);
              dJsy1(s)=dJs1(2,s);
              dJsz1(s)=dJs1(3,s);
            end
            dJx1=sum(dJsx1); dJy1=sum(dJsy1); dJz1=sum(dJsz1);
            
            % similarly to do: add calculation of dn_s to es3d version (iem=0)
            
            % dJs1 and dJs should be the same
            % dB0 and dB should be the same
            % dJ0, dJ1 and dJ should be the same
            
            fDE=Dwk*[dEx;dEy;dEz]; % D(w,k){\cdot}E, should be around zero
            DetDwk=det(Dwk); % |D(w,k)|/(max(max(abs(Dwk)))) should be around zero
            
            % to do: store this part for each (kx,kz)
          end
        end
        
        % 20-12-24 20:03 not yet, to add the group velocity relevant quantities
        Pola(jpa,jpb,jpl,22)=0.0; % vgx
        Pola(jpa,jpb,jpl,23)=0.0; % vgz
        Pola(jpa,jpb,jpl,24)=0.0; % Sx
        Pola(jpa,jpb,jpl,25)=0.0; % Sz
        Pola(jpa,jpb,jpl,26)=0.0; % S
        Pola(jpa,jpb,jpl,27)=0.0; % |n|

      else
       if(sp==1)
        if(iem==2)
          d0=eigs(M,MB,nw,wg(1)*0.999);d=double(d0);
        else
          d0=vpa(eigs(M,nw,wg(1)*0.999),16);d=double(d0);
%         d0=vpa(eigs(M,nw,'li'),16);d=double(d0); % 18-10-29 16:42
        end
       else % sp=2, 18-12-19 15:54
        d=0.*wg;
        % note: the meaning of nw here is different from sp=0,1
        for jwg=1:nw
         if(iem==2)
           d00=eigs(M,MB,1,wg(jwg)*0.999);d0=double(d00);
         else
           d00=vpa(eigs(M,1,wg(jwg)*0.999),16);d0=double(d00);
         end
         d(jwg)=d0(1);
        end
       end
      end
    end
    omega=d;
    [wi,ind]=sort(imag(omega),'descend'); % sort the solution by growth rate
    w=omega(ind);
    if(sp==2 && icalp~=1)
     wg=d; % update the initial guess for pa(jpa,jpb+1) use previous solution
    else
     wg=w(1);
    end
    
    if(jpb==1) % initial guess for pb(jpa+1,jpb)
      wg0=wg;
%       wg0=wg; % 18-10-29 15:56
    end

    wwp(jpa,jpb,1:nw)=w(1:nw);
    
    if(icalp==0) % store k, theta in the first run
      kxx(jpa,jpb)=kx;
      kzz(jpa,jpb)=kz;
      kk(jpa,jpb)=k;
      tt(jpa,jpb)=theta;
    end
  end %jpb
  
  wg=wg0; % update the new initial guess
  
end %jpa
%runtime=toc;
if(icalp==0)
  ww=wwp; % store all the solutions in the first run
end
