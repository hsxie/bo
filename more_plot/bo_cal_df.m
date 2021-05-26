% 21-04-24 06:29, Hua-sheng XIE, huashengxie@gmail.com, ENN, China.
% This file calculate f1(k,w,vx,vy,vz) for BO (PDRK) electromagnetic 3D case.
% Eq.(54) of Xie2019CPC paper

s=2;

id=55;
% w=wws(id,1)*wn;
Polatmp=squeeze(Pola(id,1,1,:));
w=Polatmp(npolaf);
kx=Polatmp(npolaf-3);
kz=Polatmp(npolaf-2);
k=Polatmp(npolaf-1);
Ex=Polatmp(1); Ey=Polatmp(2); Ez=Polatmp(3);
Jx=Polatmp(npolaf+1); Jy=Polatmp(npolaf+2); Jz=Polatmp(npolaf+3);

idps=npolaf+3+npolas*(s-1);
Jsx=Polatmp(idps+1);
Jsy=Polatmp(idps+2);
Jsz=Polatmp(idps+3);
ns=Polatmp(idps+4);

%%
% vx=1.0;vy=2; vz=5;

[vxx,vyy,vzz]=ndgrid((-3.0:0.2:3.0)*vtps(s),(-3.0:0.2:3.0)*vtps(s),...
    (-4:0.01:4)*vtzs(s));

[nvx,nvy,nvz]=size(vxx);
dfs=0.*vxx;
for  jvx=1:nvx
    jvx
    for jvy=1:nvy
        for jvz=1:nvz
            vx=vxx(jvx,jvy,jvz);
            vy=vyy(jvx,jvy,jvz);
            vz=vzz(jvx,jvy,jvz);
vparp=vz;
vperp=sqrt((vx-vdsx(s)).^2+(vy-vdsy(s)).^2);
phip=atan2((vy-vdsy(s)),(vx-vdsx(s)));
% vparp=vpar; vperp=vper;
fs0par=1/(sqrt(pi)*vtzs(s)).*exp(-((vparp-vdsz(s))/vtzs(s)).^2);
fs0pera=1/(pi*vtpsab(1,s).^2).*(...
    rsab(1,s)./(Asab(1,s)).*exp(-((vperp-vdsr(s))/vtpsab(1,s)).^2));
fs0perb=1/(pi*vtpsab(2,s).^2).*(...
    rsab(2,s)./(Asab(2,s)).*exp(-((vperp-vdsr(s))/vtpsab(2,s)).^2));
fs0a=fs0par.*fs0pera;
fs0b=fs0par.*fs0perb;
fs0=fs0a+fs0b;
fs0_b=1/(sqrt(pi^3)*vtzs(s)).*exp(-((vparp-vdsz(s))/vtzs(s)).^2).*(...
    rsab(1,s)./(Asab(1,s).*vtpsab(1,s).^2).*exp(...
     -((vperp-vdsr(s))/vtpsab(1,s)).^2)+...
    rsab(2,s)./(Asab(2,s).*vtpsab(2,s).^2).*exp(...
     -((vperp-vdsr(s))/vtpsab(2,s)).^2));
dfs0dvparp=-2*(vparp-vdsz(s)).*(fs0a./vtzs(s).^2+fs0b./vtzs(s).^2);
dfs0dvperp=-2*(vperp-vdsr(s)).*(fs0a./vtpsab(1,s).^2+fs0b./vtpsab(2,s).^2);

Us1=(w-kz*vparp).*dfs0dvperp+kz*vperp.*dfs0dvparp;
Us2=kx*vdsy(s)*dfs0dvperp;
Us3=kx*vparp*dfs0dvperp-kx*vperp*dfs0dvparp;
Us4=(w-kz*vparp-kx*vdsx(s)).*dfs0dvperp+kz*vperp.*dfs0dvparp;
Us5=kz*vdsx(s)*dfs0dvparp;
Us6=kz*vdsy(s)*dfs0dvparp;
Us7=(w-kx*vdsx(s)).*dfs0dvparp;

dfs(jvx,jvy,jvz)=0.0;
xs=-(w-kz*vparp-kx*vdsx(s)+1i*nus(s))./wcs(s); ys=kx*vperp./wcs(s)+1e-10;
N=2;
for jm=-N:N
    for jn=-N:N
        dfs(jvx,jvy,jvz)=dfs(jvx,jvy,jvz)+...
            1i*qs(s)./(ms(s)*wcs(s)*w).*(besselj(jm,ys)*exp(...
            -1i*(jn-jm)*phip)/(xs+jn))*(...
            (Us1*Ex+Us2*Ey+Us3*Ez)*jn/ys*besselj(jn,ys)+...
            Us4*Ey*1i*(-0.5)*(besselj(jn+1,ys)-besselj(jn-1,ys))+...
            (Us5*Ex+Us6*Ey+Us7*Ez)*besselj(jn,ys));
    end
end

        end
    end
end

%%
dvx=vxx(2,1,1)-vxx(1,1,1);dvy=vyy(1,2,1)-vyy(1,1,1);dvz=vzz(1,1,2)-vzz(1,1,1);
nsf=ns0(s)*sum(sum(sum(dfs)))*dvx*dvy*dvz;
Jsxf=qs(s)*ns0(s)*sum(sum(sum(vxx.*dfs)))*dvx*dvy*dvz;
Jsyf=qs(s)*ns0(s)*sum(sum(sum(vyy.*dfs)))*dvx*dvy*dvz;
Jszf=qs(s)*ns0(s)*sum(sum(sum(vzz.*dfs)))*dvx*dvy*dvz;
%%
[Jsy,Jsy,Jsz]./[Jsyf,Jsyf,Jszf]
ns/nsf
%% 
close all;
h=figure('unit','normalized','Position',[0.01 0.45 0.6 0.6],...
  'DefaultAxesFontSize',15);

subplot(131);
jvz=floor(nvz/2)+1;
contourf(squeeze(vxx(:,:,jvz))/vtps(s),squeeze(vyy(:,:,jvz))/vtps(s),...
    real(squeeze(dfs(:,:,jvz))),50);
colorbar('southoutside');
xlabel(['v_x/v_{tps}, with v_z/v_{tzs}=',num2str(vzz(1,1,jvz)/vtzs(s))]); 
ylabel('v_y/v_{tps}');axis tight;


subplot(132);
jvy=floor(nvy/2)+0;
contourf(squeeze(vxx(:,jvy,:))/vtps(s),squeeze(vzz(:,jvy,:))/vtzs(s),...
    real(squeeze(dfs(:,jvy,:))),50);
colorbar('southoutside');
xlabel(['v_x/v_{tps}, with v_y/v_{tps}=',num2str(vyy(1,jvy,1)/vtps(s))]); 
ylabel('v_z/v_{tzs}');axis tight;
title(['\delta{}f_{s=',num2str(s),'}(v_x,v_y,v_z), k/k_n=',num2str(k/kn),...
    ', \omega/\omega_n=',num2str(w/wn,3),', n_{vx}=',num2str(nvx),...
    ', n_{vy}=',num2str(nvy),', n_{vz}=',num2str(nvz)]);

subplot(133);
jvx=floor(nvx/2)+2;
contourf(squeeze(vyy(jvx,:,:))/vtps(s),squeeze(vzz(jvx,:,:))/vtzs(s),...
    real(squeeze(dfs(jvx,:,:))),50);
colorbar('southoutside');
xlabel(['v_y/v_{tps}, with v_x/v_{tps}=',num2str(vxx(jvx,1,1)/vtps(s))]); 
ylabel('v_z/v_{tzs}');axis tight;

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',[savepath,'fig_bo_',figstr,'s=',num2str(s),...
    'id=',num2str(id),',nvx=',num2str(nvx),...
    ',nvy=',num2str(nvy),',nvz=',num2str(nvz),'_plt_df.png']);

