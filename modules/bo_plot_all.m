% 18-10-05 08:00 Hua-sheng XIE, huashengxie@gmail.com, ENN-FTRC, China
% This file plot all the solution of BO's results, which can be used
% to select the initial data for separating different dispersion surface
% 20-12-23 18:58 update to include iem=3 fluid version
close all;

h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.4],...
  'DefaultAxesFontSize',15);

% normalized omega and k
wwn=ww/wn;
kkn=kk/kn;
kxxn=kxx/kn;
kzzn=kzz/kn;

if(ipa==ipb) % 1D plot

  %kkk=reshape(repmat(kkn,1,nw0),1,[]);
  papp=reshape(repmat(ppa,1,nw0),1,[]);
  www=reshape(wwn,1,[]);
  npp=length(papp);
  
  % To plot only the unstable mode
%   ind=find(imag(www)<1e-1); www(ind)=NaN+1i*NaN;

  subplot(121);
  if(iloga==0)
%     plot(papp,real(www),'g+','LineWidth',2); hold on;
    plot(papp,real(www),'.'); hold on;
    xlim([min(pa),max(pa)]);
  else
    semilogx(10.^papp,real(www),'g+','LineWidth',2); hold on;
    xlim([min(10.^pa),max(10.^pa)]);
  end  
  xlabel([strpa,', runtime=',num2str(runtime),'s']); 
  ylabel(['\omega_r/\omega_{n}, npa=',num2str(npa),',npb=',num2str(npb)]);
  title(['(a) \beta_{||}=',num2str(betasz,3),...
      ', \beta_\perp=',num2str(betasp,3)]); 
  box on;
  if(jsetplot==1 && exist('yrmin','var') && exist('yrmax','var'))
    ylim([yrmin,yrmax]);
  end
  
  subplot(122);
  if(iloga==0)
%     plot(papp,imag(www),'g+','LineWidth',2); hold on;
    plot(papp,imag(www),'.'); hold on;
    xlim([min(pa),max(pa)]);
  else
    semilogx(10.^papp,imag(www),'g+','LineWidth',2); hold on;
    xlim([min(10.^pa),max(10.^pa)]);
  end
  if(iem~=3)
    xlabel([strpa,', (S=',num2str(S),',N=',num2str(N),',J=',num2str(J),')',...
      ', iem=',num2str(iem)]); 
  else
    xlabel([strpa,', (S=',num2str(S),')',', iem=',num2str(iem)]); 
  end
  ylabel('\omega_i/\omega_{n}');
  title(['(b) v_A/c=',num2str(vA/sqrt(c2),2),', ',strpb,'=',...
    num2str(par(ipbtmp))]);
   box on; 
  
  if(jsetplot==1 && exist('yimin','var') && exist('yimax','var'))
    ylim([yimin,yimax]);
  end
  %ylim([-5,15]);%ylim([-1.0e-3,1.0e-2]);
  
else % 2D plot
  for jp=1:25 % change here to plot more surfaces
    subplot(121);
    wwjp=squeeze(wwn(:,:,jp));
    surf(ppa,ppb,real(wwjp)); hold on; box on;
    xlabel([strpa,',ilogx=',num2str(iloga)]);
    ylabel([strpb,',ilogy=',num2str(ilogb)]);
    zlabel(['\omega_r/\omega_{n},npa=',num2str(npa),',npb=',num2str(npb)]);
    title(['(a) \beta_{||}=',num2str(betasz,3),...
      ', \beta_\perp=',num2str(betasp,3)]);
    subplot(122);
    surf(ppa,ppb,imag(wwjp)); hold on; box on;
    xlabel(strpa); ylabel(strpb);
    if(iem~=3)
      zlabel(['\omega_i/\omega_{n},S=',num2str(S),',N=',num2str(N),',J=',num2str(J)]);
    else
      zlabel(['\omega_i/\omega_{n},S=',num2str(S)]);
    end
    title(['(b) runtime=',num2str(runtime),'s']);
    %%
    % zoom in the figure to find dispersion surface data for plot_select.m
    %zlim([-0.5e0,0.1]);
    %%
  end
  
end

% 20-12-23 20:35 output all solutions, to limit the output size, we cut it
% at min(nw,4*S+6)
if(iem~=3)
    fileID=fopen([savepath,'bo_kw_all_iem=',num2str(iem),'_S=',num2str(S),...
        '_J=',num2str(J),'_N=',num2str(N),'.txt'],'w');
else
    fileID=fopen([savepath,'bo_kw_all_iem=',num2str(iem),'_S=',num2str(S),...
        '_fluid.txt'],'w');
end

nkwall=min(nw,4*S+6);
kwall=zeros(npa*npb,nkwall);
for jpa = 1:npa
  for jpb = 1:npb
    jx=(jpb-1)*npa+jpa;
    kwall(jx,1)=kzzn(jpa,jpb);
    kwall(jx,2)=kxxn(jpa,jpb);
    kwall(jx,3)=kkn(jpa,jpb);
    kwall(jx,4)=tt(jpa,jpb);
    for jw=1:nkwall
      kwall(jx,4+2*jw-1)=real(wwn(jpa,jpb,jw));
      kwall(jx,4+2*jw-0)=imag(wwn(jpa,jpb,jw));
    end
  end
end
fprintf(fileID,['Below are the first nkwall=min(nw,4*S+6)=',...
    num2str(nkwall),' solutions from BO, sorted by imag part, at ',...
    datestr(datetime),'\n']);

fprintf(fileID,'%-17.16s','kz[kn]');
fprintf(fileID,'%-17.16s','kx[kn]');
fprintf(fileID,'%-17.16s','k[kn]');
fprintf(fileID,'%-17.16s','theta[deg]');
for jw=1:nkwall
  fprintf(fileID,'%-17.16s',['Re(omega',num2str(jw),'[wn])']);
  fprintf(fileID,'%-17.16s',['Im(omega',num2str(jw),'[wn])']);
end
fprintf(fileID,'\n');
fprintf(fileID,[repmat('%-17.6e',1,2*nkwall+4),'\n'],kwall.');
fclose(fileID);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

if(iem~=3)
figstr=['S=',num2str(S),'_theta=',num2str(theta),'_J=',num2str(J),'_N=',num2str(N),'_iem=',...
    num2str(iem),'_npa=',num2str(npa),'_npb=',num2str(npb),...
    '_B0=',num2str(B0),'_vdsz01=',num2str(vdsz0(1))];
else
figstr=['S=',num2str(S),'_theta=',num2str(theta),'_iem=',...
    num2str(iem),'_npa=',num2str(npa),'_npb=',num2str(npb),...
    '_B0=',num2str(B0),'_vdsz01=',num2str(vdsz0(1))];
end
print(gcf,'-dpng',[savepath,'fig_bo_',figstr,'_all.png']);
% savefig([savepath,'fig_bo_',figstr,'_all.fig']);
