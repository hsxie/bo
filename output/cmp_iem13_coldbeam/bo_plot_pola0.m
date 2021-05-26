% 20-12-24 21:29 Hua-sheng XIE, huashengxie@gmail.com, ENN-FTRC, China
% plot polarizations of select variables

for jpl=1:npl:0
    
    subplot(121); hold on;
    plot(pas,real(wws(:,1,jpl)),'--','Color',...
        pltc(jpl,:),'linewidth',2);
end

%%

% set name for each polarizations
strpola=cell(npolaf+3+npolas*S);
strpola(1:21)={ 'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz', 'U_E', 'U_B', 'U_E/(U_E+U_B)', ...
  '|Ex|^2/|E|^2', '|Ey|^2/|E|^2', '|Ez|^2/|E|^2', ...
  '|kdotE|^2/(k^2|E|^2)', '|kdotE|^2/(|kdotE|^2+|kcrossE|^2)', ...
  '|Bx|^2/|B|^2', '|By|^2/|B|^2', '|Bz|^2/|B|^2', 'Elliptic','\theta_B',...
  '|By/Bx|','angle(By/Bx)'};
strpola((npolaf-3):(npolaf+3))={'kx','kz','k','\omega/\omega_n','Jx','Jy','Jz'};
for s=1:S
    snum=['(s=',num2str(s),')'];
    strpola((npolaf+3+npolas*(s-1)+1):(npolaf+3+npolas*(s-1)+7))={['Jsx',snum], ['Jsy',snum],...
  ['Jsz',snum], ['ns',snum], ['vsx',snum], ['vsy',snum], ['vsz',snum]};
end

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

% icase=2;
jpll=[1,1]*1;
for icase=1:2;
    savepath='./';
if(icase==1)
    close all;
    load([savepath,'out_bo_S=4_B0=2e-06_fluid.mat']);
    h=figure('unit','normalized','Position',[0.01 0.05 0.8 0.8],...
        'DefaultAxesFontSize',14);
    linestr={'-','--'};
else
    load([savepath,'out_bo_S=4_J=4_N=1_B0=2e-06.mat']);
    linestr={':','-.'};
end

jpl=jpll(icase); % select which solution to plot
pkzz=Pola(:,:,jpl,npolaf-2)/kn;
% jjpola=[npolaf,1,2,3,npolaf+1,npolaf+2,npolaf+3,(npolaf+3+npolas*(s-1)+1):(npolaf+3+npolas*s)];

% select which variables to plot
jjpola=[npolaf,(npolaf+3+npolas*(1-1)+4):(npolaf+3+npolas*(1-1)+7),...
    (npolaf+3+npolas*(2-1)+4):(npolaf+3+npolas*(2-1)+7),...
    (npolaf+3+npolas*(3-1)+4):(npolaf+3+npolas*(3-1)+7),...
    (npolaf+3+npolas*(4-1)+4):(npolaf+3+npolas*(4-1)+7),1:3];

for jj=1:min(length(jjpola),20)
    
    subplot(4,5,jj);
    hold on;
    if(jjpola(jj)==npolaf)
      plot(pkzz,real(Pola(:,:,jpl,jjpola(jj))./wn),linestr{1},...
        pkzz,imag(Pola(:,:,jpl,jjpola(jj))./wn),linestr{2},'linewidth',2);
    if(icase==2)
        legend('Re fluid','Im fluid','Re kinetic','Im kinetic',...
            'location','best'); legend('boxoff');
    end
%         title(['vdsz/c=',num2str(vdsz/sqrt(c2),3)]);
    else
      plot(pkzz,real(Pola(:,:,jpl,jjpola(jj))),linestr{1},...
        pkzz,imag(Pola(:,:,jpl,jjpola(jj))),linestr{2},'linewidth',2);
    end
    grid on;
    xlabel('k/k_n');
    ylabel(strpola{jjpola(jj)});
%     xlim([0.1,5]);
end
end

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,'-dpdf',['bo_plot_pola_iem',num2str(iem),'.pdf']);
print(gcf,'-dpng',[savepath,'bo_plot_pola_iem=',num2str(iem),...
    '_jpl=',num2str(jpl),'_kmax=',num2str(k/kn),'_cmp.png']);



