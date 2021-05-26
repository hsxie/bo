% 18-10-05 16:42 Hua-sheng XIE, huashengxie@gmail.com, ENN-FTRC, China
% Try a stratege to select given dispersion surface. Steps:
% 1. select a initial point in the solutions, e.g., (k,w_j)
% 2. determined the corresponding w_j for other k, use interp or least
% squre, or nearest neighbor
% 3. store the (k,w_j) as one dispersion surface
% 4. repeat the above step for another surface
% better for smaller dk
% 18:41 Test 1D nearest neighbor can work, but not well for omega past zero
% 19:33 try interp1, works very well!
% 18-10-06 00:48 2D scan not perfect yet
% 20-01-11 1:45 fixed a bug for jselc=1
% 20-12-23 19:04 update to include iem=3 fluid version

% % Search the most close dispersion surfaces to these data.
% % Initial data for find the corresponding dispersion surfaces.
% % Please use bo_plot_all.m to visualize all the solutions, and then
% % modify here the initial point of which mode(s) you want plot/store.
% if(ipa==ipb) % for 1D scan case
% % wpdat(:,1) is pa; wpdat(:,2) is arbitrary; wpdat(:,3) is Re or Im(omega)
% wpdat=[
%     %80,0,-0.061i; % theta=75
%     83.5,0,5.7e5;
%   68.5,0,-3.23e2*1i;
%   56,5,-2.45e2*1i];
% % wpdat=[0.65,0,-0.45i;
% %   0.4,0,-0.37i;
% %   0.8,0,-1.75i]; % modify here for other dispersion surfaces
% else % for 2D scan case
%     % wpdat(:,1) is pa; wpdat(:,2) is pb; wpdat(:,3) is Re or Im(omega)
% wpdat=[70,0,-152i;
%   90.5,85,-261.5i;
%   94.5,85,-433i];
% end

% run ../input/bo_wpdat;

if(jselc==1) % alway plot the most unstable dispersion surface
  % find the most unstable/least damped solution
  maxgam=max(max(max(imag(wwn))));
%   [indww,indpa,indpb]=find(imag(wwn)==maxgam);% wrong! 2020-01-11 11:44
  [indpa,indpb,indww]=find(imag(wwn)==maxgam);
  wpdat=[pa(indpa),pb(indpb),maxgam*1i;
    wpdat];
end

% [npl,nplb]=size(wpdat);
npl=size(wpdat,1); % 18-10-27 14:27
wws=zeros(npa,npb,npl);

%% This module is to search the same dispersion surface(s) of 'wpdat'
for jpl=1:npl
% for jpl=2:2

  datstart=wpdat(jpl,:);

  if(ipa==ipb) % 1D plot
      
      % determine the start point solution
      indpa=find(abs(pa-datstart(1))==min(abs(pa-datstart(1)))); 
      indpa=indpa(1);
      pastart=pa(indpa); % pa=k*kn, etc

      if(real(datstart(3))==0) % search by imag part if real part equal 0
        indww=find(abs(imag(wwn(indpa,1,:))-imag(datstart(3)))==...
          min(abs(imag(wwn(indpa,1,:))-imag(datstart(3)))));
      else % search by real part
        indww=find(abs(real(wwn(indpa,1,:))-real(datstart(3)))==...
          min(abs(real(wwn(indpa,1,:))-real(datstart(3)))));
      end
      indww=indww(1);
      wwstart=wwn(indpa,1,indww); % wwn=ww/wn;

      % determine the solution for next parameter a
      jj=1:npa; pas=0.*jj; ws=0.*jj;
      jjpa=indpa; jjww=indww; pas(indpa)=pa(indpa); ws(indpa)=wwstart;
      while(length(jjpa)<npa)
        jjpa0=jjpa;
        if(max(jjpa)<npa) % seach the right direction of pa
          indm=max(jjpa);
          indpa=indm+1;
          jjpa=[jjpa,indpa];
        else % seach the left direction of pa
          indm=min(jjpa);
          indpa=indm-1;
          jjpa=[indpa,jjpa];
        end
        if(length(jjpa)<=2)
          wpred=ws(indm); % nearest neighbor
        else % interp1 to give the best guess
          wpred=interp1(pas(jjpa0),ws(jjpa0),pa(indpa),...
              'pchip','extrap'); % try 'cubic', 'spline' or 'pchip'
        end
        indww=find(abs(wwn(indpa,1,:)-wpred)==...
          min(abs(wwn(indpa,1,:)-wpred))); % find next solution
        indww=indww(1);
        jjww=[jjww,indww];
        pas(indpa)=pa(indpa);
        ws(indpa)=wwn(indpa,1,indww);
      end
      wws(1:npa,1,jpl)=ws(1:npa); % store the solutions

  else % 2D plot
      
      % determine the start point solution
      indpa=find(abs(pa-datstart(1))==min(abs(pa-datstart(1)))); 
      indpa=indpa(1);
      pastart=pa(indpa); % pa start point for search
      indpb=find(abs(pb-datstart(2))==min(abs(pb-datstart(2)))); 
      indpb=indpb(1);
      pbstart=pb(indpb); % pb start point for search

      if(real(datstart(3))==0) % search by imag part if real part equal 0
        indww=find(abs(imag(wwn(indpa,indpb,:))-imag(datstart(3)))==...
          min(abs(imag(wwn(indpa,indpb,:))-imag(datstart(3)))));
      else % search by real part
        indww=find(abs(real(wwn(indpa,indpb,:))-real(datstart(3)))==...
          min(abs(real(wwn(indpa,indpb,:))-real(datstart(3)))));
      end
      indww=indww(1);
      wwstart=wwn(indpa,indpb,indww); % wwn=ww/wn;
      
      % determine the solution for next parameter a, b
      ws=zeros(npa,npb); ws(indpa,indpb)=wwstart;
      jjpb=indpb; wpredb=wwstart; indpa0=indpa; indpb0=indpb;
      while(length(jjpb)<=npb)
          
          indpa=indpa0; jjpa=indpa; jjww=indww; jstart=1;
          while(length(jjpa)<npa) % for next parameter a
            jjpa0=jjpa;
               
            if(jstart==1)
               wpred=wpredb;
               jstart=0;
            else
                if(max(jjpa)<npa) % seach the right direction of pa
                  indm=max(jjpa);
                  indpa=indm+1;
                  jjpa=[jjpa,indpa];
                else % seach the left direction of pa
                  indm=min(jjpa);
                  indpa=indm-1;
                  jjpa=[indpa,jjpa];
                end
                if(length(jjpa)<=2)
                  wpred=ws(indm,indpb); % nearest neighbor
                else % interp1 to give the best guess.
                  wpred=interp1(pa(jjpa0),squeeze(ws(jjpa0,indpb)),...
                      pa(indpa),'pchip','extrap'); % 'pchip', 'cubic', 'spline'
                end
            end
            indww=find(abs(wwn(indpa,indpb,:)-wpred)==...
              min(abs(wwn(indpa,indpb,:)-wpred))); % find next solution
            indww=indww(1);
            jjww=[jjww,indww];
            ws(indpa,indpb)=wwn(indpa,indpb,indww);
          end
          
          % % for next parameter b
          jjpb0=jjpb;
          if(max(jjpb)<npb) % seach the right direction of pb
              indbm=max(jjpb);
              indpb=indbm+1;
              jjpb=[jjpb,indpb];
          else % seach the left direction of pb
              indbm=min(jjpb);
              indpb=indbm-1;
              jjpb=[indpb,jjpb];
          end
          
          if(length(jjpb)<npb)
              if(length(jjpb)<=2)
                wpredb=ws(indpa0,indbm); % nearest neighbor
              else % interp1 to give the best guess. To do: how about interp2?
                wpredb=interp1(pb(jjpb0),squeeze(ws(indpa0,jjpb0)),...
                    pb(indpb),'pchip','extrap'); % 'pchip', 'cubic', 'spline'
              end
          end
          
      end
      wws(1:npa,1:npb,jpl)=ws(1:npa,1:npb); % store the solutions
  end
  
end

% wws0=wws; ind=find(imag(wws)<0); wws(ind)=NaN+1i*NaN; % 18-10-29 19:52
%% for plot
close all;
h=figure('unit','normalized','Position',[0.01 0.45 0.6 0.4],...
  'DefaultAxesFontSize',15);
pltc=[0.0  1.0  0.0 % different colors
  1.0  0.0  0.0
  0.2  0.2  1.0
  0.8 0.8 0.0
  1.0  0.6  0.0
  0.9  0.0  0.9
  0.0  0.8  0.8
  0.0  0.0  0.0
  0.6  0.0  0.0
  0.4  0.7  0.4 
  0.0  0.0  0.5 
  0.6  0.0  0.6 
  0.0  0.5  1.0
  ];
% for jpl=1:2
for jpl=1:npl
    
  if(ipa==ipb) % 1D plot
%       subplot(221); hold on;
      subplot(121); hold on;
      if(iloga==0)
        plot(pas,real(wws(:,1,jpl)),'--','Color',...
            pltc(jpl,:),'linewidth',2);
%         xlim([0*min(pas),max(pas)]);
      else
        semilogx(10.^pas,real(wws(:,1,jpl)),'--',...
            'Color',pltc(jpl,:),'linewidth',2);
        %xlim([10^min(pas),10^max(pas)]); 
      end
      box on; grid on; %ylim([-100,600]);
      xlabel([strpa,', npa=',num2str(npa)]);
      title(['(a) \beta_{||}=',num2str(betasz,3),...
        ', \beta_\perp=',num2str(betasp,3)]);
      ylabel(['\omega_r/\omega_{n}, \alpha=',num2str(alphas,3)]);
      if(jsetplot==1 && exist('yrmin','var') && exist('yrmax','var'))
        ylim([yrmin,yrmax]);
      end
      
      subplot(122); hold on;
      if(iloga==0)
        plot(pas,imag(wws(:,1,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
%         xlim([0*min(pas),max(pas)]);
      else  
        semilogx(10.^pas,imag(wws(:,1,jpl)),'--',...
            'Color',pltc(jpl,:),'linewidth',2);
        %xlim([10^min(pas),10^max(pas)]);
      end
      box on; grid on; %ylim([-35,15]);
      xlabel([strpa,', iem=',num2str(iem)]);
      ylabel(['\omega_i/\omega_{n}, \Delta=',num2str(Deltas,3)]);
      if(iem~=3)
        title(['(b) v_A/c=',num2str(vA/sqrt(c2),2),', ',strpb,'=',...
          num2str(par(ipbtmp)),', (S=',num2str(S),',N=',num2str(N),...
          ',J=',num2str(J),')']);
      else
        title(['(b) v_A/c=',num2str(vA/sqrt(c2),2),', ',strpb,'=',...
          num2str(par(ipbtmp)),', (S=',num2str(S),')']);
      end
  
      if(jsetplot==1 && exist('yimin','var') && exist('yimax','var'))
        ylim([yimin,yimax]);
      end
    
%       subplot(2,2,3:4); hold on;
%       if(iloga==0) % 18-11-20 23:03
%         plot(pas,imag(wws(:,1,jpl))./real(wws(:,1,jpl)),'--',...
%             'Color',pltc(jpl,:),'linewidth',2);
%         xlim([min(pas),max(pas)]);
%       end
%       xlabel(strpa);
%       ylabel('\omega_i/\omega_{r}'); box on; grid on;
  else
    subplot(121);
    wwjp=squeeze(wws(:,:,jpl));
    surf(ppa,ppb,real(wwjp)); hold on; box on; set(gca,'BoxStyle','full');
    xlabel([strpa,',ilogx=',num2str(iloga)]);
    ylabel([strpb,',ilogy=',num2str(ilogb)]);axis tight;
    zlabel(['\omega_r/\omega_{n},npa=',num2str(npa),',npb=',num2str(npb)]);
    title(['(a) \beta_{||}=',num2str(betasz,3),...
      ', \beta_\perp=',num2str(betasp,3)]);
    subplot(122);
    surf(ppa,ppb,imag(wwjp)); hold on; box on;set(gca,'BoxStyle','full');
    xlabel(strpa); ylabel(strpb); axis tight;
    if(iem~=3)
      zlabel(['\omega_i/\omega_{n},S=',num2str(S),',N=',num2str(N),...
        ',J=',num2str(J)]);
    else
      zlabel(['\omega_i/\omega_{n},S=',num2str(S)]);
    end
    % zlim([-0.1,0.1]);
    title(['(b) runtime=',num2str(runtime),'s']);

  end
  
end

if(ipa==ipb) % 18-10-31 21:08 to check
  if(iloga==0)
%     kwout=[pas.',squeeze(wws(:,1,:))];
    kwout=[pas.',real(squeeze(wws(:,1,:))),imag(squeeze(wws(:,1,:)))];
    ind=1;
    for jl=1:npl
      ind=[ind,jl+1,jl+1+npl];
    end
    kwout=kwout(:,ind);
  end
  if(iem~=3)
    fileID=fopen([savepath,'bo_kwout_iem=',num2str(iem),'_S=',num2str(S),...
      '_J=',num2str(J),'_N=',num2str(N),'.txt'],'w');
  else
    fileID=fopen([savepath,'bo_kwout_iem=',num2str(iem),'_S=',num2str(S),...
      '_fluid.txt'],'w');
  end
  fprintf(fileID,'k, Re(omega1), Im(omega1), Re(omega2), Im(omega2), ... [by BO]\n');
  fprintf(fileID,[repmat('%+9.5e \t',1,2*npl+1),'\n'],kwout.');
%   fmat=['%9.5e  '];
%   for j=1:size(kwout,1)
%     fprintf(fileID,'%9.5e      ',kwout(j,:));
%     fprintf(fileID,'\n');
%   end
  fclose(fileID);
end

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',[savepath,'fig_bo_',figstr,'_select.png']);
% savefig([savepath,'fig_bo_',figstr,'_select.fig']);
