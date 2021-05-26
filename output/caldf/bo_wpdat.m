% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)

% wpdat=[0.44,0,0.05771i; % for loss cone mirror
%     0.22,0,0.05041i;
%     0.50,0,-0.04363i;
%     0.47,0,-0.1578i;
%     0.41,0,-0.2965i;
%   ];

% wpdat=[0.23,0,0.199i;
%   0.16,0,0.056i;
%   0.55,0,-0.083i;
%   0.35,0,-0.109i;
%   ];

% wpdat=[%0.34,0,0.342;
%   1.00,0,0.6488;
%   0.55,0,0.3457;
%   1.00,0,0.328i;
%   ];

wpdat=[ 0.48,0,0.40;
%     0.48,0,0.01669i;
  ];

jselc=0; % jselc=1, alway plot the most unstable dispersion surface
% wpdat=[];

%%
run ../modules/bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
