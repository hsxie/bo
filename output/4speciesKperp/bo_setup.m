% 18-10-05 08:00 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Setup for run BO (PDRK&PDRF)
% Note: set your normalization kn and omega_n in initialize.m
% 20-12-22 21:26 update to including option iem=3 for electromagnetic 
% multi-fluid version BO-F/PDRF

savepath='../output/4speciesKperp/'; % choose where to save the figures and outputs
iem=1; % =1, electromagnetic run; =0, electrostatic run; =2, Darwin model; =3 fluid model

if(iem~=3) % only needed for kinetic version
  N=4; % Number harmonics, enlarge N to make sure the results are convergent
  J=4; % J-pole, usually J=8 is sufficient; other choice: 4, 12, 16, 24
  
  % To set some species with other N, to speed up or accurate, i.e.,
  % different species can use different N: setNs(:,1) is species index;
  % setNs(:,2) is corresponding N which will be used in the computation.
  setNs=[];
  % setNs=[1, floor(N/2);
  %        3, 2*N];
else % iem=3, fluid version, 20-12-22 22:41
  % set the default polytrope exponents
  % (gammazs,gammaps)=(1.0,1.0) isothermal; =(3.0,2.0) CGL model
  gammas0=[1.0,1.0]; % default [para,perp] polytrope exponents gamma_{para,perp}
  % to set some species with other gammaps or gammazs
  % set..=[1,3.0,2.0] means set the '#1' species (gammazs,gammaps)=(3.0,2.0)
  % setgammas=[1,3.0,2.0]; % setgammas=[1,1.0,2.0;3,1.0,1.0];
  setgammas=[];
end

% sp=0, eig() for all solution; =1 sparse eigs(), 1 wg0 multi nw0;
% =2, multi wg0, 1 nw0;
sp=0;

% inital guess for sp=1&2, normalized by omega_n; not required if sp=0
if(sp<=1) % sp=1, search multi solutions use one wg0
 wg0=1.914+1i*4.455;
else % sp=2, search solution of different branches use different wg0
 wg0=[
     1.914+4.455i; % k=0.5
     ];
end

B0=100E-9; % magnetic field (Tesla)

% initialize the parameters k, theta, kz, kx
par=zeros(6,1);
par(1)=0.05; % k, =sqrt(kx^2+kz^2), normalized by kn
par(2)=90.0; % theta, the angle between k and B0, normalized by *pi/180
par(3)=cos(par(2)*pi/180)*par(1); % kz, normalized by kn
par(4)=sin(par(2)*pi/180)*par(1); % kx, normalized by kn
% par(4)=0.5;

% Choose which parameter(s) to scan, if ipa==ipb, do 1D scan; otherwise, 
% do 2D scan. % 1: k; 2: theta; 3: kz; 4: kx; 5: others.
% You can scan other parameters by modify line xxx in 'bo_si_kernel.m',
% e.g., to scan vs0(1), B0, betasz, etc
% Typical cases of (ipa,ipb): 
%   1. (1,1) scan k, fixed theta
%   2. (2,2) scan theta, fixed k
%   3. (1,2) scan 2D (k, theta)
%   4. (3,3) scan kz, fixed kx
%   5. (4,4) scan kx, fixed kz
%   6. (3,4) scan 2D (kz,kx)
%   7. (..,..) others, please modify 'bo_si_kernel.m'
ipa=3;
ipb=3;

iloga=0; % ilog=0, linear scale; =1, log scale, i.e., 10^(p1:dp:p2)
ilogb=0;
pa1=0.02; pa2=0.2; dpa=0.01; % 1st parameter a, depends on ipa
% pa1=0.2; pa2=0.23; dpa=0.01; % 1st parameter a, depends on ipa
% pa2=-0.01; pa1=-2.5; dpa=0.01; % 1st parameter a, depends on ipa
% pa1=1; pa2=89; dpa=2.5; % 1st parameter a, depends on ipa
pb1=0; pb2=89; dpb=2.5; % 2nd parameter b, depends on ipb

% wether calculate polarizations (dEx,dEy,dEz,dBx,dBy,dBz) for select omega
iout=2; % =1, only (omega,gamma); =2, also calculate (E,B); =0, skip output 

%%

% % Default SI unit
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
mu0=1/(c2*epsilon0);
kB=1.38064852e-23;
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kg

% %You can also use normalized unit, such as:
% c2=100^2; % speed of light c^2
% epsilon0=1.0;
% % mu0=1.0;
% mu0=1/(c2*epsilon0);
% kB=1.0;
% qe=1;
% mp=1836;
% me=1;

munit=mp; % set the default normalization for mass
% munit=me;

%%
% for plot ylim([ymin,ymax]), 18-12-20 19:44
jsetplot=1;
yrmin=-1.0;
yrmax=2.0;
yimin=-0.6;
yimax=0.2;
