BO (‘波’, i.e., 'wave' in Chinese) is a state-of-art tool for plasma wave 
and instability analysis. It includes two codes currently, the BO-F (PDRF, 
A general dispersion relation solver for multi-fluid plasma） and BO-K 
(PDRK, A General Kinetic Dispersion Relation Solver for both Magnetized and 
Unmagnetized Plasma).

Copyright (C) 2014-2016
 PDRF, Huasheng Xie <huashengxie@gmail.com>, 2014, IFTS-ZJU
 PDRK, Huasheng Xie <huashengxie@gmail.com> and Yong Xiao, 2014-2016, IFTS-ZJU

Copyright (C) 2018-2019
 BO, Huasheng Xie <huashengxie@gmail.com/xiehuasheng@enn.cn>, 2018-2019, CCF-ENN

This file is part of BO version 190305.

BO is free software under the BSD 3-clause License.

BO version 190305 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

----

BO-K (PDRK) is a general and powerful Kinetic Dispersion Relation solver for 
magnetized and unmagnetized Plasma, which is the first kinetic dispersion 
relation solver that can give all the important solutions at one time without 
requiring initial guess for root finding.

% 2019-01-26 23:37, Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China.
% Ackn.: Richard DENTON (Dartmouth), Xin TAO (USTC), Jin-song ZHAO (PMO),
% Zhong-wei YANG (NSSC), Chao-jie ZHANG (UCLA), Can Huang (USTC), Wen-ya 
% LI (NSSC), Liang WANG (Princeton), Kyungguk MIN (Korean), Yang Li (ENN),
% etc ...
%
% This is the most comprehensive version of BO (means 'wave' in Chinese, 
% old name PDRK/PASS-K) kinetic plasma dispersion relation solver, which
% assumes an extend non-relativitic Maxwellian distribution function. It
% includes bi-Maxwellian, parallel drift vdz, penpendicular drift across
% magnetic field (vdx, vdy), ring beam vdr, loss cone, and Krook collision.
% Assumed B0=(0,0,B0), k=(kx,0,kz).
%
% What it 'can' handle:
%  1. The species can be either magnetized (B0\neq0) or unmagnetized (B0=0),
% both EM3D and ES3D versions, and also Darwin model version;
%  2. All the EM3D parallel drift bi-Maxwellian loss cone case as in WHAMP 
% [Ronnmark1982], and also corresponding ES3D version;
%  3. The ES3D and EM3D ring beam case in [Umeda2007] and [Umeda2012].
%  4. The EM3D drift across magnetic field case in [Umeda2018], with
% extensions to including all (vdx, vdy, vdz, vdr, loss-cone) and ES3D, and
% support unmagnetized species such as in [Muschietti2017].
%  5. With Krook collision.
%  6. Support output polarizations (E,B) for electromagnetic version.
%  7. Without singularity for theta=0 & pi/2, and also support kz<0.
%
% What it 'can not' handle at this stage:
%  1. Ring beam of unmagnetized species (usually, it is not important);
%  2. Relativistic effect;
%  3. kappa-distribution, non-uniform magnetic field, etc;
%  4. Gyrokinetic model (easy, to do).
%
% In numerical aspect, for given (kx,kz), all important kinetic solutions
% of complex omega=omega_r+1i*omega_i can be obtained at one time.
% Thus:
%  1. The user can choose to solve all the solutions use eig() with sp=0,
% or several solutions around initial guess wg using eigs() with sp=1, or
% different branches from different wg with sp=2;
%  2. With module to separate dispersion surfaces;
%  3. Parallelized computation to scan k using 'parfor' (not included in
% the default version, the user need slightly modify the code) due to some
% requirement of 'parfor'.
%
% The BO (PDRK) algorithm has two steps:
%  (1) J-pole approximation for Z(zeta)~=Z_J(zeta)=sum(bj/(zeta-cj)), which
% is accurate to 10^-6 with J=8, except the less interesting strong damped
% modes;
%  (2) Transform to a equivalent matrix eigenvalue problem lambda*X=MX,
% which can be solved easily use standard library such as eig() in Matlab,
% with the solution omega=lambda, and (E,B) in X.
% 
% =======================================================================
% If you use this code for your works, please cite:
%  [Xie2016] Huasheng Xie and Yong Xiao, PDRK: A General Kinetic Dispersion
% Relation Solver for Magnetized Plasma, Plasma Science and Technology, 18,
% 2, 97 (2016). DOI: 10.1088/1009-0630/18/2/01. Update/Bugs fixed at 
% http://hsxie.me/codes/pdrk/ or https://github.com/hsxie/pdrk.
%  [Xie2019] Hua-sheng XIE, A Unified Numerically Solvable Framework for
% Complicated Kinetic Plasma Dispersion Relations, arXiv:1901.06902, 2019.
% https://arxiv.org/abs/1901.06902
% 
% Copyright@2018-2019 Center for Compact Fusion, ENN (CCF-ENN, 
% Huasheng Xie, xiehuasheng@enn.cn, huashengxie@gmail.com)
%
% Use this file to start run BO.

% 18-12-19 10:20 
% Benchmarked:
%  1. Most ES3D results, e.g., beam, loss cone, unmagnetized;
%  2. [Umeda2012] EM3D ring beam;
%  3. [Muschietti2017] vdx with mixed magnetized and unmagnetized EM3D; 
%  4. EM3D loss cone with theta \neq 0 mirror mode, to WHAMP;
%  5. Darwin model for theta=0;
%  6. kz<0;
% Not benchmarked yet:
%  1. Polarization output;
%  2. Collision;
%  3. If you need parallize the scan of k or theta, please modify the
% kernel.m to use 'parfor', by copy the matrix.m to embed to kernel.m.

-------------------------------------------------------------------------
.Readme.txt            % this file

./doc
 - *.pdf               % documents, equations, user guide

./examples             % some typical test examples
 - * 

./code
 - bo_main.m         % use this file to run the code

 ./input
  - bo.in             % input parameters for each species
  - bo_setup.m        % initial setup of how to scan the parameters
  - bo_wpdat.m        % initial data of omega for select plots
  - bo_knwn.m        % initial data of how to normalized omega and k

 ./modules
  - bo_initialize.m   % initialize the parameters for run the code
  - bo_kernel.m       % kernel, most important part of bo code
  - bo_em3d_matrix.m   % set the bo em3d matrix, used by kernel
  - bo_es3d_matrix.m   % set the bo es3d matrix, used by kernel
  - bo_plot_all.m     % output data and plot
  - bo_plot_select.m  % subtle file for separate dispersion surfaces
  - bo_output.m       % calculate the polarizations and output data
  - funAn.m           % function An for ring beam
  - funBn.m           % function Bn for ring beam
  - funCn.m           % function Cn for ring beam

 ./output               % directory for store figures and data
  - * 


-------------------------------------------------------------------------

I will long time support this code. So, please feedback to me of your
suggestions/experience.

Enjoy!

Hua-sheng XIE, huashengxie@gmail.com
CCF-ENN, China

2019-04-09 11:23


-------------------------------------------------------------------------
---------------------------- BO-K/PDRK history --------------------------
% Ref:
%  [Xie2016] H. S. Xie & Y. Xiao, PDRK: A General Kinetic Dispersion
%    Relation Solver for Magnetized Plasma, Plasma Science and Technology,
%    Vol.18, No.2, p97 (2016). (Also 2014arXiv, note: have several typos)
%
% Documents/codes/Erratum: http://hsxie.me/codes/pdrk/
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

% 18-10-03 13:14 This version should be bugs free now.
% 18-10-13 10:37 update to with loss-cone distribution, which thus can
% handle all the same input cases as WHAMP [Ronnmark1982] code.
% 18-10-20 18:14 update to support electrostatic loss-cone.

% 18-12-28 22:48
%  1. Fix a bug of p13u & p31u in em3d-u version of em3d_matrix.m;
%  2. Add more J, e.g., 16,24, and to higher precission of previous J=8,12;

% 19-01-29 12:32
%  Fix a bug of p13u in em3d-u version of em3d_matrix.m, can agree with
% C.J.Zhang's PIC simulation now.


% 20-10-25 21:36
% Fix a bug of EM3D-M p32 vdsy term.

% 20-11-14 09:26
% Fixed some known bugs (Most of them are found by Prof. R. Denton): 
% (1)
% %       Ns(js)=2*setNs(js,2)+1; % change the N for some species
%       Ns(setNs(js,1))=2*setNs(js,2)+1; % 20-11-14 09:14
% (2)
% %       plot(pas,real(Pola(:,1,jpl,2)./(1i*Pola(:,1,jpl,1))),'-','Color',...
% %           pltc(jpl,:),'linewidth',2);
%       plot(pas,imag(Pola(:,1,jpl,2)./(1i*Pola(:,1,jpl,1))),'-','Color',...
%           pltc(jpl,:),'linewidth',2); % 20-11-14 09:00 update
% (3)
% %   [indww,indpa,indpb]=find(imag(wwn)==maxgam);% wrong! 2020-01-11 11:44
%   [indpa,indpb,indww]=find(imag(wwn)==maxgam);
% (4)
%         wg=wws(jpa,jpb,jpl)*wn+eps; % add eps to avoid singular warning
%         wg=(wws(jpa,jpb,jpl)+eps)*wn; % 20-11-14 08:46 update
