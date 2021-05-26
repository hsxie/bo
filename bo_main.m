% 2019-01-26 23:37, Hua-sheng XIE, huashengxie@gmail.com, ENN-FTRC, China.
% Ackn.: Richard DENTON (Dartmouth), Xin TAO (USTC), Jin-song ZHAO (PMO),
% Zhong-wei YANG (NSSC), Chao-jie ZHANG (UCLA), Can Huang (USTC), Wen-ya 
% LI (NSSC), Liang WANG (Princeton), Kyungguk MIN (Korean), Yang Li (ENN),
% etc ...
%
% This is the most comprehensive version of BO (means 'wave' in Chinese, 
% old name PDRK&PDRF) kinetic plasma dispersion relation solver, which
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
%  6. Support output polarizations (E,B,n_s,v_s) for electromagnetic version.
%  7. Without singularity for theta=0 & pi/2, and also support kz<0.
%  8. Provide also EM3D magnetized anisotropic multi-fluid model [Xie2014]
% with drift (vdx,vdy,vdz).
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
% with the solution omega=lambda, and (E,B) in X. The density and velocity
% perturbations of each species (n_s,v_s) can also be readily calcualted.
% 
% 
% =======================================================================
% If you use this code for your works, please cite:
%  [Xie2019] Huasheng XIE, BO: A unified tool for plasma waves and
% instabilities analysis, Computer Physics Communications, 244 (2019) 343?371.
% https://doi.org/10.1016/j.cpc.2019.06.014
%  [Xie2016] Huasheng Xie and Yong Xiao, PDRK: A General Kinetic Dispersion
% Relation Solver for Magnetized Plasma, Plasma Science and Technology, 18,
% 2, 97 (2016). DOI: 10.1088/1009-0630/18/2/01. Update/Bugs fixed at 
% http://hsxie.me/codes/pdrk/ or https://github.com/hsxie/pdrk.
% 
% Copyright@2018-2021 ENN Fusion Technology R&D Center (ENN-FTRC, Huasheng 
% Xie, xiehuasheng@enn.cn, huashengxie@gmail.com)
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

% 18-12-28 22:48
%  1. Fix a bug of p13u & p31u in em3d-u version of em3d_matrix.m;
%  2. Add more J, e.g., 16,24, and to higher precission of previous J=8,12;

% 19-01-29 12:32
%  Fix a bug of p13u in em3d-u version of em3d_matrix.m, can agree with
% C.J.Zhang's PIC simulation now.

% 19-11-30 12:56 
% Change 'imethod=0' in 'em3d_matrix.m' to 'imethod=1' as default, which
% can convergent faster with smaller N.

% 20-10-25 21:40
% Fixed a bug for EM3D-M p32 vdsy term, which only affect vdsy \neq 0 cases.
% Thanks the report by Jin-song Zhao & Yu-hang Yao.
% Seems J=4 is better than J=8,12,16 when both vdsy&vdsx \neq 0, due to
% rand off error.

% 20-11-14 09:26
% Fixed some tiny known bugs found during 191130-200502 by Xie & Denton.

% 20-11-15 06:00: 
% (1) Modify the EM3D matrix em3d_matrix.m to calculate the density and
% velocity perturbations/polarizations of each species.
% (2) Write a sub-code in kernel.m to calculate the 3-by-3 tensors 
% explicitly, say D(omega,k), K, Q, sigama and Q_s, sigma_s.

% 2020-12-22 21:16
% Add the fluid version PDRF/BO-F as an option (iem=3), which uses the same
% input and output as the kinetic versions (iem=0,1,2).
% Benchmark cold beam between iem=1&3, both frequency and polarizations
% agree each other well.

% 21-04-24 06:29 add bo_cal_df.m to calculate delta f(k,w,vx,vy,vz) for 
% electromagnetic 3D case (iem=1), using Eq.(54) of Xie2019CPC paper.

close all;
clear; clc;

% % set initial parameters in 'bo.in' and 'bo_setup.m'
disp(['***** THIS IS PLASMA DISPERSION RELATION SOLVER BO (PDRK & PDRF) V210526 *****',...
   10,'******* With Most General Extended Maxwellian Version & Fluid Version ********']);
tic; runtime=0;
disp('run ./input/bo_setup.m ...');
run ./input/bo_setup;
runtime1=toc; runtime=runtime+runtime1;
disp([' use ',num2str(runtime1),' s, toal ',num2str(runtime),'s']);

tic;
% % initialize the parameters for 'bo_kernel.m'
disp('run ./modules/bo_initialize.m ...');
run ./modules/bo_initialize;
runtime2=toc; runtime=runtime+runtime2;
disp([' use ',num2str(runtime2),' s, toal ',num2str(runtime),'s']);
if(ireturn==1)
  disp(['(ipa,ipb) = (',num2str(ipa),', ',num2str(ipb),...
    '), [see ''setup.m'', 1: k; 2: theta; 3: kz; 4: kx]']);
  disp(strscan);
  return;
end

% % print some basic info
disp('---------------------------------------------------------------');
if(iem==1)
  disp(['iem = ',num2str(iem),', this is an ** electromagnetic ** run.']);
elseif(iem==2)
  disp(['iem = ',num2str(iem),', this is an ** Darwin (a.k.a.,',...
      ' magnetoinductive or magnetostatic) ** run.']);
elseif(iem==0)
  disp(['iem = ',num2str(iem),', this is an ** electrostatic ** run.']);
elseif(iem==3)
  disp(['iem = ',num2str(iem),', this is an ** electromagnetic FLUID ** run.']);
else
  disp(['Warning: iem=',num2str(iem),' is not available in BO yet.']);
  return;
end
if(iem~=3)
disp(['N = ',num2str(N),', [Number harmonics. ! Try larger N to',...
    ' make sure the results are convergent. You can also set different N',...
    ' for different species in setup.m.]']);
disp(['J = ',num2str(J),', [J-pole, usually J=8 is sufficient; ',...
    'other choice: 2,3,4 (faster), 12,16,24 (more accurate)]']);
end
disp(['sp = ',num2str(sp),', [sp=0, eig() to obtain all solutions; sp=1, ',...
    'sparse eigs() to obtain ',10,'nw0 solutions around initial guess wg; sp=2, ',...
    'sparse to obtain one solution of ',10,'different branches with several wg0]']);

disp(['nw0 = ',num2str(nw0),', [solve nw0 solutions one time]']);
disp(['iout = ',num2str(iout),', [=1, only (omega,gamma); =2, ',...
    'also calculate (E,B)]']);

disp(['(npa,npb) = (',num2str(npa),', ',num2str(npb),...
    '), [number scan for the (1st, 2nd) variables]']);
disp(['(ipa,ipb) = (',num2str(ipa),', ',num2str(ipb),...
    '), [see ''setup.m'', 1: k; 2: theta; 3: kz; 4: kx]']);
disp(['(iloga,ilogb) = (',num2str(iloga),', ',num2str(ilogb),...
    '), [=0, linear scan; =1, 10^(pa or pb) scan]']);
disp(['pa=',num2str(min(pa)),':',num2str((pa(end)-pa(1))/(npa-1+1e-10)),...
    ':',num2str(max(pa)),', [1st parameter range]']);
disp(['pb=',num2str(min(pb)),':',num2str((pb(end)-pb(1))/(npb-1+1e-10)),...
    ':',num2str(max(pb)),', [2nd parameter range]']);
disp(['This run: pa = ',strpa,', pb = ',strpb,', ',strscan]);

disp('---------------------------------------------------------------');
disp(['S [number of species] = ',num2str(S)]);
disp(['B0 [backgroud magnetic field, Tesla] = ',num2str(B0)]);

disp('---------------------------------------------------------------');

disp(['qs0 [charge, q/qe] = ',num2str(qs0)]);
disp(['ms0 [mass, m/munit] = ',num2str(ms0)]);
disp(['ns0 [desity, m^-3] = ',num2str(ns0)]);
disp(['Tzs0 [parallel temperature, eV] = ',num2str(Tzs0)]);
disp(['Tps0 [perp temperature, eV] = ',num2str(Tps0)]);
disp(['vdsz0 [para drift velocity, vdsz/c] = ',num2str(vdsz0)]);
disp(['vdsx0 [perp drift velocity 1, vdsx/c] = ',num2str(vdsx0)]);
disp(['vdsy0 [perp drift velocity 2, vdsy/c] = ',num2str(vdsy0)]);
if(iem~=3)
disp(['vdsr0 [ring beam drift velocity, vdsr/c] = ',num2str(vdsr0)]);
disp(['alphas [loss-cone size/anisotropic] = ',num2str(alphas)]);
disp(['Deltas [loss-cone depth, =0 (max) to 1 (no)] = ',num2str(Deltas)]);
disp(['nus0 [Krook collision frequency, nu_s] = ',num2str(nus0)]);
disp(['rsa [core ratio] = ',num2str(rsab(1,:))]);
disp(['rsb [loss cone ratio] = ',num2str(rsab(2,:))]);
disp(['imus [is magnetized or unmagnetized species] = ',num2str(imus)]);
disp(['Ns [# of sum_N in the computation] = ',num2str(Ns)]);
end

disp('---------------------------------------------------------------');

disp(['lambdaDs [Debye length, m] = ',num2str(lambdaDs)]);
disp(['wps [plasma frequency, Hz] = ',num2str(wps)]);
disp(['wcs [cyclotron frequency, Hz] = ',num2str(wcs)]);

disp(['rhocs [cyclotron radius, m] = ',num2str(rhocs)]);
disp(['wps [plasma frequency, Hz] = ',num2str(wps)]);
disp(['lmdT [Tpara/Tperp] = ',num2str(lmdT)]);
disp(['betasz [parallel beta] = ',num2str(betasz)]);
disp(['betasp [perp beta] = ',num2str(betasp)]);
disp(['rhoms [mass density, kg/m^3] = ',num2str(rhoms)]);
disp(['vA [=B0/sqrt(mu0*sum(ms.*ns0)), Alfven speed, m/s] = ',num2str(vA)]);
disp(['c [speed of light, m/s] = ',num2str(sqrt(c2))]);
disp(['munit [unit mass, kg] = ',num2str(munit)]);
disp(['qe [unit charge, C] = ',num2str(qe)]);
if(iem==3)% fluid version
disp('---------------below for fluid version-------------------------');
disp(['gammazs [para polytrope exponents gamma_{para}] = ',num2str(gammazs)]);
disp(['gammaps [perp polytrope exponents gamma_{perp}] = ',num2str(gammaps)]);
disp(['csz [para pressure sound speed, m/s] = ',num2str(csz)]);
disp(['csp [perp pressure sound speed, m/s] = ',num2str(csp)]);
end

disp('---------------------------------------------------------------');
disp('******* In BO plot/output, k -> k/kn, omega -> omega/wn *******');
disp(['*** ! Defautly, wn is wcs of the 1st species; kn is wps1/c. You',10,...
    'can modify them in ''knwn.m.''']);

disp(['wn [normalized omega, Hz] = ',num2str(wn)]);
disp(['kn [normalized k, m^-1] = ',num2str(kn)]);

disp('---------------------------------------------------------------');

%%
tic;
% the kernel part of bo_em3d code, donot need change for most cases
% Or, only modify 'bo_kernel.m' lines [43-65, xxx] %
icalp=0; % do not calculate polarization in the first step
disp('run ./modules/bo_kernel.m ...');
run ./modules/bo_kernel;
runtime3=toc; runtime=runtime+runtime3;
disp([' use ',num2str(runtime3),' s, toal ',num2str(runtime),'s']);

tic;
% plot all the solutions
disp('run ./modules/bo_plot_all.m ...');
run ./modules/bo_plot_all;
runtime4=toc; runtime=runtime+runtime4;
disp([' use ',num2str(runtime4),' s, toal ',num2str(runtime),'s']);
disp('---------------------------------------------------------------');

% plot select dispersion surface, very subtle and can be further improve
input(['Please update ''wpdat'' in ''./input/bo_wpdat.m'' firstly. ',10,...
    'This step is to determine which branch(es) you hope to output/store.',...
    10,'After set the ''wpdat'', press any key to continue. ',...
    'Or ctrl+c to interrupt.',10,...
    '[Note: After run kernel.m, plot_all and plot_select can also ',...
    'run separately]']);
disp('---------------------------------------------------------------');
tic;
disp('run ./modules/bo_plot_select.m use ./input/bo_wpdat ...');
% run ./modules/bo_plot_select;
run ./input/bo_wpdat; % 18-10-21 18:04
runtime5=toc; runtime=runtime+runtime5;
disp([' use ',num2str(runtime5),' s, toal ',num2str(runtime),'s']);
%%
input('Press any key to continue for the last step: run output.');
disp('---------------------------------------------------------------');
tic;
% % output the omega and/or polarization results to data file
if(iout~=0)
  disp('run ./modules/bo_output.m ...');
  run ./modules/bo_output;
else
  disp('skip the run ./modules/bo_output.m ...');
end

isave=1; % 18-12-23 23:19
if(isave==1) % save the input file to savepath
    savepath2=savepath(2:end); % % modify here
    disp(['copy the input files in ./input/ to ',savepath2,' ...']);
    copyfile('./input/bo.in',savepath2);
    copyfile('./input/bo_setup.m',savepath2);
    copyfile('./input/bo_knwn.m',savepath2);
    copyfile('./input/bo_wpdat.m',savepath2);
end

runtime6=toc; runtime=runtime+runtime6;
disp([' use ',num2str(runtime6),' s, toal ',num2str(runtime),'s']);
disp(['Finished! Figures/data and input files have been saved to ',savepath,'.']);
