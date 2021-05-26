# BO
Plasma waves and instabilities analysis tool, with both kinetic and multi-fluid dispersion relations.

This is a combination of previous several versions of PDRF/PDRK/BO, and with new features.

[Xie2019] H.S. Xie, BO: A unified tool for plasma waves and instabilities analysis, Comput. Phys. Comm. 244 (2019) 343-371. http://code.ennresearch.com/code/bo/

[Xie2016]  H.S. Xie, Y. Xiao, PDRK: A General Kinetic Dispersion Relation Solver for Magnetized Plasma, Plasma Sci. Technol. 18 (2) (2016) 97, http://dx.doi.org/10.1088/1009-0630/18/2/01, Update/bugs fixed at http://hsxie.me/codes/pdrk/ or
https://github.com/hsxie/pdrk/.

[Xie2014] H. S. Xie, PDRF: A general dispersion relation solver for magnetized multi-fluid plasma, Comput. Phys. Comm. 185 (2014) 670-675.

The March 2021 update version: 
[Xie2021] Hua-sheng Xie, Richard Denton, Jin-song Zhao and Wen Liu, BO 2.0: Plasma Wave and Instability Analysis with Enhanced Polarization Calculations, arXiv, 2021. https://arxiv.org/abs/2103.16014


Major updates:

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

