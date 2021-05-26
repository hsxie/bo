% 18-10-06 07:36 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% This file initialize the parameters for bo_kernel.m

% % read input parameters
pardat = importdata('../input/bo.in', ' ', 1);
[S, col]=size(pardat.data);
if(col~=13)
  disp('Wrong input data !!!');
end

for s=1:S
  qs0(s)=pardat.data(s,1); % charge, q/e
  ms0(s)=pardat.data(s,2); % mass, m/mp
  ns0(s)=pardat.data(s,3); % desity, m^-3
  Tzs0(s)=pardat.data(s,4); % parallel temperature, eV
  Tps0(s)=pardat.data(s,5); % perp temperature, eV
  alphas(s)=pardat.data(s,6); % loss-cone size/anisotropic
  Deltas(s)=pardat.data(s,7); % loss-cone depth, =0 (max) to 1 (no)
  vdsz0(s)=pardat.data(s,8); % para drift velocity, vds_z/c
  vdsx0(s)=pardat.data(s,9); % perp drift velocity 1, vds_x/c
  vdsy0(s)=pardat.data(s,10); % perp drift velocity 2, vds_y/c
  vdsr0(s)=pardat.data(s,11); % perp ring beam drift velocity, vds_r/c
  nus0(s)=pardat.data(s,12); % Krook collision frequency, s^{-1}
  imus(s)=pardat.data(s,13); % species magnetized (=1) or unmagnetized (=0)
  
  nrsab(s)=2; % defaultly the loss-cone has two sub-distributions
  if((alphas(s)==1) || (Deltas(s)==1)) % for no loss-cone case
    rsab(1,s)=1;
    rsab(2,s)=0;
    nrsab(s)=1;
  else
    % sigma=a, core fv ratio
    rsab(1,s)=(1-alphas(s)*Deltas(s))/(1-alphas(s));
    % sigma=b, loss cone fv ratio
    rsab(2,s)=alphas(s)*(Deltas(s)-1)/(1-alphas(s));
  end
end

Qtotal=sum(qs0.*ns0);
Jztotal=sum(qs0.*ns0.*vdsz0);
Jxtotal=sum(qs0.*ns0.*vdsx0);
Jytotal=sum(qs0.*ns0.*vdsy0);

if((Qtotal~=0) || (Jztotal~=0) || (Jxtotal~=0) || (Jytotal~=0))
  disp('Warning: Total charge or current not zero !!!');
  %input(['Warning: Total charge or current not zero!!',...
  %    'Press any key to continue.']);
end

qs=qs0*qe; % * electron charge, e -> C (coulomb)
% ms=ms0*mp; % * proton mass, m_p -> kg, 18-10-13 09:57
ms=ms0*munit; % * unit mass, -> kg, 18-12-23 22:59
Tzs=Tzs0*qe/kB; % T//, eV -> K (eV -> J * J -> K)
Tps=Tps0*qe/kB; % Tpr, eV -> K
vdsz=vdsz0*sqrt(c2); % vds_z, speed of light c -> m/s
vdsx=vdsx0*sqrt(c2); % vds_x, speed of light c -> m/s
vdsy=vdsy0*sqrt(c2); % vds_x, speed of light c -> m/s
vdsr=vdsr0*sqrt(c2); % vds_r, speed of light c -> m/s
nus=nus0; % To do, s^-1, 18-12-16 00:27

% delete these three lines if you want to see the result of magnetized
% version under B0=0
if(B0==0) % set all species to unmagnetized if B0=0
 imus=0;
end

if(iem~=3)
% find the unmagnetized species, set them vdsr=0
jds=find(imus==0);
vdsr(jds)=0;
Ns=(2*N+1)+0.*(1:S); % # of sum_N harmonic for magnetized species
if(~isempty(setNs)) % 18-12-15 23:51
  for js=1:size(setNs,1)
    if(setNs(js,1)<=S)
%       Ns(js)=2*setNs(js,2)+1; % change the N for some species
      Ns(setNs(js,1))=2*setNs(js,2)+1; % 20-11-14 09:14
    end
  end
end
Ns(jds)=nrsab(jds); % # of sum_N for unmagnetized species, inforce to nrsab
Nss=floor((Ns-1)/2); % 18-11-21 12:18
end

vtzs=sqrt(2*kB*Tzs./ms); % para thermal velocity, note the sqrt(2)
vtps=sqrt(2*kB*Tps./ms); % perp thermal velocity, note the sqrt(2)
lambdaDs=sqrt(epsilon0*kB*Tzs./(ns0.*qs.^2)); % Debye length, Tzs
kDs=1./lambdaDs;
wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
wcs=B0*qs./ms; % cyclotron frequency
% rhocs=sqrt(kB*Tps./ms)./abs(wcs); % cyclotron radius
rhocs=sqrt(kB*Tps./ms)./wcs; % cyclotron radius, 2018-06-13 21:47

wps2=wps.^2;
lmdT=Tzs./Tps;

% for sigma=a,b due to the two perp temperatures for core and loss cone fv
Tpsab(1,:)=Tps; Tpsab(2,:)=alphas.*Tps;
vtpsab(1,:)=sqrt(2*kB*Tpsab(1,:)./ms); vtpsab(2,:)=sqrt(2*kB*Tpsab(2,:)./ms);
lmdTab(1,:)=Tzs./Tpsab(1,:); lmdTab(2,:)=Tzs./Tpsab(2,:);
rhocsab(1,:)=sqrt(kB*Tpsab(1,:)./ms)./wcs;
rhocsab(2,:)=sqrt(kB*Tpsab(2,:)./ms)./wcs;
bsab(1,:)=vdsr./vtpsab(1,:); bsab(2,:)=vdsr./vtpsab(2,:);
% normalized coefficient for F(v_perp)
Asab(1,:)=exp(-bsab(1,:).^2)+sqrt(pi)*bsab(1,:).*erfc(-bsab(1,:));
Asab(2,:)=exp(-bsab(2,:).^2)+sqrt(pi)*bsab(2,:).*erfc(-bsab(2,:));

betasz=2*mu0*kB.*ns0.*Tzs./B0^2; % beta_para
betasp=2*mu0*kB.*ns0.*Tps./B0^2; % beta_perp
vA=B0/sqrt(mu0*sum(ms.*ns0)); % Alfven speed

% 20-12-22 23:13 add some new variables
rhoms=ms.*ns0;
Psz=kB.*ns0.*Tzs; % para pressure
Psp=kB.*ns0.*Tps; % perp pressure
cAs=B0./sqrt(mu0*ms.*ns0); % Alfven speed of each species
if(iem==3) % for fluid species, 20-12-22 22:56
  % set polytrope exponents for each fluid species
  gammazs=gammas0(1)+0.*ns0;
  gammaps=gammas0(2)+0.*ns0;
  if(~isempty(setgammas)) % setgammas is given in setup.m
    for js=1:size(setgammas,1)
      if(setgammas(js,1)<=S)
        gammazs(setgammas(js,1))=setgammas(js,2);
        gammaps(setgammas(js,1))=setgammas(js,3);
      end
    end
  end
  csz=sqrt(gammazs.*Psz./rhoms); % Psz=csz2.*rhos0./gammazs;
  ifluidmodel=1;
  if(ifluidmodel==2) % fluid closure model: =2 double polytrope
    csp=sqrt(Psp./rhoms); % Psp=csp2.*rhos0;
  else %  =1, Xie2014 adiabatic
    csp=sqrt(gammaps.*Psp./rhoms); % Psp=csp2.*rhos0;
  end

end

%!!- change here if you need other variables for (omega,k) normalization
% normalized by omega_c and omega_p of the first species, 18-10-18 19:25
cSs1=sqrt(2*kB*Tzs(1)/ms(1)); % sound speed
wcs1=abs(wcs(1)); % omega_{ci}
wps1=sqrt(ns0(1)*qs(1)^2/ms(1)/epsilon0);
cwp=sqrt(c2)/wps1; % c/omega_{p1}
vAwp=vA/wcs1; % v_A/omega_{c1}
wp=sqrt(sum(wps2));
lambdaD=sqrt(1./sum(1./lambdaDs.^2));


% 18-12-23 23:02, to set the wn and kn for normalization
run ../input/bo_knwn.m

%% 18-12-23 23:10 Add some more J-pole information
if(iem~=3)
% 18-12-28 22:48 calculate bj, cj to higher precission, and add J=16,24
if(J==8) % J-pole, usually J=8 is sufficient; other choice: 4, 12
  opt=3;
%   opt=2;
  if(opt==1)
    % Ronnmark1982, 8-pole for Z function, and Z', J=8, I=10
    bzj(1)=-1.734012457471826E-2-4.630639291680322E-2i;
    bzj(2)=-7.399169923225014E-1+8.395179978099844E-1i;
    bzj(3)=5.840628642184073+9.536009057643667E-1i;
    bzj(4)=-5.583371525286853-1.120854319126599E1i;
    czj(1)=2.237687789201900-1.625940856173727i;
    czj(2)=1.465234126106004-1.789620129162444i;
    czj(3)=.8392539817232638-1.891995045765206i;
    czj(4)=.2739362226285564-1.941786875844713i;
  elseif(opt==2) % J=8, I=8
    bzj(1)=  2.262756990456044 + 0.574863654552974i;
    bzj(2)=  0.007755112891899 + 0.435320280255620i;
    bzj(3)=  -0.017631167787210 - 0.001036701965463i;
    bzj(4)=  -2.752880935560734 - 4.080667522419399i;
    czj(1)=  -0.845779591068908 - 1.626106289905113i;
    czj(2)=  1.485551642694134 - 1.512332910786984i;
    czj(3)=  2.286671582218345 - 1.340185584073022i;
    czj(4)=  0.275217303800922 - 1.682797352333194i;
  else % new calculation, J=8, I=10
%     bzj(1)=  -0.0173401116032742 - 0.0463064419344598i;
%     bzj(2)=  -0.739917851897683 + 0.839518298070637i;
%     bzj(3)=  5.84063227513760 + 0.953602843950785i;
%     bzj(4)=  -5.58337431170864 - 11.2085508179677i;
%     czj(1)=   2.23768772215616 - 1.62594103256666i;
%     czj(2)=   1.46523409042510 - 1.78962030806222i;
%     czj(3)=   0.839253965702731 - 1.89199521968963i;
%     czj(4)=   0.273936217871668 - 1.94178704551807i;

    % 18-12-28 11:11
    bzj(1)=   -0.017340112270401 - 0.046306439626294i;
    bzj(2)=   -0.739917811220052 + 0.839518284620274i;
    bzj(3)=   5.840632105105495 + 0.953602751322040i;
    bzj(4)=   -5.583374181615043 -11.208550459628098i;
    czj(1)=   2.237687725134293 - 1.625941024120362i;
    czj(2)=   1.465234091939142 - 1.789620299603315i;
    czj(3)=   0.839253966367922 - 1.891995211531426i;
    czj(4)=   0.273936218055381 - 1.941787037576095i;
    
  end

  bzj(5:8)=conj(bzj(1:4));
  czj(5:8)=-conj(czj(1:4));
elseif(J==12) % from Cal_J_pole_bjcj.m
  opt=0;
  if(opt==1) % J=12; I=16; 2014;
   bzj(1)=  -0.00454786121654587 - 0.000621096230229454i;
   bzj(2)=    0.215155729087593 + 0.201505401672306i;
   bzj(3)=    0.439545042119629 + 4.16108468348292i;
   bzj(4)=  -20.2169673323552 - 12.8855035482440i;
   bzj(5)=    67.0814882450356 + 20.8463458499504i;
   bzj(6)=  -48.0146738250076 + 107.275614092570i;
  
   czj(1)=  -2.97842916245164 - 2.04969666644050i;
   czj(2)=    2.25678378396682 - 2.20861841189542i;
   czj(3)=  -1.67379985617161 - 2.32408519416336i;
   czj(4)=  -1.15903203380422 - 2.40673940954718i;
   czj(5)=    0.682287636603418 - 2.46036501461004i;
   czj(6)=  -0.225365375071350 - 2.48677941704753i;
  elseif(opt==2) % J=12; I=16; 2018-12-28 22:31
   bzj(1)= - 47.913598578418315281 - 106.98699311451399461i;
   bzj(2)= - 20.148858425809293248 + 12.874749056250453631i;
   bzj(3)= - 0.0045311004339957471789 + 0.00063311756354943215316i;
   bzj(4)= 0.2150040123642351701 + 0.20042340981056393122i;
   bzj(5)= 0.43131038679231352184 - 4.1505366661190555077i;
   bzj(6)= 66.920673705505055584 + 20.747375125403268524i;
   czj(1)= 0.2253670862838072698686 - 2.48625584284603285647646i;
   czj(2)= 1.1590491549279069691485 - 2.406192125704074076408834i;
   czj(3)= 2.9785703941315209703937799933362 - 2.0490809954949754985i;
   czj(4)= 2.25685878923092272938 - 2.208022912648570057162637i;
   czj(5)= 1.6738373878120108270775693017312 - 2.3235155478934783777i;
   czj(6)= 0.6822944098171246799968 - 2.459833442261711494628i;
  else % J=12; I=12; 18-12-28 22:35
   bzj(1)=  - 10.020983259474214017 - 14.728932929429874883i;
   bzj(2)= - 0.58878169153449514493 + 0.19067303610080007359i;
   bzj(3)= - 0.27475707659732384029 + 3.617920717493884482i;
   bzj(4)= 0.00045713742777499515344 + 0.00027155393843737098852i;
   bzj(5)=  0.017940627032508378515 - 0.036436053276701248142i;
   bzj(6)=  10.366124263145749629 - 2.5069048649816145967i;
   czj(1)= 0.22660012611958088507627 - 2.0716877594897791206264i;
   czj(2)= - 1.70029215163003500750575 - 1.8822474221612724460388i;
   czj(3)= 1.17139325085601178534269 - 1.97725033192085410977458i;
   czj(4)= 3.0666201126826972102007 - 1.59002082593259971758095i;
   czj(5)= 2.307327490410578276422 - 1.7546732543728200653674i;
   czj(6)= 0.687200524906019065672977 - 2.040288525975844018682i;
  end

  bzj(7:12)=conj(bzj(1:6));
  czj(7:12)=-conj(czj(1:6));
elseif(J==16)
  opt=1;
  if(opt==1) % J=16; I=18;  
   bzj(1)= - 86.416592794839804566 - 147.57960545984972964i;
   bzj(2)= - 22.962540986214500398 + 46.211318219085729914i;
   bzj(3)= - 8.8757833558787660662 - 11.561957978688249474i;
   bzj(4)= - 0.025134802434111256483 + 0.19730442150379382482i;
   bzj(5)= - 0.0056462830661756538039 - 0.0027884991898011769583i;
   bzj(6)= 0.000028262945845046458372 + 0.000026335348714810255537i;
   bzj(7)= 2.3290098166119338312 - 0.57238325918028725167i;
   bzj(8)= 115.45666014287557906 - 2.8617578808752183449i; 
   czj(1)= 0.1966439744113664608458045976392 - 2.5854046363167904820930238267552i;
   czj(2)= 1.000427687089304511157923736374 - 2.5277610669350594581215470753126i;
   czj(3)= 1.4263380087098663428834704281261 - 2.4694803409658086505344546718783i;
   czj(4)= 2.382753075769737513956410751299 - 2.2903917960623787648467068236658i;
   czj(5)= 2.9566517643704010426779572638885 - 2.1658992556376956216621262314559i;
   czj(6)= - 3.6699741330155866185481740497527 - 2.008727613312046260114172119472i;
   czj(7)= 1.8818356204685089975461092960437 - 2.3907395820644127767937911780402i;
   czj(8)= 0.5933003629474285223202174828712 - 2.5662607006180515205167080595386i;
  end
  bzj(9:16)=conj(bzj(1:8));
  czj(9:16)=-conj(czj(1:8));
elseif(J==24) 
  opt=1;
  if(opt==1)% J=24,I=24
    bzj(1)= - 579.77656932346560644 - 844.01436313629880827i;
    bzj(2)= - 179.52530851977905732 - 86.660002027244731382i;
    bzj(3)= - 52.107235029274485215 + 453.3246806707749413i;
    bzj(4)= - 2.1607927691932962178 + 0.63681255371973499384i;
    bzj(5)= - 0.018283386874895507814 - 0.21941582055233427677i;
    bzj(6)= - 0.00006819511737162705016 + 0.00032026091897256872621i;
    bzj(7)= - 0.0000028986123310445793648 - 0.00000099510625011385493369i;
    bzj(8)= 0.0000000023382228949223867744 - 0.0000000040404517369565098657i;
    bzj(9)= 0.01221466589423530596 + 0.00097890737323377354166i;
    bzj(10)= 7.3718296773233126912 - 12.575687057120635407i;
    bzj(11)= 44.078424019374375065 - 46.322124026599601416i;
    bzj(12)= 761.62579175738689742 + 185.11797721443392707i;
    czj(1)= 0.16167711630587375808393823760988 - 2.9424665391729649010502939606152i;
    czj(2)= 1.1509135876493567244599398043479 - 2.8745542965490153159866506667543i;
    czj(3)= 0.81513635269214329286824152984179 - 2.9085569383176322446978082849749i;
    czj(4)= 2.2362950589041724110736073820844 - 2.7033607074680388479084431872604i;
    czj(5)= 2.6403561313404041541230494846625 - 2.6228400297078984516779261304916i;
    czj(6)= 3.5620497451197056657834990483967 - 2.4245607245823420555878190731282i;
    czj(7)= 4.116925125710675393072860873751 - 2.3036541720854573608940600179944i;
    czj(8)= 4.8034117493360317933109830717707 - 2.1592490859689535412501218722927i;
    czj(9)= 3.0778922349246567316482750461458 - 2.5301774598854448463007864644617i;
    czj(10)= - 1.8572088635240765003561090479193 - 2.7720571884094886583775397071469i;
    czj(11)= 1.4969881322466893380396663902149 - 2.8290855580900544693059801078858i;
    czj(12)= - 0.48636891219330428093331493852099 - 2.9311741817223824196339069754696i;
  end
  bzj(13:24)=conj(bzj(1:12));
  czj(13:24)=-conj(czj(1:12));
  
elseif(J==4)
  opt=2;
  if(opt==1)
  % Martin1980, 4-pole, J=4, I=5
    bzj(1)=0.5468-0.0372i;
    bzj(2)=-1.0468+2.1018i;
    czj(1)=-1.2359-1.2150i;
    czj(2)=-0.3786-1.3509i;
  else % new calculation, J=4, I=5
    bzj(1)=0.546796859834032 + 0.037196505239277i;
    bzj(2)=-1.046796859834027 + 2.101852568038518i;
    czj(1)=1.23588765343592 - 1.21498213255731i;
    czj(2)=-0.378611612386277 - 1.350943585432730i;
  end
  
  bzj(3:4)=conj(bzj(1:2));
  czj(3:4)=-conj(czj(1:2));
elseif(J==3)
  % Martin1980, 3-pole, J=3, I=3
  bzj(1)=0.1822+0.5756i;
  bzj(2)=-1.3643;
  czj(1)=-0.9217-0.9091i;
  czj(2)=-1.0204i;
  
  bzj(3)=conj(bzj(1));
  czj(3)=-conj(czj(1));
elseif(J==2)
  % Huba2009, 2-pole
  bzj(1)=-(0.5+0.81i);
  czj(1)=0.51-0.81i;
  bzj(2)=conj(bzj(1));
  czj(2)=-conj(czj(1));
else
  % Martin1979, 2-pole, J=2, I=3
  bzj(1)=-(0.5+1.2891i);
  czj(1)=0.5138-1.0324i;
  bzj(2)=conj(bzj(1));
  czj(2)=-conj(czj(1));
end
% J=length(bzj);
% sum(bzj) % should be -1
% sum(bzj.*czj) % should be 0
% sum(bzj.*czj.^2) % should be -1/2

% SNJ=S*(2*N+1)*J;
% SNJ=S*(sum(Ns))*J; % 18-11-10 03:58, wrong
SNJ=(sum(Ns))*J; % 18-12-16 18:06, should be this one
% SNJ1=SNJ+1;
% if(iem==1 || iem==2) % electromagnetic case
%   SNJ3=3*SNJ1;
%   NN=SNJ3+6;
% else % electrostatic case
%   NN=SNJ1;
% end
SNJ1=SNJ+S; % 20-11-14 10:16 update for EM3D case to can calculate dJ_s 
end

if(iem==1 || iem==2) % electromagnetic case
  SNJ3=3*SNJ1;
  NN=SNJ3+6;

  % 20-11-14 23:40 to obtain the index of each species in the full matrix
  % 20-11-15 10:56 fixed bug
  ids0=zeros(1,S); ids1=0.*ids0; ids2=0.*ids0;
  for  s=1:S
      if(s==1)
          ids0(s)=1;
      else
          ids0(s)=sum(Ns(1:(s-1)))*J+1;
      end
      ids1(s)=ids0(s)-1+Ns(s)*J;
      ids2(s)=sum(Ns)*J+s;
      %   indjs=[ids0(s):1:ids1(s),ids2(s)]; % index [v_snj, js] for species #s
  end

elseif(iem==0) % electrostatic case
  NN=SNJ+1;
elseif(iem==3) % multi-fluid electromagnetic case, 20-12-22 21:31
  SJ=4*S;
  NN=SJ+6;
else
  NN=1;
  disp(['Warning: iem=',num2str(iem),' is not available in BO yet.']);
%   return;
end

if(sp==0)
  nw=NN; % number of roots to obtain
elseif(sp==1) % using sparse matrix, size(wg0)=1
%   nw=1; % !! only nw solutions around the initial guesses are given
  nw=10; % Modify here if you do not wish nw=10
  wg=wg0*wn;
elseif(sp==2) % using sparse matrix, size(wg0)>=1, 18-12-19 15:30
% thus you can set multi initial guess value, which is particular useful
% for scan roots when N is large.
  wg=wg0*wn;
  nw=size(wg0,1); % here, the meaning of nw is different from sp=0,1
end
nw0=nw;
sp0=sp;

npa=round((pa2-pa1)/dpa)+1; pa=pa1+(0:npa-1)*dpa;  
if(ipa==ipb) % if ipa==ipb, do only 1D scan of pa
  npb=1;
  pb=pa(1); % 18-10-06 01:00
else % do 2D scan (pa,pb)
  npb=round((pb2-pb1)/dpb)+1; pb=pb1+(0:npb-1)'*dpb;
end

% Typical cases of (ipa,ipb): 
%   1. (1,1) scan k, fixed theta
%   2. (2,2) scan theta, fixed k
%   3. (1,2) scan 2D (k, theta)
%   4. (3,3) scan kz, fixed kx
%   5. (4,4) scan kx, fixed kz
%   6. (3,4) scan 2D (kz,kx)
%   7. (..,..) others, please modify 'bo_si_kernel.m'
ireturn=0;
if(ipa==1 && ipb==1)
  strpa='k/k_n'; strpb='\theta^\circ';
  strscan='1. (1,1) scan k, fixed theta';
  ipbtmp=2;
elseif(ipa==2 && ipb==2)
  strpa='\theta^\circ'; strpb='k/k_n';
  strscan='2. (2,2) scan theta, fixed k';
  ipbtmp=1;
elseif(ipa==3 && ipb==3)
  strpa='k_z/k_n'; strpb='k_x/k_n';
  strscan='4. (3,3) scan kz, fixed kx';
  ipbtmp=4;
elseif(ipa==4 && ipb==4)
  strpa='k_x/k_n'; strpb='k_z/k_n';
  strscan='5. (4,4) scan kx, fixed kz';
  ipbtmp=3;
elseif(ipa==1 && ipb==2)
  strpa='k/k_n'; strpb='\theta';
  strscan='3. (1,2) scan 2D (k, theta)';
elseif(ipa==3 && ipb==4)
  strpa='k_z/k_n'; strpb='k_x/k_n';
  strscan='6. (3,4) scan 2D (kz,kx)';
else
  strpa='--'; strpb='--';
  strscan=['7. (..,..), not support this scan (ipa,ipb) yet, ',...
      'please modify ''kernel.m''!!!'];
  ireturn=1;
end

if(~exist(savepath,'dir')) % in case savepath not exist
  mkdir(savepath);
end
