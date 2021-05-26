% WHAMP data, 2018-10-14 12:05
% kx=0.5, scan kz
whpdat=[0.5000000      0.0250000    3.6828034E-02 -1.15E-03    
    0.5000000      0.0500000    7.3684177E-02 -2.29E-03    
    0.5000000      0.0750000    1.1060079E-01 -3.40E-03    
    0.5000000      0.1000000    1.4762028E-01 -4.47E-03    
    0.5000000      0.1250000    1.8480585E-01 -5.50E-03    
    0.5000000      0.1500000    2.2226146E-01 -6.48E-03    
    0.5000000      0.1750000    2.6016991E-01 -7.41E-03    
    0.5000000      0.2000000    2.9897149E-01 -8.31E-03    
    0.5000000      0.2250000    3.3998528E-01 -8.78E-03    
    0.5000000      0.2500000    3.8463804E-01 -5.26E-03    
    0.5000000      0.2750000    4.2735500E-01  7.15E-03    
    0.5000000      0.3000000    4.6121091E-01  2.31E-02    
    0.5000000      0.3250000    4.8821067E-01  3.67E-02    
    0.5000000      0.3500000    5.1125847E-01  4.69E-02    
    0.5000000      0.3750000    5.3189965E-01  5.38E-02    
    0.5000000      0.4000000    5.5091639E-01  5.76E-02    
    0.5000000      0.4250000    5.6872861E-01  5.87E-02    
    0.5000000      0.4500000    5.8557270E-01  5.70E-02    
    0.5000000      0.4750000    6.0158109E-01  5.28E-02    
    0.5000000      0.5000000    6.1681851E-01  4.60E-02    
    0.5000000      0.5250000    6.3129744E-01  3.65E-02    
    0.5000000      0.5500000    6.4498180E-01  2.44E-02    
    0.5000000      0.5750000    6.5778182E-01  9.44E-03    
    0.5000000      0.6000000    6.6953883E-01 -8.52E-03    
    0.5000000      0.6250000    6.7999380E-01 -2.98E-02    
    0.5000000      0.6500000    6.8872267E-01 -5.47E-02    
    0.5000000      0.6750000    6.9499203E-01 -8.40E-02    
    0.5000000      0.7000000    6.9739061E-01 -1.18E-01    
    0.5000000      0.7250000    6.9269029E-01 -1.60E-01    
    0.5000000      0.7500000    6.7176398E-01 -2.07E-01    
    0.5000000      0.7750000    6.2237387E-01 -2.39E-01    
    0.5000000      0.8000000    5.8044501E-01 -2.36E-01    
    0.5000000      0.8250000    5.5785062E-01 -2.28E-01    
    0.5000000      0.8499999    5.4536556E-01 -2.21E-01    
    0.5000000      0.8749999    5.3818703E-01 -2.18E-01    
    0.5000000      0.8999999    5.3407403E-01 -2.18E-01    
    0.5000000      0.9249999    5.3187209E-01 -2.20E-01    
    0.5000000      0.9499999    5.3092862E-01 -2.24E-01    
    0.5000000      0.9749999    5.3084659E-01 -2.29E-01    
    0.5000000      0.9999999    5.3137013E-01 -2.36E-01 ];
kxwhp=whpdat(:,1);
kzwhp=whpdat(:,2);
wrwhp=whpdat(:,3);
wiwhp=whpdat(:,4);
subplot(121); hold on; grid on;
plot(kzwhp(1:2:end),wrwhp(1:2:end),'kx','Linewidth',2);
% xlabel('k_zc/\omega_{pi}');
subplot(122); hold on; grid on;
plot(kzwhp(1:2:end),wiwhp(1:2:end),'kx','Linewidth',2);
ylim([-1.0,0.1]);
% xlabel('k_zc/\omega_{pi}');
legend('bo-mode 1','bo-mode 2','bo-mode 3','bo-mode 4',...
    'whamp','location','southwest');
legend('boxoff');
