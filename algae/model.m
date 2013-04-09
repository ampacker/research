function err = model(params,m)
%MODEL Solves model, returns MSE, & produces figures for Packer et. al. 2011
%<http://dx.doi.org/10.1016/j.biortech.2010.06.029>
%   err = MODEL uses parameters from paper
%
%   err = MODEL(params) uses values in the array params and m = 2 for
%   exponent of p_m function
%
%         k = params(1);
%         gamma = params(2);
%         q = params(3);
%         q_m = params(4);
%         v_m = params(5);
%         v_h = params(6);
%         muq_m= params(7);
%         rho_m = params(8);
%         p_m0 = params(9);
%         c = params(10);
%         Q0 = params(11);
%         Ch0 = params(12);
%
%   err = MODEL(params, n) uses values in the array params and m as exponent
%   for p_m function
%
%(c) 2013 Aaron Packer.
%Licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License
%<http://creativecommons.org/licenses/by-sa/3.0/deed.en_US>

global VARS


I0          = 600/100;                  % umol photons / dm^2 / s (fixed irradiance) (really is 5.99537037037037)
z           = 0.03*10;                  % dm (fixed depth of reactor)

if nargin == 0
    % these are unrounded parameter values from paper. Note conversion
    % factors below. In this script, we are working in L (=dm^3), seconds,
    % and micromoles. Therefore the values from paper need to be converted
    % accordingly.
    
    %1 day = 86400 seconds
    %1 umol = 1e-6 mole
    %1 mole = 1000000 umoles
    %1 dm^2 = 0.01 m^2
    %1 m^2 = 100 dm^2
    
    k = 4.82*100;                   % dm^2 / g chl
    gamma = (9.84e-2)*12*(1e-6);       % g C / umol photons
    q = 0.0278;
    q_m = 0.0935;
    v_m = 0.596/86400;
    v_h = (1.03e-5)*0.001;
    muq_m = 3.26/86400;
    rho_m = 0.283;
    p_m0 = 90.1/86400;              % g C / g chl / s
    c = 0.610;                      % g C / g dw
    Q0 = 0.06;                      
    Ch0 = 0.0052;
    m = 2;
else
    k = params(1);
    gamma = params(2);
    q = params(3);
    q_m = params(4);
    v_m = params(5);
    v_h = params(6);
    muq_m= params(7);
    rho_m = params(8);
    p_m0 = params(9);
    c = params(10);
    Q0 = params(11);
    Ch0 = params(12);
    if nargin == 1, m = 2; end
end

VARS(1) = I0;
VARS(2) = z;
VARS(3) = k;
VARS(4) = gamma;
VARS(5) = q;
VARS(6) = v_m;
VARS(7) = v_h;
VARS(8) = q_m;
VARS(9) = muq_m;
VARS(10) = rho_m;
VARS(11) = p_m0;
VARS(12) = c;
VARS(13) = m;

t0 = 0;
tf = 288*3600;      % seconds (3600 s in 1 hour)

A0 = 0.150;
N0 = 0.240;
L0 = 0;

Y0 = [A0 Q0 N0 Ch0 L0];

options = odeset('RelTol',1e-8,'AbsTol',1e-9);

Y0(3) = 0.244;
sol1 = ode15s(@m2,[t0 tf],Y0,options);
[sol,deriv] = deval(sol1,[t0:3600:tf]);
A1 = sol(1,:);
Q1 = sol(2,:);
Q1b = Q1.*A1;
N1 = sol(3,:);
Ch1 = sol(4,:);
Ch1b = Ch1.*A1;
L1 = sol(5,:);

A1d = deriv(1,:);
L1d = deriv(5,:);
t1 = [t0:3600:tf]/3600/24;
% A1 = A1 + L1;
% A1d = (A1d+L1d);
Q1b = Q1b./(A1+L1);
Ch1b = Ch1b./(A1+L1);


Y0(3) = 0.060;
sol2 = ode15s(@m2,[t0 tf],Y0,options);
[sol,deriv] = deval(sol2,[t0:3600:tf]);
A2 = sol(1,:);
Q2 = sol(2,:);
Q2b = Q2.*A2;
N2 = sol(3,:);
Ch2 = sol(4,:);
Ch2b = Ch2.*A2;
L2 = sol(5,:);

A2d = deriv(1,:);
L2d = deriv(5,:);
t2 = [t0:3600:tf]/3600/24;
% A2 = A2 + L2;
% A2d = (A2d+L2d);
Q2b = Q2b./(A2+L2);
Ch2b = Ch2b./(A2+L2);

Y0(3) = 0;
sol3 = ode15s(@m2,[t0 tf],Y0,options);
[sol,deriv] = deval(sol3,[t0:3600:tf]);
A3 = sol(1,:);
Q3 = sol(2,:);
Q3b = Q3.*A3;
N3 = sol(3,:);
Ch3 = sol(4,:);
Ch3b = Ch3.*A3;
L3 = sol(5,:);

A3d = deriv(1,:);
L3d = deriv(5,:);
t3 = [t0:3600:tf]/3600/24;
% A3 = A3+L3;
% A3d = (A3d+L3d);
Q3b = Q3b./(A3+L3);
Ch3b = Ch3b./(A3+L3);

t0 = t0/24/3600; tf = tf/24/3600;

HLdays=[0 1 2 3 4 6 8 10 12];

HL000 = [0.15 0.29 0.31 0.33 0.41 0.44 0.56 0.55 0.56];
HL000E =[0.02 0.01 0.04 0.02 0.04 0.05 0.05 0.06 0.03];
HL000Yield(2:9) = (HL000(2:9) - HL000(1:8))./(HLdays(2:9)-HLdays(1:8));
HL000YieldE(2:9) = (((HL000(2:9)+HL000E(2:9))-(HL000(1:8)+HL000E(1:8)))-HL000Yield(2:9))./(HLdays(2:9)-HLdays(1:8));

HL060 = [0.15 0.96 2.41 2.93 3.52 4.71 5.28 5.75 5.91];
HL060E =[0.02 0.11 0.12 0.19 0.29 0.45 0.34 0.27 0.33];
HL060Yield(2:9) = (HL060(2:9) - HL060(1:8))./(HLdays(2:9)-HLdays(1:8));
HL060YieldE(2:9) = (((HL060(2:9)+HL060E(2:9))-(HL060(1:8)+HL060E(1:8)))-HL060Yield(2:9))./(HLdays(2:9)-HLdays(1:8));

HL240 = [0.15 0.78 2.38 3.87 4.65 5.85 7.01 7.60 8.02];
HL240E =[0.02 0.06 0.16 0.24 0.33 0.31 0.42 0.38 0.29];
HL240Yield(2:9) = (HL240(2:9) - HL240(1:8))./(HLdays(2:9)-HLdays(1:8));
HL240YieldE(2:9) = (((HL240(2:9)+HL240E(2:9))-(HL240(1:8)+HL240E(1:8)))-HL240Yield(2:9))./(HLdays(2:9)-HLdays(1:8));

LHL000 = [0.00 0.00 2.80 8.10 18.1 29.6 36.2 43.6 44.1].*HL000/100;
LHL000E =[0.00 0.00 0.40 1.10 2.10 3.25 3.36 1.50 2.30].*HL000/100;
LHL000Yield(2:9) = (LHL000(2:9) - LHL000(1:8))./(HLdays(2:9)-HLdays(1:8));
LHL000YieldE(2:9) = (((LHL000(2:9)+LHL000E(2:9))-(LHL000(1:8)+LHL000E(1:8)))-LHL000Yield(2:9))./(HLdays(2:9)-HLdays(1:8));
LHL000 = [0.00 0.00 2.80 8.10 18.1 29.6 36.2 43.6 44.1];
LHL000E =[0.00 0.00 0.40 1.10 2.10 3.25 3.36 1.50 2.30];

LHL060 = [0.00 0.00 4.60 22.8 32.8 42.0 46.1 52.5 55.1].*HL060/100;
LHL060E =[0.00 0.00 1.01 2.00 2.10 0.60 1.59 2.30 3.29].*HL060/100;
LHL060Yield(2:9) = (LHL060(2:9) - LHL060(1:8))./(HLdays(2:9)-HLdays(1:8));
LHL060YieldE(2:9) = (((LHL060(2:9)+LHL060E(2:9))-(LHL060(1:8)+LHL060E(1:8)))-LHL060Yield(2:9))./(HLdays(2:9)-HLdays(1:8));
LHL060 = [0.00 0.00 4.60 22.8 32.8 42.0 46.1 52.5 55.1];
LHL060E =[0.00 0.00 1.01 2.00 2.10 0.60 1.59 2.30 3.29];

LHL240 = [0.00 0.00 0.00 4.60 11.5 22.7 30.9 38.3 42.1].*HL240/100;
LHL240E =[0.00 0.00 0.00 0.70 0.40 0.96 2.80 0.86 2.39].*HL240/100;
LHL240Yield(2:9) = (LHL240(2:9) - LHL240(1:8))./(HLdays(2:9)-HLdays(1:8));
LHL240YieldE(2:9) = (((LHL240(2:9)+LHL240E(2:9))-(LHL240(1:8)+LHL240E(1:8)))-LHL240Yield(2:9))./(HLdays(2:9)-HLdays(1:8));
LHL240 = [0.00 0.00 0.00 4.60 11.5 22.7 30.9 38.3 42.1];
LHL240E =[0.00 0.00 0.00 0.70 0.40 0.96 2.80 0.86 2.39];

HLN000 = [0 0 0 0 0 0 0 0]/1000;
HLN000E = [0 0 0 0 0 0 0 0]/1000;
HLN060 = [59.4 0 0 0 0 0 0 0]/1000;
HLN060E = [2.35 0 0 0 0 0 0 0]/1000;
HLN240 = [244.1 211.5 59.5 0.486 0 0 0 0]/1000;
HLN240E = [6.2 5.9 2.2 0.04 0 0 0 0]/1000;
HLNdays = [0 1 2 3 4 6 8 10];

days = [0:24*3600:tf*24*3600];
s1 = deval(sol1,days);
s2 = deval(sol2,days);
s3 = deval(sol3,days);

for i = 2:length(days)
    yield1(i) = s1(1,i)+s1(5,i)-s1(1,i-1)-s1(5,i-1);
    Lyield1(i) = s1(5,i)-s1(5,i-1);
end
for i = 2:length(days)
    yield2(i) = s2(1,i)+s2(5,i)-s2(1,i-1)-s2(5,i-1);
    Lyield2(i) = s2(5,i)-s2(5,i-1);
end
for i = 2:length(days)
    yield3(i) = s3(1,i)+s3(5,i)-s3(1,i-1)-s3(5,i-1);
    Lyield3(i) = s3(5,i)-s3(5,i-1);
end
yield1(1) = 0;  Lyield1(1) = 0;
yield2(1) = 0;  Lyield2(1) = 0;
yield3(1) = 0;  Lyield3(1) = 0;

figure(10)
clf
hold on;
plot(t2, A2+L2, '--', t3, A3+L3, '-')
errorbar(HLdays, HL060, HL060E, 'ks', 'MarkerSize',4);
errorbar(HLdays, HL000, HL000E, 'kx', 'MarkerSize',4);
legend('25% (Model)','0% (Model)','25% (Data)','0% (Data)','Location','Best')
title('Dry weight');
ylabel('g/L')
xlabel('Days')
axis([0 12 0 7]);

% Uncomment for DW figure (100%)
% figure(11)
% clf
% hold on;
% plot(t1, A1+L1, '-')
% errorbar(HLdays, HL240, HL240E, 'kd', 'MarkerSize',4);
% legend('100% (Model)','100% (Data)','Location','Best')
% title('Dry weight');
% ylabel('g/L')
% xlabel('Days')
% axis ([0 12 0 19]);

figure(12)
clf;
hold on;
errorbar(HLdays, HL240, HL240E, 'kd-', 'MarkerSize',6);
errorbar(HLdays, HL060, HL060E, 'ks-', 'MarkerSize',6);
errorbar(HLdays, HL000, HL000E, 'kx-', 'MarkerSize',6);
title('Dry Weght (Data only)')
ylabel('g dw/L')
xlabel('Days')
legend('100%','25%','0%','Location','Best')
axis([0 12 0 9]);

figure(20)
clf;
hold on;
plot(days/24/3600,yield2,'--',days/24/3600,yield3,'-');
errorbar(HLdays, HL060Yield, HL060YieldE, 'ks', 'MarkerSize',4);
errorbar(HLdays, HL000Yield, HL000YieldE, 'kx', 'MarkerSize',4);
title('Biomass yield')
ylabel('g dw/L/day')
xlabel('Days')
legend('25% (Model)','0% (Model)','25% (Data)','0% (Data)','Location','Best')
axis([0 12 0 2]);

% Uncomment for DW yield figure (100%)
% figure(21)
% clf;
% hold on;
% plot(days/24/3600,yield1,'-');
% errorbar(HLdays, HL240Yield, HL240YieldE, 'kd', 'MarkerSize',4);
% title('Biomass yield')
% ylabel('g dw/L/day')
% xlabel('Days')
% legend('100% (Model)','100% (Model)','Location','Best')
% axis([0 12 0 3]);

figure(30)
clf;
hold on;
plot(t2,100*L2./(A2+L2),'--',t3,100*L3./(A3+L3),'-');
errorbar(HLdays, LHL060, LHL060E, 'ks', 'MarkerSize',4);
errorbar(HLdays, LHL000, LHL000E, 'kx', 'MarkerSize',4);
title('NL % of biomass');
legend('25% (Model)','0% (Model)','25% (Data)','0% (Data)','Location','Best')
ylabel('%')
xlabel('Days')
axis([t0 tf 0 80]);

figure(32)
clf;
hold on;
errorbar(HLdays, LHL240, LHL240E, 'kd-', 'MarkerSize',6);
errorbar(HLdays, LHL060, LHL060E, 'ks-', 'MarkerSize',6);
errorbar(HLdays, LHL000, LHL000E, 'kx-', 'MarkerSize',6);
title('NL % of biomass (data only)');
ylabel('%')
xlabel('Days')
legend('100%','25%','0%','Location','Best')
axis([0 12 0 60]);

figure(40)
clf;
hold on;
plot(t2,L2,'--',t3,L3,'-');
errorbar(HLdays, LHL060.*HL060/100, LHL060E.*HL060/100, 'ks', 'MarkerSize',4);
errorbar(HLdays, LHL000.*HL000/100, LHL000E.*HL000/100, 'kx', 'MarkerSize',4);
title('NL density');
legend('25% (Model)','0% (Model)','25% (Data)','0% (Data)','Location','Best')
ylabel('g NL / L')
xlabel('Days')
axis([0 12 0 4]);

% Uncomment for NL yield figure (0% and 25%)
% figure(50)
% clf;
% hold on;
% plot(days/24/3600,Lyield2,'--',days/24/3600,Lyield3,'-');
% errorbar(HLdays, LHL060Yield, LHL060YieldE, 'ks', 'MarkerSize',4);
% errorbar(HLdays, LHL000Yield, LHL000YieldE, 'kx', 'MarkerSize',4);
% title('NL yield')
% ylabel('g dw/L/day')
% xlabel('Days')
% legend('25% (Model)','0% (Model)','25% (Data)','0% (Data)','Location','Best')
% axis([0 12 0 1]);

% Uncomment for NL yield figure (100%)
% figure(51)
% clf;
% hold on;
% plot(days/24/3600,Lyield1,'-');
% errorbar(HLdays, LHL240Yield, LHL060YieldE, 'kd', 'MarkerSize',4);
% title('NL yield')
% ylabel('g dw/L/day')
% xlabel('Days')
% legend('100% (Model)','100% (Model)','Location','Best')
% axis([0 12 0 1.7]);

figure(60)
clf
hold on
plot(t1,N1,':',t2,N2,'--',t3,N3,'-');
errorbar(HLNdays, HLN240, HLN240E, 'kd', 'MarkerSize',4);
errorbar(HLNdays, HLN060, HLN060E, 'ks', 'MarkerSize',4);
errorbar(HLNdays, HLN000, HLN000E, 'kx', 'MarkerSize',4);
title('Extracellular nitrogen concentration');
ylabel('g N/L')
xlabel('Days')
legend('100% (Model)','25% (Model)','0% (Model)','100% (Data)','25% (Data)','0% (Data)','Location','Best')
axis([t0 max(HLNdays) 0 0.250]);

% TODO: assume Q(0) is true then calculate Q(t) at each day using the fact
% that the system is closed for N. Then compare to model predictions.
figure(70)
clf;
hold on;
plot(t1,Q1,':',t2,Q2,'--',t3,Q3,'-');
title('Q(t)');
ylabel('g N/g dw')
xlabel('Days')
legend('100% (Model)','25% (Model)','0% (Model)','Location','Best')
axis([t0 tf 0 1.05*max(Q1)])

figure(80)
clf;
hold on;
plot(t1,Ch1,':',t2,Ch2,'--',t3,Ch3,'-');
title('H(t)');
ylabel('g chl a / g dw')
xlabel('Days')
legend('100% (Model)','25% (Model)','0% (Model)','Location','Best')
axis xy

%calculate error
tt = [0 1 2 3 4 6 8 10 12]*24 +1;       % indices for DW and NL values
ttN = [0 1 2 3 4 6 8 10]*24 +1;         % indices for N values

A1e = A1(tt);
L1e = L1(tt);
N1e = N1(ttN);
A2e = A2(tt);
L2e = L2(tt);
N2e = N2(ttN);
A3e = A3(tt);
L3e = L3(tt);
N3e = N3(ttN);

err100 = [norm(A1e+L1e-HL240), norm(L1e-(LHL240.*HL240/100)), norm(N1e-HLN240)];
err25 = [norm(A2e+L2e-HL060), norm(L2e-(LHL060.*HL060/100)), norm(N2e-HLN060)];
err0 = [norm(A3e+L3e-HL000), norm(L3e-(LHL000.*HL000/100)), 0];

err = [err0; err25; err100];

%-------------------------------------------------
function dydt = m2(t,y)
global VARS

I0 = VARS(1);
z = VARS(2);
k = VARS(3);
gamma = VARS(4);
q = VARS(5);
v_m = VARS(6);
v_h = VARS(7);
q_m = VARS(8);
muq_m = VARS(9);
rho_m = VARS(10);
p_m0 = VARS(11);
c = VARS(12);
m = VARS(13);

A = y(1);
Q = y(2);
N = y(3);
H = y(4);
L = y(5);
Q_tilde = Q*A/(A+L);

v = v_m*((q_m-Q)/(q_m-q))*(N/(v_h+N));

p_m = p_m0*Q_tilde^m/(Q_tilde^m+q^m);
l = I0*(1-exp(-k*H*A*z))/(k*A*H*z);
fl = (1-exp(-k*gamma*l/p_m));
p = H*p_m*fl;

fq = 1-q/Q;
muq = muq_m*fq;
muq = max([0,muq]);
mu = min([muq,p/c]);

rho = rho_m*c*mu/p;
mul = p-c*mu;

dydt(1,1) = mu*A;
dydt(2,1) = v-Q*mu;
dydt(3,1) = -v*A;
dydt(4,1) = rho*v-H*mu;
dydt(5,1) = mul*A;
%-------------------------------------------------