function err = model(params)
%MODEL Solves model, returns MSE, & produces figures for Packer et. al. 2011
%<http://dx.doi.org/10.1016/j.biortech.2010.06.029>
%   err = MODEL uses parameters from paper
%
%   err = MODEL(params) uses values in the array params
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
%(c) 2013 Aaron Packer.

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
    Ch0 = 0.0061;
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
end

t0 = 0;
tf = 12*24*3600;      % seconds (3600 s in 1 hour)
tsol = t0:3600:tf;
tplot = tsol/24/3600;

A0 = 0.150;
N0 = 0.240;
L0 = 0;

Y0 = [A0 Q0 N0 Ch0 L0];

options = odeset('RelTol',1e-8,'AbsTol',1e-9);

Y0(3) = 0.244;
sol1 = ode15s(@m2,[t0 tf],Y0,options);
[sol,deriv] = deval(sol1,tsol);
A1 = sol(1,:);
Q1 = sol(2,:);
Q1b = Q1.*A1;
N1 = sol(3,:);
Ch1 = sol(4,:);
Ch1b = Ch1.*A1;
L1 = sol(5,:);
DW1 = A1+L1;
Qtil1 = Q1b./DW1;

Y0(3) = 0.060;
sol2 = ode15s(@m2,[t0 tf],Y0,options);
[sol,deriv] = deval(sol2,tsol);
A2 = sol(1,:);
Q2 = sol(2,:);
Q2b = Q2.*A2;
N2 = sol(3,:);
Ch2 = sol(4,:);
Ch2b = Ch2.*A2;
L2 = sol(5,:);
DW2 = A2+L2;
Qtil2 = Q2b./DW2;

Y0(3) = 0;
sol3 = ode15s(@m2,[t0 tf],Y0,options);
[sol,deriv] = deval(sol3,tsol);
A3 = sol(1,:);
Q3 = sol(2,:);
Q3b = Q3.*A3;
N3 = sol(3,:);
Ch3 = sol(4,:);
Ch3b = Ch3.*A3;
L3 = sol(5,:);
DW3 = A3+L3;
Qtil3 = Q3b./DW3;


%% data
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

LHLf000 = [0.00 0.00 2.80 8.10 18.1 29.6 36.2 43.6 44.1];
LHLf000E =[0.00 0.00 0.40 1.10 2.10 3.25 3.36 1.50 2.30];
LHL000 = LHLf000.*HL000/100;
LHL000E = LHLf000E.*HL000/100;
LHL000Yield(2:9) = (LHL000(2:9) - LHL000(1:8))./(HLdays(2:9)-HLdays(1:8));
LHL000YieldE(2:9) = (((LHL000(2:9)+LHL000E(2:9))-(LHL000(1:8)+LHL000E(1:8)))-LHL000Yield(2:9))./(HLdays(2:9)-HLdays(1:8));

LHLf060 = [0.00 0.00 4.60 22.8 32.8 42.0 46.1 52.5 55.1];
LHLf060E =[0.00 0.00 1.01 2.00 2.10 0.60 1.59 2.30 3.29];
LHL060 = LHLf060.*HL060/100;
LHL060E = LHLf060E.*HL060/100;
LHL060Yield(2:9) = (LHL060(2:9) - LHL060(1:8))./(HLdays(2:9)-HLdays(1:8));
LHL060YieldE(2:9) = (((LHL060(2:9)+LHL060E(2:9))-(LHL060(1:8)+LHL060E(1:8)))-LHL060Yield(2:9))./(HLdays(2:9)-HLdays(1:8));

LHLf240 = [0.00 0.00 0.00 4.60 11.5 22.7 30.9 38.3 42.1];
LHLf240E =[0.00 0.00 0.00 0.70 0.40 0.96 2.80 0.86 2.39];
LHL240 = LHLf240.*HL240/100;
LHL240E =LHLf240E.*HL240/100;
LHL240Yield(2:9) = (LHL240(2:9) - LHL240(1:8))./(HLdays(2:9)-HLdays(1:8));
LHL240YieldE(2:9) = (((LHL240(2:9)+LHL240E(2:9))-(LHL240(1:8)+LHL240E(1:8)))-LHL240Yield(2:9))./(HLdays(2:9)-HLdays(1:8));

AHL240 = HL240-LHL240;
AHL240E = HL240E-LHL240E;
AHL060 = HL060-LHL060;
AHL060E = HL060E-LHL060E;
AHL000 = HL000-LHL000;
AHL000E = HL000E-LHL000E;

HLN000 = [0 0 0 0 0 0 0 0]/1000;
HLN000E = [0 0 0 0 0 0 0 0]/1000;
HLN060 = [59.4 0 0 0 0 0 0 0]/1000;
HLN060E = [2.35 0 0 0 0 0 0 0]/1000;
HLN240 = [244.1 211.5 59.5 0.486 0 0 0 0]/1000;
HLN240E = [6.2 5.9 2.2 0.04 0 0 0 0]/1000;
HLNdays = [0 1 2 3 4 6 8 10];

HLQA240 = Q0*A0+HLN240(1)-HLN240;
HLQA060 = Q0*A0+HLN060(1)-HLN060;
HLQA000 = Q0*A0+HLN000(1)-HLN000;

HLQtil240 = HLQA240./HL240(1:end-1);
HLQtil240E = (HLQA240./HL240(1:end-1)).*HL240E(1:end-1)./(HL240(1:end-1)-HL240E(1:end-1));
HLQtil060 = HLQA060./HL060(1:end-1);
HLQtil060E = (HLQA060./HL060(1:end-1)).*HL060E(1:end-1)./(HL060(1:end-1)-HL060E(1:end-1));
HLQtil000 = HLQA000./HL000(1:end-1);
HLQtil000E = (HLQA000./HL000(1:end-1)).*HL000E(1:end-1)./(HL000(1:end-1)-HL000E(1:end-1));

HLQ240 = HLQA240./AHL240(1:end-1);
HLQ240E = (HLQA240./AHL240(1:end-1)).*AHL240E(1:end-1)./(AHL240(1:end-1)-AHL240E(1:end-1));
HLQ060 = HLQA060./AHL060(1:end-1);
HLQ060E = (HLQA060./AHL060(1:end-1)).*AHL060E(1:end-1)./(AHL060(1:end-1)-AHL060E(1:end-1));
HLQ000 = HLQA000./AHL000(1:end-1);
HLQ000E = (HLQA000./AHL000(1:end-1)).*AHL000E(1:end-1)./(AHL000(1:end-1)-AHL000E(1:end-1));

%% calculate DW and lipid yield
days = (0:12)*24;
daysA = days(2:end)+1; daysB = days(1:end-1)+1;
yield1 = zeros(length(days),1); Lyield1=yield1; yield2=yield1; Lyield2=Lyield1; yield3=yield1; Lyield3=Lyield2;
yield1(2:end) = DW1(daysA)-DW1(daysB);
Lyield1(2:end) = L1(daysA)-L1(daysB);
yield2(2:end) = DW2(daysA)-DW2(daysB);
Lyield2(2:end) = L2(daysA)-L2(daysB);
yield3(2:end) = DW3(daysA)-DW3(daysB);
Lyield3(2:end) = L3(daysA)-L3(daysB);

%% plots
figure(10)
clf
hold on;
plot(tplot, DW2, 'b-');
plot(tplot, DW3, 'g-', 'Color',[0 .5 0]);
errorbar(HLdays, HL060, HL060E, 'bs', 'MarkerSize',6);
errorbar(HLdays, HL000, HL000E, 'gx', 'MarkerEdgeColor',[0 .5 0], 'Color',[0 .5 0], 'MarkerSize',6);
legend('25%','0%','Location','Best')
% title('Dry weight');
ylabel('g/L')
xlabel('d')
axis([0 12 0 7]);

% Uncomment for DW figure (100%)
figure(11)
clf
hold on;
plot(tplot, DW1, '-')
errorbar(HLdays, HL240, HL240E, 'kd', 'MarkerSize',4);
legend('100% (Model)','100% (Data)','Location','Best')
title('Dry weight');
ylabel('g/L')
xlabel('d')
axis ([0 12 0 19]);

figure(30)
clf;
hold on;
plot(tplot,100*L2./DW2,'b-');
plot(tplot,100*L3./DW3,'g-','Color',[0 .5 0]);
errorbar(HLdays, LHLf060, LHLf060E, 'bs', 'MarkerSize',6);
errorbar(HLdays, LHLf000, LHLf000E, 'gx', 'MarkerEdgeColor',[0 .5 0], 'Color',[0 .5 0], 'MarkerSize',6);
% title('NL % of biomass');
legend('25%','0%','Location','Best')
ylabel('%')
xlabel('d')
axis([t0 12 0 80]);

figure(40)
clf;
hold on;
plot(tplot,L2,'b-');
plot(tplot,L3,'g-','Color',[0 .5 0]);
errorbar(HLdays, LHL060, LHL060E, 'bs', 'MarkerSize',6);
errorbar(HLdays, LHL000, LHL000E, 'gx', 'MarkerEdgeColor',[0 .5 0], 'Color',[0 .5 0], 'MarkerSize',6);
% title('NL concentration');
legend('25%','0%','Location','Best')
ylabel('g NL/L')
xlabel('d')
axis([0 12 0 4]);
% axis tight;
% xlim([0 12]);

% N data + model
figure(60)
clf
hold on
plot(tplot,N1,'k-',tplot,N2,'b-');
plot(tplot,N3,'g-','Color',[0 .5 0]);
errorbar(HLNdays, HLN240, HLN240E, 'kd', 'MarkerSize',6);
errorbar(HLNdays, HLN060, HLN060E, 'bs', 'MarkerSize',6);
errorbar(HLNdays, HLN000, HLN000E, 'gx', 'MarkerEdgeColor',[0 .5 0], 'Color',[0 .5 0], 'MarkerSize',6);
hold off
% title('Extracellular nitrogen concentration');
ylabel('g N/L')
xlabel('d')
legend('100%','25%','0%','Location','Best');
% legend('100% (model)','25%','0%','100% (data)','25%','0%','Location','Best');
axis([t0 max(HLNdays) 0 0.250]);

figure(70)
clf;
hold on;
% plot(tplot,Qtil1,'k-');
plot(tplot,Qtil2,'b-');
plot(tplot,Qtil3,'g-','Color',[0 .5 0]);
% plot(tplot,Q1,'k--');
plot(tplot,Q2,'b--');
plot(tplot,Q3,'g--', 'Color',[0 .5 0]);
% plot(HLNdays,HLQ240,'kd');
plot(HLNdays,HLQtil060,'bs')
plot(HLNdays,HLQtil000,'gd','Color',[0 .5 0]);
% plot(HLNdays,HLQtil240,'kd');
% plot(HLNdays,HLQ060,'bs','MarkerFaceColor','b');
% plot(HLNdays,HLQ000,'gd','Color',[0 .5 0],'MarkerFaceColor',[0 .5 0]);
hold off;
% title('\textsf{$Q$ and $\widetilde{Q}$}','interpreter','latex');
ylabel('g N/g dw')
xlabel('d')
% legend('100% (model)','25%','0%','Location','Best')
legend('25% Q','0% Q','Location','Best')
axis([0 12 0 1.05*max(Q1)])

figure(80)
clf;
hold on;
plot(tplot,Ch1,'k-');
plot(tplot,Ch2,'b-');
plot(tplot,Ch3,'g-','Color',[0 .5 0]);
hold off
% title('H(t)');
ylabel('g Chl/g dw')
xlabel('d')
legend('100%','25%','0%','Location','Best')
axis xy;
xlim([0 12]);

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

err100 = [norm(A1e+L1e-HL240), norm(L1e-LHL240), norm(N1e-HLN240)];
err25 = [norm(A2e+L2e-HL060), norm(L2e-LHL060), norm(N2e-HLN060)];
err0 = [norm(A3e+L3e-HL000), norm(L3e-LHL000), 0];

err = [err0; err25; err100];

    function dydt = m2(~,y)
        A = y(1);
        Q = y(2);
        N = y(3);
        H = y(4);
        L = y(5);
        Q_tilde2 = (Q*A/(A+L))^2;
        
        v = v_m*((q_m-Q)/(q_m-q))*(N/(v_h+N));
        
        p_m = p_m0*Q_tilde2/(Q_tilde2+q^2);
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
    end

end
