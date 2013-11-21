function [tsol,xsol] = quota_model(trange,x0,p,A_fit,calcQ0)
%QUOTA_MODEL solves quota model with two subpopulations using ODE15s
%   

[um, q, dm, dk, delta0, sigma, sigmak, sigma0, Rt, K, k1, Kxy, um2, q2, dm2, dk2, delta02, sigma2, sigmak2, sigma02, Rt2, K2, k2, Kyx, m, n, a, b] = deal(p{:});

qm = q^m;
q2m = q2^m;

dkn = dk^n;
dk2n = dk2^n;

Kxyb = Kxy^b;
Kyxb = Kyx^b;
% k1 = k1*3;

%% PSA vars
sigmaka = sigmak^a;       % nM^{m}
betaP = log(2)/0.12/24; %0.08;       % 1/day
sigmak2a = sigmak2^a;       % nM^{m}

%% AR kinetic vars
%Rt = 50;
kaT = 0.14*24;         % 1/nM/hr
kdT = 0.069*24;        % 1/hr
kaD = 0.053*24;        % 1/nM/hr
kdD = 0.018*24;        % 1/hr
ka = 0.8*kaD+0.2*kaT;       % 1/nM/hr;
kd = 0.8*kdD+0.2*kdT;       % 1/hr;

kaT2 = kaT;         % 1/nM/hr
kdT2 = kdT;        % 1/hr
kaD2 = kaD;       % 1/nM/hr
kdD2 = kdD;        % 1/hr
ka2 = 0.8*kaD2+0.2*kaT2;       % 1/nM/hr;
kd2 = 0.8*kdD2+0.2*kdT2;       % 1/hr;

betaR = log(2)/3*24;   % 1/hr
betaT = log(2)/3*24;   % 1/hr
betaD = log(2)/9*24;   % 1/hr
betaA = 0.8*betaD+0.2*betaT;    % 1/hr (weighted sum of betaT and betaD assuming steady state proportion of free DHT and T relative to total free DHT+T.

betaR2 = betaR;   % 1/hr
betaT2 = betaT;   % 1/hr
betaD2 = betaD;   % 1/hr
betaA2 = (0.8*betaD2+0.2*betaT2)/2;    % 1/hr (weighted sum of betaT and betaD assuming steady state proportion of free DHT and T relative to total free DHT+T.

alpha = 5;        % mg/L
kcat = 8.276*24;          % nmol/hr/mg (+/- 15)
dhtM = alpha*kcat;
dhtK = 70;       % from Eikenberry (Jain11: 0.104?!)
dhtM2 = dhtM;
dhtK2 = dhtK;

if calcQ0
    dhtM=2*8.276*24;
    Ts = interp_fit(A_fit,0);
    Ti = 0.5*((K*Ts-dhtM)./(K+betaT)-dhtK+sqrt( ((K*Ts-dhtM)./(K+betaT)-dhtK).^2 + 4*dhtK*K*Ts./(K+betaT) ));
    Di = dhtM/betaD/2*Ti/(Ti+dhtK);
    theta = kaT/kdT*Ti+kaD/kdD*Di;
    Qi = theta/(1+theta)*Rt;
    
    Ti2 = 0.5*((K2*Ts-dhtM2)./(K2+betaT2)-dhtK2+sqrt( ((K2*Ts-dhtM2)./(K2+betaT2)-dhtK2).^2 + 4*dhtK2*K2*Ts./(K2+betaT2) ));
    Di2 = dhtM2/betaD2/2*Ti2/(Ti2+dhtK2);
    theta2 = kaT2/kdT2*Ti2+kaD2/kdD2*Di2;
    Qi2 = theta2/(1+theta2)*Rt2;
%     
%     fx2 = 1-14.95/15;
%     fp = betaP/(sigma*Qi.^a./(Qi.^a + sigmak^a)+sigma0);
%     Xi = (1-fp*psa0);    
    
    %x0 = [Xi; Qi; Xi2; Qi2; Pi];
    x0([2 4]) = [Qi Qi2];
end
dhtM = alpha*kcat;

[tsol,xsol]=ode15s(@(t,x)f(t,x),trange,x0);
% X = xsol(:,1);
% Q = xsol(:,2);
% X2 = xsol(:,3);
% Q2 = xsol(:,4);
% P = xsol(:,5);
% 
% h = K+betaT;
% h2 = K2+betaT2;
% Ts = interp_fit(A_fit,tsol)';
% v = (kaT*K*Ts/h+kaD*dhtM/betaD*K*Ts./(K*Ts+dhtK*h));%.*(Rt-Q);
% v2 = (kaT2*K2*Ts/h2+kaD2/betaD2*dhtM2*K2*Ts./(K2*Ts+dhtK2*h2));%.*(Rt2-Q2);
% 
% vb = (kaT*K*Ts/h.*(K*Ts+dhtK*h)./(K*Ts+dhtK*h+dhtM)+kaD/betaD*dhtM*K*Ts.*(K*Ts+dhtK*h)./((K*Ts+dhtK*h).^2+dhtK*h*dhtM));%.*(Rt-Q);
% vb2 =(kaT2*K2*Ts/h2.*(K2*Ts+dhtK2*h2)./(K2*Ts+dhtK2*h2+dhtM2)+kaD2/betaD2*dhtM2*K2*Ts.*(K2*Ts+dhtK2*h2)./((K2*Ts+dhtK2*h2).^2+dhtK2*h2*dhtM2));%.*(Rt2-Q2);

% figure(10865)
% clf;
% plot(tsol,v,':',tsol,v2,':',tsol,vb,'-.',tsol,vb2,'-.');
% legend('v','v_2','v_b','v_b2');
% % ylim([0 100]);


%% model ODE function
    function dxdt = f(t,x)
        X = x(1); Q = x(2); X2 = x(3); Q2 = x(4); P = x(5);
        Ts = interp_fit(A_fit,t);
        
        %Droop model
%         mu = um*(1 - q/Q);
%         mu2 = um2*(1 - q2/Q2);
        %
        mu = um*Q^m/(Q^m+qm);
        mu2 = um2*Q2^m/(Q2^m+q2m);
        
        D = dm*dkn/(Q^n+dkn)+delta0;
        D2 = dm2*dk2n/(Q2^n+dk2n)+delta02;
        
        mxy = k1*Kxyb/(Q^b + Kxyb);
        myx = k2*Q2^b/(Q2^b + Kyxb);
        
        h = K+betaT;
        h2 = K2+betaT2;
        %First "order" approximations
%         v = (kaT*K*Ts/h+kaD*dhtM/betaD*K*Ts/(K*Ts+dhtK*h))*(Rt-Q);
%         v2 = (kaT2*K2*Ts/h2+kaD2/betaD2*dhtM2*K2*Ts/(K2*Ts+dhtK2*h2))*(Rt2-Q2);
        
        %2nd "order" approximations
        v = (kaT*K*Ts/h*(K*Ts+dhtK*h)/(K*Ts+dhtK*h+dhtM)+kaD/betaD*dhtM*K*Ts*(K*Ts+dhtK*h)/((K*Ts+dhtK*h)^2+dhtK*h*dhtM))*(Rt-Q);
        v2 =(kaT2*K2*Ts/h2*(K2*Ts+dhtK2*h2)/(K2*Ts+dhtK2*h2+dhtM2)+kaD2/betaD2*dhtM2*K2*Ts*(K2*Ts+dhtK2*h2)/((K2*Ts+dhtK2*h2)^2+dhtK2*h2*dhtM2))*(Rt2-Q2);
        
        dX = (mu - D - mxy)*X + myx*X2;
        dX2 = (mu2 - D2 - myx)*X2 + mxy*X;
        dQ = v - mu*Q - betaA*Q;
        dQ2 = v2 - mu2*Q2 - betaA2*Q2;
        dP = sigma*X*Q^a/(Q^a + sigmaka) + sigma2*X2*Q2^a/(Q2^a + sigmak2a) + sigma0*(X + X2) - betaP*P;
        
        dxdt = [dX; dQ; dX2; dQ2; dP];
    end
end
