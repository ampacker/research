function [tsol,xsol] = quota_model_onepop(trange,x0,p,A_fit,calcQ0)
%quota_model_onepop solves AR quota model with one subpopulation using
%ode15s
%  

[um, q, dm, dk, delta0, sigma, sigmak, sigma0, Rt, K, m, n, a] = deal(p{:});

qm = q^m;
dkn = dk^n;

%% PSA vars
sigmaka = sigmak^a;       % nM^{m}
betaP = log(2)/0.12/24; %0.08;       % 1/day

%% AR kinetic vars
kaT = 0.14*24;         % 1/nM/hr
kdT = 0.069*24;        % 1/hr
kaD = 0.053*24;        % 1/nM/hr
kdD = 0.018*24;        % 1/hr
ka = 0.8*kaD+0.2*kaT;       % 1/nM/hr;
kd = 0.8*kdD+0.2*kdT;       % 1/hr;

betaR = log(2)/3*24;   % 1/hr
betaT = log(2)/3*24;   % 1/hr
betaD = log(2)/9*24;   % 1/hr
betaA = 0.8*betaD+0.2*betaT;    % 1/hr (weighted sum of betaT and betaD assuming steady state proportion of free DHT and T relative to total free DHT+T.

alpha = 5;        % mg/L
kcat = 8.276*24;          % nmol/hr/mg (+/- 15)
dhtM = alpha*kcat;
dhtK = 70;       % from Eikenberry (Jain11: 0.104?!)

h = K+betaT;    % for ``uptake'' function

if calcQ0
    dhtM = 2*kcat;
    Ts = interp_fit(A_fit,0);
    Ti = 0.5*((K*Ts-dhtM)./(K+betaT)-dhtK+sqrt( ((K*Ts-dhtM)./(K+betaT)-dhtK).^2 + 4*dhtK*K*Ts./(K+betaT) ));
    Di = dhtM/betaD*Ti/(Ti+dhtK);
    theta = kaT/kdT*Ti+kaD/kdD*Di;
    Qi = theta/(1+theta)*Rt;
    x0(2) = Qi;
end
dhtM = alpha*kcat;
[tsol,xsol]=ode15s(@(t,x)f(t,x),trange,x0);

%% model ODE function
    function dxdt = f(t,x)
        X = x(1); Q = x(2);  P = x(3);
        Ts = interp_fit(A_fit,t);
        
%         %Droop model
%         mu = um*(1 - q/Q);
        
        mu = um*Q^m/(Q^m+qm);
        D = dm*dkn/(Q^n+dkn)+delta0;

        %First "order" approximations
%         v = (kaT*K*Ts/h+kaD*dhtM/betaD*K*Ts/(K*Ts+dhtK*h))*(Rt-Q);
        
        %2nd "order" approximations
        v = (kaT*K*Ts/h*(K*Ts+dhtK*h)/(K*Ts+dhtK*h+dhtM)+kaD/betaD*dhtM*K*Ts*(K*Ts+dhtK*h)/((K*Ts+dhtK*h)^2+dhtK*h*dhtM))*(Rt-Q);
        
        dX = (mu - D)*X;
        dQ = v - mu*Q - betaA*Q;
        dP = sigma*X*Q^a/(Q^a + sigmaka) + sigma0*X - betaP*P;
        
        dxdt = [dX; dQ; dP];
    end
end
