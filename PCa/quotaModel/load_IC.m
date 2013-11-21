function [params, pIndices, paramNames] = load_IC(patient, onepop)
%LOAD_IC return initial guesses for specified patient
%  

% fixed     MSE 1.81 3.18  x.xx 3.38   0.62 0.98 1.62
% AR quota  MSE 1.74 2.77  x.xx 3.38   0.62 0.98 1.92
% CDR model MSE 2.36 2.02 45.05 3.74  11.51 0.38 3.20
% PKN model MSE 2.96 2.79 46.42 3.924 41.81 0.40 3.71
if nargin<2
    onepop = false;
end

% model parameters
paramNames2 = {'um', 'q', 'dm', 'dk', 'delta0', 'sigma', 'sigmak', 'sigma0', 'Rt', 'K', 'k1', 'Kxy', 'um2', 'q2', 'dm2', 'dk2', 'delta02', 'sigma2', 'sigmak2', 'sigma02', 'Rt2', 'K2', 'k2', 'Kyx', 'm', 'n', 'a', 'b'};
paramNames1 = {'um', 'q', 'dm', 'dk', 'delta0', 'sigma', 'sigmak', 'sigma0', 'Rt', 'K', 'm', 'n', 'a'}; %for one population model
if onepop
    paramNames = paramNames1; %for one population model
else
    paramNames = paramNames2; 
end

% selected free parameters (can change for specific patients)
% pNames = {'um', 'q', 'dm','sigma', 'sigmak', 'k1', 'um2', 'dm2', 'q2','sigma2','sigmak2','Rt2'};
pNames = {'um', 'q', 'dm','delta0','sigma', 'sigmak'};
pNames = {'um', 'q', 'dm','sigma','um2','q2','dm2','sigma2','Rt2','k1'};
% pNames = {'um','dm','sigma','um2','dm2','sigma2','k1'};
% pNames = {'um2','q2','dm2','dk2','sigma2','sigmak2','Rt2','k1'};

%"Fixed" parameter values
m = 4;
n = 4;
a = 4;
b = 2;

switch patient
    case 1.1
        um = 0.013;
        q = 81;

        um2 = 0.0094;
        q2 = q;
        
        dm = 0.019;
        dk = 85.35;
        delta0 = 0.000127;
        
        dm2 = 0.019;
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 180;
        Rt2 = 220;
        
        sigma = 0.36;
        sigmak = 143;
        sigma0 = 0.00142;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1042*24;         % 1/hr
        K2 = 1*K;

        k1 = 0.00036;
        Kxy = 48;
        k2 = 0.000212;
        Kyx = 150;
        
    case 1   
         um = 0.013;
        q = 81;

        um2 = 0.0094;
        q2 = q;
        
        dm = 0.019;
        dk = 85.35;
        delta0 = 0.000127;
        
        dm2 = 0.019;
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 180;
        Rt2 = 220;
        
        sigma = 0.36;
        sigmak = 143;
        sigma0 = 0.00142;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1042*24;         % 1/hr
        K2 = 1*K;

        k1 = 0.00036;
        Kxy = 48;
        k2 = 0.000212;
        Kyx = 150;

    case 3
%         um = 0.0416;%0.013;
%         q =  71.4;%81;
% 
%         um2 = 0.0108;%0.0094;
%         q2 = 69;
%         
%         dm = 0.0136;%0.019;
%         dk = 85.35;
%         delta0 = 0.0011;%0.00127;
%         
%         dm2 = 0.00868;%0.019;
%         dk2 = 87.5;%dk;
%         delta02 = delta0;
%         
%         Rt = 150;
%         Rt2 = 200;
%         
%         sigma =0.443;%0.36;
%         sigmak =  80;
%         sigma0 = 0.00142;
%         
%         sigma2 = sigma;
%         sigmak2 = sigmak;
%         sigma02 = sigma0;
                
        um = 0.025;%0.0416;%0.013;
        q =  75.4;%71.4;%81;

        um2 = 0.012;
        q2 = q;
        
        dm = 0.008;%0.0136;%0.019;
        dk = 85.35;
        delta0 = 0.00113;%0.0011;%0.00127;
        
        dm2 = 0.010;
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 150;
        Rt2 = 200;
        
        sigma = 0.342;%0.443;%0.36;
        sigmak = 70;% 39.7;%59;
        sigma0 = 0.00142;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1242*24;         % 1/hr
        K2 = 1*K;

        k1 = 0.00051;
        Kxy = 40;
        k2 = 0.000112;
        Kyx = 150;

    case 2
        um = 0.019;
        q = 100;

        um2 = 0.016;
        q2 = q;
        
        dm = 0.019; %0.15;
        dk = 100.35;
        delta0 = 0.000127;
        
        dm2 = 0.021; %0.010;
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 180;
        Rt2 = 230;
        
        sigma = 0.30;
        sigmak = 130;
        sigma0 = 0.00142;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1042*24;         % 1/hr
        K2 = 1*K;

        k1 = 0.00036;
        Kxy = 48;
        k2 = 0.000212;
        Kyx = 150;
          
    case 3
        um = 0.019;
        q = 140;

        um2 = 0.009;
        q2 = q;
        
        dm = 0.008;
        dk = 50.35;
        delta0 = 0.000127;
        
        dm2 = 0.013;
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 190;
        Rt2 = 210;
        
        sigma = 0.33;
        sigmak = 130;
        sigma0 = 0.00142;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1542*24;         % 1/hr
        K2 = 1*K;

        k1 = 0.00036;
        Kxy = 48;
        k2 = 0.000212;
        Kyx = 150;
        
    case 4
        um = 0.005;
        q = 150;

        um2 = 0.011;
        q2 = q;
        
        dm = 0.015; %0.15;
        dk = 40.35;
        delta0 = 0.000127;
        
        dm2 = 0.017; %0.010;
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 180;
        Rt2 = 195;
        
        sigma = 0.2;
        sigmak = 60;
        sigma0 = 0.00142;%0.001;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1042*24;         % 1/hr
        K2 = 1*K;

        k1 = 0.00015;
        Kxy = 48;
        k2 = 0.000112;
        Kyx = 180;
        
    case 5  %not what was used for the actual case 5 parameters used in paper. I lost those initial guesses      
        um = 0.014;
        q = 90;

        um2 = 0.011;
        q2 = q;
        
        dm = 0.014; %0.15;
        dk = 85.35;
        delta0 = 0.000127;
        
        dm2 = 0.017; %0.010;
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 180;
        Rt2 = 195;
        
        sigma = 0.30;
        sigmak = 60;
        sigma0 = 0.00142;%0.001;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1042*24;         % 1/hr
        K2 = 1*K;

        k1 = 0.00015;
        Kxy = 48;
        k2 = 0.000112;
        Kyx = 180;
        
    case 6
        um = 0.013;
        q = 81;

        um2 = 0.0094;
        q2 = q;
        
        dm = 0.019;
        dk = 85.35;
        delta0 = 0.000127;
        
        dm2 = 0.019; 
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 180;
        Rt2 = 220;
        
        sigma = 0.36;
        sigmak = 143;
        sigma0 = 0.00142;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1042*24;         % 1/hr
        K2 = 1*K;

        k1 = 0.00036;
        Kxy = 48;
        k2 = 0.000212;
        Kyx = 150;
        
    case 7
        pNames = {'um','um2','dm','dm2','sigma','sigma2','sigmak','sigmak2'};
        um = 0.86*0.0133;
        q = 7.06;

        um2 = 0.94*0.0121;
        q2 = 7.06;
        
        dm =  0.0192;
        dk = 85.35;
        delta0 = 8.79e-05;
        
        dm2 = 0.78*0.0242;
        dk2 = dk;
        delta02 = delta0;
        
        Rt = 180;
        Rt2 = 183;
        
        sigma = 0.797/2;
        sigmak = 233/2;
        sigma0 = 0.00142;
        
        sigma2 = sigma;
        sigmak2 = sigmak;
        sigma02 = sigma0;
        
        K = 0.1042*24;
        K2 = 1*K;

        k1 = 0.000265*3;
        Kxy = 48;
        k2 = 0.000212;
        Kyx = 150;
end

% params = [um, q, dm, dk, delta0, sigma, sigmak, sigma0, Rt, K, k1, Kxy, um2, q2, dm2, dk2, delta02, sigma2, sigmak2, sigma02, Rt2, K2, k2, Kyx, m, n, a, b]
% params = [um, q, dm, dk, delta0, sigma, sigmak, sigma0, Rt, K, m, n, a]
params = eval(['[' strjoin(paramNames,',') ']'])
pIndices = getIndices(pNames);
paramNames(pIndices)

    %%
    function pIndices = getIndices(pNames)
        N = length(pNames);
        pIndices = zeros(1,N);
        for k=1:N
            pindex = find(strcmp(pNames{k},paramNames));
            if ~isempty(pindex)
                pIndices(k) = pindex;                
            end
        end
        if any(pIndices==0)
            invalid = pIndices==0;
            % don't throw error if param is from two subpopulation model.
            if onepop
                %logical indices of paramNames where paramNames~=paramNames1
                check = [~strcmp(paramNames2(1:length(paramNames1)),paramNames1), true(1,length(paramNames2)-length(paramNames1))];
                for k = find(pIndices==0)
                    if any(strcmp(pNames{k},paramNames2(check)))
                        invalid(k) = false; %it's a valid two pop. model parameter
                    end
                end
            end
            if any(invalid)
                errstr = ['Invalid parameter names: ' strjoin(pNames(invalid),',')];
                error(errstr);
            else
                pIndices = pIndices(pIndices>0);
            end
        end
        pIndices = sort(pIndices);  %sort because I'm OCD
    end
end
