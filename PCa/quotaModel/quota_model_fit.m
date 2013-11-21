function [params,err] = quota_model_fit(patient, fixdata, subpopulations, test)
%QUOTA_MODEL_FIT use fminsearch to find a local minimum of MSE for quota
%model with one or two subpopulations
%   

if nargin<4
    test = false;
    if nargin<3
        subpopulations = 2;
        if nargin<2
            fixdata = false;
            if nargin<1
                error('Need to specity case # as first parameter.');
            end
        end
    end
end

if subpopulations > 2 || subpopulations < 0
  error('invalid number of subpopulations');
end

onepop = subpopulations==1;

[allparams, pIndices, allparamnames] = load_IC(patient,onepop);
[A_fit,psa_data,A_data,trange,x0] = load_data(patient,fixdata);
param0 = allparams(pIndices);

if onepop
    x0 = x0([1 2 5]);
    funhandle = @runOnePop;
    saveFile = sprintf('onepop_quota_params_%d',patient);
else
    funhandle = @run;
    saveFile = sprintf('quota_params_%d',patient);
end

if ~test
    [params,err] = fminsearch(funhandle,param0);
else
    saveFile = ['test_' saveFile]
    err = funhandle(param0);
    params = param0;
end
mse = err/size(psa_data,1);
    
allparams(pIndices) = abs(params);
% allparams = num2cell(allparams);
save(saveFile,'allparams','allparamnames','mse');

    function err = run(params)
%         p = {um, q, dm, dk, delta0, sigma, sigmak, sigma0, Rt, K, k1, Kxy, um2, q2, dm2, dk2, delta02, sigma2, sigmak2, sigma02, Rt2, K2, k2, Kyx, m, n, a, b};
        allparams(pIndices) = abs(params);
        p = num2cell(allparams);
        [Tsol,Xsol]=quota_model(trange,x0,p,A_fit,true);
        Psol = Xsol(:,5);
        err = norm(interp1(Tsol,Psol,psa_data(:,1),'linear','extrap')-psa_data(:,2))^2;
        fprintf('%d:\terr = %f\n',patient,err);
    end

    function err = runOnePop(params)
%         p = {um, q, dm, dk, delta0, sigma, sigmak, sigma0, Rt, K, m, n, a};
        allparams(pIndices) = abs(params);
        p = num2cell(allparams);
        [Tsol,Xsol]=quota_model_onepop(trange,x0,p,A_fit,true);        
        Psol = Xsol(:,3);        
        err = norm(interp1(Tsol,Psol,psa_data(:,1),'linear','extrap')-psa_data(:,2))^2;
        fprintf('%d:\terr = %f\n',patient,err);
    end
end
