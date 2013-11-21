function quotaModelTest(patient, fixdata, subpopulations)

if nargin<3
    subpopulations = 2;
    if nargin<2
        fixdata = false;
        if nargin<1
            error('Need to specify case # as first parameter');
        end
    end
end

quota_model_fit(patient,fixdata,subpopulations,1);
plotPsa(patient,0,1,1,fixdata,subpopulations);
