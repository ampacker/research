function fitCase(patient,fixdata)
if nargin == 1
    fixdata = 0;
end

[params,err] = quota_model_fit(patient,fixdata,2,0);

try
plotPsa(patient,0,1,0,fixdata);
catch
   disp('fuck'); 
end

end
