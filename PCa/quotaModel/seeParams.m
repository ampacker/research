function seeParams(cases, subpopulations)
if nargin==0
    cases = 1:7;
elseif nargin<2
   subpopulations = 2; 
end
load(getfile(cases(1),subpopulations==1),'allparams','allparamnames','mse');
dispFormat = '%8s:\t%8.3g';
if subpopulations==1
    dispCell = [[allparamnames 'mse']; num2cell([allparams mse])];
else
    sortindex = zeros(1,length(allparamnames));
    sortindex(1:2:end-4) = 1:(length(allparamnames)-4)/2;
    sortindex(2:2:end-4) = (length(allparamnames)-4)/2+1:length(allparamnames)-4;
    sortindex(end-3:end) = length(allparamnames)-3:length(allparamnames);
    dispCell = [[allparamnames(sortindex) 'mse']; num2cell([allparams(sortindex) mse])];
end

for j = 2:length(cases)
    load(getfile(cases(j)), 'allparams','allparamnames','mse');
    dispCell(j+1,:) = num2cell([allparams(sortindex) mse]);
    dispFormat = [dispFormat '\t%8.3g'];   
end
dispFormat = [dispFormat '\n'];
cases = num2cell(cases);
fprintf('%8s ',' ');
fprintf('\t    %d    ',cases{:});
fprintf('\n');
fprintf(dispFormat,dispCell{:});
end

function filename = getfile(patient,onepop)
if ~onepop
    filename = sprintf('quota_params_%d',patient);
else
    filename = sprintf('onepop_quota_params_%d',patient);
end
end
