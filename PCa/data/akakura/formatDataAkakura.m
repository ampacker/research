function formatDataAkakura
%FORMATDATAAKAKURA Take Akakura (1993) .mat data as used by Portz (2012) and 
%re-save in the format I'm using for the large 100 patient data set.

for caseNum = 1:7
   formatCaseData(caseNum); 
end
end

function formatCaseData(caseNum)
% periodDuration: column1=duration of treatment periods; col2=duration of
% off-treatment periods. values (months) are those reported in original paper
switch caseNum
    case 1
        periodDuration = [6 7; 7 7; 7 6; 7 0];
    case 2
        periodDuration = [8 5; 10 2; 5 0];
    case 3
        periodDuration = [8.5 2.5; 17 4; 5.5 0];
    case 4
        periodDuration = [13 6; 7 5; 6 0];
    case 5
        periodDuration = [7 10; 4 0];
    case 6
        periodDuration = [10 7; 4.5 0];
    case 7
        periodDuration = [17 11; 6 0];
end
periodDuration = 30*periodDuration; %approx given month values into days

cycleDuration = sum(periodDuration, 2);
periodInt = cumsum(periodDuration(:));
n = length(periodInt);
periodInt = periodInt([(1:2:n)' (2:2:n)']);

load(sprintf('data/akakura/originalFormat/case%d_data.mat',caseNum));
psa_data(:,1) = 30*psa_data(:,1);
A_data(:,1) = 30*A_data(:,1);

%% create ordered list of unique data times
T = union(psa_data(:,1), A_data(:,1));

%% add blank entries in psa and andr data sets to prepare to write to file
psaData = NaN(size(T));
andrData = NaN(size(T));
[inT,psaLocInT] = ismember(psa_data(:,1),T);
[inT,andrLocInT] = ismember(A_data(:,1),T);
psaData(psaLocInT) = psa_data(:,2);
andrData(andrLocInT) = A_data(:,2);

%% cycle # and treatment period (binary on or off treatment)
cycleList = zeros(size(T));
period = zeros(size(T));

indexOn = find(T < periodInt(1,1));
indexOff = find(T < periodInt(1,2) & T >= periodInt(1,1));
period(indexOn) = 1;
period(indexOff) = 0;
cycleList(indexOn) = 1;
cycleList(indexOff) = 1;
for j = 2:length(cycleDuration)    
    indexOn = find(T < periodInt(j,1) & T >= periodInt(j-1,2));
    indexOff = find(T < periodInt(j,2) & T >= periodInt(j,1));
    period(indexOn) = 1;
    period(indexOff) = 0;
    cycleList(indexOn) = j;
    cycleList(indexOff) = j;
end

dataTable = [caseNum*ones(length(T),1), T, psaData, andrData, cycleList, period]';
%% replace NaN's with empty cells so they don't appear in the file
nanIndex = isnan(dataTable);
dataTable = num2cell(dataTable);
dataTable(nanIndex) = cell(1);

% open a file for writing
file = sprintf('data/akakura/case%d.txt',caseNum);
fid = fopen(file, 'w+');

% print col names in quotes, followed by a blank line
cols = {'Patient','Day','PSA','Ts','Cycle','Treatment'};
fprintf(fid, '\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\r\n',cols{:});

% print values in column order
fprintf(fid, '%d,%f,%f,%f,%d,%d\r\n', dataTable{:});
fclose(fid);

end
