function plotPsa(patient, plotAndr, plotQ, testresult, fixdata, subpopulations)
%PLOTPSA plot PSA data vs simulation results
%   
if nargin<6
    subpopulations = 2;
    if nargin<5
        fixdata = 0;
        if nargin<4
            testresult = 0;
            if nargin<3
               plotQ = 1;
               if nargin<2
                  plotAndr = 1; 
               end
            end
        end
    end
end

onepop = subpopulations==1;

[A_fit,psa_data,A_data,trange,x0,psa_data_full,psa_data_bad] = load_data(patient,fixdata);
if onepop
    filename = sprintf('onepop_quota_params_%d',patient);
else
    filename = sprintf('quota_params_%d',patient);
end
if testresult
    filename = ['test_' filename];
end

load(filename, 'allparams','allparamnames','mse');
fprintf('Case %d MSE=%f\n',patient,mse);
for j = 1:length(allparamnames)
   if strcmp(allparamnames(j),'q')
       q = allparams(j);
   elseif strcmp(allparamnames(j),'q2')
       q2 = allparams(j);
   elseif strcmp(allparamnames(j),'Rt')
       Rt = allparams(j);
   elseif strcmp(allparamnames(j),'Rt2')
       Rt2 = allparams(j);
    elseif strcmp(allparamnames(j),'dk')
       dk = allparams(j);
   elseif strcmp(allparamnames(j),'dk2')
       dk2 = allparams(j);
    elseif strcmp(allparamnames(j),'sigmak')
       sigmak = allparams(j);
   elseif strcmp(allparamnames(j),'sigmak2')
       sigmak2 = allparams(j);
    elseif strcmp(allparamnames(j),'K')
       K = allparams(j);
    elseif strcmp(allparamnames(j),'K2')
       K2 = allparams(j);
   end
end
          
if onepop
    fhandle = @quota_model_onepop;
    x0 = x0([1 2 5]);
else
    fhandle = @quota_model;
end

[T,x]=fhandle(trange,x0,num2cell(allparams),A_fit,true);
if onepop
    X = x(:,1);
    Q = x(:,2);
    P = x(:,3);
else
    X = x(:,1);
    Q = x(:,2);
    X2 = x(:,3);
    Q2 = x(:,4);
    P = x(:,5);
end

mse = norm(interp1(T,P,psa_data(:,1),'linear','extrap')-psa_data(:,2))^2/length(psa_data(:,1));
fprintf('%d:\tMSE = %f\n',patient,mse);
rse = norm((interp1(T,P,psa_data(:,1),'linear','extrap')-psa_data(:,2))./psa_data(:,2))^2/length(psa_data(:,1));
fprintf('%d:\tRSE = %f\n',patient,rse);

titlestr = sprintf('Patient %d %s',patient);

figure(10*patient)
clf
plot(psa_data_full(:,1),psa_data_full(:,2),'-.sr','MarkerEdgeColor','r');
% plot(psa_data(:,1),psa_data(:,2),'-.sr','MarkerEdgeColor','r');
hold on
plot(T,P,'-k');
plot(T,X,'b','LineWidth',1.0)
if ~onepop
    plot(T,X2,'--g','LineWidth',1.0)
end

plot(psa_data_bad(:,1),psa_data_bad(:,2),'xr','MarkerEdgeColor','k');
hold off
% axis tight;
xlim(trange);
if ~onepop
    hl = legend('PSA data','PSA model','CS cells','CR cells','Location','NorthWest');
else
    hl = legend('PSA data','PSA model','cells','Location','NorthWest');
end
% hl = columnlegend(2, {'PSA data','PSA model','CS cells','CR cells'},'boxon')
% set(hl, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
xlabel('days')
yunits = sprintf('Serum PSA (ng/mL) and Cells (millions)');
ylabel(yunits);
% title(titlestr);

if plotQ
    figure(10*patient+1)
    clf
    trange = T([1 end]);
    plot(T,Q,'-k',trange,[Rt; Rt],'--b',trange,[dk; dk],'-.r',trange,[q; q],'--b');%,trange,[sigmak, sigmak],':');
    legend('Q','R_t, q','d_k')
        hold on
        if ~onepop
            plot(T,Q2,':k',trange,[Rt2; Rt2],':b',trange,[dk2; dk2],':r',trange,[q2; q2],':b');%,trange,[sigmak, sigmak],':');
        end
    hold off
%     axis tight;
%     ylim([0 max(Rt,Rt2)]);
%     ylim([0 Rt]);
    xlim(trange);
%     hl = legend('Q','R_t','q','Q_2','{R_t}_2','q_2');
    hl = legend('Q','R_t, q','d_k');%,'\sigma_k');
%     set(hl, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'Color', 'w','Orientation','Horizontal');
    xlabel('days')
    yunits = sprintf('nM');
    ylabel(yunits);
%     title(titlestr);
end

end
