function [A_fit,psa_data,A_data,trangedata,x0,psa_data_full,psa_data_bad] = load_data(patient,fix)
%LOAD_DATA load and format data for specified patient
%
if nargin == 1
    fix = 0;
    psa_data_bad = zeros(0,2);
end

data_file = sprintf('data/case%d_data.mat',patient);
load(data_file);
psa_data(:,1) = 30*psa_data(:,1);
psa_data_full = psa_data;
psa_data_bad = zeros(0,2);
t_data = 30*A_data(:,1)'; A_data = A_data(:,2)';
gamma = 2.7726;
trangedata = [psa_data(1,1) max(psa_data(end,1),t_data(end))];

% fit androgen data
switch patient
    case 1
        Ai1 = 10.7; ti1 = 0; Af1 = 0.52;
        p1 = pchip([5 6 t_data(2:10)], [decay([5 6],Ai1,ti1,Af1,gamma) A_data(2:10)]);
        Ai2 = ppval(p1,t_data(10)); ti2 = t_data(10); Af2 = 0.6;
        p2 = pchip([400 405 t_data(11:22)], [decay([400 405],Ai2,ti2,Af2,gamma) A_data(11:22)]);
        Ai3 = ppval(p2,t_data(22)); ti3 = t_data(22); Af3 = 0.67;
        p3 = pchip([815 820 t_data(23:33)], [decay([815 820],Ai3,ti3,Af3,gamma) A_data(23:33)]);
        Ai4 = ppval(p3,t_data(33)); ti4 = t_data(33); Af4 = 0.57;
        p4 = pchip([1215 1220 t_data(34:39)], [decay([1215 1220],Ai4,ti4,Af4,gamma) A_data(34:39)]);
        Ai = [Ai1 Ai2 Ai3 Ai4]; ti = [ti1 ti2 ti3 ti4]; Af = [Af1 Af2 Af3 Af4];
        splines = [p1 p2 p3 p4];
        regions = [5 t_data(10) 400 t_data(22) 815 t_data(33) 1215];
        trange = [0 1400];
        
        if fix
            badindex = find(psa_data(:,2)>2 & psa_data(:,2)<2.09);
            psa_data_bad = psa_data(badindex,:);
            keep = find(psa_data(badindex-1,2)>3);
            badindex(keep) = [];
%             psa_data(badindex,2) = mean(psa_data(psa_data(:,2)<1, 2));
            psa_data(badindex,:) = [];
        end
        
        fp = 1;
    case 2
        Ai1 = 11.57; ti1 = 0; Af1 = 0.62;
        p1 = pchip([10 15 t_data(2:9)], [decay([10 15],Ai1,ti1,Af1,gamma) A_data(2:9)]);
        Ai2 = ppval(p1,t_data(9)); ti2 = t_data(9); Af2 = 0.58;
        p2 = pchip([405 410 t_data(10:16)], [decay([405 410],Ai2,ti2,Af2,gamma) A_data(10:16)]);
        Ai3 = ppval(p2,t_data(16)); ti3 = t_data(16); Af3 = 0.91;
        p3 = pchip([770 775 t_data(17:20)], [decay([770 775],Ai3,ti3,Af3,gamma) A_data(17:20)]);
        Ai = [Ai1 Ai2 Ai3]; ti = [ti1 ti2 ti3]; Af = [Af1 Af2 Af3];
        splines = [p1 p2 p3];
        regions = [10 t_data(9) 405 t_data(16) 770];
        trange = [0 900];
        
        if fix
            badindex = find(psa_data(:,2)>2.1 & psa_data(:,2)<2.21);
            psa_data_bad = psa_data(badindex,:);
            keep = find(psa_data(badindex-1,2)>3);
            badindex(keep) = [];
            psa_data(badindex,:) =  [];%mean(psa_data(psa_data(:,2)<1, 2));
        end
        
        fp = 1;
    case 3
        Ai1 = 10.48; ti1 = 0; Af1 = 0.468;
        p1 = pchip([20 25 t_data(2:9)], [decay([20 25],Ai1,ti1,Af1,gamma) A_data(2:9)]);
        Ai2 = ppval(p1,t_data(9)); ti2 = t_data(9); Af2 = 0.50;
        p2 = pchip([380 385 t_data(10:26)], [decay([380 385],Ai2,ti2,Af2,gamma) A_data(10:26)]);
        Ai3 = ppval(p2,t_data(26)); ti3 = t_data(26); Af3 = 0.487;
        p3 = pchip([1010 1015 t_data(27:30)], [decay([1010 1015],Ai3,ti3,Af3,gamma) A_data(27:30)]);
        Ai = [Ai1 Ai2 Ai3]; ti = [ti1 ti2 ti3]; Af = [Af1 Af2 Af3];
        splines = [p1 p2 p3];
        regions = [20 t_data(9) 380 t_data(26) 1010];
        trange = [0 1200];
        
        %no psa readings below detection?
        badindex = [];
        
        fp = 1;
    case 4
        Ai1 = 7.3; ti1 = 0; Af1 = 0.54;
        p1 = pchip([20 25 t_data(2:14)], [decay([20 25],Ai1,ti1,Af1,gamma) A_data(2:14)]);
        Ai2 = ppval(p1,t_data(14)); ti2 = t_data(14); Af2 = 0.65;
        p2 = pchip([620 625 t_data(15:23)], [decay([620 625],Ai2,ti2,Af2,gamma) A_data(15:23)]);
        Ai3 = ppval(p2,t_data(23)); ti3 = t_data(23); Af3 = 0.47;
        p3 = pchip([980 985 t_data(24:27)], [decay([980 985],Ai3,ti3,Af3,gamma) A_data(24:27)]);
        Ai = [Ai1 Ai2 Ai3]; ti = [ti1 ti2 ti3]; Af = [Af1 Af2 Af3];
        splines = [p1 p2 p3];
        regions = [20 t_data(14) 620 t_data(23) 980];
        trange = [0 1150];
        
        badindex = [];
        
        fp = 0.7;
    case 5
        Ai1 = 20; ti1 = 0; Af1 = 1.4;
        p1 = pchip([20 25 t_data(2:11)], [decay([20 25],Ai1,ti1,Af1,gamma) A_data(2:11)]);
        Ai2 = ppval(p1,t_data(11)); ti2 = t_data(11); Af2 = 1.8;
        p2 = pchip([530 535 t_data(12:14)], [decay([530 535],Ai2,ti2,Af2,gamma) A_data(12:14)]);
        Ai = [Ai1 Ai2]; ti = [ti1 ti2]; Af = [Af1 Af2];
        splines = [p1 p2];
        regions = [20 t_data(11) 530];
        trange = [0 650];
        
        badindex = [];

        fp = 0.8;
    case 6
        Ai1 = 14; ti1 = 0; Af1 = 0.54;
        p1 = pchip([20 25 t_data(2:16)], [decay([20 25],Ai1,ti1,Af1,gamma) A_data(2:16)]);
        Ai2 = ppval(p1,t_data(16)); ti2 = t_data(16); Af2 = 0.55;
        p2 = pchip([535 540 t_data(17:21)], [decay([535 540],Ai2,ti2,Af2,gamma) A_data(17:21)]);
        Ai = [Ai1 Ai2]; ti = [ti1 ti2]; Af = [Af1 Af2];
        splines = [p1 p2];
        regions = [20 t_data(16) 535];
        trange = [0 650];
        
        badindex = [];
        
        fp = 0.6;
    case 7
        Ai1 = 12.1; ti1 = 0; Af1 = 0.51;
        p1 = pchip([20 25 t_data(2:20)], [decay([20 25],Ai1,ti1,Af1,gamma) A_data(2:20)]);
        Ai2 = ppval(p1,t_data(20)); ti2 = t_data(20); Af2 = 0.50;
        p2 = pchip([890 895 t_data(21:25)], [decay([890 895],Ai2,ti2,Af2,gamma) A_data(21:25)]);
        Ai = [Ai1 Ai2]; ti = [ti1 ti2]; Af = [Af1 Af2];
        splines = [p1 p2];
        regions = [20 t_data(20) 890];
        trange = [0 1025];
        
        if fix
            badindex = find(psa_data(:,2)>2.0295 & psa_data(:,2)<2.055);
            psa_data_bad = psa_data(badindex,:);
            keep = find(psa_data(badindex-1,2)>3);
            badindex(keep) = [];
%             psa_data(badindex,2) = mean(psa_data(psa_data(:,2)<1, 2));
            psa_data(badindex,:) = [];
        end
        
        fp = 1.1;
end
A_fit = {gamma Ai ti Af splines regions};
A_data = [t_data', A_data'];
xy = 14.9/15;
psa = interp1(psa_data_full(:,1),psa_data_full(:,2),trangedata(1),'cubic','extrap');
Xi = xy*psa*fp;
Qi = 35;%1.2*q;

Xi2 = (1-xy)*Xi;%(1 - xy)*psa*fp;
Qi2 = 35;%1.2*q2;

Pi = psa;
x0 = [Xi; Qi; Xi2; Qi2; Pi];

end

%% exponential decay function for androgen data  fitting
function A = decay(t,Ai,ti,Af,gamma)
A = Af + (Ai - Af)*exp(-gamma*(t-ti));
end
