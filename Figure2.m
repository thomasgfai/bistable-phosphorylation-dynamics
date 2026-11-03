% Figure2.m
% Generates Figure 2 in the main text.

% X. B. 2026-10-01 Initial work.

%% Initialization
Ncurves = 3;
tcc = cell(Ncurves,1);
cnoc = cell(Ncurves,1);
cnopc = cell(Ncurves,1);
cnonullxc = cell(Ncurves,1);
cnonullyc = cell(Ncurves,1);
cnopnullxc = cell(Ncurves,1);
cnopnullyc = cell(Ncurves,1);

%% Parameters
% Time
ts = 0;
% te = 1e1;
% te = 5e1;
% te = 1e2;
% te = 2e2;
% te = 5e2;
te = 1e3;
% te = 2e3;
% te = 5e3;
tspan = [ts, te];
% Nominal values of parameters
% Gene fractions
phin = 0.08; % Nuclear proteins
phic = 1-phin; % Cytoplasmic proteins
% Volume coefficients
% Acytoc = 10; % (mmol)^(-1)*L
Acytoc = 1e-2; % (μmol)^(-1)*L
% Translation rate
kt = 0.005; % s^(-1); Rough estimate
% kt = 0.004; % s^(-1); Rough estimate
% kt = 0.001; % s^(-1); Rough estimate
kn = kt*phin;
kc = kt*phic;
kvn = kn/Acytoc;
kvc = kc/Acytoc;
% Bistable phosphorylation threshold
kkc = [1 2.75 4.2]; % μM
% Nucleocytoplasmic transport rate
knt = 1e-1; % Nominal value
% Phosphorylation and dephosphorylation rate constant
% tau = 1;
tau = 5e-1;
% tau = 1e-1;
% tau = 1e-2;

%% Initial conditions
% %{
cno0 = 0; % μM
cnop0 = 0; % μM
% %}
%{
cno0 = 2; % μM
cnop0 = 5; % μM
%}
%{
cno0 = 6; % μM
cnop0 = 4; % μM
%}

% Vectors of initial conditions
xc0 = [cno0,cnop0];

%% Solve the ODE and record the max and min of solutions
% options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');

options = odeset('RelTol',1e-10,'AbsTol',1e-8);
[tc1,xc1] = ode23s(@(t,xc) ODEFUN_BistableCS2D1(t,xc,kvn,kvc,knt,kkc(1),tau,Acytoc),tspan,xc0);
[tc2,xc2] = ode23s(@(t,xc) ODEFUN_BistableCS2D1(t,xc,kvn,kvc,knt,kkc(2),tau,Acytoc),tspan,xc0);
[tc3,xc3] = ode23s(@(t,xc) ODEFUN_BistableCS2D1(t,xc,kvn,kvc,knt,kkc(3),tau,Acytoc),tspan,xc0);
tcc{1,1} = tc1;
cnoc{1,1} = xc1(:,1);
cnopc{1,1} = xc1(:,2);
tcc{2,1} = tc2;
cnoc{2,1} = xc2(:,1);
cnopc{2,1} = xc2(:,2);
tcc{3,1} = tc3;
cnoc{3,1} = xc3(:,1);
cnopc{3,1} = xc3(:,2);
clear tc1 xc1 tc2 xc2 tc3 xc3

%% Calculate nullclines 
% Plots are not really useful but cannot be turned off for now
cnomargin = 2;
cnopmargin = 2;
figure(1)
for k = 1:Ncurves
   cnomax = max(cnoc{k,1})+cnomargin;
   cnopmax = max(cnopc{k,1})+cnopmargin;
   cnonull = fimplicit(@(cno,cnop) dcnodt(cno,cnop,kvn,kvc,knt,Acytoc,kkc(k),tau),[0 cnomax 0 cnopmax],'Linewidth',1.5);
   hold on
   cnopnull = fimplicit(@(cno,cnop) dcnopdt(cno,cnop,kvn,kvc,knt,Acytoc,kkc(k),tau),[0 cnomax 0 cnopmax],'Linewidth',1.5);
   hold off
   cnonullxc{k,1} = cnonull.XData;
   cnonullyc{k,1} = cnonull.YData;
   cnopnullxc{k,1} = cnopnull.XData;
   cnopnullyc{k,1} = cnopnull.YData;
   clear cnonull cnopnull
end

%% Plot cnop/(cno+cnop) vs cno+cnop
figure(2)
k = 1;
plot(cnonullxc{k,1}+cnonullyc{k,1},cnonullyc{k,1}./(cnonullxc{k,1}+cnonullyc{k,1}),'--','Color','blue','Linewidth',2)
hold on
for k = 2:Ncurves
   plot(cnonullxc{k,1}+cnonullyc{k,1},cnonullyc{k,1}./(cnonullxc{k,1}+cnonullyc{k,1}),'--','Color','blue','Linewidth',2)
end
for k = 1:Ncurves
   plot(cnopnullxc{k,1}+cnopnullyc{k,1},cnopnullyc{k,1}./(cnopnullxc{k,1}+cnopnullyc{k,1}),':','Color','magenta','Linewidth',2)
end
plot(cnoc{1,1}+cnopc{1,1},cnopc{1,1}./(cnoc{1,1}+cnopc{1,1}),'Color','#EF9B00','Linewidth',2)
plot(cnoc{2,1}+cnopc{2,1},cnopc{2,1}./(cnoc{2,1}+cnopc{2,1}),'Color','#09961F','Linewidth',2)
plot(cnoc{3,1}+cnopc{3,1},cnopc{3,1}./(cnoc{3,1}+cnopc{3,1}),'Color','#0690B1','Linewidth',2)
hold off
axis([0 3.5 0 1])
xlabel('c_{no}+c_{nop}')
ylabel('Fraction of c_{nop}')

%% Plot individual trajectories
figure(3)
plot(tcc{1,1},cnopc{1,1}./(cnoc{1,1}+cnopc{1,1}),'Color','#EF9B00','Linewidth',2)
hold on
plot(tcc{2,1},cnopc{2,1}./(cnoc{2,1}+cnopc{2,1}),'Color','#09961F','Linewidth',2)
plot(tcc{3,1},cnopc{3,1}./(cnoc{3,1}+cnopc{3,1}),'Color','#0690B1','Linewidth',2)
hold off
xlabel('Time')
ylabel('Fraction of c_{nop}')

%% Functions
function propensity = dcnodt(cno,cnop,kvn,kvc,knt,Acytoc,kkc,tau)
   propensity = -phos(cno,cnop,kkc,tau,1)+dephos(cno,cnop,tau)+(kvn-Acytoc.*cno.*(kvn+kvc-knt.*cnop)); % Vectorized
end

function propensity = dcnopdt(cno,cnop,kvn,kvc,knt,Acytoc,kkc,tau)
   propensity = phos(cno,cnop,kkc,tau,1)-dephos(cno,cnop,tau)+(-knt.*cnop-Acytoc.*cnop.*(kvn+kvc-knt.*cnop)); % Vectorized
end

function propensity = phos(cno,cnop,K,tau,scale)
% n = 30;
% n = 10;
n = 5;
T = cno+cnop;
if T == 0
    propensity = 0;
else
    % propensity = scale*(1/tau)*T^(n+1)/(T^n+(linear_threshold(K,cnop))^n);
    % propensity = scale.*(1/tau)*T.^(n+1)./(T.^n+(linear_threshold(K,cnop)).^n); % Vectorized
    % propensity = scale*(1/tau)*T^(n+1)/(T^n+(threshold(K,cnop))^n);
    propensity = scale.*(1./tau)*T.^(n+1)./(T.^n+(threshold(K,cnop)).^n); % Vectorized
end
end

function propensity = dephos(cno,cnop,tau)
T = cno+cnop;
if T == 0
    propensity = 0;
else
    propensity = (1./tau).*cnop;
end
end

%{
function thr = linear_threshold(K,cnop)
thr = max(0.125*K, 1.25*K-4*cnop);
% That is, thr = K*max(0.125, 1.25-4*cnop/K);
% thr=K;
end
%}

function thr = threshold(K,cnop)
% Make a smooth approximation of linear_threshold (?)
% nthr = 4;
nthr = 3;
% thr = K*(5/4-(5/4-1/8)*(cnop/K)^nthr/((cnop/K)^nthr+(9/64)^nthr));
thr = K.*(5/4-(5/4-1/8).*(cnop./K).^nthr./((cnop./K).^nthr+(9/64).^nthr)); % Vectorized
% thr = K*(5/4-9/8*(cnop/K)^nthr/((cnop/K)^nthr+(9/64)^nthr));
% thr = K*(5/4-9/8*(cnop/K).^nthr/((cnop/K).^nthr+(9/64)^nthr)); % Vectorized
end
