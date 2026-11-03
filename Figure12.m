% Calculate the coefficient of variation and mean period of the stochastic solutions.

% J. T. 2025-12-23
% X. B. 2025-06-17 Added data export.
% X. B. 2025-06-16 Added various plots.
% X. B. 2025-06-16 Bistability of $c_{nop}/c_{no}+c_{nop}$ may be the key to
% explaining the (qualitative) mechanism behind bifurcation diagrams.
% X. B. 2025-01-13 Initial work.

%% Parameters

PLOTS = 0;
SAVE = 1;

% Time
ts = 0;
% te = 1e1;
% te = 5e1;
% te = 1e2;
% te = 2e2;
% te = 5e2;
% te = 2e3;
% te = 5e3;
te = 5e5;
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
% Nucleocytoplasmic transport rate
% knt = 1;
% knt = 3e-1;
knt = 1e-1; % Nominal value
% knt = 5e-2;
% knt = 2e-2;
% knt = 1e-2;
% knt = 5e-3;
% knt = 3.97e-3;
% Phosphorylation and dephosphorylation rate constant
taumin = 0.5;
% taumin = 1;
% taumin = 10;
taumax = 60;
tau_list = [taumin, taumax];
ntau = length(tau_list);
% Bistable phosphorylation threshold
% ns = 5;
% ns = 10;
% ns = 50;
ns = 60; % It took about 9800s or 2.72h to run.
Means = zeros(ntau,ns);
Stds = zeros(ntau,ns);
npks = zeros(ntau,ns);
kkcmin = 0.5;
% kkcmin = 1;
% kkcmin = 1.5;
% kkcmax = 4.5;
kkcmax = 5;
KKCCs = linspace(kkcmin,kkcmax,ns);
Sig=0.01;
ktau=0;
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
for tau=tau_list
ktau=ktau+1;

% options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
options = odeset('RelTol',1e-10,'AbsTol',1e-8);

% Stochastic
Multiplicative=0;
ProtectNegative=1;
dt=0.01;
time_vect=ts:dt:te;
ntimes=length(time_vect);
X_stoch=zeros(2,ntimes);
X_stoch(:,1)=xc0;
for ks=1:ns
% Sig=Sigmas(ks)*[1;1];
kkc=KKCCs(ks);
[tc,xc] = ode23s(@(t,xc) ODEFUN_BistableCS2D1(t,xc,kvn,kvc,knt,kkc,tau,Acytoc),tspan,xc0);
if (max(xc(ceil(end/2:end),1))-min(xc(ceil(end/2:end),1)))>0.1
for k=1:ntimes-1
    t=time_vect(k);
    xs=X_stoch(:,k);
    if Multiplicative ~=0
        Noise=Sig.*sqrt(abs(Noise_FUN_BistableCS2D1(t,xs,kvn,kvc,knt,kkc,tau,Acytoc))).*randn(2,1);
    else
        Noise=Sig.*randn(2,1);
    end
    X_stoch(:,k+1)=xs+dt*ODEFUN_BistableCS2D1(t,xs,kvn,kvc,knt,kkc,tau,Acytoc)+...
        sqrt(dt)*Noise;

    if ProtectNegative
        if X_stoch(1,k+1)<0
            X_stoch(1,k+1)=X_stoch(1,k+1)-2*sqrt(dt)*Noise(1);
        end
        if X_stoch(2,k+1)<0
            X_stoch(2,k+1)=X_stoch(2,k+1)-2*sqrt(dt)*Noise(2);
        end
    end

end
if PLOTS ~=0
figure;plot(tc,xc(:,1))
hold on
plot(time_vect,X_stoch(1,:))
figure;findpeaks(-smooth(X_stoch(1,:),5),'MinPeakProminence',0.2*(max(smooth(X_stoch(1,ceil(end/2):end),5))-min(smooth(X_stoch(1,ceil(end/2):end),5))));
end


if max(smooth(X_stoch(1,ceil(end/2):end),5))-min(smooth(X_stoch(1,ceil(end/2):end),5))>0.3
    [pks,locs]=findpeaks(-smooth(X_stoch(1,:),5),'MinPeakProminence',0.2*(max(smooth(X_stoch(1,ceil(end/2):end),5))-min(smooth(X_stoch(1,ceil(end/2):end),5))));
% figure;
% findpeaks(-smooth(X_stoch(1,:),5),'MinPeakProminence',0.6*(max(smooth(X_stoch(1,:),5))-min(smooth(X_stoch(1,:),5))));
% pause()
Means(ktau,ks)=mean(diff(time_vect(locs)));
Stds(ktau,ks)=std(diff(time_vect(locs)));
npks(ktau,ks)=length(locs);
end
end
xc0=xc(end,:);
end
end
%%
close all;
for i=1:ntau
    figure(100);
    hold on
    errorbar(KKCCs, Means(i,:),Stds(i,:)./sqrt(npks(i,:)));
    figure(101)
    hold on
    plot(KKCCs,Stds(i,:)./Means(i,:),'o-');
end

%% Save data
if SAVE ~=0
save(['BistableC2D2_CV_kkc=' num2str(kkcmin) '_to_' num2str(kkcmax) '_tau=' num2str(taumin) '_to_' num2str(taumax) '_sigma=' num2str(Sig) '.mat'])
end

%% Functions
function propensity = phos(cno,cnop,K,tau,scale)
% n = 30;
% n = 10;
n = 5;
T = cno+cnop;
if T == 0
    propensity = 0;
else
    propensity = scale*(1/tau)*T.^(n+1)./(T.^n+(threshold(K,cnop)).^n);
end
end

function pswitch = phos_switch(cno,cnop,K)
% n = 30;
% n = 10;
n = 5;
T = cno+cnop;
if T == 0
    pswitch = 0;
else
    pswitch = T.^(n+1)./(T.^n+(threshold(K,cnop)).^n);
end
end

function thr = threshold(K,cnop)
% Make a smooth approximation of linear_threshold (?)
% nthr = 4;
nthr = 3;
thr = K*(5/4-(5/4-1/8)*(cnop./K).^nthr./((cnop./K).^nthr+(9/64)^nthr));
% thr = K*(5/4-9/8*(cnop./K).^nthr./((cnop./K).^nthr+(9/64)^nthr));
end
