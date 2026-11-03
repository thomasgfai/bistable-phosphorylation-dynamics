% Figure10c_11.m
% % Generates Figure 10(c1) (Sigmas = 0.001), (c2) (Sigmas = 0.01), (c3) (Sigmas = 0.08), and Figure 11 (Sigmas = 0.01) in the main text.

% X. B. 2026-10-01 Initial work.

%% Parameters

PLOTS = 1;
SAVE = 1;

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
% Step length for Euler-Maruyama
dt = 0.001;
time_vect = ts:dt:te;
ntimes = length(time_vect);

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
kkc = 4.2; % μM;
% Nucleocytoplasmic transport rate
knt = 1e-1; % Nominal value
% Phosphorylation and dephosphorylation rate constant
tau = 0.5;

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
[tc,xc] = ode23s(@(t,xc) ODEFUN_BistableCS2D1(t,xc,kvn,kvc,knt,kkc,tau,Acytoc),time_vect,xc0);

% Stochastic
Multiplicative=0;
ProtectNegative=1;
X_stoch=zeros(2,ntimes);
X_stoch(:,1)=xc0;
% Sigmas = 0.001;
Sigmas = 0.01;
% Sigmas = 0.08;
ns=length(Sigmas);
Means=zeros(1,ns);
Stds=zeros(1,ns);
npks=zeros(1,ns);

for ks=1:ns
    Sig=Sigmas(ks)*[1;1];
    for k=1:ntimes-1
        t=time_vect(k);
        xs=X_stoch(:,k);
        if Multiplicative ~= 0
            Noise=Sig.*sqrt(abs(Noise_FUN_BistableCS2D1(t,xs,kvn,kvc,knt,kkc,tau,Acytoc))).*randn(2,1);
        else
            Noise=Sig.*randn(2,1);
        end
        X_stoch(:,k+1)=xs+dt*ODEFUN_BistableCS2D1(t,xs,kvn,kvc,knt,kkc,tau,Acytoc)+...
            sqrt(dt)*Noise;
        
        if ProtectNegative ~= 0
            X_stoch(:,k+1)=abs(X_stoch(:,k+1));
            % if X_stoch(1,k+1)<0
            %     X_stoch(1,k+1)=X_stoch(1,k+1)-2*sqrt(dt)*Noise(1);
            %
            % end
            % if X_stoch(2,k+1)<0
            %     X_stoch(2,k+1)=X_stoch(2,k+1)-2*sqrt(dt)*Noise(2);
            % end
        end
    end
    
    if PLOTS ~= 0
        times_plotted=1:ntimes;
        % times_plotted=ceil(ntimes/2):ntimes;
        figure(1)
        plot(time_vect(times_plotted),X_stoch(1,times_plotted),'Linewidth',2)
        hold on
        plot(tc(times_plotted),xc(times_plotted,1),'Linewidth',2)
        hold off
        figure(2)
        plot(X_stoch(2,times_plotted)./(X_stoch(1,times_plotted)+X_stoch(2,times_plotted)),(X_stoch(1,times_plotted)+X_stoch(2,times_plotted)),'Linewidth',2);
        hold on
        plot(xc(times_plotted,2)./(xc(times_plotted,1)+xc(times_plotted,2)),(xc(times_plotted,1)+xc(times_plotted,2)),'Linewidth',2);
        scatter(X_stoch(2,end)./(X_stoch(1,end)+X_stoch(2,end)),(X_stoch(1,end)+X_stoch(2,end)),60,'filled')
        scatter(xc(end,2)./(xc(end,1)+xc(end,2)),(xc(end,1)+xc(end,2)),60,'k','filled')
        hold off
    end
end

%% Calculate nullclines 
% Plots are not really useful but still cannot be turned off
cnomargin = abs(kkc-1.2);
cnopmargin = abs(kkc-1.2);
figure(1)
cnomax = max(xc(:,1))+cnomargin;
cnopmax = max(xc(:,2))+cnopmargin;
cnonull = fimplicit(@(cno,cnop) dcnodt(cno,cnop,kvn,kvc,knt,Acytoc,kkc,tau),[0 cnomax 0 cnopmax],'Linewidth',1.5);
hold on
cnopnull = fimplicit(@(cno,cnop) dcnopdt(cno,cnop,kvn,kvc,knt,Acytoc,kkc,tau),[0 cnomax 0 cnopmax],'Linewidth',1.5);
hold off
cnonullx = cnonull.XData;
cnonully = cnonull.YData;
cnopnullx = cnopnull.XData;
cnopnully = cnopnull.YData;
clear cnonull cnopnull

%% Plots
ntimes = length(time_vect);
times_plotted=1:ntimes;
% times_plotted=ceil(ntimes/2):ntimes;
figure(2)
plot((X_stoch(1,times_plotted)+X_stoch(2,times_plotted)),X_stoch(2,times_plotted)./(X_stoch(1,times_plotted)+X_stoch(2,times_plotted)),'Color','green','Linewidth',2);
hold on
plot(cnonullx+cnonully,cnonully./(cnonullx+cnonully),'--','Color','blue','Linewidth',3)
plot(cnopnullx+cnopnully,cnopnully./(cnopnullx+cnopnully),':','Color','magenta','Linewidth',3)
plot((xc(times_plotted,1)+xc(times_plotted,2)),xc(times_plotted,2)./(xc(times_plotted,1)+xc(times_plotted,2)),'Color','red','Linewidth',2);
scatter((X_stoch(1,end)+X_stoch(2,end)),X_stoch(2,end)./(X_stoch(1,end)+X_stoch(2,end)),60,'filled')
scatter((xc(end,1)+xc(end,2)),xc(end,2)./(xc(end,1)+xc(end,2)),60,'k','filled')
hold off
xlabel('c_{no}+c_{nop}')
ylabel('Fraction of c_{nop}')

% %{
figure(3)
plot(time_vect(times_plotted),X_stoch(1,times_plotted),'Color','green','Linewidth',2)
hold on
plot(tc(times_plotted),xc(times_plotted,1),'Color','red','Linewidth',2)
hold off
xlabel('Time')
ylabel('C_{no}')

figure(4)
plot(time_vect(times_plotted),X_stoch(2,times_plotted),'Color','green','Linewidth',2)
hold on
plot(tc(times_plotted),xc(times_plotted,2),'Color','red','Linewidth',2)
hold off
xlabel('Time')
ylabel('C_{nop}')
% %}

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
