% Figure10b.m
% % Generates Figure 10(b1) (kkc = 1.25) and (b2) (kkc = 4.2) in the main text.
% Stochastic resonance curve of the 2D bistable phosphorylation model using FFT.
% Parallel computation.
% Based on the codes by Prof. Jonathan Touboul.

% X. B. 2026-10-01 Initial work.

%% Parameters
PLOTS1 = 1;
LOGPLOT = 1;
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
te = 1e4;
tspan = [ts, te];
dt=0.01;
time_vect=ts:dt:te;
% ntimes=length(time_vect);
% Exclude relaxation from FFT
relaxtime = 5e2;
relaxind = relaxtime/dt;
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
% Bistable phosphorylation threshold; change Sigmamin and Sigmamax accordingly.
kkc = 4.2; % μM
% kkc = 1.25; % μM
% Nucleocytoplasmic transport rate
% knt = 1;
% knt = 3e-1;
knt = 1e-1; % Nominal value
% Phosphorylation and dephosphorylation rate constant
tau = 0.5;
% Number of noise levels
% n_sigma = 1;
% n_sigma = 6;
% n_sigma = 15;
% n_sigma = 30;
% n_sigma = 40;
% n_sigma = 50;
n_sigma = 60;
% %{
Sigmamin = 5e-4; % For kkc = 4.2
Sigmamax = 8e-2; % For kkc = 4.2
% %}
%{
Sigmamin = 1e-3; % For kkc = 1.25
Sigmamax = 3e-1; % For kkc = 1.25
%}
Sigmas = exp(linspace(log(Sigmamin),log(Sigmamax),n_sigma));
% Sigmas = 0.0005;
% Sigmas = 0.0013;
% Sigmas = 0.01;
% Sigmas = 0.05;
% Number of samples for each sigma value
% MC = 1;
% MC = 2;
% MC = 5;
% MC = 20;
% MC = 30;
% MC = 40;
MC = 80;
% Sigma vector for parallel computation
SigmaMat = repmat(Sigmas,MC,1);
SigmaVec = reshape(SigmaMat,MC*n_sigma,1);

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

%% Parallel computation
Multiplicative = 0;
ProtectNegative = 1;
FFTVec = zeros(MC*n_sigma,1);
parfor k = 1:MC*n_sigma
% for k = 1:MC*n_sigma
    Sigma = SigmaVec(k);
    val = BistableCS2DStochasticFFT(xc0,dt,time_vect,Sigma,kvn,kvc,knt,kkc,tau,Acytoc,Multiplicative,ProtectNegative,relaxind);
    FFTVec(k) = val;
end
FFTMat = reshape(FFTVec,MC,n_sigma);
if PLOTS1 ~= 0
    figure;
    if LOGPLOT ~= 0
        errorbar(log10(Sigmas), mean(FFTMat),std(FFTMat)/sqrt(MC),'Linewidth',2);
    else
        errorbar(Sigmas, mean(FFTMat),std(FFTMat)/sqrt(MC),'Linewidth',2);
    end
end

%% Save data
if SAVE ~= 0
    save(['BistableC2D2_SRFFT_kkc=' num2str(kkc) '_tau=' num2str(tau) '_sigma=' num2str(Sigmamin) '_to_' num2str(Sigmamax) '.mat'])
end
