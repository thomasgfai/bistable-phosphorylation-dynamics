% Function for the Euler-Maruyama steps and FFT of the solution of the
% 2D bistable phosphorylation model.
% Based on the codes by Prof. Jonathan Touboul.

% X. B. 2026-01-03 Initial work.

% function [X_stoch,val] = BistableCS2DStochasticFFT(xc0,dt,time_vect,Sigma,kvn,kvc,knt,kkc,tau,Acytoc,Multiplicative,ProtectNegative,relaxind)
function val = BistableCS2DStochasticFFT(xc0,dt,time_vect,Sigma,kvn,kvc,knt,kkc,tau,Acytoc,Multiplicative,ProtectNegative,relaxind)
ntimes = length(time_vect);
X_stoch = zeros(2,ntimes);
X_stoch(:,1) = xc0;
Sig = Sigma*[1;1];
for k = 1:ntimes-1
    t = time_vect(k);
    xs = X_stoch(:,k);
    if Multiplicative ~= 0
        Noise = Sig.*sqrt(abs(Noise_FUN_BistableCS2D1(t,xs,kvn,kvc,knt,kkc,tau,Acytoc))).*randn(2,1);
    else
        Noise = Sig.*randn(2,1);
    end
    X_stoch(:,k+1) = xs+dt*ODEFUN_BistableCS2D1(t,xs,kvn,kvc,knt,kkc,tau,Acytoc)+sqrt(dt)*Noise;
    if ProtectNegative ~= 0
        X_stoch(:,k+1) = abs(X_stoch(:,k+1));
    end
end
fourier = fft(X_stoch(1,relaxind:20:end)-mean(X_stoch(1,relaxind:20:end)));
af = abs(fourier);
[a,b] = max(af);
uu = 2;
% uu = 20;
if b>uu
    val = mean(af(b-uu:b+uu));
else
    val = af(b);
end
end

function dxcdt = ODEFUN_BistableCS2D1(t,xc,kvn,kvc,knt,kkc,tau,Acytoc)
% cno = xc(1);
% cnop = xc(2);
dxcdt = [kvn-phos(xc(1),xc(2),kkc,tau,1)+dephos(xc(1),xc(2),tau)-xc(1)*Acytoc*(kvn+kvc-knt*xc(2)); % cno
         phos(xc(1),xc(2),kkc,tau,1)-dephos(xc(1),xc(2),tau)-knt*xc(2)-xc(2)*Acytoc*(kvn+kvc-knt*xc(2))]; % cnop
end

function dxdt = NoiseFUN_CircadianPS4varND0(t,x,lambda_c,lambda_ni,A,Kd,Kb,tau)
dxdt = [ft(x(4),A,Kd)+fc(x(1),lambda_c)+x(1);
        fc(x(1),lambda_c)+fp(x(2),x(3),Kb,tau)+fdp(x(2),x(3),tau)+x(2);
        fp(x(2),x(3),Kb,tau)+fdp(x(2),x(3),tau)+lambda_ni*x(3)+x(3);
        lambda_ni*x(3)+x(4)];
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
    propensity = scale*(1/tau)*T^(n+1)/(T^n+(threshold(K,cnop))^n);
end
end

function propensity = dephos(cno,cnop,tau)
T = cno+cnop;
if T == 0
    propensity = 0;
else
    propensity = (1/tau)*cnop;
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
thr = K*(5/4-(5/4-1/8)*(cnop/K)^nthr/((cnop/K)^nthr+(9/64)^nthr));
% thr = K*(5/4-9/8*(cnop/K)^nthr/((cnop/K)^nthr+(9/64)^nthr));
end
