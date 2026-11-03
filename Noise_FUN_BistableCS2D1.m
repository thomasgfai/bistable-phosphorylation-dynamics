% ODE function for the simplified translation model
% with bistable phosphorylation of the nuclear cargo proteins.
% "Protein concentration" viewpoint
% 2 variables: cno, cnop
% For bifurcation diagrams
% The translation rates are set to be proportional to the cytoplasmic volume.

% X. B. 2025-01-13 Initial work.

function dxcdt = Noise_FUN_BistableCS2D1(t,xc,kvn,kvc,knt,kkc,tau,Acytoc)

% cno = xc(1);
% cnop = xc(2);

dxcdt = [kvn+phos(xc(1),xc(2),kkc,tau,1)+dephos(xc(1),xc(2),tau)+xc(1)*Acytoc*(kvn+kvc-knt*xc(2)); % cno
         phos(xc(1),xc(2),kkc,tau,1)+dephos(xc(1),xc(2),tau)+knt*xc(2)+xc(2)*Acytoc*(kvn+kvc+knt*xc(2))]; % cnop
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

function thr = linear_threshold(K,cnop)
thr = max(0.125*K, 1.25*K-4*cnop);
% That is, thr = K*max(0.125, 1.25-4*cnop/K);
% thr=K;
end

function thr = threshold(K,cnop)
% Make a smooth approximation of linear_threshold (?)
% nthr = 4;
nthr = 3;
thr = K*(5/4-(5/4-1/8)*(cnop/K)^nthr/((cnop/K)^nthr+(9/64)^nthr));
% thr = K*(5/4-9/8*(cnop/K)^nthr/((cnop/K)^nthr+(9/64)^nthr));
end

