clearvars; %close all;

load('Data/Realisations_Uniform_1D_N5_Pe20_Da10.mat')
psi=(2*Pe*phi + Pe^2 + 2*Da)*exp(phi) + (2*Pe*phi - Pe^2 - 2*Da)*exp(-phi);

%% Functions
CH = @(x) (Pe/psi).*exp((Pe/2).*x).*( (2*phi - Pe).*exp(phi.*(x - 1)) + (2*phi + Pe).*exp(phi.*(1 - x)) );
G = @(x,y) -(1/(2*phi)).*exp((Pe/2).*(x - y) - phi*abs(x-y));
F = @(x,y,varsigma) normpdf(x,y,varsigma);

Cov_F = @(x,y,varsigma) rho*(F(x,y,sqrt(2)*varsigma)-1);
I_C2 = @(x) integralN(@(xx,xxx) G(x,xx).*G(xx,xxx).*CH(xxx).*Cov_F(xx,xxx,varsigma),0,1,0,1,'AbsTol',1e-10);
I_C2_delta_1 = @(x) integralN(@(xx) G(x,xx).*G(xx,xx).*CH(xx),0,1,'AbsTol',1e-10);
I_C2_delta_2 = @(x,xx) integralN(@(xx,xxx) G(x,xx).*G(xx,xxx).*CH(xxx),0,1,0,1,'AbsTol',1e-10);

%% E[C_1] & E[C_2]
Exp_C1=zeros(1,M);
Exp_C2=zeros(1,M);
tic
for i=1:M %x
    Exp_C2(i) = I_C2(x(i));
end
toc

%% Expectation - Simplified
Exp_C1=zeros(1,M);
Exp_C2_simplified=zeros(1,M);
tic
for i=1:M %x
    Exp_C2_simplified(i) = - (rho/(2*phi)).*integralN(@(xx) G(x(i),xx).*CH(xx),0,1,'AbsTol',1e-10) - rho.*integralN(@(xx,xxx) G(x(i),xx).*G(xx,xxx).*CH(xxx),0,1,0,1,'AbsTol',1e-10);
end
toc

%% Expectation - delta
Exp_C1_delta=zeros(1,M);
Exp_C2_delta=zeros(1,M);
for i=1:M %x
    Exp_C2_delta(i) = rho*(I_C2_delta_1(x(i)) - I_C2_delta_2(x(i)));
end

%% Saving/Loading
% save('Data/Uniform_1D_Exp_N5_Pe20_Da10.mat','N','rho','h','M','Num','x','Da','Pe','phi','psi','varsigma','g_mean','g_var','C_mean','C_var','C_mean_delta','C_var_delta','Exp_C1','Exp_C2','Exp_C2_simplified','Exp_C1_delta','Exp_C2_delta')