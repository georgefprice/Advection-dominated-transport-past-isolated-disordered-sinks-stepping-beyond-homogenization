clearvars; %close all;

load('Data/Realisations_Uniform_2D_N5_Pe20_Da10_L3.mat')

%% Functions
CH = @(x) (Pe/psi).*exp((Pe/2).*x).*( (2*phi - Pe).*exp(phi.*(x - 1)) + (2*phi + Pe).*exp(phi.*(1 - x)) );
G_1D = @(x,y) -(1/(2*phi)).*exp((Pe/2).*(x - y) - phi*abs(x-y));
G = @(x,y,xx,yy) -(1/(2*pi)).*exp((Pe/2).*(x - xx)).*besselk(0,phi.*sqrt((x-xx).^2+(y-yy).^2));
F = @(x,y,xx,yy,varsigma) normpdf(x,xx,varsigma).*normpdf(y,yy,varsigma);

%% Leading order expression for Exp[C_2]
Exp_C2=zeros(1,M_x);
tic
for i=1:M_x %x
    Exp_C2(i) = -(rho^2/(4*pi))*(eulergamma - 2*log(2*varsigma*phi)).*integralN(@(xx) G_1D(x(i),xx).*CH(xx),0,1,'AbsTol',1e-8);
end
toc

%% Exp[C_2] including the term involving the integral of GGCH (This correction is small when N_y is large (i.e. L_s is large))
Exp_C2_2=zeros(1,M_x);
N_y=floor(L_s*N);
tic
for i=1:M_x %x
    Exp_C2_2(i) = Exp_C2(i) - (rho/(2*N_y+1)).*integralN(@(xx,xxx) G_1D(x(i),xx).*G_1D(xx,xxx).*CH(xxx),0,1,0,1,'AbsTol',1e-8);
end
toc

%% Leading order expression for Exp[C_2] in 3D
Exp_C2_3D=zeros(1,M_x);
tic
for i=1:M_x %x
    Exp_C2_3D(i) = -(rho^3/(4*pi^(3/2)))*(1/varsigma).*integralN(@(xx) G_1D(x(i),xx).*CH(xx),0,1,'AbsTol',1e-8);
end
toc

%% Load or save data
% save('Data/Uniform_2D_Exp_N5_Pe_20_Da_10_L3.mat','N','rho','h','M_x','M_y','L','Pe','Da','phi','psi','x','y','Exp_C2','Exp_C2_2','Exp_C2_3D')