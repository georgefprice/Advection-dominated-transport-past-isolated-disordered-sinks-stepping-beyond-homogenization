clearvars; %close all;

N=5; %N is the number of sinks per unit length
rho=1/N; %Sink density per unit length
h=0.02/N; %Distance between each step of the iterative method
M_x=(1)/h + 1; %Nodes required in the x and y-direction
L=3; %Domain length in y-direction
L_s=2.5; %Domain length in y-direction
M_y=(2*L)/h + 1; %Nodes required in the y-direction
x=linspace(0,1,M_x); y=linspace(-L,L,M_y);
M_y_strip=(1/N)/h + 1; y_strip=linspace(-1/(2*N),1/(2*N),M_y_strip);
M_y_Ls=2*L_s/h + 1; y_Ls=linspace(-L_s,L_s,M_y_Ls);
Pe=20; Da=10; phi=sqrt(Pe^2/4 + Da);
psi=(2*Pe*phi + Pe^2 + 2*Da)*exp(phi) + (2*Pe*phi - Pe^2 - 2*Da)*exp(-phi);
varsigma = 0.01; %Used to generate a Gaussian profile about each sink location, where varsigma is the width

%% Functions
CH = @(x) (Pe/psi).*exp((Pe/2).*x).*( (2*phi - Pe).*exp(phi.*(x - 1)) + (2*phi + Pe).*exp(phi.*(1 - x)) );
G_1D = @(x,y) -(1/(2*phi)).*exp((Pe/2).*(x - y) - phi*abs(x-y));
G = @(x,y,z,xx,yy,zz) -(1/((4*pi).*sqrt((x-xx).^2+(y-yy).^2+(z-zz).^2))).*exp((Pe/2).*(x - xx) - phi.*sqrt((x-xx).^2+(y-yy).^2+(z-zz).^2));
F = @(x,y,z,xx,yy,zz,varsigma) normpdf(x,xx,varsigma).*normpdf(y,yy,varsigma).*normpdf(z,zz,varsigma);

%% Leading order expression for Exp[C_2]
Exp_C2=zeros(1,M_x);
tic
for i=1:M_x %x
    Exp_C2(i) = -(rho^3/(4*pi^(3/2)))*(1/varsigma).*integralN(@(xx) G_1D(x(i),xx).*CH(xx),0,1,'AbsTol',1e-8);
end
toc

%% Exp[C_2] including the term involving the integral of GGCH (This correction is small when N_y is large (i.e. L_s is large))
Exp_C2_2=zeros(1,M_x);
N_y=floor(L_s*N);
tic
for i=1:M_x %x
    Exp_C2_2(i) = Exp_C2(i) - (rho/((2*N_y+1)^2)).*integralN(@(xx,xxx) G_1D(x(i),xx).*G_1D(xx,xxx).*CH(xxx),0,1,0,1,'AbsTol',1e-8);
end
toc

%% Load or save data
% save('Data/Uniform_3D_Exp_N5_Pe_20_Da_10.mat','N','rho','h','M_x','M_y','L','Pe','Da','phi','psi','x','y','Exp_C2','Exp_C2_2')