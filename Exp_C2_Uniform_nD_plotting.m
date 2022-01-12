clearvars; %close all;

N=5; %N_hat is the number of sinks per unit length
rho=1/N; %Sink density per unit length
h=0.02/N; %h is the distance between each step of the iterative method
M=(1)/h + 1; %N is how many nodes required to make this many steps in the x-direction
x=linspace(0,1,M);
Pe=20; Da=10; phi=sqrt(Pe^2/4 + Da);
varsigma = 0.01; %Used to generate a Gaussian profile about each sink location, where varsigma is the width
psi=(2*Pe*phi + Pe^2 + 2*Da)*exp(phi) + (2*Pe*phi - Pe^2 - 2*Da)*exp(-phi);

%% Functions
CH = @(x) (Pe/psi).*exp((Pe/2).*x).*( (2*phi - Pe).*exp(phi.*(x - 1)) + (2*phi + Pe).*exp(phi.*(1 - x)) );

%% Useful integrals
Int_GCH = -(Pe/(4*phi^2*psi))*exp((Pe/2)*x).*((Pe+2*phi)*(1+2*phi*x).*exp(phi*(1-x)) - 2*Pe*exp(phi*(x-1)));
Int_GGCH = (Pe/(8*phi^4*psi))*exp(Pe*x/2).*((- 5*Pe/2 - phi*(1 + Pe + 2*phi)).*exp(phi*(x-1)) + (2*phi+Pe)*(1+phi*x).^2.*exp(phi*(1-x)));
Int_GCHGCH = (Pe^2/(16*phi^3*psi^2))*exp(Pe*x).*((2*phi+Pe)^2*(4*phi*x+1).*exp(2*phi*(1-x)) + 4*(4*phi^2-Pe^2)*(2-exp(-2*phi*x)) + 4*(Pe^2 - 2*phi*Pe-4*phi^2).*exp(2*phi*(x-1)) );

%% Plotting
figure();
box on
hold on
p3=plot(x,CH(x),'--','LineWidth',1.5);

%% 1D
load('Data/Uniform_1D_Exp_N5_Pe20_Da10.mat')
% Calculate errors
error1D_CH = C_mean - CH(x);
error1D_CHC2 = C_mean - CH(x)-(Da^2).*Exp_C2;
error1D_CHC2_simplified = C_mean - CH(x)-(Da^2).*Exp_C2_simplified;

% Plotting
p1=plot(x,C_mean,'-o','LineWidth',1.5,'MarkerSize',7,'MarkerIndices',[round(M/30),round(7*M/30),round(13*M/30),round(19*M/30),round(25*M/30)]);
% p4=plot(x,CH(x)+(Da^2).*Exp_C2,'-.o','LineWidth',1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)]);
p4=plot(x,CH(x)+(Da^2).*Exp_C2_simplified,':o','LineWidth',1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)]);

%% 2D
load('Data/Realisations_Uniform_2D_N5_Pe20_Da10_L3.mat')
load('Data/Uniform_2D_Exp_N5_Pe_20_Da_10_L3.mat')
% Reduced moments
y_reduction=0.5; % This will find the moments between -(L_s-y_reduction) and L_s-y_reduction 
reduce=2*round(y_reduction*(M-1)); % this calculates how many y-rows need to be removed from matricies
C_mean_Ls_x_reduced = mean(C_mean_Ls(:,1+reduce:2*L_s/h+1-reduce),2); % This removes the rows from the matrix
% Calculate errors
error2D_CH = C_mean_Ls_x_reduced' - CH(x);
error2D_CHC2 = C_mean_Ls_x_reduced' - CH(x)-(Da^2).*Exp_C2;
% Plotting
p2=plot(x,C_mean_Ls_x_reduced,'-s','LineWidth',1.5,'MarkerSize',7,'MarkerIndices',[round(M/30),round(7*M/30),round(13*M/30),round(19*M/30),round(25*M/30)]);
p5=plot(x,CH(x)+(Da^2).*Exp_C2,':s','LineWidth',1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)]);

%% 3D
load('Data/Uniform_3D_Exp_N5_Pe_20_Da_10.mat')
% Plotting
p6=plot(x,CH(x)+(Da^2).*Exp_C2,':d','LineWidth',1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)]);

%% Axis and legend
xlabel('$x_1$','FontSize',22,'Interpreter','latex');
leg=legend([p1,p2,p3,p4,p5,p6],'$E [C(x_1;\omega)]$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2}$','$C_H(x_1)$','$C_H(x_1)$ + Da$^2 E[\widehat{C}_2(x_1;\omega)]$','$C_H(x_1)$ + Da$^2 E[\widehat{C}_2(x_1,0;\omega)]$','$C_H(x_1)$ + Da$^2 E[\widehat{C}_2(x_1,0,0;\omega)]$');
set(leg,'FontSize',12,'Interpreter','latex','Location','SouthWest');

%% Plotting - errors
figure()
box on
xlabel('$x_1$','FontSize',22,'Interpreter','latex');

hold on
set(gca,'ColorOrderIndex',1)
plot(x,error1D_CH,'--o','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
set(gca,'ColorOrderIndex',1)
plot(x,error2D_CH,'--s','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
set(gca,'ColorOrderIndex',3)
% plot(x,error1D_CHC2,'-.o','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
plot(x,error1D_CHC2_simplified,':o','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
set(gca,'ColorOrderIndex',5)
plot(x,error2D_CHC2,':s','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
leg=legend('$E [C(x_1;\omega)] - C_H(x_1)$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2} - C_H(x_1)$','$E [C(x_1;\omega)] - C_H(x_1) -$ Da$^2 E[\widehat{C}_2(x_1;\omega)]$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2} - C_H(x_1) -$ Da$^2 E[\widehat{C}_2(x_1,0;\omega)]$');
set(leg,'FontSize',12,'Interpreter','latex','Location','NorthWest');