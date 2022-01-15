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
p3=plot(x, CH(x), '--', 'color', [0, 0.4470, 0.7410], 'LineWidth', 1.5);

%% 1D
load('Data/Uniform_1D_Exp_N5_Pe20_Da10.mat')
% Calculate errors
error1D_CH = C_mean - CH(x);
error1D_CHC2 = C_mean - CH(x)-(Da^2).*Exp_C2;
error1D_CHC2_simplified = C_mean - CH(x)-(Da^2).*Exp_C2_simplified;

%% 2D
load('Data/Realisations_Uniform_2D_N5_Pe20_Da10_L3.mat')
load('Data/Uniform_2D_Exp_N5_Pe_20_Da_10_L3.mat')
% Reduced moments
y_reduction=0.5; % This will find the moments between -(L_s-y_reduction) and L_s-y_reduction 
reduce=2*round(y_reduction*(M-1)); % this calculates how many y-rows need to be removed from matricies
C_mean_Ls_x_reduced = mean(C_mean_Ls(:,1+reduce:2*L_s/h+1-reduce),2); % This removes the rows from the matrix
% Effective Uptake
Da_eff_2D = Da*(1-(rho^2/(4*pi))*Da*(eulergamma - 2*log(varsigma*sqrt(Pe^2 + 4*Da))));
% Da_eff_2D = Da*(1+(rho^2/(2*pi))*Da*log(varsigma*sqrt(Pe^2 + 4*Da))); % Assumes eulergamma << 2*log(varsigma*sqrt(Pe^2 + 4*Da)), i.e. 0.5772 << 3.1236
phi_eff_2D = sqrt(Pe^2/4 + Da_eff_2D);
psi_eff_2D = (2*Pe*phi_eff_2D + Pe^2 + 2*Da_eff_2D)*exp(phi_eff_2D) + (2*Pe*phi_eff_2D - Pe^2 - 2*Da_eff_2D)*exp(-phi_eff_2D);
CH_eff_2D = @(x) (Pe/psi_eff_2D).*exp((Pe/2).*x).*( (2*phi_eff_2D - Pe).*exp(phi_eff_2D.*(x - 1)) + (2*phi_eff_2D + Pe).*exp(phi_eff_2D.*(1 - x)) );
% Calculate errors
error2D_CH = C_mean_Ls_x_reduced' - CH(x);
error2D_CHC2 = C_mean_Ls_x_reduced' - CH(x)-(Da^2).*Exp_C2;
error2D_CHCeff = C_mean_Ls_x_reduced' - CH_eff_2D(x);
% Plotting
p2=plot(x,C_mean_Ls_x_reduced, '-s', 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerIndices', [round(M/30),round(7*M/30),round(13*M/30),round(19*M/30),round(25*M/30)]);
p5=plot(x,CH(x)+(Da^2).*Exp_C2, ':s', 'color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerIndices', [round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)]);
p7=plot(x,CH_eff_2D(x), '-.s', 'color',	[0, 0.75, 0.75], 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerIndices', [round(5*M/30),round(11*M/30),round(17*M/30),round(23*M/30),round(29*M/30)]);

%% 3D
load('Data/Uniform_3D_Exp_N5_Pe_20_Da_10.mat')
% Effective Uptake
% Da_eff_3D = Da*(1 - rho^3/(4*pi^(3/2)*varsigma)*Da*(1-(varsigma*sqrt(pi)/2)*sqrt(Pe^2+4*Da)));
Da_eff_3D = Da*(1 - rho^3/(4*pi^(3/2)*varsigma)*Da); % Assumes 1 >> (varsigma*sqrt(pi)/2)*sqrt(Pe^2+4*Da), i.e. 1 >> 0.1859
phi_eff_3D = sqrt(Pe^2/4 + Da_eff_3D);
psi_eff_3D = (2*Pe*phi_eff_3D + Pe^2 + 2*Da_eff_3D)*exp(phi_eff_3D) + (2*Pe*phi_eff_3D - Pe^2 - 2*Da_eff_3D)*exp(-phi_eff_3D);
CH_eff_3D = @(x) (Pe/psi_eff_3D).*exp((Pe/2).*x).*( (2*phi_eff_3D - Pe).*exp(phi_eff_3D.*(x - 1)) + (2*phi_eff_3D + Pe).*exp(phi_eff_3D.*(1 - x)) );
% Plotting
p6=plot(x,CH(x)+(Da^2).*Exp_C2, ':d', 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerIndices', [round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)]);
p8=plot(x,CH_eff_3D(x), '-.d', 'color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerIndices', [round(5*M/30),round(11*M/30),round(17*M/30),round(23*M/30),round(29*M/30)]);

%% Axis and legend
xlabel('$x_1$','FontSize',22,'Interpreter','latex');
% leg=legend([p1,p2,p3,p4,p5,p6,p7,p8],'$E [C(x_1;\omega)]$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2}$','$C_H(x_1)$','$C_H(x_1)$ + Da$^2 E[\widehat{C}_2(x_1;\omega)]$','$C_H(x_1)$ + Da$^2 E[\widehat{C}_2(x_1,0;\omega)]$','$C_H(x_1)$ + Da$^2 E[\widehat{C}_2(x_1,0,0;\omega)]$','$C_H^{UR}(x_1) \quad$ (2D)','$C_H^{UR}(x_1) \quad$ (3D)');
leg=legend([p2,p3,p5,p6,p7,p8],'$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2}$','$C_H(x_1)$','$C_H(x_1)$ + Da$^2 E[\widehat{C}_2(x_1,0;\omega)]$','$C_H(x_1)$ + Da$^2 E[\widehat{C}_2(x_1,0,0;\omega)]$','$C_H^{UR}(x_1) \quad$ (2D)','$C_H^{UR}(x_1) \quad$ (3D)');
set(leg,'FontSize',12,'Interpreter','latex','Location','SouthWest');

%% Plotting - errors
figure()
box on
xlabel('$x_1$','FontSize',22,'Interpreter','latex');

hold on
set(gca,'ColorOrderIndex',1)
% plot(x,error1D_CH,'--o','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
set(gca,'ColorOrderIndex',1)
plot(x,error2D_CH,'--s','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
set(gca,'ColorOrderIndex',3)
% % plot(x,error1D_CHC2,'-.o','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
% plot(x,error1D_CHC2_simplified,':o','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
set(gca,'ColorOrderIndex',5)
plot(x,error2D_CHC2,':s','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
plot(x,error2D_CHCeff,'m-.s','linewidth', 1.5,'MarkerSize',7,'MarkerIndices',[round((3*M-15)/30),round((9*M-15)/30),round((15*M-15)/30),round((21*M-15)/30),round((27*M-15)/30)])
% leg=legend('$E [C(x_1;\omega)] - C_H(x_1)$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2} - C_H(x_1)$','$E [C(x_1;\omega)] - C_H(x_1) -$ Da$^2 E[\widehat{C}_2(x_1;\omega)]$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2} - C_H(x_1) -$ Da$^2 E[\widehat{C}_2(x_1,0;\omega)]$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2} - C_H^{UR}(x_1)$');
leg=legend('$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2} - C_H(x_1)$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2} - C_H(x_1) -$ Da$^2 E[\widehat{C}_2(x_1,0;\omega)]$','$\langle E [C(x_1,x_2;\omega)] \rangle_{x_2} - C_H^{UR}(x_1)$');
set(leg,'FontSize',12,'Interpreter','latex','Location','NorthWest');