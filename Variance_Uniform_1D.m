clearvars; %close all;

load('Data/Realisations_Uniform_1D_N5_Pe20_Da10.mat')

%% Functions
CH = @(x) (Pe/psi).*exp((Pe/2).*x).*( (2*phi - Pe).*exp(phi.*(x - 1)) + (2*phi + Pe).*exp(phi.*(1 - x)) );
G = @(x,y) -(1/(2*phi)).*exp((Pe/2).*(x - y) - phi.*sqrt((x - y).^2));
GCH = @(x,y) G(x,y).*CH(y);
F = @(x,y,varsigma) (1/(2*pi*varsigma^2).^(1/2)).*exp(-(1/(2*varsigma^2)).*(x-y).^2);

%% Delta function
Var_C1_delta=zeros(1,M);
for i=1:M %x
    Var_C1_delta(i)=rho.*(integralN(@(xx) GCH(x(i),xx).^2,0,1,'AbsTol',1e-12) - (integralN(@(xx) GCH(x(i),xx),0,1,'AbsTol',1e-12)).^2);
end

%% Different varsigma
Var_C1=zeros(1,M);
for i=1:M %x
    Var_C1(i)=integralN(@(xx,yy) GCH(x(i),xx).*GCH(x(i),yy).*rho.*(F(xx,yy,sqrt(2)*varsigma)-1),0,1,0,1,'AbsTol',1e-12);
end

%% Plotting
figure()
hold on
plot(x,C_var,'linewidth', 1.5)
plot(x,(Da^2)*Var_C1,':','linewidth', 1.5)
plot(x,(Da^2)*Var_C1_delta,'-.','linewidth', 1.5)
xlabel('$x_1$','FontSize',22,'Interpreter','latex');
leg=legend('Var$[C(x_1;\omega)]$','Da$^2$Var$[\widehat{C}_1(x_1;\omega)]$','Da$^2$Var$_{\varsigma \rightarrow 0}[\widehat{C}_1(x_1;\omega)]$');
set(leg,'FontSize',14,'Interpreter','latex');
ylim([-inf inf])