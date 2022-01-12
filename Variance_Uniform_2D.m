clearvars; %close all;

%% Load Numerics
load('Data/Realisations_Uniform_2D_N5_Pe20_Da10_L3.mat')
psi=(2*Pe*phi + Pe^2 + 2*Da)*exp(phi) + (2*Pe*phi - Pe^2 - 2*Da)*exp(-phi);

%% Functions
CH = @(x) (Pe/psi).*exp((Pe/2).*x).*( (2*phi - Pe).*exp(phi.*(x - 1)) + (2*phi + Pe).*exp(phi.*(1 - x)) );
G = @(x,y,xx,yy) -(1/(2*pi)).*exp((Pe/2).*(x - xx)).*besselk(0,phi.*sqrt((x-xx).^2+(y-yy).^2));
G_1D = @(x,y) -(1/(2*phi)).*exp((Pe/2).*(x - y) - phi*abs(x-y));

%% Analytic Variance
delta=1e-5;
Var_C1=zeros(1,M_x);
N_y = 2*floor(L_s/rho)+1;
tic
for i=1:M_x %x
%     Var_C1(i)=(rho^2).*integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,0,1,-L,L,'AbsTol',1e-10) - (rho/N_y).*(integralN(@(xx) G_1D(x(i),xx).*CH(xx),0,1,'AbsTol',1e-10)).^2;
    if x(i)- delta < 0
        Var_C1(i) = (rho^2).*(integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,x(i)+delta,1,-L,- delta,'AbsTol',1e-8) + integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,x(i)+delta,1,delta,L,'AbsTol',1e-8));
    elseif x(i)+delta > 1
        Var_C1(i) = (rho^2).*(integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,0,x(i)-delta,-L,- delta,'AbsTol',1e-8) + integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,0,x(i)-delta,delta,L,'AbsTol',1e-8));
    else
        Var_C1(i) = (rho^2).*(integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,0,x(i)-delta,-L,- delta,'AbsTol',1e-8) + integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,x(i)+delta,1,-L,- delta,'AbsTol',1e-8) + integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,0,x(i)-delta,delta,L,'AbsTol',1e-8) + integralN(@(xx,yy) (G(x(i),0,xx,yy).*CH(xx)).^2,x(i)+delta,1,delta,L,'AbsTol',1e-8));
    end
    Var_C1(i) = Var_C1(i) - (rho/N_y).*(integralN(@(xx) G_1D(x(i),xx).*CH(xx),0,1,'AbsTol',1e-10)).^2;
end
toc

%% Reductions
C_var_Ls_x=mean(C_var_Ls,2);
y_reduction=0.5; % This will find the moments between -(L_s-y_reduction) and L_s-y_reduction 
reduce=2*round(y_reduction*(M_x-1)); % this calculates how many y-rows need to be removed from matricies
C_var_Ls_x_reduced = mean(C_var_Ls(:,1+reduce:2*L_s/h+1-reduce),2); % This removes the rows from the matrix

%% Plotting - All variances
figure()
hold on
for i=1+reduce:1:2*L_s/h+1-reduce
    p(i)=plot(x,C_var_Ls(:,i),'k','linewidth', 1.5);
    p(i).Color(4) = 0.02;
end
set(gca,'ColorOrderIndex',1)
p1=plot(x,C_var_Ls_x_reduced,'linewidth', 1.5);
% plot(x,C_var_Ls(:,L_s/h),'m','linewidth', 1);
set(gca,'ColorOrderIndex',3)
p2=plot(x,(Da^2).*Var_C1,'-.','linewidth', 1.5);
xlabel('$x_1$','FontSize',22,'Interpreter','latex');
leg=legend([p1,p2],'$\langle \mbox{Var}[C(\mathbf{x};\omega)] \rangle_{x_2}$','Da$^2$Var$_{\varsigma \rightarrow 0}[\widehat{C}_1(x_1,0;\omega)]$');
set(leg,'FontSize',14,'Interpreter','latex','Location','SouthEast');