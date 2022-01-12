clearvars; %close all;

global M Pe Da phi;
Pe=20; Da=10; phi = sqrt(Pe^2/4 + Da);
M=1000;
x=linspace(0,1,M+1);
xx=linspace(0,1,M+1);
G=zeros(M+1,M+1);
GG=zeros(M+1,M+1);
psi = (2*Pe*phi + Pe^2 + 2*Da)*exp(phi) + (2*Pe*phi - Pe^2 - 2*Da)*exp(-phi);
Psi = 2*Pe^2 - 2*Da*exp(-Pe);
for i=1:M+1
    for j=1:M+1
        if x(i)<xx(j)
            G(i,j) = -(1/(4*phi*psi))*exp((Pe/2)*(x(i)-xx(j)))*((2*phi + Pe)^2*exp(phi*(x(i)-xx(j)+1)) + (2*phi - Pe)^2*exp(-phi*(x(i)-xx(j)+1)) + 4*Da*(exp(phi*(x(i)+xx(j)-1)) + exp(-phi*(x(i)+xx(j)-1))));
            GG(i,j) = -(1/(2*phi))*exp((Pe/2)*(x(i) - xx(j)) + phi*(x(i) - xx(j)));
        else
            G(i,j) = -(1/(4*phi*psi))*exp((Pe/2)*(x(i)-xx(j)))*((2*phi + Pe)^2*exp(phi*(xx(j)-x(i)+1)) + (2*phi - Pe)^2*exp(-phi*(xx(j)-x(i)+1)) + 4*Da*(exp(phi*(x(i)+xx(j)-1)) + exp(-phi*(x(i)+xx(j)-1))));
            GG(i,j) = -(1/(2*phi))*exp((Pe/2)*(x(i) - xx(j)) - phi*(x(i) - xx(j)));
        end
    end
end

%% Plotting - Exact Green's function
figure()
contour(x,xx,G',25,'linewidth', 1.5);
colormap(gca,'autumn');
c=colorbar;
shading interp
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_1''$','FontSize',20,'Interpreter','latex');
daspect([1 1 1])
caxis([-inf 0]);

%% Plotting - free-space Green's function
figure()
contour(x,xx,GG',25,'linewidth', 1.5);
colormap(gca,'autumn');
c=colorbar;
shading interp
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_1''$','FontSize',20,'Interpreter','latex');
daspect([1 1 1])
caxis([-inf 0]);

%% Plotting - Error
figure()
contour(x,xx,G'-GG',20,'linewidth', 1.5);
colormap(gca,'autumn');
c=colorbar;
shading interp
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_1''$','FontSize',20,'Interpreter','latex');
daspect([1 1 1])
caxis([-inf 0]);

figure()
plot(x,G(:,round((M+1)/2)));
hold on
plot(x,GG(:,round((M+1)/2)));
xlabel('$x$','FontSize',20,'Interpreter','latex');
leg=legend('$(x,0.5)$ exact','$G(x,0.5)$ free-space');
set(leg,'FontSize',18,'Interpreter','latex');
ylim([-inf 0]);