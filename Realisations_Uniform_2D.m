clearvars; %close all;

global h M_x M_y M_y_Ls N rho L L_s Pe Da phi;
N=5; %N is the number of sinks per unit length
rho=1/N; %Sink density per unit length
h=0.02/N; %Distance between each step of the iterative method
M_x=(1)/h + 1; %Nodes required in the x and y-direction
L=3; %Domain length in y-direction
L_s=2.5; %Domain length in y-direction
M_y=(2*L)/h + 1; %Nodes required in the y-direction
x=linspace(0,1,M_x); y=linspace(-L,L,M_y);
M_y_Ls=2*L_s/h + 1; y_Ls=linspace(-L_s,L_s,M_y_Ls);
Pe=20; Da=10; phi=sqrt(Pe^2/4 + Da);
psi=(2*Pe*phi + Pe^2 + 2*Da)*exp(phi) + (2*Pe*phi - Pe^2 - 2*Da)*exp(-phi);
varsigma = 0.01; %Used to generate a Gaussian profile about each sink location, where varsigma is the width

%% Produce g(x,y) for uniformly random sink locations
g=zeros(M_x,M_y);
varsigmaX = varsigma; varsigmaY = varsigma; %Used to generate a Gaussian profile about each sink location, where varsigmaX and varsigmaY are the s.d in the x- and y-directions
for k=1:N
    for kk=-floor(L_s*N):1:floor(L_s*N)
        g = g + rho^2*(1/(2*pi*(varsigmaX*varsigmaY)))*exp(-((x'-rand).^2/(2*varsigmaX^2)+(y-2*L_s*(rand-1/2)).^2/(2*varsigmaY^2))); % generate a Gaussian function
    end
end

%% Calculate C(x,y) using a finite difference solver
g_vector=g(:); %Puts g(i,j) into a vector of length N*M
m=1; %Used when storing the matrix as vectors (good for sparce matrices)
i=zeros(6*(M_x)+6*(M_y-2)+5*(M_x-2)*(M_y-2),1); j=i; v=i;
for k=1:M_x*M_y
    if mod(k-1,M_x)==0 % x=0 BC
        i(m)=k; j(m)=k;     v(m)=1+3/(2*Pe*h);
        m=m+1;
        i(m)=k; j(m)=k+1;	v(m)=-2/(Pe*h);
        m=m+1;
        i(m)=k; j(m)=k+2;	v(m)=1/(2*Pe*h);
        m=m+1;
    elseif mod(k,M_x)==0 % x=1 BC
        i(m)=k; j(m)=k-2;   v(m)=1;
        m=m+1;
        i(m)=k; j(m)=k-1;	v(m)=-4;
        m=m+1;
        i(m)=k; j(m)=k;  	v(m)=3;
        m=m+1;
    elseif k<M_x+1 % y=0 BC
        i(m)=k; j(m)=k;         v(m)=-3;
        m=m+1;
        i(m)=k; j(m)=k+M_x;     v(m)=4;
        m=m+1;
        i(m)=k; j(m)=k+2*M_x;   v(m)=-1;
        m=m+1;
    elseif k > M_x*(M_y-1) % y=1 BC
        i(m)=k; j(m)=k-2*M_x;   v(m)=1;
        m=m+1;
        i(m)=k; j(m)=k-M_x;     v(m)=-4;
        m=m+1;
        i(m)=k; j(m)=k;         v(m)=3;
        m=m+1;
    else
        i(m)=k; j(m)=k-1;	v(m)=1/(h^2) + Pe/(2*h);
        m=m+1;
        i(m)=k; j(m)=k+1;	v(m)=1/(h^2) - Pe/(2*h);
        m=m+1;
        i(m)=k; j(m)=k;     v(m)=-4/(h^2)-Da*g_vector(k);
        m=m+1;
        i(m)=k; j(m)=k-M_x;	v(m)=1/(h^2);
        m=m+1;
        i(m)=k; j(m)=k+M_x;	v(m)=1/(h^2);
        m=m+1;
    end
end
A = sparse(i,j,v,M_x*M_y,M_x*M_y);
b=zeros(M_x*M_y,1);
for k=1:M_x*M_y
    if mod(k-1,M_x)==0
        b(k)=1;
    end
end
Z=A\b;
C=reshape(Z,M_x,M_y);

%% Calculate some (crude) moments of g and C
g_Ls=g(:,L/h-L_s/h+1:L/h+L_s/h+1);
C_Ls=C(:,L/h-L_s/h+1:L/h+L_s/h+1);
g_bar = simpson2Dy(g_Ls,-L_s,L_s)/(2*L_s);
g_var = simpson2Dy((g_Ls-g_bar).^2,-L_s,L_s)/(2*L_s);
C_bar = simpson2Dy(C_Ls,-L_s,L_s)/(2*L_s);
C_var = simpson2Dy((C_Ls-C_bar).^2,-L_s,L_s)/(2*L_s);
g_bar_bar = simpson1D(g_bar,-L_s,L_s)/(2*L_s);
g_var_var = simpson1D((g_bar-g_bar_bar).^2,-L_s,L_s)/(2*L_s);
C_bar_bar = simpson1D(C_bar,-L_s,L_s)/(2*L_s);
C_var_var = simpson1D((C_bar-C_bar_bar).^2,-L_s,L_s)/(2*L_s);

%% Plotting - g(x,y)
figure();
p1=pcolor(x,y,g');
daspect([1 1 1]);
set(p1,'edgecolor','none');
colormap jet;
c=colorbar;
% caxis([0,inf]);
shading interp
xlabel('$x_1$','FontSize',16,'Interpreter','latex');
ylabel('$x_2$','FontSize',16,'Interpreter','latex');
% ylabel(c,'$g(x,y)$','FontSize',16,'Interpreter','latex');

%% Plotting - C(x,y)
figure();
p1=pcolor(x,y,C');
daspect([1 1 1]);
set(p1,'edgecolor','none');
colormap jet;
c=colorbar;
caxis([0.6,1]);
shading interp
xlabel('$x_1$','FontSize',16,'Interpreter','latex');
ylabel('$x_2$','FontSize',16,'Interpreter','latex');
% ylabel(c,'$C(x,y)$','FontSize',16,'Interpreter','latex');