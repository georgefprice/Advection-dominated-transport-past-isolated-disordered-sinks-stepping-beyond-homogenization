clearvars; %close all;
rng shuffle

% global h M_x M_y M_y_strip M_y_Ls N rho L L_s Pe Da phi;
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

%% Realsiations
Num=10000;
g_total   	= zeros(M_x,M_y,Num);
C_total 	= zeros(M_x,M_y,Num);
parpool(8);
parfor eta=1:Num
    %Produce g(x,y) for uniformly random sink locations
    g=zeros(M_x,M_y);
    for k=1:N
        for kk=-floor(L_s*N):1:floor(L_s*N)
            g = g + rho^2*(1/(2*pi*varsigma^2))*exp(-((x'-rand).^2/(2*varsigma^2)+(y-2*L_s*(rand-1/2)).^2/(2*varsigma^2))); % generate a Gaussian function
        end
    end
    
    %Calculate C(x,y) using a finite difference solver
    g_vector=g(:); %Puts g(i,j) into a vector of length N*M
    m=1; %Used when storing the matrix as vectors (good for sparce matrices)
    i=zeros(6*(M_x)+6*(M_y-2)+5*(M_x-2)*(M_y-2),1); j=i; v=i;
    for k=1:M_x*M_y
        if mod(k-1,M_x)==0 % x=0 BC
            i(m)=k; j(m)=k;     v(m)=1+3/(2*Pe*h);
            m=m+1;
            i(m)=k; j(m)=k+1;   v(m)=-2/(Pe*h);
            m=m+1;
            i(m)=k; j(m)=k+2;   v(m)=1/(2*Pe*h);
            m=m+1;
        elseif mod(k,M_x)==0 % x=1 BC
            i(m)=k; j(m)=k-2;       v(m)=1;
            m=m+1;
            i(m)=k; j(m)=k-1;       v(m)=-4;
            m=m+1;
            i(m)=k; j(m)=k;         v(m)=3;
            m=m+1;
        elseif k<M_x+1 % y=0 BC
            i(m)=k; j(m)=k;         v(m)=-3;
            m=m+1;
            i(m)=k; j(m)=k+M_x;     v(m)=4;
            m=m+1;
            i(m)=k; j(m)=k+2*M_x;	v(m)=-1;
            m=m+1;
        elseif k > M_x*(M_y-1) % y=1 BC
            i(m)=k; j(m)=k-2*M_x;   v(m)=1;
            m=m+1;
            i(m)=k; j(m)=k-M_x;     v(m)=-4;
            m=m+1;
            i(m)=k; j(m)=k;         v(m)=3;
            m=m+1;
        else % middle of domain
            i(m)=k; j(m)=k-1;       v(m)=1/(h^2) + Pe/(2*h);
            m=m+1;
            i(m)=k; j(m)=k+1;       v(m)=1/(h^2) - Pe/(2*h);
            m=m+1;
            i(m)=k; j(m)=k;         v(m)=-4/(h^2) - Da*g_vector(k);
            m=m+1;
            i(m)=k; j(m)=k-M_x;     v(m)=1/(h^2);
            m=m+1;
            i(m)=k; j(m)=k+M_x;     v(m)=1/(h^2);
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
    
    %Storing C and g
    g_total(:,:,eta) = g;
    C_total(:,:,eta) = C;
end
g_mean     	= mean(g_total,3);
g_var     	= var(g_total,0,3);
g_mean_Ls   = mean(g_total(:,L/h-L_s/h+1:L/h+L_s/h+1,:),3);
g_var_Ls    = var(g_total(:,L/h-L_s/h+1:L/h+L_s/h+1,:),0,3);
g_mean_strip= mean(g_total(:,(M_y-M_y_strip+2)/2:(M_y+M_y_strip)/2,:),3);
g_var_strip = var(g_total(:,(M_y-M_y_strip+2)/2:(M_y+M_y_strip)/2,:),0,3);
C_mean     	= mean(C_total,3);
C_var     	= var(C_total,0,3);
C_mean_Ls   = mean(C_total(:,L/h-L_s/h+1:L/h+L_s/h+1,:),3);
C_var_Ls    = var(C_total(:,L/h-L_s/h+1:L/h+L_s/h+1,:),0,3);
C_mean_strip= mean(C_total(:,(M_y-M_y_strip+2)/2:(M_y+M_y_strip)/2,:),3);
C_var_strip = var(C_total(:,(M_y-M_y_strip+2)/2:(M_y+M_y_strip)/2,:),0,3);

delete(gcp('nocreate'))

save('Data/Realisations_Uniform_2D_N5_Pe20_Da10_L3.mat','M_x','M_y','M_y_strip','h','N','rho','L','L_s','Pe','Da','phi','psi','x','y','y_strip','y_Ls','varsigma','Num','g_mean','g_var','g_mean_Ls','g_var_Ls','g_mean_strip','g_var_strip','C_mean','C_var','C_mean_Ls','C_var_Ls','C_mean_strip','C_var_strip')
