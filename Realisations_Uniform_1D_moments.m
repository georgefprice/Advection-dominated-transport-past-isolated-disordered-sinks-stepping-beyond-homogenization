clearvars; %close all;
rng shuffle

global h M N rho Pe Da phi;
N=5; %N_hat is the number of sinks per unit length
rho=1/N; %Sink density per unit length
h=0.02/N; %h is the distance between each step of the iterative method
M=(1)/h + 1; %N is how many nodes required to make this many steps in the x-direction
x=linspace(0,1,M);
Pe=20; Da=10; phi=sqrt(Pe^2/4 + Da);
psi=(2*Pe*phi + Pe^2 + 2*Da)*exp(phi) + (2*Pe*phi - Pe^2 - 2*Da)*exp(-phi);
varsigma = 0.01; %Used to generate a Gaussian profile about each sink location, where varsigma is the width

%% Numeric - varsigma
Num=10^6; %Do 10^6 for figures ;
g_total=zeros(Num,M);
C_total=zeros(Num,M);
parpool(8);
parfor eta=1:Num
    %Produce g(x) for uniformly random sink locations
    g=zeros(1,M);
    xi=rand(N,1);
    for k=1:N
        g = g + rho*normpdf(x,xi(k),varsigma); % generate g(x)
    end
    g_total(eta,:)=g;

    % Calculate C(x) using a finite difference solver
    A=zeros(M,M);
    % Inlet BC
    A(1,1) = 1  + 3/(2*Pe*h);
    A(1,2) =    - 2/(Pe*h)  ;
    A(1,3) =    + 1/(2*Pe*h);
    % Governing ODE
    for i=2:M-1
        A(i,i-1) =          + Pe/(2*h)  + 1/(h^2);
        A(i,i)   = - Da*g(i)           	- 2/(h^2);
        A(i,i+1) =          - Pe/(2*h)  + 1/(h^2);
    end
    % Outlet BC
    A(M,M-2)    = + 1/(2*h);
    A(M,M-1)    = - 2/h    ;
    A(M,M)      = + 3/(2*h);
    
    bb=zeros(M,1);
    bb(1)=1;
    C_total(eta,:) = A\bb;
end
g_mean = mean(g_total);
g_var = var(g_total);
C_mean = mean(C_total);
C_var = var(C_total);
delete(gcp('nocreate'))

%% Numerics - delta
C_total = zeros(Num,M);
parpool(6);
parfor eta = 1:Num
    bb=zeros(2*(N+1),1); bb(1)=1;
    C = zeros(1,M);
    A=zeros(2*(N+1),2*(N+1));
    A(1,N+2)=1;
    for i=2:N+1
        A(i,i) = rho*Da/Pe;
        A(i,N+i) = - 1;
        A(i,N+1+i) = 1 + rho*Da/Pe;
    end
    for i=1:N
        A(N+1+i,i+1)=rho*Da/Pe - 1;
        A(N+1+i,N+1+i+1)=rho*Da/Pe;
    end

    r = rand(N,1);
    r = sort(r)';
    %     r = linspace(1/(2*N),1 - 1/(2*N),N);
    xi = [0 r 1];
    
    for i=1:N+1
        A(N+1+i,i)=exp(Pe*(xi(i+1)-xi(i)));
    end
    Z=linsolve(A,bb);
    for i=1:M
        for j=1:N+1
            if (xi(j) <= x(i)) && (x(i) <= xi(j+1))
                C(i)=Z(j)*exp(Pe*(x(i)-xi(j)))+Z(j+N+1);
                break
            end
        end
    end
    C_total(eta,:) = C;
end
C_mean_delta = mean(C_total);
C_var_delta = var(C_total);
delete(gcp('nocreate'))

%% Save data
save('Data/Realisations_Uniform_1D_N5_Pe20_Da10.mat','N','rho','h','M','Num','x','Da','Pe','phi','psi','varsigma','g_mean','g_var','C_mean','C_var','C_mean_delta','C_var_delta')