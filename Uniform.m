clearvars; %close all;

M=500;
N=19; epsilon=1/(N+1);
S=epsilon;
x=linspace(0,N+1,M*(N+1)); X=epsilon*x;
Num=10000;

%% Numerics
C_total = zeros(Num,M*(N+1));
xi=zeros(1,N+2);
xi(1)=0;
xi(N+2)=epsilon^(-1);
for eta = 1:Num
    r = epsilon^(-1)*rand(N,1);
    s = sortrows(r);
    for l=1:N
        xi(l+1) = s(l);
    end
    
    H=zeros(1,N+1);
    C=ones(1,M*(N+1));
    for i=1:M*(N+1)
        for j=1:N
            if x(i)>xi(j+1)
                H(j)=1;
            else
                break
            end
        end
        for j=1:N
            C(i)=C(i)-S*H(j)*(1+S)^(-j);
        end
    end
    C_total(eta,:) = C;
end
C_min=min(C)-0.01;

C_mean = mean(C_total);
C_var = var(C_total);
C_std = std(C_total);
C_minus = C_mean - 1.96*C_std;
C_plus =  C_mean + 1.96*C_std;
C_median = median(C_total);
% C_lower = credible_interval(C_total,0.975);
% C_upper = credible_interval(C_total,0.025);

%% Mean

Mean=ones(1,M*(N+1));
for i=1:M*(N+1)
    for j=1:1:N
        Mean(i) = Mean(i) - S*(1+S)^(-j)*betainc(epsilon*x(i),j,N-j+1);
    end
end

%% Variance
Var=zeros(1,M*(N+1));
for i=1:M*(N+1)
    for j=1:N
        for k=1:N
            Var(i) = Var(i) - S^2*(1+S)^(-j)*betainc(epsilon*x(i),j,N-j+1)*(1+S)^(-k)*betainc(epsilon*x(i),k,N-k+1);
        end
        Var(i) = Var(i) + S*(2 - (2+S)*(1+S)^(-j))*(1+S)^(-j)*betainc(epsilon*x(i),j,N-j+1);
    end
end

Var_leading = epsilon*(1-epsilon)*S^2*(x.*(epsilon^(-1)-x)).*exp(-2*S*x);

%% CDF
CDF=zeros(M*(N+1),N);
for i=1:M*(N+1)
    for j=1:N+1
        CDF(i,j)=1 - betacdf(epsilon*x(i),j,N-j+1);
    end
    CDF(i,N+1)=1;
end

C_0025=ones(M*(N+1),1);
r=0.025;
for i=1:M*(N+1)
    for j=1:N
        C_0025(i)=C_0025(i)-S*heaviside(x(i) - epsilon^(-1)*betaincinv(r,j,N-j+1))*(1+S)^(-j);
    end
end
C_05=ones(M*(N+1),1);
r=0.5;
for i=1:M*(N+1)
    for j=1:N
        C_05(i)=C_05(i)-S*heaviside(x(i) - epsilon^(-1)*betaincinv(r,j,N-j+1))*(1+S)^(-j);
    end
end
C_0975=ones(M*(N+1),1);
r=0.975;
for i=1:M*(N+1)
    for j=1:N
        C_0975(i)=C_0975(i)-S*heaviside(x(i) - epsilon^(-1)*betaincinv(r,j,N-j+1))*(1+S)^(-j);
    end
end

%% Realisations plot
figure()
hold on
% for i=1:Num
for i=1:Num
    p(i)=plot(epsilon*x,C_total(i,:),'k');
    p(i).Color(4) = 0.02;
end

p3=plot(epsilon*x,Mean-1.96.*sqrt(Var),'-bd','LineWidth',1.25,'MarkerSize',6,'MarkerIndices',[4*M,8*M,12*M,16*M]);
p3=plot(epsilon*x,Mean+1.96.*sqrt(Var),'-bd','LineWidth',1.25,'MarkerSize',6,'MarkerIndices',[4*M,8*M,12*M,16*M]);
p2=plot(epsilon*x,Mean,':rd','LineWidth',1.25,'MarkerSize',6,'MarkerIndices',[2*M,6*M,10*M,14*M,18*M]);

p5=plot(epsilon*x,C_0025,'-c*','LineWidth',1.25,'MarkerSize',7,'MarkerIndices',[2*M,6*M,10*M,14*M,18*M]);
p5=plot(epsilon*x,C_0975,'-c*','LineWidth',1.25,'MarkerSize',7,'MarkerIndices',[2*M,6*M,10*M,14*M,18*M]);
p4=plot(epsilon*x,C_05,':g*','LineWidth',1.25,'MarkerSize',7,'MarkerIndices',[4*M,8*M,12*M,16*M]);

plot(epsilon*x,C_total(Num,:),'m','LineWidth',1.25)

xlabel('$x_1$','FontSize',22,'Interpreter','latex');
ylabel('$C(x_1)$','FontSize',22,'Interpreter','latex');
leg=legend([p2,p3,p4,p5],'$E[C(x_1)]$','$E[C(x_1)] \pm 1.96 \sqrt{\mathrm{Var}[C(x_1)]}$','$CI(x_1;0.5)$','$CI(x_1;0.5\pm0.475)$');
set(leg,'FontSize',13,'Interpreter','latex','location','southwest');
