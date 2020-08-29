n=100; % number of individuals
rho=8;% dist variable, change here
N=1000; %number of time steps
lambda=1; % infection rate
mu=8;
time=1.5; % time for deterministic models to run for

SIM=100; % number of simulations
Y=zeros(SIM,N+1); % matrix for number of infectives
T=NaN(SIM,N+1); % matrix of times
T(:,1)=0;

for sim=1:SIM
    
Degrees = poissrnd(rho,1,n); %change here
L=sum(Degrees);
A=cumsum(Degrees);
HE=zeros(1,L);
HE(1:A(1))=1;
for k=(2:n)
    HE(A(k-1)+1:A(k))=k;
end
FE=HE(randperm(L));
a=FE(1:floor(L/2));
b=FE((floor(L/2)+1):2*floor(L/2));
c=ones(1,floor(L/2));
M=sparse(a,b,c,n,n);
U = M + transpose(M);
U(U > 1) = 1;
D=diag(diag(U));
U=U-D;
%G=plot(graph(U));



    initinf=datasample([1:n],10); % number of initial infectious
    Y(:,1)=length(initinf);
    y=zeros(1,n); % vector of infectives
    y(initinf)=1; % location of initial infected
   for t=1:N
    noI=sum(y); % number of infected
    S=find(~y); % indices of subseptibles 
    I=find(y); % indices of infected
    Unew=U(I,S); % matrix of subseptible nodes attached to an infected node
    SI=sum(sum(Unew)); % Number of SI pairs
    q=lambda*SI+mu*noI; % total rate
    X=rand; % random U~(0,1) sample
    if q*X > lambda*SI
        idy = datasample([1:noI],1); 
        y(I(idy))=0;
    else
        [k,l]=find(Unew);
        idx = datasample([1:length(k)],1); % random non zero element of Unew
        % transform back to indices of U
        z=S(l(idx));
        % Now we have a random SI pair
        y(z)=1;
    end
    R = exprnd(1/q);
    y;
    Y(sim,t+1)=sum(y);
    T(sim,t+1:N+1)= T(sim,t)+R;
    
    if sum(y)==0
        break
    end
   end
end 

%for sim=1:SIM
%plot(T(sim,:),Y(sim,:),'color',[.6 .6 .6])
%hold on
%end

xlabel('Time')
ylabel('I(t)')

Tnew=T;
Tnew(isnan(Tnew))=0;
Tnewer=zeros(1,N);
Ymean=zeros(1,N+1);
for i=1:N+1
    Ymean(i)=sum(Y(:,i))/SIM;
end

for i=1:N+1
    Tnewer(i)=sum(Tnew(:,i))/SIM;
end

[Tdet,~,Idet]=lindquistEtalODE(15,n,'poi',rho,length(initinf)/n,lambda,mu,time);
[tp,~,Ip]=PastorODE(15,'poi',rho,length(initinf)/n,lambda,mu,time);
plot(Tnewer,Ymean,'LineWidth', 3,'Color', 'blue') % Configuration model
hold on
%plot(tp,Ip*n,'--','LineWidth', 3 ,'Color', [.5 1 0]) % Pastor-Satoras model
%plot(Tdet,Idet,'--','LineWidth', 3,'Color', [1 .5 0]) % Lindquist model


ylim([0 100])
xlabel('Time')
ylabel('I(t)')
hold off
%clc
%clear