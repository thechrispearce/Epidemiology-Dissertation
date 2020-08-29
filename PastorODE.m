function [t,X,I]=PastorODE(M,typeD,paramD,epsilon,beta,gamma,tMax)

pmf = @(x)probsD(typeD,paramD,x,0);

%PastorODE(M,N,typeD,paramD,epsilon,beta,gamma,tMax)

%Here M is the max degree in Pastor's ODEs, N is the population size, 
%type D and paramD specify the degree distribution (so use 'poi' and 5 
%for a Poisson(5) degree distribution), epsilon is the proportion 
%initially infected, beta is the infection rate, gamma the recovery rate 
%and tMax the time to integrate until.

%numerically integrate the Pastor ODEs from time 0 to tMax, keeping
%track of degrees 1,2,...,M, for a system of N individuals, degree
%distribution specified by typeD and paramD, a proportion epsilon initially
%infected (chosen uniformly) infection rate beta, recovery rate gamma.

tspan = [0,tMax];

X0 = PastorInitODE(M,typeD,paramD,epsilon);


[t,X] = ode45(@(t,y)PastorRHS(y,beta,gamma,M,typeD,paramD),tspan,X0);

I=zeros(length(t),1);

for k=0:M
    I=I+X(:,k+1)*pmf(k);
end



%figure
%plot(t,I,'-k')
%xlabel('t')
%ylabel('I(t)')
