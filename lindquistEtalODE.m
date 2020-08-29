function [t,X,I]=lindquistEtalODE(M,N,typeD,paramD,epsilon,beta,gamma,tMax)



%lindquistEtalODE(M,N,typeD,paramD,epsilon,beta,gamma,tMax)

%Here M is the max degree in Lindquist's ODEs, N is the population size, 
%type D and paramD specify the degree distribution (so use 'poi' and 5 
%for a Poisson(5) degree distribution), epsilon is the proportion 
%initially infected, beta is the infection rate, gamma the recovery rate 
%and tMax the time to integrate until.

%numerically integrate the Lindquist etal ODEs from time 0 to tMax, keeping
%track of degrees 1,2,...,M, for a system of N individuals, degree
%distribution specified by typeD and paramD, a proportion epsilon initially
%infected (chosen uniformly) infection rate beta, recovery rate gamma.

tspan = [0,tMax];

X0 = lindquistEtalInitODE(M,typeD,paramD,epsilon)*N;


[t,X] = ode45(@(t,y)lindquistEtalRHS(y,beta,gamma,M),tspan,X0);


I = zeros(length(t),1);
I(:,1) = sum(X(:,((M+1)*(M+2)/2)+1:end),2);

%figure
%plot(t,I,'-k')
%xlabel('t')
%ylabel('I(t)')
