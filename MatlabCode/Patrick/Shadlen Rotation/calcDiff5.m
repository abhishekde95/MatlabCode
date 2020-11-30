% calcDiff5 implements Fokker-Plank equation.
% Same bound hight A = B

function [t1,t2,p] = calcDiff5(cohs, theta)

% input> theta=[k,A,t1_res,t2_res]

k      = theta(1);
A      = theta(2);
t1_res = theta(3);
t2_res = theta(4);

for i = 1:length(cohs)
C = cohs(i);

if C ~= 0
    t1(i,1)=A./(k.*C).*tanh(k.*C.*A)+t1_res;
    t2(i,1)=A./(k.*C).*tanh(k.*C.*A)+t2_res;
else
    t1(i,1)=A.^2+t1_res;
    t2(i,1)=A.^2+t2_res;
end

p(i,1)=1./(1+exp(-2.*k.*C.*A));
end