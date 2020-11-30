% Email: jpweller@u.washington.edu
% Homework #7

clear all
close all

%% Question 4

A = [1 1; 1 2; 1 3; 1 5];
b = [1 0 -2 -3]';

x = A\b;

t = 0:.01:6;
y = x(1) + x(2)*t;

figure(1); clf; hold on; grid on;
plot(A(:,2),b,'r*')
plot(t,y,'k--')
title('Question 4')
xlabel('t')
ylabel('y')

