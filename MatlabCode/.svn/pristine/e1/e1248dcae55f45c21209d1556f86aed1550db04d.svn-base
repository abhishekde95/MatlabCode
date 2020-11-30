% Email: jpweller@u.washington.edu
clear all

% Homework 1
% Exercies 3-7

%% Exercise 3

% Independent variable
x = 0:.05:1;

% Dependent variables
f = cos(2*pi*x);
g = (exp(x) + exp(-x))./2;
h = (1-x)./(1+x);

% Plot
figure(1); clf; grid on; hold on;
plot(x,f,'--')
plot(x,g,'g*--')
plot(x,h,'ro--')
legend('f(x)','g(x)','h(x)','Location','NorthWest')
xlabel('x value (independent variable)')
ylabel('y value (dependent variables)')
title('Exercise #3')


%% Exercise 4

% a)
coeffs = [6 1; 3 -2];
answers = [5; 5];

x_y_vals_4a = coeffs \ answers

% For thoroughness, I will type out the solution to this first problem
% long-hand:

% 6x + y  = 5
% 3x - 2y = 5; <- Multiply this by 2, then subtract.

%  6x + y  = 5
% -6x + 4y = -10
% ---------------
%       5y = -5
%        y = -1 <- Substitute into original equation

% 6x + (-1) = 5
% 6x = 6
% x = 1

% b)
coeffs = [2 -1 2; -1 -1 3; 3 0 -2];
answers = [2; 1; 1];

x_y_vals_4b = coeffs \ answers

% c)
coeffs = [1 1 -1; 2 -1 3; -1 -1 0];
answers = [0; 3; 6];

x_y_vals_4c = coeffs \ answers


%% Exercise 5

% a)
coeffs = [1 2; -3 -4];
answers = [1; 4];

c_d_vals_5a = coeffs \ answers

% b)
coeffs = [1 -1; 1 2];
answers = [7; 3];

c_d_vals_5b = coeffs \ answers

% c)
coeffs = [2 1 2; -1 3 3; 4 -3 0];
answers = [3; -2; 7];

c_d_e_vals_5c = coeffs \ answers


%% Exercise 6

% a)
a_6a = (2*-1 + -3*-2)/3

% b)
a_6b = (1*1 + 1*-1)/1


%% Exercise 7

v1 = [1 0 -1]';
v2 = [-1 1 0]';

onenorm_v1 = norm(v1,1)
twonorm_v1 = norm(v1,2)
infnorm_v1 = norm(v1,Inf)

onenorm_v2 = norm(v2,1)
twonorm_v2 = norm(v2,2)
infnorm_v2 = norm(v2,Inf)

