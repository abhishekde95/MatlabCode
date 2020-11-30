% Email: jpweller@u.washington.edu
% Homework #2

clear all
close all


%% Exercise 3

% Give x1 a range (min and max values that x1 can take)
x1 = -2:.1:2;

% Compute x2 values
x2 = abs(3 * ((.5 * abs(x1)) -1));

% Sanity Check: Confirm that all norms = 1
unitcircle = .5 * abs(x1) + 1/3 * abs(x2)

% Plot unit circle
figure(1); clf; hold on; grid on;
plot(x1,x2)
plot(x1,x2.*-1)
axis equal
xlabel('x1 values')
ylabel('x2 values')
title('Homework 1 Exercise 3')
