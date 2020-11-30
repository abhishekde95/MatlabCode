% Email: jpweller@u.washington.edu
% Homework #5

clear all
close all

%% Question 2
% Part a

a = 6;
b = 5;
c = -4;

xplus_1 = (-b + sqrt(b.^2 - 4*a*c)) / (2 * a)
xminus_1 = (-b - sqrt(b.^2 - 4*a*c)) / (2 * a)

save xplus_1.dat xplus_1 -ascii
save xminus_1.dat xminus_1 -ascii


%% Part b
a = 6e154;
b = 5e154;
c = -4e154;

xplus_2 = (-b + b * sqrt(1 - 4*(a/b)*(c/b))) / (2 * a)
xminus_2 = (-b - b * sqrt(1 - 4*(a/b)*(c/b))) / (2 * a)

save xplus_2.dat xplus_2 -ascii
save xminus_2.dat xminus_2 -ascii


%% Part c
a = 1;
b = -4;
c = 3.999999;

xplus_3 = (-b + sqrt(b.^2 - 4*a*c)) / (2 * a)
xminus_3 = (-b - sqrt(b.^2 - 4*a*c)) / (2 * a)

save xplus_3.dat xplus_3 -ascii
save xminus_3.dat xminus_3 -ascii

