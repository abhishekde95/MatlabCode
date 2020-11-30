% Email: jpweller@u.washington.edu
% Homework #2

clear all
close all

%% Exercise 1a
x = [1; 2; 3; 5; 7];

onenorm_1 = norm(x,1);
twonorm_1 = norm(x,2);
fournorm_1 = norm(x,4);

save Norm_X_1.dat onenorm_1 -ascii
save Norm_X_2.dat twonorm_1 -ascii
save Norm_X_4.dat fournorm_1 -ascii


%% Exercise 1b

y = [1; 2; 4; 6; 8];

dot_xy = x' * y;

angle_xy = subspace(x,y)/pi*180;

save dot_xy.dat dot_xy -ascii
save angle_xy.dat angle_xy -ascii