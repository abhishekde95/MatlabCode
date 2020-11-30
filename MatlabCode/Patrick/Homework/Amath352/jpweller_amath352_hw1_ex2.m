% Email: jpweller@u.washington.edu
clear all

% Homework 1
% Exercise 2

TC = -50:5:50;
TF = nan(size(TC));

for n = 1:numel(-50:5:50)
    TF(n) = (9/5) * TC(n) + 32;
end

TC = TC';
TF = TF';

save TC.dat TC -ascii
save TF.dat TF -ascii