% Email: jpweller@u.washington.edu
% Homework #7

clear all
close all

%% Question 2
% Part A

y(1) = 1;
y(2) = 1;

for n = 2:49
    y(n+1) = y(n-1) + y(n);
end

y50 = y(50);

save y50.dat y50 -ascii


%% Part B

z = log(y)';

A = cat(2,ones(size(z)),[1:50]');

c = A\z;
c1 = c(1);
c2 = c(2);

save c1.dat c1 -ascii
save c2.dat c2 -ascii


%% Part c


beta = exp(c2);

save beta.dat beta -ascii

