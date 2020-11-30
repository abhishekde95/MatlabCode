% Email: jpweller@u.washington.edu
clear all

% Exercise 1

A1.dat = 1;
for n = 3:4:1003
    A1.dat = A1.dat - 1/n + 1/(n+2);
end
save A1.dat


A2.dat = 0;
for n = 1:2:1000
    A2.dat = A2.dat + 1/(n.^2 + (n+2).^2);
end
save A2.dat