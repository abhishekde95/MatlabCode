% Email: jpweller@u.washington.edu
% Homework #2

clear all
close all

%% Exercise 2a

s20 = 0;
for n = 1:20
    s20 = s20 + 1/n^2;
end
%s20 = sqrt(s20 * 6);
save s20.dat s20 -ascii


%% Exercise 2b

e20 = abs(pi - sqrt(s20 * 6));
save e20.dat e20 -ascii


%% Exersise 2c

n04 = 0; n = 0; error = pi;
while error > 1e-4
    n = n+1;
    n04 = n04 + 1/n^2;
    error = abs(pi - sqrt(6*n04));
    if n == 100000
        disp('Danger! Potential infinite loop! Call an adult!')
        return
    end
end

save n04.dat n -ascii

