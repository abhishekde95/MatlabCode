% Email: jpweller@u.washington.edu
% Homework #4
% Question 2

clear all
close all


%% Part 1

% Setting some variables
A = [6 5 -5; 2 6 -2; 2 5 -1];

% Save the matrix
save A.dat A -ascii

g = 0;
for i = 1:size(A,1)
    for j = 1:size(A,2)
        g = g + (-1).^(i+j) * A(i,j);
    end
end

save gA.dat g -ascii


%% Part 2

% Set random state
rand('state',1)

% Draw random numbers
m = randi(100,1);
n = randi(100,1);
B = rand(m,n);

% save files
save B.dat B -ascii

% Perform function 'g' on B
g = 0;
for i = 1:size(B,1)
    for j = 1:size(B,2)
        g = g + (-1).^(i+j) * B(i,j);
    end
end

% Save results
save gB.dat g -ascii


