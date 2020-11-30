% Email: jpweller@u.washington.edu
% Homework #5

clear all
close all

%% Question 1
% Parts a and b

A = [2 -1 0; 1 6 -2; 4 -3 8];

R = zeros(size(A));
Q = zeros(size(A));

% Gram-Schmidt
for n = 1:size(A,2)
    tempRQ = zeros(size(A,1),1);
    for k = 1:n
        if k < n
            R(k,n) = A(:,n)' * Q(:,k);
            tempRQ = tempRQ - R(k,n) * Q(:,k);
        elseif k == n
            R(n,n) = norm(A(:,n) + tempRQ,2);
        end
    end
    Q(:,n) = (1/R(n,n)) * (A(:,n) + tempRQ);
end        


save A1.dat A -ascii
save QA.dat Q -ascii
save RA.dat R -ascii


%% Part c

rand('state',3)
%B = rand(111,23);
B = [1 1 7; 2 1 8; 3 1 9]

R = zeros(size(B));
Q = zeros(size(B));

% Gram-Schmidt
for n = 1:size(B,2)
    tempRQ = zeros(size(B,1),1);
    for k = 1:n
        if k < n
            R(k,n) = B(:,n)' * Q(:,k);
            tempRQ = tempRQ - R(k,n) * Q(:,k);
        elseif k == n
            R(n,n) = norm(B(:,n) + tempRQ,2);
        end
    end
    Q(:,n) = (1/R(n,n)) * (B(:,n) + tempRQ);
end

R = R(1:size(R,2),:);

save B.dat B -ascii
save QB.dat Q -ascii
save RB.dat R -ascii

