% Email: jpweller@u.washington.edu
% Homework #6

clear all
close all

%% Question 2

%Part 1
e = 1e-8;
A = [2 -1 -3; -1 6 5; -3 5 8];
itera = 0;

while max(max(max(abs(triu(A,1))),max(max(abs(tril(A,-1)))))) > e*max(diag(abs(A)))
    
    itera = itera+1;
    [Q,R] = qr(A);
    A = R*Q;

    if itera == 1000
        disp('Whoa there, Nelly! Possible infinite loop occurring in Part A...')
        break
    end
    
end

diaga = diag(A);

save diaga.dat diaga -ascii
save itera.dat itera -ascii


% Part 2
rand('state',1)
B = rand(6,6);
B = B+B';
save b.dat B -ascii

iterb = 0;
while max(max(max(abs(triu(B,1))),max(max(abs(tril(B,-1)))))) > e*max(diag(abs(B)))
    
    iterb = iterb+1;
    [Q,R] = qr(B);
    B = R*Q;

    if iterb == 1000
        disp('Whoa there, Nelly! Possible infinite loop occurring in Part B...')
        break
    end
    
end

diagb = diag(B);

save diagb.dat diagb -ascii
save iterb.dat iterb -ascii
    