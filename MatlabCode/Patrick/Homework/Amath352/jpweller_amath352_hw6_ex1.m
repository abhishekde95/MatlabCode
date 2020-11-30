% Email: jpweller@u.washington.edu
% Homework #6

clear all
close all

%% Question 1

rand('state',2)
A = rand(5,5);

L = A;
U = eye(5,5);

save l1.dat L -ascii


for k = 1:size(A,2)
    for i = k+1:size(A,2)
        U(k,i) = L(k,i)/L(k,k);
        for j = k:size(A,2)
            L(j,i) = L(j,i) - U(k,i)*L(j,k);
        end
    end
    if k == 1
        save l1.dat L -ascii
    end
end


save lfinal.dat L -ascii
save ufinal.dat U -ascii

