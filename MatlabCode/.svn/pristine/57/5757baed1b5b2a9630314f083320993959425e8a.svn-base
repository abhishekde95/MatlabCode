% Email: jpweller@u.washington.edu
% Homework #4

clear all
close all


%% Question 1

% Set some parameters
A = [2 -1 0; 1 6 -2; 4 -3 8];
b = [2; -4; 5];
x = zeros(1,size(A,1));
eta = 1e-8;
r(1,:) = b-A*x';

% Gauss-Seidel loop
for k = 1:1000
    for i = 1:size(A,1)
        j = 1;
        tempsum1 = 0;
        tempsum2 = 0;
        if j<i
            for j = 1:(i-1)
                tempsum1 = tempsum1 + A(i,j) * x(k+1,j);
            end
        end
        for j = (i+1):size(A,2)
            tempsum2 = tempsum2 + A(i,j) * x(k,j);
        end
        x(k+1,i) = (1/A(i,i)) * (b(i) - tempsum1 - tempsum2);
    end
    r(k+1,:) = b-A*x(k+1,:)';
    if sqrt(sum(r(k+1,:).^2)) < (eta * sqrt(sum(b.^2)))
        disp('Did it!')
        break
    end
end

% Redefining variables for saving ease
x1 = x(2,:); % Using 2nd row b/c x(1,:) = x0 (indexing in Matlab...)
xfinal = x(end,:);
countR = size(r,1);

% Save relevant variables
save A1.dat A -ascii
save x1.dat x1 -ascii
save xfinal.dat xfinal -ascii
save countR.dat countR -ascii
            
