% Email: jpweller@u.washington.edu
% Homework #5

clear all
close all

%% Question 3

x = 1.92:.001:2.08;

p1 = (x-2).^9;
p2 = x.^9 - 18*x.^8 + 144*x.^7 - 672*x.^6 + 2016*x.^5 ...
    - 4032*x.^4 + 5376*x.^3 - 4608*x.^2 + 2304*x - 512;
q1 = (x-2).^4;
q2 = x.^4 - 8*x.^3 + 24*x.^2 - 32*x + 16;


% Plot results
figure(1); clf; grid on; hold on;
title('P Expressions')
xlabel('x-values')
ylabel('p(x)')
plot(x,p1,'b--')
plot(x,p2,'g')
legend('Developed Expression','Reduced Expression')

figure(2); clf; grid on; hold on;
title('Q Expressions')
xlabel('x-values')
ylabel('q(x)')
plot(x,q1,'b-*')
plot(x,q2,'g-o')
legend('Developed Expression','Reduced Expression')




