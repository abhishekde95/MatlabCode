% Email: jpweller@u.washington.edu
% Homework #7

clear all
close all

%% Question 1
% One degree precision

t = [1900:10:2010]';

b = [75.995 91.972 105.711 123.203 131.669 150.697 179.323...
    203.212 226.505 249.633 281.422 308.745]';

A = cat(2,ones(size(t)),t);
[Q,R] = qr(A);
x = R\(Q\b);

y = x(1) + x(2)*t;

% Plot Results
% figure(1); hold on; grid on;
% plot(t,b,'r*')
% plot(t,y,'k--')
% xlabel('Time (years)')
% ylabel('Total US Population')
% title('Model of Population Growth')

pop2020_deg1 = x(1) + x(2)*2020;

save pop2020_deg1.dat pop2020_deg1 -ascii


%% Two degree precision

A = cat(2,ones(size(t)),t,t.^2);
[Q,R] = qr(A);
x = R\(Q\b);

y = x(1) + x(2)*t + x(3)*t.^2;

% Plot Results
% figure(2); hold on; grid on;
% plot(t,b,'r*')
% plot(t,y,'k--')
% xlabel('Time (years)')
% ylabel('Total US Population')
% title('Model of Population Growth')

pop2020_deg2 = x(1) + x(2)*2020 + x(3) * 2020.^2;


save pop2020_deg2.dat pop2020_deg2 -ascii