% For first TCM

close all
clearvars

% for subunits
A = 0.5 * ones(100,100,3);
vec1 = A; vec2 = A; vec1inv = A; vec2inv = A;
vec1(20:50,30:70,2) = 0.3; vec1inv(20:50,30:70,2) = 0.7;
vec2(50:80,30:70,2) = 0.7; vec2inv(50:80,30:70,2) = 0.3;
figure(1),image(vec1), set(gca,'XTick',[],'YTick',[]);
figure(2),image(vec2), set(gca,'XTick',[],'YTick',[]);
figure(3),image(vec1inv), set(gca,'XTick',[],'YTick',[]);
figure(4),image(vec2inv), set(gca,'XTick',[],'YTick',[]);

%for multi - mechanisms
B = 0.5 * ones(100,100,3);
vec1 = B; vec2 = B; vec1inv = B; vec2inv = B;
vec1(20:50,30:70,3) = 0.3; vec1(50:80,30:70,3) = 0.7;
vec2(20:50,30:70,2) = 0.3; vec2(50:80,30:70,2) = 0.7;
vec1inv(20:50,30:70,3) = 0.7; vec1inv(50:80,30:70,3) = 0.3;
vec2inv(20:50,30:70,2) = 0.7; vec2inv(50:80,30:70,2) = 0.3;

figure(1),image(vec1), set(gca,'XTick',[],'YTick',[]);
figure(2),image(vec2), set(gca,'XTick',[],'YTick',[]);
figure(3),image(vec1inv), set(gca,'XTick',[],'YTick',[]);
figure(4),image(vec2inv), set(gca,'XTick',[],'YTick',[]);

