% For first TCM

close all
clearvars

% for subunits
A = 0.5 * ones(100,100,3);
vec1 = A; vec2 = A; vec1inv = A; vec2inv = A;
vec1(20:50,30:70,2) = 0.3; vec1inv(20:50,30:70,2) = 0.7;
vec2(50:80,30:70,2) = 0.7; vec2inv(50:80,30:70,2) = 0.3;
% figure(1),image(vec1), set(gca,'XTick',[],'YTick',[]);
% figure(2),image(vec2), set(gca,'XTick',[],'YTick',[]);
% figure(3),image(vec1inv), set(gca,'XTick',[],'YTick',[]);
% figure(4),image(vec2inv), set(gca,'XTick',[],'YTick',[]);

% for multi - mechanisms
B = 0.5 * ones(100,100,3);
vec3 = B; vec3inv = B;

vec3(20:50,30:70,2) = 0.45; vec3(50:80,30:70,2) = 0.7;
vec3inv(20:50,30:70,2) = 0.3; vec3inv(50:80,30:70,2) = 0.55;
% figure(5),image(vec3), set(gca,'XTick',[],'YTick',[]);
% figure(6),image(vec3inv), set(gca,'XTick',[],'YTick',[]);


% For isoresponse plot
vec41 = B; vec41(20:50,30:70,2) = 0.45;
vec42 = B; vec42(20:50,30:70,2) = 0.40;
vec43 = B; vec43(20:50,30:70,2) = 0.35;
vec44 = B; vec44(20:50,30:70,2) = 0.30;
% figure(7),image(vec41), set(gca,'XTick',[],'YTick',[]);
% figure(8),image(vec42), set(gca,'XTick',[],'YTick',[]);
% figure(9),image(vec43), set(gca,'XTick',[],'YTick',[]);
% figure(10),image(vec44), set(gca,'XTick',[],'YTick',[]);

vec51 = B; vec51(50:80,30:70,2) = 0.55;
vec52 = B; vec52(50:80,30:70,2) = 0.60;
vec53 = B; vec53(50:80,30:70,2) = 0.65;
vec54 = B; vec54(50:80,30:70,2) = 0.70;
% figure(11),image(vec51), set(gca,'XTick',[],'YTick',[]);
% figure(12),image(vec52), set(gca,'XTick',[],'YTick',[]);
% figure(13),image(vec53), set(gca,'XTick',[],'YTick',[]);
% figure(14),image(vec54), set(gca,'XTick',[],'YTick',[]);

vec61 = B; vec61(20:50,30:70,2) = 0.45; vec61(50:80,30:70,2) = 0.55;
vec62 = B; vec62(20:50,30:70,2) = 0.40; vec62(50:80,30:70,2) = 0.60;
vec63 = B; vec63(20:50,30:70,2) = 0.35; vec63(50:80,30:70,2) = 0.65;
vec64 = B; vec64(20:50,30:70,2) = 0.30; vec64(50:80,30:70,2) = 0.70;
figure(15),image(vec61), set(gca,'XTick',[],'YTick',[]);
figure(16),image(vec62), set(gca,'XTick',[],'YTick',[]);
figure(17),image(vec63), set(gca,'XTick',[],'YTick',[]);
figure(18),image(vec64), set(gca,'XTick',[],'YTick',[]);