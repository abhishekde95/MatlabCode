% Author - Abhishek De, 5/18
% for Thesis committee meeting 3
close all; clearvars;
x = linspace(0,2,51);
z = linspace(-1,2,101);
[X,Y,Z] = meshgrid(x,x,z);
Z1 = X.^2 + Y.^2 + Z.^2;
Z2 = X.^2 + Y.^2 - 5*Z.^2;
val1 = prctile(Z1(:),[25 50 75]);
val2 = prctile(Z2(:),[25 50 75]);
plot_counter = 1;
for ii = 1:3
    figure(plot_counter), isosurface(X,Y,Z,Z1,val1(ii)); hold on;
    figure(plot_counter+1), isosurface(X,Y,Z,Z2,val1(ii)); hold on;
end
figure(plot_counter), xlabel('R'), ylabel('G'), zlabel('B'), hold off;
figure(plot_counter+1), xlabel('R'), ylabel('G'), zlabel('B'), hold off;
plot_counter = plot_counter + 2;