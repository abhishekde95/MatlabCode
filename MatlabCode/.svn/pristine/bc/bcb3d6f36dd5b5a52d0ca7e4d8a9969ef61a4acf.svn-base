close all; clearvars;
%% Model expectations 2-D (More applicable to spatial integration problems betweens 2 different parts of the space
x = linspace(-10,10,101);
y = linspace(-10,10,101);
[X,Y] = meshgrid(x,y);
Z1 = X+Y;
Z2 = X.^2 + Y.^2;
Z3 = (X + Y).^2;
Z4 = X.^2 - Y.^2;
figure(1);
subplot(221),contour(X,Y,Z1); xlabel('X'),ylabel('Y'); title('X+Y');
subplot(222),contour(X,Y,Z2); xlabel('X'),ylabel('Y'); title('X^2+Y^2');
subplot(223),contour(X,Y,Z3); xlabel('X'),ylabel('Y'); title('(X+Y)^2')
subplot(224),contour(X,Y,Z4); xlabel('X'),ylabel('Y'); title('X^2-Y^2')