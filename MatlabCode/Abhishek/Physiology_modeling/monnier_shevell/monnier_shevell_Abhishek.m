%% Author - Abhishek
% Hacking around with the Monnier and Shevell 2004 model, antagonistic center surround receptive field

%****************************************
% Defining Parameters 
%****************************************
close all;
clearvars;
pixperhalfcycle = 10;
pixmargin = 50;
ncycles = 5;
sigma_c = 15;
sigma_s = 2*sigma_c;

%****************************************
% Creating center-surround RF using DOG
%****************************************
z = [-3*sigma_s:1:3*sigma_s];
y_center = normpdf(z,0,sigma_c);
y_surround = normpdf(z,0,sigma_s);
RF = y_center - y_surround;
figure(2), subplot(339); plot(z,y_center,'g','Linewidth',2); hold on;
plot(z,-1*y_surround,'r','Linewidth',2);
plot(z,RF,'k','Linewidth',2); title('DOG RF'); ylabel('Weights'); 
xlabel('Positions in pixels'); legend('center','surround','RF');hold off;

%****************************************
% Patterned background #1
%****************************************
x1 = [1 -1];
x1 = repmat(x1, pixperhalfcycle, ncycles);
x1 = x1(:);
x1 = [zeros(pixmargin,1); x1; zeros(pixmargin,1)];
tmp_wotest1 = conv(x1,RF,'same');
x1(91:100) = 0;
tmp_wtest1 = conv(x1,RF,'same');
diff1 = tmp_wtest1 - tmp_wotest1;
figure(2), subplot(331); plot(x1,'b-','Linewidth',2); hold on;
% plot(tmp_wotest1,'g-','Linewidth',2);
plot(tmp_wtest1,'k-','Linewidth',2); title('Patterened background #1'); ylim([-1 1]); hold off;

%****************************************
% Patterned background #2
%****************************************
x2 = [-1 1];
x2 = repmat(x2, pixperhalfcycle, ncycles);
x2 = x2(:);
x2 = [zeros(pixmargin,1); x2; zeros(pixmargin,1)];
tmp_wotest2 = conv(x2,RF,'same');
x2(91:100) = 0;
tmp_wtest2 = conv(x2,RF,'same');
diff2 = tmp_wtest2 - tmp_wotest2;
figure(2), subplot(332); plot(x2,'b-','Linewidth',2); hold on;
% plot(tmp_wotest2,'g-','Linewidth',2);
plot(tmp_wtest2,'k-','Linewidth',2); title('Patterned background #2'); ylim([-1 1]);hold off;

%****************************************
% Uniform background #1
%****************************************
x3 = [1 1];
x3 = repmat(x3, pixperhalfcycle, ncycles);
x3 = x3(:);
x3 = [zeros(pixmargin,1); x3; zeros(pixmargin,1)];
tmp_wotest3 = conv(x3,RF,'same');
x3(91:100) = 0;
tmp_wtest3 = conv(x3,RF,'same');
diff3 = tmp_wtest3 - tmp_wotest3;
figure(2), subplot(334); plot(x3,'b-','Linewidth',2); hold on;
% plot(tmp_wotest3,'g-','Linewidth',2);
plot(tmp_wtest3,'k-','Linewidth',2); title('Uniform background #1'); ylim([-1 1]);hold off;

%****************************************
% Uniform background #2
%****************************************
x4 = [-1 -1];
x4 = repmat(x4, pixperhalfcycle, ncycles);
x4 = x4(:);
x4 = [zeros(pixmargin,1); x4; zeros(pixmargin,1)];
tmp_wotest4 = conv(x4,RF,'same');
x4(91:100) = 0;
tmp_wtest4 = conv(x4,RF,'same');
diff4 = tmp_wtest4 - tmp_wotest4;
figure(2), subplot(335); plot(x4,'b-','Linewidth',2); hold on;
% plot(tmp_wotest4,'g-','Linewidth',2);
plot(tmp_wtest4,'k-','Linewidth',2); title('Uniform background #2'); ylim([-1 1]);hold off;

%****************************************
% Plotting the color shift
%****************************************
figure(2),subplot(337); plot(tmp_wtest1-tmp_wtest3,'r','Linewidth',2); title('Patterned #1 - Uniform #1');
subplot(338); plot(tmp_wtest2-tmp_wtest4,'r','Linewidth',2);  title('Patterned #2 - Uniform #2');

figure(2),subplot(333); plot(tmp_wtest1-tmp_wtest2,'m','Linewidth',2); title('Patterned #1 - Patterened #2');
subplot(336); plot(tmp_wtest3-tmp_wtest4,'m','Linewidth',2);  title('Uniform #1 - Uniform #2');