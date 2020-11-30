%% Author - Greg
% Hacking around with the Monnier and Shevell 2004 model, antagonistic
% center surround receptive field
 
pixperhalfcycle = 10;
pixmargin = 50;
ncycles = 5;
sigma_c = 12;
sigma_s = 2*sigma_c;
 
figure;
% Patterned background #1
x = [1 -1];
x = repmat(x, pixperhalfcycle, ncycles);
x = x(:);
x = [zeros(pixmargin,1); x; zeros(pixmargin,1)];
x(91:100) = 0;
 
y = normpdf([-3*sigma_s:1:3*sigma_s],0,sigma_c)-...
            normpdf([-3*sigma_s:1:3*sigma_s],0,sigma_s)
subplot(2,2,1); hold on;
plot(x,'b-')
plot(conv(x,y,'same'),'k-');
title('Patterned background #1')
 
% Now a uniform lime background
xx = x;
xx(pixmargin+1:91) = -1;
xx(110:end-pixmargin) = -1;
subplot(2,2,2); hold on;
plot(xx,'b-');
plot(conv(xx,y,'same'),'k-');
title('Lime background');
 
% Now a uniform purple background
xx = x;
xx(pixmargin+1:91) = 1;
xx(100:end-pixmargin) = 1;
subplot(2,2,3); hold on;
plot(xx,'b-')
plot(conv(xx,y,'same'),'k-')
title('Purple background');
 
% Patterned background #2
x = [-1 1];
x = repmat(x, pixperhalfcycle, ncycles);
x = x(:);
x = [zeros(pixmargin,1); x; zeros(pixmargin,1)];
x(91:100) = 0;
 
y = normpdf([-3*sigma_s:1:3*sigma_s],0,sigma_c)-...
            normpdf([-3*sigma_s:1:3*sigma_s],0,sigma_s)
subplot(2,2,4); hold on;
plot(x,'b-')
plot(conv(x,y,'same'),'k-')
title('Patterned background #2')