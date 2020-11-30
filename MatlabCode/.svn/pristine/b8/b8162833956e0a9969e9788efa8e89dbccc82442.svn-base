function [plot_counter] = Selectpixel(plot_counter,STAs)
global maskidx spikeidx nstixperside seedidx nframesidx stimonidx muidxs sigmaidxs
global msperframe maxT xx yy M

figure(plot_counter-1)
[a,b] = ginput(1); % y-rows, x-columns
a = round(a);
b = round(b);
fprintf('You have selected Row-%1d and Column-%1d \n',b,a);
RGB_temp_plot = zeros(3,maxT);
LMS_temp_plot = zeros(3,maxT);
for m=1:maxT
    tmp_var = reshape(STAs(:,m),[nstixperside nstixperside 3]);
    RGB_temp_plot(:,m) = tmp_var(b,a,:);
end
LMS_temp_plot = M*RGB_temp_plot;

figure(plot_counter); subplot(2,1,1);
plot(RGB_temp_plot(1,:),'r','LineWidth',2); hold on;
plot(RGB_temp_plot(2,:),'g','LineWidth',2); hold on;
plot(RGB_temp_plot(3,:),'LineWidth',2); hold on;
xlabel('frames'),ylabel('Intensity');legend('R','G','B'); hold off;

subplot(2,1,2);
plot(LMS_temp_plot(1,:),'r','LineWidth',2); hold on;
plot(LMS_temp_plot(2,:),'g','LineWidth',2); hold on;
plot(LMS_temp_plot(3,:),'LineWidth',2); hold on;
xlabel('frames'),ylabel('Intensity');legend('L','M','S'); hold off;

prompt = 'Do u want to analyze another pixel? (Y/N)';
str = input(prompt,'s');
if ((strcmp('Y',str) || strcmp('y',str)))
    eval('Selectpixel(plot_counter,STAs)');
end

end





