% This code is for Neubeh 528 homework due 5/20/11

clear all
close all

% Question 1

v = linspace(-3,3,20);
w = linspace(-3,3,20);
I = .1;
epsilon = 0.1; 
b0 = 0.825; 
b1 = .5;
dv = nan(length(v),length(w));
dw = dv;

for i = 1:length(v)
    for j = 1:length(w)
        dv(i,j) = v(i) - (v(i)^3)/3 -w(j) + I;
        dw(i,j) = epsilon*(b0 + b1*v(i) - w(i));
    end
end

figure(1);clf;hold on;
axis([-3 3 -3 3])
quiver(v,w,dv,dw)
plot(v, v - v.^3/3,'k');
plot(v, b0+ b1*v,'k');

% Try to get the voltage trace
nv(1) = 0;
nw(1) = 0;
for t = 2:250
    dv = nv(t-1) - (nv(t-1)^3)/3 -nw(t-1) + I;
    dw = epsilon*(b0 + b1*nv(t-1) - nw(t-1));
    nv(t) = (nv(t-1) + dv);
    nw(t) = (nw(t-1) + dw);
    figure(1);hold on;drawnow;
    plot(nv,nw)
    %keyboard
end

figure(2);clf;hold on;grid on;
xlabel('Time (ms)')
ylabel('Voltage (mv)')
plot(1:t,nv)

% fr_I = [0 4 12 32 36 36 36 44 40 44 48]; *These values were determined by
% %varying I and counting the voltage spikes in a 250ms window, then multiplying by 4.
% figure(3);clf,hold on; grid on; 
% title('Firing Rate vs Injected Current')
% xlabel('Injected Current')
% ylabel('Firing Rate (sp/s')
% plot(0:.1:1,fr_I)

%% Question 2

% Try to get the v,t line
ve(1) = 0;
vi(1) = 0;
Mee = 1.25;
Mie = 1;
Mii = 0;
Mei = -1;
gammae = -10;
gammai = 10;
te = 10;
%ti = 50;

for ti = [30 35 36 37 38 39 40 50]
    for t = 2:1000
        dve = (-ve(t-1) + (Mee*ve(t-1) + Mei*vi(t-1) - gammae))/te
        dvi = (-vi(t-1) + (Mii*vi(t-1) + Mie*ve(t-1) - gammai))/ti
        ve(t) = (ve(t-1) + dve);
        vi(t) = (vi(t-1) + dvi);
    end
    figure(ti)
    subplot(2,2,1:2);hold on;grid on;
    title(['Ti = ',num2str(ti),'(ms)'])
    xlabel('Ve (Hz)')
    ylabel('Vi (Hz)')
    plot(ve,vi)
    subplot(2,2,3);hold on;grid on;
    xlabel('Time (ms)')
    ylabel('Ve (Hz)')
    plot(1:t,ve)
    subplot(2,2,4);hold on;grid on;
    xlabel('Time (ms)')
    ylabel('Vi (Hz)')
    plot(1:t,vi)
end