% This code is for Steve Perlmutter's neubeh 503 hw due 5/26
% Modified May_2011 JPW

clear all
close all

%% Quesiton 1a
disp('Qustion 1a')

% Set variables
k  = [1.2 -1.1 -1.5];
A  = [25 50 40];
mu = [0 -20 10];
sig = 10;
meridian = -40:40;
eyePos   = [-20:10:20];
f = nan(numel(meridian),numel(eyePos),numel(k));

% Loop through 3 cells with different x's and e's
for n = 1:length(k)
    for e = 1:numel(eyePos)
        for x = 1:numel(meridian)
            f(x,e,n) = (k(n)*eyePos(e)+A(n))*exp(-(meridian(x)-mu(n))^2/(2*sig^2));
        end
    end
end

% Plot results
figure(1);clf;hold on;
for n = 1:3
    subplot(3,1,n); hold on;grid on;
    axis([-40 40 0 max(max(max(f)))*1.2])
    title(['Neuron ',num2str(n)])
    if n == 3
        xlabel('Position of Stimulus Along Retina (degrees)')
    end
    ylabel('Firing Rate (sp/s)')
    plot(meridian,f(:,:,n))
    legend('-20','-10','0','10','20','Location','EastOutside');
    if n ==1
        hl = findobj(gcf,'Tag','legend');
        hl_t = get(hl,'Title');
        set(hl_t,'String','Eye Position (degrees)');
    end
end


%% Quesiton 1c
disp('Question 1c')

meridian_h = -60:60;
f_h = nan(numel(meridian_h),numel(eyePos),numel(k));

for n = 1:length(k)
    for e = 1:numel(eyePos)
        meridian_e = meridian_h-eyePos(e); 
        for x = 1:numel(meridian_e)
            f_h(x,e,n) = (k(n)*eyePos(e)+A(n))*exp(-(meridian_e(x)-mu(n))^2/(2*sig^2));
        end
    end
end


figure(2);clf;hold on;grid on;
for n = 1:3
        subplot(3,1,n);hold on;grid on;
        axis([min(meridian_h) max(meridian_h) 0 max(max(max(f)))*1.2])
        title(['Neuron ',num2str(n)])
        if n == 3
            xlabel('Position of Stimulus in Head-Centered Coordinates (degrees)')
        end
        ylabel('Firing Rate (sp/s)')
        plot(meridian_h,f_h(:,:,n))
        legend('-20','-10','0','10','20','Location','EastOutside');
        if n ==1
            hl = findobj(gcf,'Tag','legend');
            hl_t = get(hl,'Title');
            set(hl_t,'String','Eye Position (degrees)');
        end
end



%% Question 1d
disp('Question 1d')

fRates1010 = f(meridian==10,eyePos==10,:)


%% Question 1e
disp('Question 1e')

n1 = f_h(:,:,1)<1; 
n2 = f_h(:,:,2)<50 & f_h(:,:,2)>40; 
n3 = f_h(:,:,3)<50 & f_h(:,:,3)>40;

possLoc = n1==1 & n2==1 & n3==1    



