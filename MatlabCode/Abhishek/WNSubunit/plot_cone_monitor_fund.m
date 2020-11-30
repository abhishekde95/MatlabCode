function plot_counter =  plot_cone_monitor_fund(stro,plot_counter)
% Plot the cone and the monitor color fundamentals
cone = reshape(stro.sum.exptParams.fundamentals,[81 3]);
mon = reshape(stro.sum.exptParams.mon_spd,[101 3]);
mon = spline([380:4:780], mon', [380:5:780]);

figure(plot_counter);
subplot(211), plot([380:5:780],fliplr(cone)); xlabel('Wavelength'); ylabel('Energy'); title('Absorption Spectra'); legend('S','M','L');
subplot(212), plot([380:5:780],fliplr(mon')); xlabel('Wavelength'); ylabel('Energy'); title('Emission Spectra'); legend('blue','green','red');
plot_counter = plot_counter + 1;
end

