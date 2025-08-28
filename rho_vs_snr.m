clear;
close;
home;
pkg load statistics
fontsize = 16;
markersize = 2;
signal_linewidth = 4;
noise_linewidth = 4;
results = load('all_results_2.dat');
snr = results(:,7);
rho = results(:,5);
rho_noise = results(:,8);
x_bin_delta = 0.02;
x_bin_min = 0.02;
x_bin_max = 5.0;
y_bin_delta = 0.02;
y_bin_min = 0.02;
y_bin_max = 1.0;
xbins = x_bin_min:x_bin_delta:x_bin_max;
ybins = y_bin_min:y_bin_delta:y_bin_max;
contour_min = 0.0;
contour_max = 0.002;
contour_delta = 0.0004;
contours = contour_min:contour_delta:contour_max;
[h,c] = hist3([snr,rho],{xbins,ybins});
[h2,c2] = hist3([snr,rho_noise],{xbins,ybins});
h/=sum(sum(h));
h2/=sum(sum(h2));

figure(1);
hold on;
contour(20*log10(c{1}), c{2}, h',contours,'linewidth',signal_linewidth);
contour(20*log10(c2{1}), c2{2}, h2',contours,'linewidth',noise_linewidth);
axis([-5 15 -0.05 1.05]);

colormap('gray');
%cb = colorbar();
caxis([contour_min, contour_max]);
%set(cb,'YTick',contours,'fontname','courier','fontsize',fontsize,'fontweight','bold');

set(gca(),'fontname','courier','fontsize',fontsize,'box','on');
xlabel('SNR (dB)','fontname','courier','fontsize',fontsize,'position',[5 -0.18]);
ylabel('Correlation Coefficient','fontname','courier','fontsize',fontsize,'position',[-8 0.5]);
print('Aug15_plot2.pdf','-dpdf','-r300');