clear;
close;
home;
pkg load statistics
fontsize = 12;
markersize = 2;
signal_linewidth = 3.0;
noise_linewidth = 1.5;
results = load('all_results.dat');
snr = results(:,7);
rho = results(:,5);
x_bin_delta = 0.02;
x_bin_min = 0.0;
x_bin_max = 5.0;
y_bin_delta = 0.02;
y_bin_min = 0.0;
y_bin_max = 1.0;
xbins = x_bin_min:x_bin_delta:x_bin_max;
ybins = y_bin_min:y_bin_delta:y_bin_max;
contour_min = 0.0;
contour_max = 0.003;
contour_delta = 0.0003;
contours = contour_min:contour_delta:contour_max;
[h,c] = hist3([snr,rho],{xbins,ybins});
h/=sum(sum(h));

figure(1);
contour(c{1}, c{2}, h',contours,'linewidth',signal_linewidth);
set(gca, 'XScale', 'log');
axis([1e-1 10 0 1]);

colormap('gray');
cb = colorbar();
caxis([contour_min, contour_max]);
set(cb,'YTick',contours,'fontname','courier','fontsize',fontsize);

set(gca(),'fontname','courier','fontsize',fontsize,'box','on');
xlabel('Signal to Noise Ratio','fontname','courier','fontsize',fontsize);
ylabel('Correlation Coefficient','fontname','courier','fontsize',fontsize);
print('Aug15_plot2.png','-dpng');