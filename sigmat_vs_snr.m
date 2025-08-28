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
sigma = results(:,6);
x_bin_delta = 0.02;
x_bin_min = 0.02;
x_bin_max = 5.0;
y_bin_delta = 0.25;
y_bin_min = 0.0;
y_bin_max = 5.0;
xbins = x_bin_min:x_bin_delta:x_bin_max;
ybins = y_bin_min:y_bin_delta:y_bin_max;
contour_min = 0.0;
contour_max = 0.005;
contour_delta = 0.001;
contours = contour_min:contour_delta:contour_max;
[h,c] = hist3([snr,sigma],{xbins,ybins});
h/=sum(sum(h));

figure(1);
hold on;
contour(20*log10(c{1}), c{2}, h',contours,'linewidth',signal_linewidth);
axis([-5 15 0.0 5.0]);

colormap('gray');
%cb = colorbar();
caxis([contour_min, contour_max]);
%set(cb,'YTick',contours,'fontname','courier','fontsize',fontsize,'fontweight','bold');

set(gca(),'fontname','courier','fontsize',fontsize,'box','on');
xlabel('SNR (dB)','fontname','courier','fontsize',fontsize,'position',[5 -0.75]);
ylabel('\sigma_{t} (ns)','fontname','courier','fontsize',fontsize,'position',[-6.75 2.5]);
print('Aug27_plot1.pdf','-dpdf','-r300');