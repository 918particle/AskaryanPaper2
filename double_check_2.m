clear;
close;
home;
pkg load statistics
fontsize = 16;
markersize = 2;
signal_linewidth = 4;
noise_linewidth = 4;
results = load('all_results_3.dat');
sig = results(find(results(:,9)>2),6);
rho = results(find(results(:,9)>2),5);
x_bin_delta = 1.0;
x_bin_min = 0.0;
x_bin_max = 20.0;
y_bin_delta = 0.01;
y_bin_min = 0.0;
y_bin_max = 1.0;
xbins = x_bin_min:x_bin_delta:x_bin_max;
ybins = y_bin_min:y_bin_delta:y_bin_max;
contour_min = 0.0;
contour_max = 0.005;
contour_delta = 0.001;
contours = contour_min:contour_delta:contour_max;
[h,c] = hist3([sig,rho],{xbins,ybins});
h/=sum(sum(h));

figure(1);
hold on;
contour(c{1}, c{2}, h',contours,'linewidth',signal_linewidth);
axis([0.0 20.0 0.0 1.0]);

colormap('gray');
caxis([contour_min, contour_max]);

set(gca(),'fontname','courier','fontsize',fontsize,'box','on');
xlabel('\sigma_{t} (ns)','fontname','courier','fontsize',fontsize,'position',[2.5 -0.15]);
ylabel('Correlation Coefficient','fontname','courier','fontsize',fontsize,'position',[-0.75 0.5]);
print('Sept2_plot1.pdf','-dpdf','-r300');