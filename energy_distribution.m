clear;
close;
home;
fontsize = 16;
markersize = 6;
linewidth = 1.5;
results = load('all_results_2.dat');
c = 0.3/1.78;
sinC = sin(55.8*pi/180.0);
dTheta = 5.77*pi/180.0;
sigma = results(:,6);
snr = results(:,7);
cem = 0.825;
log10E = (c*sigma./cem/dTheta/sinC).^2/log(10.0)+8.0;
mean_log10E = mean(log10E);
std_log10E = std(log10E);
x_bin_delta = 1.0;
x_bin_min = 8.0;
x_bin_max = 30.0;
xbins = x_bin_min:x_bin_delta:x_bin_max;
[y,x] = hist(log10E,xbins,1);

%Prediction from NuRadioMC, 10, 20, ..., 100 PeV
m = mean(log10((10:10:100)*10^(15)));
s = 2;
shade_delta = x_bin_delta/100.0;
shade_bins = x_bin_min:shade_delta:x_bin_max;
region = zeros(size(shade_bins));
region(find(and(shade_bins>=(m-s),shade_bins<=(m+s)))) = 1;

figure(1);
hold on;
bar(x,y,'facecolor','white','edgecolor','black','linewidth',linewidth);
axis([7 31 0.0 0.25]);
h = area(shade_bins,region);
set(h,'facecolor',[0.5 0.5 0.5]);
set(h,'edgecolor',[0.5 0.5 0.5]);
set(h,'facealpha',0.2);
set(h,'edgealpha',0.9);
plot([m m],[0 1],'-','color','black','linewidth',linewidth);
labelString = ['mean: ',num2str(mean_log10E,3),' std: ',num2str(std_log10E,3)];
text(16.0,0.26,labelString,'fontname','courier','fontsize',fontsize)
set(gca(),'fontname','courier','fontsize',fontsize,'box','on');
xlabel('log_{10} E_{C}','fontname','courier','fontsize',fontsize);
ylabel('Probability (Normalized)','fontname','courier','fontsize',fontsize);
print('Aug19_plot1.pdf','-dpdf');