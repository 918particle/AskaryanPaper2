clear;
close;
home;
pkg load optim
fontsize = 16;
markersize = 6;
linewidth = 1.5;
results = load('all_results.dat');
theta = results(:,9)*180.0/pi;
x_bin_delta = 2;
x_bin_min = 0.0;
x_bin_max = 90.0;
xbins = x_bin_min:x_bin_delta:x_bin_max;
[y1,x1] = hist(theta,xbins,1);

p0 = [55.8,10,0.33];
[p1,~] = lsqcurvefit(@fit_function_3,p0,x1,y1);
xfit = linspace(x_bin_min,x_bin_max,length(xbins)*100);
yfit = fit_function_3(p1,xfit);
cAngle = acos(1.0/1.78)*180.0/pi;

figure(1);
hold on;
semilogy(x1,y1,'o','color','black','markersize',markersize);
semilogy(xfit,yfit,'-','color','black','linewidth',linewidth);
axis([0 90 1e-4 10]);
plot([cAngle cAngle],[1e-4 10],'--','color','#444444','linewidth',linewidth)
legendString = ['\mu: ' num2str(p1(1),3) ' \sigma: ' num2str(p1(2),3)];
h = legend('NuRadioMC output',legendString,'Cherenkov Angle','location','northwest');
set(h,'fontname','courier','fontsize',fontsize)
set(gca(),'fontname','courier','fontsize',fontsize,'box','on');
xlabel('Viewing Angle (deg)','fontname','courier','fontsize',fontsize);
ylabel('Probability (Normalized)','fontname','courier','fontsize',fontsize);
print('Aug18_plot1.pdf','-dpdf');