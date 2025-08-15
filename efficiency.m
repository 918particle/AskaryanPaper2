clear;
close;
home;
pkg load optim

results = load('all_results.dat');
bin_min = 0.100;
bin_max = 1.0;
bin_delta = 0.01;
bins = bin_min:bin_delta:bin_max;
n = length(bins);
[y1,x1] = hist(results(:,8),bins,1);
[y2,x2] = hist(results(:,5),bins,1);

p0 = [0.08,0.05,0.05];
[p1,~] = lsqcurvefit(@fit_function,p0,x1,y1);
xfit = linspace(bin_min,bin_max,length(bins)*100);
yfit = fit_function(p1,xfit);
figure(1,'position',[0 500 0 500]);
hold on;
semilogy(x1,y1,'o','color','black','markersize',8);
semilogy(xfit,yfit,'-','color','black','linewidth',0.5)
semilogy(x2,y2,'-','color','black','linewidth',2)
axis([0 1 3.01e-8 0.301]);
set(gca(),'fontname','courier','fontsize',20,'box','on');
xlabel('Correlation coefficient','fontname','courier','fontsize',20,'offset',[0 -0.1]);
ylabel('Probability density','fontname','courier','fontsize',20,'offset',[-0.1 0]);