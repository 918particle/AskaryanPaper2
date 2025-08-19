clear;
close;
home;
pkg load optim
fontsize = 16;
markersize = 8;
signal_linewidth = 3.0;
noise_linewidth = 1.5;
results = load('all_results.dat');
bin_delta = 0.02;
bin_min = 0.1;
bin_max = 1.0;
bins = bin_min:bin_delta:bin_max;
n = length(bins);
[y1,x1] = hist(results(:,8),bins,1);
[y2,x2] = hist(results(:,5),bins,1);
[~,pk_index] = max(y2);
pk_corr = x2(pk_index);

[~,exp_limit_index] = max(y1);
exp_limit = x1(exp_limit_index)+7*bin_delta;
exp_limit_index = find(x1>=exp_limit,1);
p0 = [0.08,0.05,0.05];
p1 = [exp_limit,0.012];
[p2,~] = lsqcurvefit(@fit_function,p0,x1,y1);
[p3,~] = lsqcurvefit(@fit_function_2,p1,x1(exp_limit_index:end),y1(exp_limit_index:end));
xfit = linspace(bin_min,bin_max,length(bins)*100);
yfit1 = fit_function(p2,xfit);
yfit2 = fit_function_2(p3,xfit);
noise_events = 100;
cut_limit = 0.4;
signal_efficiency = 0;
noise_rate = 1;
while(noise_events>1)
	cut_limit += 0.001;
	sig_ll = find(x2>=cut_limit,1);
	noise_ll = find(xfit>=cut_limit,1);
	signal_efficiency = sum(y2(sig_ll:end));
	noise_events = 5*365*24*60*60*noise_rate*sum(yfit2(noise_ll:end));
endwhile
display([exp_limit pk_corr cut_limit signal_efficiency noise_events])

figure(1);
hold on;
semilogy(x1,y1,'o','color','black','markersize',markersize);
semilogy(xfit,yfit1,'--','color','#999999','linewidth',noise_linewidth)
semilogy(xfit,yfit2,'-','color','black','linewidth',noise_linewidth)
semilogy(x2,y2,'-','color','#999999','linewidth',signal_linewidth)
plot([cut_limit cut_limit],[1e-9 1],'--','color','#444444','linewidth',signal_linewidth)
axis([0 1 1e-7 1]);
set(gca(),'fontname','courier','fontsize',fontsize,'box','on');
xlabel('Correlation coefficient','fontname','courier','fontsize',fontsize);
ylabel('Probability (Normalized)','fontname','courier','fontsize',fontsize);
print('Aug15_plot1.pdf','-dpdf');