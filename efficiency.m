clear;
close;
home;
results = load('all_results.dat');
bin_min = 0.0;
bin_max = 1.0;
bin_delta = 0.005;
bins = bin_min:bin_delta:bin_max;
n = length(bins);
[y1,x1] = hist(results(:,8),bins);
[y2,x2] = hist(results(:,5),bins);
eff = zeros(size(bins));
for i=1:n
	if(y1(i)>0)
		eff(i) = y2(i)/y1(i);
	endif
endfor
plot(eff)
hold on
plot(y1)
plot(y2)