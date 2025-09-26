clear;
close;
home;
results = load('all_results_100PeV.dat');
c = 0.3/1.78;
sinC = sin(55.8*pi/180.0);
jj = find(and(abs(results(:,9)*180/pi-55.8)>0,results(:,6)<3.0));
sigma = mean(results(jj,6));
dTheta = std(results(jj,9)-55.8*pi/180.0);
cem = 0.86;
log10E = (c*sigma./cem./dTheta/sinC).^2/log(10.0)+8.0
