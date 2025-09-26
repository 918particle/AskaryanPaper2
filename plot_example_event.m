clear;
close;
home;

dt = 1; %ns

results = load('example_event_bad.dat');
[n m] = size(results);
times = [0:1:n-1]*dt;
csw = results(:,1);
model = results(:,2);
model = circshift(model,[44,0]);

plot(csw)
hold on
plot(model)

results = [times' csw model]