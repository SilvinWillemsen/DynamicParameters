clear all;
close all;

x = rand(1, 21) - 0.5;
% x = -1.05:0.1:1;
x = sort(x);
h = max(diff(x));
% h = 1
% y = rand(length(x), 1);
freq = 1;
y = sin(2 * pi * freq * x');
xI = 0:0.0001:1;
xI = 0;
bmax = pi / h;
sinc = sin(bmax * (xI' - x)) ./ (bmax * (xI' - x));
sinc (isnan(sinc)) = 1;
yInterp = sinc * y;
scatter(x,y)
hold on; scatter(xI, yInterp)


%% plot sinc
figure('Position', [0, 0, 700, 200])
sincRange = (-5:0.01:5);
sincToPlot = sin(pi * sincRange)./(pi * sincRange);
sincToPlot(isnan(sincToPlot)) = 1;
plot(sincRange, sincToPlot, 'k', 'Linewidth', 2)
grid on
set(gca, 'Linewidth', 2, 'Fontsize', 16, ...
'Position',  [0.0485714285714286 0.121 0.938571428571429 0.834])
set(gcf, 'color', 'w')