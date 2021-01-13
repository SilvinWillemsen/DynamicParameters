clear all;
close all;

figure('Position', [0, 0, 500, 500])
% x = rand(1, 21) - 0.5;
x = 0:1:9;
x = sort(x);
h = max(diff(x));
xI = min(x):0.001:max(x);
% h = 1
% y = rand(length(x), 1);
freq = 1/max(x);
y = sin(2 * pi * freq * x');
yTrue = sin(2 * pi * freq * xI');
% xI = 0;
bmax = pi / h;
dist = (xI' - x);
sincRes = sin(bmax * dist) ./ (bmax * dist);
sincRes (isnan(sincRes)) = 1;
yInterp = sincRes * y;

subplot(311)
plot(xI, yTrue, 'k', 'Linewidth', 2);
hold on;
scatter(x,y, 'r', 'Linewidth', 2)
set(gca, 'Position', [0.01 0.69 0.97 0.3])
grid on;
set(gca, 'Linewidth', 2)
set(gca,'xticklabel',{[]})
text((x(end-1) + x(end)) * 0.5, 0.8, "a)", 'Fontsize', 16, 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')

subplot(312)
plot(xI, sincRes .* y', '--', 'Linewidth', 0.75, 'color', [0.5, 0.5, 0.5])
grid on;
hold on;
plot(xI, sincRes(:,3) .* y(3), 'Linewidth', 1.5, 'color', 'k')
scatter(x,y, 'r', 'Linewidth', 2)

set(gca, 'Linewidth', 2)
set(gca, 'Position', [0.01 0.38 0.97 0.3])
set(gca,'xticklabel',{[]})
text((x(end-1) + x(end)) * 0.5, 0.8, "b)", 'Fontsize', 16, 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')

subplot(313)
if length(xI) == 1
    scatter(xI, yInterp)
else
    plot(xI, yInterp, 'color', 'b', 'Linewidth', 2);
end
grid on;
hold on;
scatter(x,y, 'r', 'Linewidth', 2)
set(gca, 'Linewidth', 2)
set(gca, 'Position', [0.01 0.07 0.97 0.3], 'Fontsize', 14)
xlabel("$l$", 'interpreter', 'latex', 'Position', [mean(xlim), -1.2])
text((x(end-1) + x(end)) * 0.5, 0.8, "c)", 'Fontsize', 16, 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')

set(gcf, 'color', 'w')
%% plot sinc
% figure('Position', [0, 0, 700, 200])
% sincRange = (-5:0.01:5);
% sincToPlot = sin(pi * sincRange)./(pi * sincRange);
% sincToPlot(isnan(sincToPlot)) = 1;
% plot(sincRange, sincToPlot, 'k', 'Linewidth', 2)
% grid on
% set(gca, 'Linewidth', 2, 'Fontsize', 16, ...
% 'Position',  [0.0485714285714286 0.121 0.938571428571429 0.834])
% set(gcf, 'color', 'w')