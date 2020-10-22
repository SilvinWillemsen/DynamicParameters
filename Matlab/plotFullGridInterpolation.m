clear all;
close all;

figure1 = figure('Color',[1 1 1], 'Position', [0, 800, 500, 300]);
axes1 = axes('Parent',figure1,...
    'Position',[0.048 0.162629757785467 0.94 0.829036908881207]);


N = 30;
u = zeros(N+1, 1);
u(N/4 - 2 :  N/4 + 2) = hann(5);
u(3*N/5 - 2 :  3*N/5 + 2) = hann(5);

diff = 1;
% interp1((1:Nprev)*hPrev,u,(1:N)*h,'cubic', 'extrap')';
w = interp1((0:N) / N, u, (0:N+diff) / (N+diff), 'linear', 'extrap')';

hold off
for i = 0.2:0.2:0.8
    plot([0.002, 1-0.002], [i, i], 'Color', [0, 0, 0, 0.1], 'LineWidth', 2);
    hold on;
end

unm1 = plot((0:N)/N, u, 'b-o', 'Linewidth', 1.5, 'Markersize', 8);
un = plot((0:N+diff)/(N+diff), w, 'r--x', 'Linewidth', 1.5, 'Markersize', 9);
legend([unm1, un], ["$N^{n-1} = " + num2str(N) + "$", "$N^n = " + num2str((N+diff)) + "$"], 'Fontsize', 16, 'interpreter', 'latex')
yticks([]);
xlabel("$N^n/N^n$", 'interpreter', 'latex')
ylabel("$u$", 'interpreter', 'latex', 'Position', [-0.01, 0.5])
ylim([0, 1.05])
set(gca, 'FontSize', 16, 'Linewidth', 2);
set(gcf,'color','w');
grid on
