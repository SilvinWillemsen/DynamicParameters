%to be used with modalAnalysis.m
modesSave(modesSave==0) = nan;
h = plot(real(modesSave));

colours = [];
for colLoop = 1:floor(maxNumberOfPoints)
    if mod(colLoop,2) == 0
        colours = [colours; 0,0,1];
    else
        colours = [colours; 1,0,0];
    end
end
%     figure
set(h, {'color', 'Linewidth'}, [num2cell(colours, 2), num2cell(2 * ones(floor(maxNumberOfPoints), 1))])
title ("Modal Analysis $N = " + loopNStart + " \rightarrow" + loopNend + "$", 'interpreter', 'latex');
xlabelsave = num2cell(NinitSave:Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', loopStart, 'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex')
xlabel("$N$", 'interpreter', 'latex')
ylabel("Frequency (Hz)", 'interpreter', 'latex')
ylim([0, fs / 2])
grid on