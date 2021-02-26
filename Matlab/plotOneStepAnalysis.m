%to be used with modalAnalysis.m
modesSave(modesSave==0) = nan;
if (Nend-NinitSave) < 0
    modesSaveRange = 2:size(modesSave, 1);
    loopStartRange1 = 2:length(loopStart);
else
    modesSaveRange = 1:size(modesSave,1);
    loopStartRange1 = 1:length(loopStart);
end

colorMap = zeros(length(modesSaveRange), 3);
imagesc(exp(sigmaSave))
normSigmaSave = abs(sigmaSave) / max(max(abs(sigmaSave)));
for scatLoop = 1:size(modesSave, 2)
    damp = exp(sigmaSave(:, scatLoop));
%     colorMap = repmat(1-damp, 1, 3);
    colorMap = [(1-damp) * 0.75,(1-damp) * 0.75, 1-damp];
    sz = 10 * damp + 0.5;
%     sz = 5;
    scatter(modesSaveRange, real(modesSave(modesSaveRange, scatLoop)), sz, colorMap);
    hold on;
end
hold off
title ("Modal Analysis $N = " + loopNStart + " \rightarrow" + loopNend + "$", 'interpreter', 'latex');
xlabelsave = num2cell(NinitSave:sign(Nend-NinitSave):Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', loopStart(loopStartRange1), 'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex')
xlabel("$N$", 'interpreter', 'latex')
ylabel("Frequency (Hz)", 'interpreter', 'latex')
ylim([0, fs / 2])
grid on