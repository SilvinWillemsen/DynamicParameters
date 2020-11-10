% close all;
clear all;

fs = 44100;                 % sample rate
k = 1/fs;                   % time step

Ninit = 15;                 % number of intervals (!)
Nend = 20;
h = 1/Ninit;                % grid spacing
c = h/k;                    % wave speed
N = floor(1/h);             % recalculate (should be exactly equal to Ninit)
NinitSave = Ninit;

cInit = h/k;
cEnd = 1/(Nend * k);
loopAmount = 10000;
cVec = linspace (cInit, cEnd, loopAmount + 1);

% lambdaSq = c^2 * k^2 / h^2; % is 1

% Create B-matrix
BFull = zeros(N, N);

M = ceil(N/2);
Mw = floor(N/2);

Dxxu = (sparse(2:M, 1:M-1, ones(1, M-1), M, M) + ...
            sparse(1:M, 1:M, -2 * ones(1, M), M, M) + ...
            sparse(1:M-1, 2:M, ones(1, M-1), M, M));
BFull(1:M, 1:M) = 2 * eye(M) + Dxxu;

Dxxw = (sparse(2:Mw, 1:Mw-1, ones(1, Mw-1), Mw, Mw) + ...
            sparse(1:Mw, 1:Mw, -2 * ones(1, Mw), Mw, Mw) + ...
            sparse(1:Mw-1, 2:Mw, ones(1, Mw-1), Mw, Mw));
        
BFull((M+1):end, (M+1):end) = 2 * eye(Mw) + Dxxw;

BFullInit = BFull;

interpolation = "cubic";
maxNumberOfPoints = max(Ninit, Nend);
modesSave = zeros(loopAmount, maxNumberOfPoints);
NPrev = N;
j = 1;
loopStart = zeros(abs(Nend-Ninit), 1);
loopStart(j) = 1;
addLastPoint = true;

for i = 1:loopAmount
    h = cVec(i) * k;
    Ninit = 1/h;
    N = floor(Ninit);
    alf  = Ninit - N;

    if N ~= NPrev
        BFull = zeros(N, N);
        j = j + 1;
        loopStart(j) = i;
        if i == loopAmount
            addLastPoint = false;
        end
        
        M = ceil(N/2);
        Mw = floor(N/2);

        Dxxu = (sparse(2:M, 1:M-1, ones(1, M-1), M, M) + ...
                    sparse(1:M, 1:M, -2 * ones(1, M), M, M) + ...
                    sparse(1:M-1, 2:M, ones(1, M-1), M, M));
        BFull(1:M, 1:M) = 2 * eye(M) + Dxxu;

        Dxxw = (sparse(2:Mw, 1:Mw-1, ones(1, Mw-1), Mw, Mw) + ...
                    sparse(1:Mw, 1:Mw, -2 * ones(1, Mw), Mw, Mw) + ...
                    sparse(1:Mw-1, 2:Mw, ones(1, Mw-1), Mw, Mw));

        BFull((M+1):end, (M+1):end) = 2 * eye(Mw) + Dxxw;

        BFullInit = BFull;
    end
    NPrev = N;

    if interpolation == "linear"
        BFull(M, (M+1):(length(u) + 2)) = [alf, (1-alf)];
        BFull(M + 1, (M-1) : M) = [(1-alf), alf];
    else
        ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
        Ainv = inv([1, -ip(4); -ip(4), 1]);
        BFull(M, (M + 1):(M + 3)) = ip(3:-1:1) * Ainv(1, 1);
        BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
        BFull(M+1, (M + 1):(M + 3)) = BFullInit(M+1, (M + 1):(M + 3)) + ip(3:-1:1) * Ainv(2, 1);
        BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);

    end
    
    modesSave(i, 1:N) = sort(1/(2 * pi * k) * acos (1/2 * eig(BFull)));
end
if addLastPoint
    loopStart(j+1) = loopAmount;
end
figure
hold off;
modesSavePlot = modesSave;
modesSavePlot(modesSavePlot==0) = nan;
h = plot(modesSavePlot);
% colours = [hsv(floor(maxNumberOfPoints/4)); ...
%            hsv(ceil(maxNumberOfPoints/4)); ...
%            hsv(floor(maxNumberOfPoints/4)); ...
%            hsv(ceil(maxNumberOfPoints/4))];
% colours = parula(floor(maxNumberOfPoints));
% colours = [parula(floor(maxNumberOfPoints / 2)); ...
%     parula(ceil(maxNumberOfPoints/2))];
colours = [1, 0, 0; 0, 0, 1];
colours = repmat(colours, Nend/2, 1);
set(h, {'color', 'Linewidth'}, [num2cell(colours, 2), num2cell(2 * ones(maxNumberOfPoints, 1))])
title ("Modal Analysis $N = 15 \rightarrow 20$", 'interpreter', 'latex');
xlabelsave = num2cell(NinitSave:Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', loopStart, 'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex')
xlabel("$N$", 'interpreter', 'latex')
ylabel("Frequency (Hz)", 'interpreter', 'latex')
grid on
xData = 1:loopAmount;
goodness = zeros(N, 1);
for i = 1:N
    [~, got] = fit (xData', modesSave(:,i), 'poly1');
    goodness(i) = got.adjrsquare;
%     eval(['mode', num2str(i), ' = modesSave(:,', num2str(i), ');']);
end
plot(goodness)
