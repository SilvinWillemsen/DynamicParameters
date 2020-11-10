% close all;
clear all;
loopingN = false;
loopNStart = 20; % also use for Ninit an Nend
loopNend = 25;

if loopingN
    range = loopNStart:0.01:loopNend;
else
    range = 1;
end

changeNumFromBound = true;

for Nloop = range 
    fs = 44100;                 % sample rate
    k = 1/fs;                   % time step

    if loopingN
        Ninit = Nloop;          % number of intervals (!)
        Nend = Nloop;
    else
        Ninit = loopNStart;
        Nend = loopNend;
    end
    
    h = 1/Ninit;                % grid spacing
    c = h/k;                    % wave speed
    N = floor(1/h);             % recalculate (should be exactly equal to Ninit)
    NinitSave = Ninit;

    cInit = h/k;
    cEnd = 1/(Nend * k);
    loopAmount = 1000;
    cVec = linspace (cInit, cEnd, loopAmount + 1);

    % Create B-matrix
    BFull = zeros(N, N); 

    if changeNumFromBound
        if mod(N, 2) == 0
            numFromBound = 2;
        else
            numFromBound = 3;
        end
    else
        numFromBound = 2;
    end
    
    M = N-numFromBound;
    Mw = numFromBound;

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
    modesSave = zeros(loopAmount, floor(maxNumberOfPoints));
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
            if i~= 1
                j = j + 1;
                loopStart(j) = i;
            end
            if i == loopAmount
                addLastPoint = false;
            end
       
            if changeNumFromBound 
                if mod(N, 2) == 1
                    numFromBound = numFromBound + 1;
                else
                    numFromBound = numFromBound - 1;
                end
            end
            M = N-numFromBound;
            Mw = numFromBound;

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
            BFull(M, (M + 1):(M + 2)) = ip(3:-1:2) * Ainv(1, 1);
            BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
            BFull(M+1, (M + 1):(M + 2)) = BFullInit(M+1, (M + 1):(M + 2)) + ip(3:-1:2) * Ainv(2, 1);
            BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);

        end

        modesSave(i, 1:N) = sort(1/(2 * pi * k) * acos (1/2 * eig(BFull)));
    end
    if addLastPoint
        loopStart(j+1) = loopAmount;
    end
%     figure
    hold off;
    modesSave = modesSave;
    modesSave(modesSave==0) = nan;
    h = plot(real(modesSave));
    
    colours = [];
    for j = 1:floor(Nend)
        if mod(j,2) == 0
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
    drawnow;
end
xData = 1:loopAmount;
goodness = zeros(N, 1);
figure;
for i = 1:N
    [~, got] = fit (xData(~isnan(modesSave(:,i)))', modesSave(~isnan(modesSave(:,i)),i), 'poly1');
    goodness(i) = got.adjrsquare;
    got
%     eval(['mode', num2str(i), ' = modesSave(:,', num2str(i), ');']);
end
plot(goodness)

