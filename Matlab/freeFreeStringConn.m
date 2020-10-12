clear all;
close all;
% clc;

drawSpeed = 10;
drawStart = 1000000;

fs = 44100;             % Sample rate
k = 1/fs;               % Time step
lengthSound = fs*5;     % Length of the simulation

Ninit = 140;           % edit how many points you want
h = 1/Ninit;

rho = 7850;
r = 0.0005;
A = pi * r^2;
kappa = 1;

cInit = sqrt((h^4 - 4 * kappa^2 * k^2) / (k^2 * h^2));            % calculate wave speed
c = cInit;

T = c^2 * rho * A;
h = sqrt((c^2 * k^2 + sqrt((c^4 * k^4) + 16 * kappa^2 * k^2)) / 2);                % calculate h from wavespeed
N = floor(1/h);         % calculate points from h
% h = 1/N;

lambdaSq = (c*k/h)^2
muSq = kappa^2 * k^2 / h^4

lambdaSq + 4 * muSq                 % should always be 1 as h is not recalculated

alf = Ninit - N;        % fractional remainder for the grid point

%% initialise states
uNext = zeros(ceil(N/2), 1);
u = zeros(ceil(N/2), 1);

u(floor(N/5)-4:floor(N/5)+4) = 1 * hann (9); % use hann window for excitation
uPrev = u;

wNext = zeros(floor(N/2), 1);
w = zeros(floor(N/2), 1);
wPrev = w;

zNext = zeros(N - 1, 1);
z = zeros(N - 1, 1);
z(floor(N/5)-4:floor(N/5)+4) = 1 * hann (9); % use hann window for excitation
zPrev = z;

% initialise laplacians
eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
Dxxxxu = spdiags([eu -4*eu 6*eu -4*eu eu], -2:2, length(u),length(u));
Dxxxxu (1,1) = 5;

ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
Dxxxxw = spdiags([ew -4*ew 6*ew -4*ew ew], -2:2, length(w),length(w));
Dxxxxw (end,end) = 5;

ez = ones(length(z), 1);
Dxxz = spdiags([ez -2*ez ez], -1:1, length(z),length(z));
% Dxxz(end, end) = 0;
% Dxxz(end, end-1) = 0;

Dxxxxz = spdiags([ez -4*ez 6*ez -4*ez ez], -2:2, length(z),length(z));
Dxxxxz (1, 1) = 5;
Dxxxxz (end, end) = 5;
% Dxxxxz (end,end) = 2;
% Dxxxxz (end,end-1) =  -4;
% Dxxxxz (end,end-2) = 2;
% Dxxxxz (end-1, end) = -2;
% Dxxxxz (end-1, end-1) = 5;

interpol = "linear";
outFree = zeros(lengthSound, 1);

changeC = true; % set to true for dynamic changes in wavespeed
interpolatedPoints = [0; 0];

for n = 1:lengthSound  
 
    % change wave speed
    if changeC
        c = cInit * (1 + sin(2 * pi * n/fs));
    else
        c = c;
    end
    cSave(n) = c;
    
    % save previous state for comparison later
    NPrev = N;

    % recalculate gridspacing, points lambda^2 and alpha from new wave speed
    h = sqrt((c^2 * k^2 + sqrt((c^4 * k^4) + 16 * kappa^2 * k^2)) / 2);                % calculate h from wavespeed

    Ninit = 1/h;
    N = floor(1/h);         % calculate points from h
%     h = 1/N;
    
    lambdaSq = (c*k/h)^2;
    muSq = kappa^2 * k^2 / h^4;

    Nsave(n) = N;
    hSave(n) = h;
    
    alf = (Ninit - N);
  
    % calculate interpolator
    if interpol == "cubic"
        ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
    else
        ip = [0, (1-alf), alf, 0];
    end

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    
    % add point if N^n > N^{n-1}
    if N > NPrev
        if mod(N,2) == 1
            uNext = [uNext; (ip(4) * uNext(end-1) + ip(3) * uNext(end) + ip(2) * wNext(1) + ip(1) * wNext(2))];
            u = [u; (ip(4) * u(end-1) + ip(3) * u(end) + ip(2) * w(1) + ip(1) * w(2))];
            uPrev = [uPrev; (ip(4) * uPrev(end-1) + ip(3) * uPrev(end) + ip(2) * wPrev(1) + ip(1) * wPrev(2))];
        else 
            wNext = [(ip(1) * uNext(end-1) + ip(2) * uNext(end) + ip(3) * wNext(1) + ip(4) * wNext(2)); wNext];
            w = [(ip(1) * u(end-1) + ip(2) * u(end) + ip(3) * w(1) + ip(4) * w(2)); w];
            wPrev = [(ip(1) * uPrev(end-1) + ip(2) * uPrev(end) + ip(3) * wPrev(1) + ip(4) * wPrev(2)); wPrev];
        end
        eu = ones(length(u), 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
        Dxxxxu = spdiags([eu -4*eu 6*eu -4*eu eu], -2:2, length(u),length(u));
        Dxxxxu(1,1) = 5;
        
        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
        Dxxxxw = spdiags([ew -4*ew 6*ew -4*ew ew], -2:2, length(w),length(w));
        Dxxxxw(end,end) = 5;

    end
    
    % remove point if N^n < N^{n-1}
    if N < NPrev
        if mod(N,2) == 0
            uNext = uNext(1:end-1);
            u = u(1:end-1);
            uPrev = uPrev(1:end-1);
        else 
            wNext = wNext(2:end);
            w = w(2:end);
            wPrev = wPrev(2:end);
        end
        eu = ones(length(u), 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
        Dxxxxu = spdiags([eu -4*eu 6*eu -4*eu eu], -2:2, length(u),length(u));
        Dxxxxu(1,1) = 5;

        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
        Dxxxxw = spdiags([ew -4*ew 6*ew -4*ew ew], -2:2, length(w),length(w));
        Dxxxxw(end,end) = 5;
        
    end
     
    % calculate interpolated points
    interpolatedPointsPrev = interpolatedPoints;
    interpolatedPoints = [1, -ip(4); -ip(4), 1] \ [ip(3) * w(1) + ip(2) * w(2) + ip(1) * w(3);
                                        ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];
    uMP2 = ip(1) * w(4) + ip(2) * w(3) + ip(3) * w(2) + ip(4) * w(1);
    w0M2 = ip(1) * u(end-3) + ip(2) * u(end-2) + ip(3) * u(end-1) + ip(4) * u(end);

    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - muSq * Dxxxxu * u - uPrev;
    uNext(end-1) = uNext(end-1) - muSq * interpolatedPoints(1);
    uNext(end) = uNext(end) + lambdaSq * interpolatedPoints(1) + 4 * muSq * interpolatedPoints(1) - muSq * uMP2;
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - muSq * Dxxxxw * w - wPrev;
    wNext(2) = wNext(2) - muSq * interpolatedPoints(2);
    wNext(1) = wNext(1) + lambdaSq * interpolatedPoints(2) + 4 * muSq * interpolatedPoints(2) - muSq * w0M2;

    %% orig string
    zNext = 2 * z + lambdaSq * Dxxz * z - muSq * Dxxxxz * z - zPrev;
    
    %% save output
    outFree(n) = w(end - 5);
    outZ(n) = z(end-5);
    %% draw stuff
    if mod(n, drawSpeed) == 0 && n > drawStart
        hLocsLeft = (1:(length(u))) * h;
        hLocsRight = flip(1 - ((1:(length(w))) * h));
   
%         subplot(311)
        hold off;
        plot(hLocsLeft, u, 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 

        hold on;
        wOffset = 0.05;
        plot(hLocsRight, w, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')
        plot((1:length(z)) * h, z, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 5, 'Color', 'g')


%         xlim([0.4, 0.6])

        ylim([-1.5, 1.5])
        grid on;
        title("Sample = " + num2str(n) + "   N = " + num2str(floor(Ninit * 10) / 10))
        set(gca, 'Fontsize', 16, 'Linewidth', 2)
        if n == 16
            disp("wait")
        end

%         subplot(2,1,2)
% 
%         window = 1024;
%         if n>window
%             outfft = fft(outFree(n-window+1 : n));
%             data = abs(outfft);
%             semilogy([0:window-1]'*fs/window, data, 'r');
%             if ~changeC % not accurate if frequency is changing on the fly
%                 title("Peak should be at " + num2str(c/2))
%             end
%             xlim([0 3*cInit])
%         end
        pause(0.2)
        drawnow;
%         if frame <= length(M) && filmFlag == true
%             M(frame) = getframe(gcf);
%             frame = frame + 1;
%         else
%             v = VideoWriter('dynamicGridCenter.mp4', 'MPEG-4');
%             v.FrameRate = 15;
%             open(v)
%             writeVideo(v, M);
%             close(v)
%             filmFlag = false;
%             break;
%         end
        
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext; 
    
    zPrev = z;
    z = zNext;
    
end

subplot(2,1,1)
hold off
plot(outFree)

subplot(2,1,2)
outfft = fft(outFree);
semilogy([0:lengthSound-1]'*fs/lengthSound, abs(outfft), 'r');
xlim([0 3*c])

c/2