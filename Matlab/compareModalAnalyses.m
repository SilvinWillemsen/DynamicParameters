% clear all;
close all;

fs = 44100;
k = 1/fs;
N = 15;

BFull = zeros(N, N);
numFromBound = 1;
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

alf = 0.5;
% Ainv = inv([1, 0; 0, 1]);
%connect the two
ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                    (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                    alf * (alf + 1) * (alf - 2) / -2, ...
                    alf * (alf + 1) * (alf - 1) / 6];
            Ainv = inv([1, -ip(4); -ip(4), 1]);
            
if numFromBound == 2
    BFull(M, (M + 1):(M + 2)) = ip(3:-1:2) * Ainv(1, 1);
    BFull(M, (M-2 : M)) =  BFull(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
    BFull(M+1, (M + 1):(M + 2)) = BFull(M+1, (M + 1):(M + 2)) + ip(3:-1:2) * Ainv(2, 1);
    BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
elseif numFromBound == 1
    BFull(M, M+1) = (ip(3) - ip(1)) * Ainv(1, 1);
    BFull(M, (M-2 : M)) =  BFull(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
    BFull(M+1, M+1) = BFull(M+1, M+1) + (ip(3) - ip(1)) * Ainv(2, 1);
    BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
else
    BFull(M, (M + 1):(M + 3)) = ip(3:-1:1) * Ainv(1, 1);
    BFull(M, (M-2 : M)) =  BFull(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
    BFull(M+1, (M + 1):(M + 3)) = BFull(M+1, (M + 1):(M + 3)) + ip(3:-1:1) * Ainv(2, 1);
    BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
end

u = zeros(N, 1);
% u = rand(N, 1) - 0.5; % without boundaries, plus overlap
% u(M) = u(M+1);
% u(3:7) = hann(5);
u(3:4) = 0.75;
uPrev = u;
uNext = zeros(N, 1);
v = zeros(N-1, 1);
% v = rand(N-1, 1)-0.5;
% v(3:7) = hann(5);
v(3:4) = 0.75;

vPrev = v;
vNext = zeros(N-1, 1);
N = N-1;
Dxxv = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, ones(1, N-1), N, N));
Bnormal = 2 * eye(N) + Dxxv;

N = N + 1;
lengthSound = 10000;
outu = zeros(lengthSound, 1);
outw = zeros(lengthSound, 1);

drawstuff = true;
% plot(sort(1/(2 * pi * k) * acos (1/2 * eig(Bnormal))));
% hold on;
% plot(real(sort(1/(2 * pi * k) * acos (1/2 * eig(BFull)))));
% figure
for n = 1:lengthSound
    uNext = BFull * u - uPrev;
    vNext = Bnormal * v - vPrev;

    outu(n) = uNext(3);
    outv(n) = vNext(3);
    
    uPrev = u;
    u = uNext;
    
    vPrev = v;
    v = vNext;
    if drawstuff && mod(n,100) == 0
        hold off
        plot(1:M, u(1:M), 'o--');
        hold on;
        plot((M:N-1) + alf, u(M+1:end), 'x--');
%         plot(1:N-1, v);
        pause(0.5);
        drawnow;
        if n == 200
            disp("wait")
        end
    end
end
spectrogram(outu, 512,64,512, fs)
view([90, -90])