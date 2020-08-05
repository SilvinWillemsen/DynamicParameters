clear all;
close all;

fs = 44100;

%% String
k = 1/fs;
L1 = 1; % Length of the string [m]
L2 = 1; % Length of the string [m]
c1 = 300;
c2 = 300;
[B1, C1, N1, h1, D1] = create1Dwave (c1, L1, k);
[B2, C2, N2, h2, D2] = create1Dwave (c2, L2, k);


%% Initialise state vectors
u = zeros(N1, 1);
w = zeros(N2, 1);

exciteString1 = true;
exciteString2 = false;

%% Excite
if exciteString1
    exciterPos = 0.25;
    rcW = 6;
    u(1 + floor(exciterPos*N1-rcW/2):1 + floor(exciterPos*N1+rcW/2)) = (1-cos(2*pi*(0:rcW)/rcW)) * 0.5;

end
if exciteString2
    exciterPos = 0.25;
    rcW = 6;
    w(N1 + 1 + floor(exciterPos*N2-rcW/2):N1 + 1 + floor(exciterPos*N2+rcW/2)) = (1-cos(2*pi*(0:rcW)/rcW)) * 0.5;
end
uPrev = u;
uNext = u;

wPrev = w;
wNext = w;

%% Connection points
% connS = 0.5;
% connP = 0.5;

Iu = zeros(N1, 1);
Ju = zeros(N1, 1);

Iw = zeros(N2, 1);
Jw = zeros(N2, 1);

alf = 0.25;
conn1 = 0.5;
Iu(1 + floor(conn1*N1):1 + floor(conn1*N1+1)) = [1-alf; alf];
Ju(1 + floor(conn1*N1):1 + floor(conn1*N1+1)) = [1-alf; alf] / h1;

conn2 = 0.5;
Iw(1 + floor(conn2*N2):1 + floor(conn2*N2+1)) = -[1-alf; alf];
Jw(1 + floor(conn2*N2):1 + floor(conn2*N2+1)) = -[1-alf; alf] / h2;

%% Length of the sound
lengthSound = round(fs);

potEnergyString1 = zeros(lengthSound, 1);
kinEnergyString1 = zeros(lengthSound, 1);
totEnergyString1 = zeros(lengthSound, 1);
totEnergyString1Plot = zeros(lengthSound, 1);
potEnergyString2 = zeros(lengthSound, 1);
kinEnergyString2 = zeros(lengthSound, 1);
totEnergyString2 = zeros(lengthSound, 1);
totEnergyString2Plot = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);
totEnergyPlot = zeros(lengthSound, 1);
energyConnection = zeros(lengthSound, 1);
totEnergyPlusConn = zeros(lengthSound, 1);
totEnergyPlusConnPlot = zeros(lengthSound, 1);

out = zeros(lengthSound, 1);
out2 = zeros(lengthSound, 1);

stringVec1 = 2:N1-1;
stringVec2 = 2:N2-1;

connected = true;

JuFc = zeros(N1, 1);
JwFc = zeros(N2, 1);

drawState = false;
drawEnergy = true;

w0 = 100000000;
w1 = 0;
sx = 0;

% L = [I(1:N1) * h1; I(N1+1:end) * h2]';

etaPrev = 0;
eta = 0;
for n = 1 : lengthSound
    eta = Iu' * u + Iw' * w;
    uNext = B1 * u + C1 * uPrev;
    wNext = B2 * w + C2 * wPrev;

    if connected
        varPhi = w0/4 + (w1 * eta^2)/2;
        Fc = (-(Iu' * uNext + Iw' * wNext) * varPhi - (w0 * eta) / 2 - etaPrev * varPhi) / ((Iu' * Ju + Iw' * Jw) * k^2 * varPhi + 1);        
        JuFc = Ju*Fc;
        JwFc = Jw*Fc;

    else 
        JuFc = 0;
        JwFc = 0;
    end
    uNext = uNext + JuFc * k^2; 
    wNext = wNext + JwFc * k^2; 

    %     I' * u
    potEnergyString1(n) = c1^2 / 2 * sum (1 / h1 * (u(stringVec1 + 1) - u(stringVec1)) .* (uPrev(stringVec1 + 1) - uPrev(stringVec1)));
    kinEnergyString1(n) = 1 / 2 * sum (h1 * ((1 / k * (u(stringVec1) - uPrev(stringVec1))).^2));
    
    potEnergyString2(n) = c2^2 / 2 * sum (1 / h2 * (w(stringVec2 + 1) - w(stringVec2)) .* (wPrev(stringVec2 + 1) - wPrev(stringVec2)));
    kinEnergyString2(n) = 1 / 2 * sum (h2 * ((1 / k * (w(stringVec2) - wPrev(stringVec2))).^2));
    
%     energyConnection(n) = w0/4 * (eta^2 + etaPrev^2); % when using center time difference 
%     energyConnection(n) = w0/2 * (1/2 * (eta + etaPrev))^2; % when using forward * backward time difference (mu_tt}
%     energyConnection(n) = w1 / 2 * (1/2 * eta^2 * etaPrev^2); 
    energyConnection(n) = w0/2 * (1/2 * (eta + etaPrev))^2 + w1 / 2 * (1/2 * eta^2 * etaPrev^2);
    etaPrev = eta;
    
    totEnergyString1(n) = kinEnergyString1(n) + potEnergyString1(n);
    totEnergyString1Plot(n) = (totEnergyString1(n)-totEnergyString1(1))/totEnergyString1(1);
    totEnergyString2(n) = kinEnergyString2(n) + potEnergyString2(n);
    totEnergyString2Plot(n) = (totEnergyString2(n)-totEnergyString2(1))/totEnergyString2(1);
    totEnergy(n) = totEnergyString1(n) + totEnergyString2(n);
    totEnergyPlot(n) = (totEnergy(n)-totEnergy(1))/totEnergy(1);
    totEnergyPlusConn(n) = totEnergyString1(n) + totEnergyString2(n) + energyConnection(n);
    totEnergyPlusConnPlot(n) = (totEnergyPlusConn(n)-totEnergyPlusConn(1))/totEnergyPlusConn(1);
        
        
    if drawEnergy == true && mod(n, 10) == 0

        clf
        subplot(3,2,1)
%         plot(kinEnergyString1(1:n))
%         hold on;
%         plot(potEnergyString1(1:n));
%         title("Kinetic and Potential Energy String ")
%         subplot(3,2,2)
        plot(totEnergyString1(1:n))
        title("Total String1 Energy")

        
        subplot(3,2,2)
%         plot(kinEnergyString2(1:n))
%         hold on;
%         plot(potEnergyString2(1:n));
%         title("Kinetic and Potential Energy String 2")
%         subplot(3,2,4)
        plot(totEnergyString2(1:n))
        title("Total String2 Energy")

        subplot(3,2,3)
        plot(energyConnection(1:n))
        title("Connection Energy")
        
        subplot(3,2,4);
        plot(totEnergy(1:n))
        title("TotEnergy Without Connection")
        
        subplot(3,2,5)
        plot(totEnergyPlusConnPlot(1:n))
        title("Total Energy")
        drawnow;
       
%        subplot(3, 1, 1)
%        plot(energyConnection(1:n))
%        title("Connection Energy");
%        subplot(3, 1, 2);
%        plot(totEnergy(1:n))
%        title("Tot Energy Min Connection");
%        subplot(3, 1, 3);
%        plot(totEnergyPlusConn(1:n))
%        title("Tot Energy Plus Connection");
%        drawnow;
    end
    
    uPrev = u;
    u = uNext;
        
    wPrev = w;
    w = wNext;

    if drawState == true % && mod(n, 10) == 0
        clf
        subplot(2,1,1)
        plot (u, 'o-')
         subplot(2,1,2)
        plot(w, 'o-')
        pause(0.2);
        drawnow;
    end
    
    out(n) = uNext(floor(N1/pi));
    out2(n) = w(round(N2 - 15));
%     clf;
%     scatter(uNext(idx(1)), 1)
%     hold on;
%     scatter(uNext(idx(2)), 2)
%     drawnow;
end
% subplot(2,1,1)
% plot(out)
% subplot(2,1,2)
% plot(out2)


figure;
totEnergyString1 = kinEnergyString1 + potEnergyString1;
totEnergyString1Plot = (totEnergyString1-totEnergyString1(1))/totEnergyString1(1);
plot(totEnergyString1Plot)

hold on;
totEnergyString2 = kinEnergyString2 + potEnergyString2;
totEnergyString2Plot = (totEnergyString2-totEnergyString2(1))/totEnergyString2(1);
plot(totEnergyString2Plot)

figure;
totEnergy = totEnergyString1+totEnergyString2+energyConnection;
totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
maxTotEnergyDiff = (max(totEnergy) + abs(min(totEnergy)));
plot(totEnergy)
drawnow;
        