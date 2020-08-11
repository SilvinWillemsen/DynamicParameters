% S. Bilbao
% 23/07/2020
% AAG, Edinburgh

clear all
close all


plot_on = 0;
interp_on = 1; % turn interpolation on or off
stefansMethod = true;

L = 1;
SR = 44100;
k = 1/SR;

N = 35.9;
h = 1/N;
c = h/k;

Tf = 1;

h = c*k;
lambda = c*k/h;

xctr = 0.5;
wid = 0.1;

N_seg = L/h;

% alf here is the fractional remainder for the grid point...

if stefansMethod
    alf = (interp_on)*((N_seg)-floor(N_seg))
    N = floor(N_seg) - 1;
else
    alf = (interp_on)*((N_seg - 0.5)-floor(N_seg - 0.5));
    N = floor(N_seg);
    if (alf >= 0.5)
        N = N - 1;
    end
end


e = ones(N,1);
Dxx = spdiags([e -2*e e], -1:1,N,N);

% the key adjustments to the Laplacian:
% % 
% Dxx(N, N-1) = 2/((1+alf)^2 + 1);
% Dxx(N, N) = -4 / ((1+alf)^2 + 1);
if stefansMethod
    Dxx(N,N-1) = 2/(2+alf);
    Dxx(N,N) = -2/(2+alf)-2/((2+alf)*(1+alf));
else
    Dxx(N, N-1) = 4 / (2*alf + 3);
    Dxx(N, N) = -4 / (2*alf + 1);
end

Nf = floor(Tf*SR);
out = zeros(Nf,1);

u0 = zeros(N,1);
for ll=1:N
    x = ll*h;
    dist = abs(x-xctr);
    if(dist<=wid)
        u0(ll) = 0.5*(1+cos(pi*dist/wid));
    end
end
    
u  = zeros(N,1);
u1 = u0;
u2 = u0;

P = eye(N);
if stefansMethod
    P(N,N) = (2+alf)/2;
else
    P(N,N) = (2*alf + 3)/4;
end
for n=1:Nf
%     c = c-0.005;
    h = c*k;
    lambda = c*k/h;
    N_seg = L/h;

    % alf here is the fractional remainder for the grid point...

    
    if stefansMethod
        alf = (interp_on)*((N_seg)-floor(N_seg));
        N = floor(N_seg) - 1;
    else
        alf = (interp_on)*((N_seg - 0.5)-floor(N_seg - 0.5));
        N = floor(N_seg);
        if (alf >= 0.5)
            N = N - 1;
        end
    end


    NPrev = N;
    N = floor(N_seg)-1;
    
    if NPrev < N
        u = [u; 0];
        u1 = [u1; (u1(end)) / 2];
        u2 = [u2; (u2(end)) / 2];
        e = ones(N,1);
        Dxx = spdiags([e -2*e e], -1:1,N,N);

        % the key adjustments to the Laplacian:
    
    end
    if stefansMethod
        Dxx(N,N-1) = 2/(2+alf);
        Dxx(N,N) = -2/(2+alf)-2/((2+alf)*(1+alf));
    else
        Dxx(N, N-1) = 4 / (2*alf + 3);
        Dxx(N, N) = -4 / (2*alf + 1);
    end

    P = eye(N);
    if stefansMethod
        P(N,N) = (2+alf)/2;
    else
        P(N,N) = (2*alf + 3)/4;
    end

    u = 2*u1-u2+lambda^2*Dxx*u1;
    
    scaling = ones(N,1);
%     scaling(N) = 1/(alf+1);
    scaling(N) = (2+alf)/2;
    range = 1:N-1;
    H(n) = 0.5*(u-u1)'*P*(u-u1)/k^2-0.5*(c^2/h^2)*u'*P*Dxx*u1;
    
%     hold off;
%     plot(kinEnergy(1:n))
%     hold on;
%     plot(potEnergy(1:n))
%     plot(H(1:n))
%     drawnow;
    u2 = u1;
    u1 = u;
    out(n) = u(floor(1*N/6));
    if((plot_on==1)&&(mod(n,10000)==0))
        plot([0:N, N_seg]'*h, [0; u; 0], 'k', 'Marker', '.', 'MarkerSize', 20);
        pause(0.3);
        drawnow
    end
end

figure(1)
outfft = fft(out);
semilogy([0:Nf-1]'*SR/Nf, abs(outfft), 'r');
xlim([0 3*c/L])

figure(2)

Herr = (H(2:end)-H(1:end-1))/H(1);

plot(Herr)
% soundsc(out,SR)
c/2
