clc;

loadFiles = true;
if loadFiles
%     clear all;
    mode = "Release";
    loadCppFiles;
end

if curFs == 44100
    close all;
end

plotState = false;
k = 1/curFs;
drawSpeed = 1; 

plotSpectrogram = true;

% if origVersion == "I"
%     dynamicString = 0;
%   %     if curFs == 44100
% % %         stateAtI = reshape (stateAt, NSave(1)-1, 1000);
% %         outputI = output;
% %     else
% %         outputINot44100 = output;
% %     end
% elseif origVersion == "D"
%     dynamicString = 1;
% 
% %     if curFs == 44100
% % %         stateAtD = reshape (stateAt, NSave(1), 1000);
% % %         stateAtD(floor(NSave(1) * 0.5), :) = [];
% %         outputD = output;
% %     else
% %         outputDNot44100 = output;
% %     end
% elseif origVersion == "B"
%     
% end

NChangedIdx = 0;
startIdxU = 1;
startIdxW = 1;
startSample = 0; %should come from freeFree1DWaveConn

% for i = 1:length(stateAtI)
%     hold off;
%     plot(stateAtI(:,i)); 
%     hold on;
%     plot(stateAtD(:,i)); 
%     pause(0.1)
%     drawnow;
% end
startSample = 0;
if plotState    
    if origVersion == "D" || origVersion == "B"
        for n = 1:4*drawSpeed:length(plotIdxD)
            if n >= startSample * 4
        %     hold off;
                Mu = plotIdxD(n+1)-plotIdxD(n) + 1;
                Mw = plotIdxD(n+3)-plotIdxD(n+2) + 1;
                if n > length(plotIdxD) * 0.9
                    disp ("wait");
                end
                cSaveIdx = ((n-1)/4 + 1);
%                 h = cSaveD(round(((n-1)/4) * length(outputD) / (length(plotIdxD)/4))+1) * k
                h = cSaveD(cSaveIdx) * k;
                hLocsLeft = (1:Mu) * h;
                hLocsRight = 1 - fliplr((1:Mw) * h);
                hold off;
                plot(hLocsLeft, stateAtD(plotIdxD(n):plotIdxD(n+1)) ,'-o');
                hold on;
                plot(hLocsRight, stateAtD(plotIdxD(n+2):plotIdxD(n+3)), '-x');
                
%                 endIdxU = startIdxU + MuSave((n-1)/4 - startSample + 1)-1;
%                 endIdxW = startIdxW + MwSave((n-1)/4 - startSample + 1)-1;
% 
%                 plot(hLocsLeft, uSave(startIdxU:endIdxU))
%                 plot(hLocsRight, wSave(startIdxW:endIdxW))
%                 startIdxU = endIdxU + 1;
%                 startIdxW = endIdxW + 1;
                title(100 * n/(length(plotIdxD)))
                ylim([-1, 1])
%                 xlim([0.9, 1])
%                 pause(0.5);
                if cSaveIdx > 4 && (NSaveD(cSaveIdx) ~= NSaveD(cSaveIdx-4))
%                     disp("wait");
                end
                drawnow;
                if n == 200 * 4 + 1
                    disp("wait");
                end
                    
            end
        end 
    end
    if origVersion == "I" || origVersion == "B"
        hold off;
        Nprev = plotIdxI(2)-plotIdxI(1) + 2;
        for n = 1:2*drawSpeed:length(plotIdxI)/2
        %     hold off;
            N = plotIdxI(n+1)-plotIdxI(n) + 2;

        %     semilogy((start1:NPrev-1)/NPrev, abs(stateAt((plotIdx(n)+start1-1):plotIdx(n+1))));
        %     subplot(211)
            plot((1:N-1)/N, stateAtI(plotIdxI(n):plotIdxI(n+1)));
            
            if N ~= Nprev
                NChangedIdx = NChangedIdx + 1;
                disp("NChanged" + NChangedIdx)
            end
            
            Nprev = N;
        %     plot((1:1500)/NPrev, stateAt((plotIdx(n)):(plotIdx(n)+1499)));
    %         plot((1:NPrev-1)/ NPrev, stateAt((plotIdx(n)):plotIdx(n+1)))
    %         xlim([0, 1]);
        %     hold on;
        %     subplot(212)
        %     plot((start1+1:NPrev-2)/NPrev, stateAt((plotIdx(n)+start1 + 1):plotIdx(n+1)) -2 * stateAt((plotIdx(n)+start1):(plotIdx(n+1)-1)) + stateAt((plotIdx(n)+start1-1):(plotIdx(n+1)-2)));
        %     hold on;
        %     N = plotIdx(n+3)-plotIdx(n+2) + 2;
        %     start2 = max(1, (N-stateStartFromRightBound));

        %     semilogy((start2: zN-1)/N, abs(stateAt((plotIdx(n+2)+start2-1):plotIdx(n+3))));
        %     plot((start1+1:NPrev-2)/NPrev, diff(diff(stateAt((plotIdx(n)+start1-1):plotIdx(n+1)))));
        %     title(NChange((n-1)/4 + 1, 1))
        %     title(NSave((n-1)/2 + 1, 1))
            title(100 * n/(length(plotIdxI)/2))
            ylim([-1, 1])
    %     xlim([0.45, 0.55])
            pause(0.1);
            drawnow;
        end 
    end
    
end

if plotSpectrogram
    figure;
    % plot(output((21900:21960) * curFs / 44100))
    % plot(outputI)
    % hold on;
    % plot(outputD)

    % plot(abs(fft(outputI)))
    % hold on;
    % plot(abs(fft(outputD)))

    % if curFs == 44100
    %     outputIToPlot = outputI; 
    %     outputDToPlot = outputD; 
    % else
    %     outputIToPlot = outputI;
    %     outputDToPlot = outputDNot44100;
    % end
    if origVersion == "B"
        subplot(211)
    end
    if origVersion == "I" || origVersion == "B"
        fsRatio = curFs / 44100;
        spectrogram(outputI,512*fsRatio,64*fsRatio,512*fsRatio, curFs)
        view(90, -90)
        xlim([0, 22.050])
        title("Interpolation")
    end

    if origVersion == "B"
        subplot(212)
    end
    if origVersion == "D" || origVersion == "B"
%         subplot(211)
        spectrogram(outputD,512,64,512, curFs, 'yaxis');
%         subplot(212)
%         plot(outputD)
        % Create figure
    %     figure1 = figure;
    % 
    %     % Create axes
    %     axes1 = axes('Parent',figure1);
    %     hold(axes1,'on');
    % 
    %     % Create surf
    %     surf(20*log10(abs(S)),'Parent',axes1,'LineStyle','none',...
    %         'EdgeColor','none');
    % 
    %     % Create xlabel
    %     xlabel('Time (s)');
    % 
    %     % Create ylabel
    %     ylabel('Frequency (kHz)');
    % 
    %     % Uncomment the following line to preserve the X-limits of the axes
    %     % xlim(axes1,[0 22.05]);
    %     % Uncomment the following line to preserve the Y-limits of the axes
    %     % ylim(axes1,[0.00580498866213152 9.99183673469388]);
    %     view(axes1,[0 90]);
    %     grid(axes1,'on');
    %     hold(axes1,'off');
    %     % Set the remaining axes properties
    %     set(axes1,'CLim',[-156.53559774527 -18.2268523216544],'TickDir','out');
    %     % Create colorbar
    %     colorbar(axes1,'TickLabels',{'-140','-120','-100','-80','-60','-40','-20'});
    %     title("Dynamic Grid")
    end
else
    if curFs == 44100
        figure('Position', [200, 200, 600, 350])
    end
    if origVersion == "I"
        outfft = fft(outputI);
        lengthSound = length(outputI);

    end
    if origVersion == "D" || origVersion == "B"
        outfft = fft(outputD);
        lengthSound = length(outputD);
    end
    outfftInDb = 20 * log10(abs(outfft));
    [pks,locs] = findpeaks(outfftInDb, 'Threshold',0.1, 'MinPeakHeight',-10);
    locs = locs - 1;
    triWidth = 500;
    triHeight = 4;
    if curFs == 44100
        rangeFFT = 0:lengthSound-1;
        
        fftPlot44100 = plot(rangeFFT(~isinf(outfftInDb))'*curFs/lengthSound, outfftInDb(~isinf(outfftInDb)), 'k', 'Linewidth', 2);
        locsSave = locs;
        hold on;
        for i = 1:length(locs)
            xTri = [locs(i) locs(i)-triWidth locs(i)-triWidth];%x coordinates of vertices
            yTri = [pks(i) pks(i)-triHeight*0.5 pks(i)+triHeight*0.5];%y coordinates of vertices
            patch(xTri,yTri,'k', 'EdgeColor', 'none') %plotting triangle in white color
        end
    else
        greyVal = 0.6;
        fftPlot88200 =plot([0:lengthSound-1]'*curFs/lengthSound, 20 * log10(abs(outfft)), ':', 'color', [greyVal, greyVal, greyVal], 'Linewidth', 2);
        for i = 1:length(locs)
            xTri = [locs(i) locs(i)+triWidth locs(i)+triWidth];%x coordinates of vertices
            yTri = [pks(i) pks(i)-triHeight*0.5 pks(i)+triHeight*0.5];%y coordinates of vertices
            patch(xTri,yTri, [greyVal, greyVal, greyVal],  'EdgeColor', 'none') %plotting triangle in white color
        end

    end
    % plot([0:lengthSound-1]'*curFs/lengthSound, 20 * log10(abs(outfftI)), 'Linewidth', 2);
    xlim([0 0.5 * 44100])
    ylim([-25, 100])
    if curFs == 88200
        legend([fftPlot44100, fftPlot88200], ["$f_\textrm{\fontsize{7}{0}\selectfont s} = 44100 \quad N = 15.5$", ...
        "$f_\textrm{\fontsize{7}{0}\selectfont s} = 88200 \quad N = 31$"], ...
        'interpreter', 'latex')
    end
    grid on
    xlabel("$f$ (in Hz)", 'interpreter', 'latex')
    ylabel("Magnitude (in dB)", 'interpreter', 'latex')
    title("Spectra of systems with different sample rates")
    set(gca, 'Fontsize', 16, 'Linewidth', 2)
    set(gcf, 'Color', 'w')
    % 
    % c/2
    % hold on;
    % plot(outFree)
    % ylim([-0.2, 0.2])

    % soundsc(output(1:(curFs / 44100):end), 44100);
end
