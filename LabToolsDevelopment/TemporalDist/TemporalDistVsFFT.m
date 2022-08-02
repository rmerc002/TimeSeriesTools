function TemporalDistVsFFT(ts, mm, Fs, TD)
    if nargin < 4
        TD = TemporalDist(ts, mm);
    end

    L = length(ts);
%     f = Fs./linspace(1,ceil(L/2)); %Fs*(0:(L/2))/L;
    Y = fft(ts);
%     P2 = abs(Y/L);
%     P1 = P2(1:length(f));
%     P1(2:end-1) = 2*P1(2:end-1);
% 
%     X = Fs./f;

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;

    figure;
    tiledlayout(2,1);
    ax1 = nexttile();
    plot(TD);

    ax2 = nexttile();
%     plot(X(end:-1:1),P1(end:-1:1)) 
    plot(f,P1(end:-1:1)) 

    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')

%     linkaxes([ax1, ax2],'x');

end