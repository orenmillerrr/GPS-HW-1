clear;clc;close all

%% Problem 5
freq = 1;
t = -1 : 0.0001 : 1;
func = sin(2*pi*freq*t);

[aCorr,shift] = autoCorr(func);
aCorr_norm = aCorr/max(aCorr);

figure
subplot(3,1,1)
plot(t,func)
title("1 Hz Sine Wave")
grid on
xlabel("Time (sec)")
xlim([-1 1])
ylabel("Amplitude")
subplot(3,1,2)
plot(shift,aCorr)
title("Autocorrelation of 1 Hz Sine Wave (Unbiased)")
grid on
xlabel("Shift")
xlim([shift(1) shift(end)])
ylabel("Correlation")
subplot(3,1,3)
plot(shift,aCorr_norm)
title("Normalized Autocorrelation of 1 Hz Sine Wave")
grid on
xlabel("Shift")
xlim([shift(1) shift(end)])
ylabel("Correlation")

function [aCorr,shift] = autoCorr(seq)

    n = length(seq);
    m = 2*n-1;
    for i = 1 : n
        aCorr(i) = sum(seq(n-i+1:n) .* seq(1:i));
        aCorr(m+1-i) = aCorr(i); % autocorrelation is symmetric 
    end

    aCorr = aCorr;  % unbiased correlation

    shift = [-n+1:n-1];
end