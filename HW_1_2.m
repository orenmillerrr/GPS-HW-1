clear;clc;close all

%% Problem 2

seq1 = 2*ceil((rand(100,1)-.5))-1;
seq2 = 2*ceil((rand(100,1)-.5))-1;

% Part A - Plot Histogram of Each Sequence

figure
subplot(1,2,1)
histogram(seq1)
title("Histogram of Sequence 1")
xlabel("Value")
ylabel("# of Occurances")
subplot(1,2,2)
histogram(seq2)
title("Histogram of Sequence 2")
xlabel("Value")
ylabel("# of Occurances")

% Part B - Plot Spectral Analysis of Each Sequence

figure
subplot(2,1,1)
plot(abs(fft(seq1)))
title("Spectral Analysis of Sequence 1")
xlabel("Frequency (hz)")
ylabel("Amplitude")
subplot(2,1,2)
plot(abs(fft(seq2)))
title("Spectral Analysis of Sequence 2")
xlabel("Frequency (hz)")
ylabel("Amplitude")

% Part C - Plot Autocorrelation for Each Sequence
[aCorr1,shift1] = autoCorrelation(seq1);
[aCorr2,shift2] = autoCorrelation(seq2);

figure
subplot(2,1,1)
plot(shift1,aCorr1)
title("Autocorrelation of Sequence 1 (Unbiased)")
xlabel("Shift")
ylabel("Correlation")
subplot(2,1,2)
plot(shift2,aCorr2)
title("Autocorrelation of Sequence 2 (Unbiased)")
xlabel("Shift")
ylabel("Correlation")

function [aCorr,shift] = autoCorrelation(seq)

    n = length(seq);
    m = 2*n-1;
    for i = 1 : n
        aCorr(i) = sum(seq(n-i+1:n) .* seq(1:i));
        aCorr(m+1-i) = aCorr(i); % autocorrelation is symmetric 
    end

    aCorr = aCorr/n;

    shift = [-n+1:n-1];
end

% Part D - Plot Cross Correlation Between the Two Sequences

xCorr = crossCorr(seq1,seq2);
figure
plot(xCorr)
title("Cross Correlation of Sequence 1 & 2")
xlabel("Shift")
ylabel("Correlation")

function xCorr = crossCorr(seq1,seq2)

    n = length(seq1);
    for i = -n : n
        if i < 0
            xCorr(i+n+1) = sum(seq1(1:end+i) .* seq2(-i+1:end));
        elseif i > 0 
            xCorr(i+n+1) = sum(seq1(1:end-i) .* seq2(i+1:end));
        else
            xCorr(i+n+1) = sum(seq1 .* seq2);
        end
    end
end

%% Problem 2 Bonus
seq1 = 2*ceil((rand(1000,1)-.5))-1;
seq2 = 2*ceil((rand(1000,1)-.5))-1;

% Part A - Plot Histogram of Each Sequence

subplot(1,2,1)
histogram(seq1)
title("Histogram of Sequence 1")
xlabel("Value")
ylabel("# of Occurances")
subplot(1,2,2)
histogram(seq2)
title("Histogram of Sequence 2")
xlabel("Value")
ylabel("# of Occurances")
sgtitle('2a Bonus')

% Part B - Plot Spectral Analysis of Each Sequence

figure
subplot(2,1,1)
periodogram(seq1)
title("Spectral Analysis of Sequence 1")
xlabel("Frequency (hz)")
ylabel("Amplitude")
subplot(2,1,2)
periodogram(seq1)
title("Spectral Analysis of Sequence 2")
xlabel("Frequency (hz)")
ylabel("Amplitude")
sgtitle('2b Bonus')

% Part C - Plot Autocorrelation for Each Sequence

[aCorr1,shift1] = autoCorrelation(seq1);
[aCorr2,shift2] = autoCorrelation(seq2);

figure
subplot(2,1,1)
plot(shift1,aCorr1)
title("Autocorrelation of Sequence 1 (Unbiased)")
xlabel("Shift")
ylabel("Correlation")
subplot(2,1,2)
plot(shift2,aCorr2)
title("Autocorrelation of Sequence 2 (Unbiased)")
xlabel("Shift")
ylabel("Correlation")
sgtitle('2c Bonus')

% Part D - Plot Cross Correaltion Between the Two Sequences

xCorr = xcorr(seq1,seq2);
figure
plot(xCorr)
title("Bonus - Cross Correlation of Sequence 1 & 2")
xlabel("Shift")
ylabel("Correlation")