clear;clc;close all
%% Problem 2
seq1 = 2*ceil((rand(100,1)-.5))-1;
seq2 = 2*ceil((rand(100,1)-.5))-1;

% Part A - Plot Histogram of Each Sequence
subplot(1,2,1)
histogram(seq1)
grid on
title("Histogram of Sequence 1")
xlabel("Value")
ylabel("# of Occurances")
subplot(1,2,2)
histogram(seq2)
grid on
title("Histogram of Sequence 2")
xlabel("Value")
ylabel("# of Occurances")
snapnow

% Part B - Plot Spectral Analysis of Each Sequence
figure
subplot(2,1,1)
plot(abs(fft(seq1)))
grid on
title("Spectral Analysis of Sequence 1")
xlabel("Frequency (hz)")
ylabel("Amplitude")
subplot(2,1,2)
plot(abs(fft(seq2)))
grid on
title("Spectral Analysis of Sequence 2")
xlabel("Frequency (hz)")
ylabel("Amplitude")
snapnow

% Part C - Plot Autocorrelation for Each Sequence
[aCorr1,shift1] = autoCorrelation(seq1,"normalized");
[aCorr2,shift2] = autoCorrelation(seq2,"normalized");

figure
subplot(2,1,1)
plot(shift1,aCorr1)
grid on
title("Autocorrelation of Sequence 1 (Normalized)")
xlabel("Shift")
ylabel("Correlation")
subplot(2,1,2)
plot(shift2,aCorr2)
grid on
title("Autocorrelation of Sequence 2 (Normalized)")
xlabel("Shift")
ylabel("Correlation")
snapnow

% Part D - Plot Cross Correlation Between the Two Sequences
[xCorr,shift] = crossCorr(seq1,seq2,"normalized");

figure
plot(shift,xCorr)
grid on
title("Cross Correlation of Sequence 1 & 2 (Normalized)")
xlabel("Shift")
ylabel("Correlation")
snapnow

%% Problem 2 Bonus 1
seq1 = 2*ceil((rand(1000,1)-.5))-1;
seq2 = 2*ceil((rand(1000,1)-.5))-1;

% Part A - Plot Histogram of Each Sequence
subplot(1,2,1)
histogram(seq1)
grid on
title("Histogram of Sequence 1")
xlabel("Value")
ylabel("# of Occurances")
subplot(1,2,2)
histogram(seq2)
grid on
title("Histogram of Sequence 2")
xlabel("Value")
ylabel("# of Occurances")
sgtitle('2a Bonus')
snapnow

% Part B - Plot Spectral Analysis of Each Sequence
figure
subplot(2,1,1)
plot(abs(fft(seq1)))
grid on
title("Spectral Analysis of Sequence 1")
xlabel("Frequency (hz)")
ylabel("Amplitude")
subplot(2,1,2)
plot(abs(fft(seq2)))
grid on
title("Spectral Analysis of Sequence 2")
xlabel("Frequency (hz)")
ylabel("Amplitude")
sgtitle('2b Bonus')
snapnow

% Part C - Plot Autocorrelation for Each Sequence
[aCorr1,shift1] = autoCorrelation(seq1,"normalized");
[aCorr2,shift2] = autoCorrelation(seq2,"normalized");

figure
subplot(2,1,1)
plot(shift1,aCorr1)
grid on
title("Autocorrelation of Sequence 1 (Normalized)")
xlabel("Shift")
ylabel("Correlation")
subplot(2,1,2)
plot(shift2,aCorr2)
grid on
title("Autocorrelation of Sequence 2 (Normalized)")
xlabel("Shift")
ylabel("Correlation")
sgtitle('2c Bonus')
snapnow

% Part D - Plot Cross Correaltion Between the Two Sequences
[xCorr,shift] = crossCorr(seq1,seq2,"normalized");

figure
plot(shift,xCorr)
grid on
title("Bonus - Cross Correlation of Sequence 1 & 2 (Normalized)")
xlabel("Shift")
ylabel("Correlation")
snapnow

%% Question 1 Bonus 2

seq1 = 2*ceil((rand(100,1)-.5))-1;
for i = 1 : 100
    
    seq2 = 2*ceil((rand(100,1)-.5))-1;

    [xCorr,shift] = crossCorr(seq1,seq2);

    x_bar(i) = mean(xCorr);
    sigma(i) = std(xCorr);

end

figure
plot(1:100,x_bar)
hold on
grid on 
plot(1:100,sigma)
title("Monte Carlo Simulation Bonus")
legend(["Mean","Standard Deviation"])
xlabel("Sequence #")
ylabel("Value")

function [aCorr,shift] = autoCorrelation(seq,type)

    if nargin < 2
        type = "";
    end

    n = length(seq);
    m = 2*n-1;
    for i = 1 : n
        aCorr(i) = sum(seq(n-i+1:n) .* seq(1:i));
        aCorr(m+1-i) = aCorr(i); % autocorrelation is symmetric 
    end

    if type == "normalized"
        aCorr = aCorr/n;
    end

    shift = [-n+1:n-1];
end

function [xCorr,shift] = crossCorr(seq1,seq2,type)

    if nargin < 3
        type = "";
    end

    n = length(seq1);
    for i = -n+1 : n-1
        if i < 0
            xCorr(i+n) = sum(seq1(1:end+i) .* seq2(-i+1:end));
        elseif i > 0 
            xCorr(i+n) = sum(seq1(1:end-i) .* seq2(i+1:end));
        else
            xCorr(i+n) = sum(seq1 .* seq2);
        end
    end

    if type == "normalized"
        xCorr = xCorr/n;
    end
    
    shift = [-n+1:n-1];
end
