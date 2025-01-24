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

aCorr1 = autoCorrelation(seq1);
aCorr2 = autoCorrelation(seq2);

figure
subplot(2,1,1)
plot(aCorr1)
title("Autocorrelation of Sequence 1")
xlabel("Shift")
ylabel("Correlation")
subplot(2,1,2)
plot(aCorr2)
title("Autocorrelation of Sequence 2")
xlabel("Shift")
ylabel("Correlation")

function aCorr = autoCorrelation(seq)

    n = length(seq);
    m = 2*n-1;
    for i = 1 : n
        aCorr(i) = sum(seq(n-i+1:n) .* seq(1:i));
        aCorr(m+1-i) = aCorr(i); % autocorrelation is symmetric 
    end
end

% Part D - Plot Cross Correaltion Between the Two Sequences

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

aCorr1 = xcorr(seq1,seq1);
aCorr2 = xcorr(seq2,seq2);

figure
subplot(2,1,1)
plot(aCorr1)
title("Autocorrelation of Sequence 1")
xlabel("Shift")
ylabel("Correlation")
subplot(2,1,2)
plot(aCorr2)
title("Autocorrelation of Sequence 2")
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

%% Problem 3
for i = 1 : 3
    A = 3 + 3*randn(1000,1);
    B = 5 + 5*randn(1000,1);
    C = A + B;
    D = 3*A + 4*B;
    E = 3*A - 4*B;

    DATA = [A B C D E];
    
    % Part A - Find Mean and Variance for A, B, C, D, and E
    
    meanA(i) = mean(A);
    meanB(i) = mean(B);
    meanC(i) = mean(C);
    meanD(i) = mean(D);
    meanE(i) = mean(E);
    
    varA(i) = var(A);
    varB(i) = var(B);
    varC(i) = var(C);
    varD(i) = var(D);
    varE(i) = var(E);
    
    meanDATA(i,:) = mean(DATA);
    
    covDATA(:,:,i) = cov(DATA);

    fprintf("\nSequence %g\n",i)
    fprintf("--------------------------\n")
    fprintf("Mean A = %.4f\n",meanA(i))
    fprintf("Mean B = %.4f\n",meanB(i))
    fprintf("Mean C = %.4f\n",meanC(i))
    fprintf("Mean D = %.4f\n",meanD(i))
    fprintf("Mean E = %.4f\n\n",meanE(i))

    fprintf("Var A = %.4f\n",varA(i))
    fprintf("Var B = %.4f\n",varB(i))
    fprintf("Var C = %.4f\n",varC(i))
    fprintf("Var D = %.4f\n",varD(i))
    fprintf("Var E = %.4f\n\n",varE(i))

    fprintf("Mean DATA = \n")
    fprintf("%.4f %.4f %.4f %.4f %.4f\n\n",meanDATA(i,:))

    fprintf("Cov DATA = \n")
    fprintf("%10.4f %10.4f %10.4f %10.4f %10.4f\n",covDATA(:,:,i))

end

subplot(2,1,1)
hold on 
grid on
plot(meanA,'.')
plot(meanB,'.')
plot(meanC,'.')
plot(meanD,'.')
plot(meanE,'.')
title("Mean for Each Sequence")
xlabel("Sequence Number")
xlim([0,4])
ylabel("Mean")
legend('A','B','C','D','E')
subplot(2,1,2)
hold on 
grid on
plot(varA,'.')
plot(varB,'.')
plot(varC,'.')
plot(varD,'.')
plot(varE,'.')
title("Variance for Each Sequence")
xlabel("Sequence Number")
xlim([0,4])
ylabel("Variance")
legend('A','B','C','D','E')

%% Problem 5
freq = 1;
t = -1 : 0.01 : 1;
func = sin(2*pi*freq*t);

aCorr = autoCorr(func);
aCorr_norm = aCorr/length(t);

figure
subplot(2,1,1)
plot(aCorr)
title("Autocorrelation of 1 Hz Sine Wave")
grid on
xlabel("Shift")
xlim([0 2*length(func)-1])
ylabel("Correlation")
subplot(2,1,2)
plot(aCorr_norm)
title("Normalized Autocorrelation of 1 Hz Sine Wave")
grid on
xlabel("Shift")
xlim([0 2*length(func)-1])
ylabel("Correlation")

function aCorr = autoCorr(seq)

    n = length(seq);
    m = 2*n-1;
    for i = 1 : n
        aCorr(i) = sum(seq(n-i+1:n) .* seq(1:i));
        aCorr(m+1-i) = aCorr(i); % autocorrelation is symmetric 
    end
end

