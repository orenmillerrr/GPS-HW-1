clear;clc;close all
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

