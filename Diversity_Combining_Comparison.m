%% Diversity Performance Comparision
%% 
% * Compare *BER versus SNR* between *SC, EGC, MRC*
% * Each Combining applied for *two-path of Rayleigh Fading* with same average 
% branch SNR 
% Rayleigh Realization

clear; clc;
realization = 1000000;
rayleigh1 = random('Rayleigh', 1/sqrt(2), realization, 1);  
rayleigh2 = random('Rayleigh', 1/sqrt(2), realization, 1);

% System Parameter

x = [0: 30];
SNR = 10 .^ (x/10);

SC_SER = zeros([30,1]); EGC_SER = zeros([30,1]);   MRC_SER = zeros([30,1]);
MRC_dBGain_EGC = zeros([30,1]);
for i = 1:30
  
branchSNR_1 = rayleigh1.^2 * SNR(i); 
branchSNR_2 = rayleigh2.^2 * SNR(i);
% Selection Combining

SC_SNR = max(branchSNR_1, branchSNR_2);
SC_SERinc = qfunc(sqrt(2 * SC_SNR));
SC_SER(i) = mean(SC_SERinc);
% Equal Gain Combining

EGC_SNR = 0.5 * branchSNR_1 + 0.5 * branchSNR_2 + sqrt(branchSNR_1 .* branchSNR_2);
EGC_SERinc = qfunc(sqrt(2 * EGC_SNR));
EGC_SER(i) = mean(EGC_SERinc);
% Maximal Ratio Combining

MRC_SNR = branchSNR_1 + branchSNR_2;
MRC_SERinc = qfunc(sqrt(2 * MRC_SNR));
MRC_SER(i) = mean(MRC_SERinc);
MRC_dBGain_EGC(i) = 10*log10(EGC_SER(i) / MRC_SER(i))  ;
end
%% Performance Plot

plt = semilogy(0:29, [SC_SER, EGC_SER, MRC_SER]);
legend({'SC', 'EGC', 'MRC'})
markers = {'o', 'diamond', 'pentagram'};
set(plt, {'Marker'}, markers(:))
xlabel('SNR in dB')
ylabel('Symbol Error Rate')

mean(MRC_dBGain_EGC) % We can find about 1dB gain by using MRC than EGC