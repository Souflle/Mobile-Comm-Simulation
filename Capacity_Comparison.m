%% 
% *Parameter Setting*

clear;  clc;    
SNRdB = 1:30;

M = 10000;                                  % sample #
N = 3;                                      % channel #
N0 = 1;                                     % Noise Power 
AWGN = zeros(N, length(SNRdB)); 
CSIR_only = zeros(N, length(SNRdB)); 
CSIT = zeros(N, length(SNRdB));
truncated_inv = zeros(N, length(SNRdB));
zero_outage = zeros(N, length(SNRdB));
%% Channel Model
% $E\left\lbrace |h|^2 \right\rbrace =1$ for comparison
%% 
% # Rayleigh Fading with $\sigma {\;}^2 =\frac{1}{2}$
% # Nakagami Fading with m =  2, $\sigma {\;}^2 =\frac{1}{2}$  ( ~ 2 antenna 
% diversity)
% # Log Normal Shadowing

h = zeros(N, M);
h(1, :) = [random('Rayleigh', sqrt(1/2), [1, M-1]), 0];
h(2, :) = random('Nakagami', 2, 1, [1, M]);
h(3, :) = random('Lognormal', -0.65, 0.8, [1, M]);

channel_gain = mean(abs(h).^2,2)
%% 
% 
%% Ergodic Capacity of Channel
% $$C=\int_0^{\infty } B\;\log_2 \left(1+\gamma \;\right)\Pr \left(\gamma \;\right)d\gamma 
% \;=E\left\lbrace {\textrm{Blog}}_2 \left(1+\gamma \;\right)\right\rbrace$$
% 
% $$\gamma =\textrm{SNR}\cdot \frac{\left|h{\left|\right.}^2 \right.}{\;E\left\lbrace 
% \left|h{\left|\right.}^2 \right.\right\rbrace }$$
% 
% 

for n = 1:3
for snrdb = 1:30
snr = 10 ^ (snrdb/10);
    
AWGN(n, snrdb) = log2(1 + snr);  

gamma =  abs(h(n, :)).^2 * snr / mean(abs(h(n, :)).^ 2);

% Waterfilling
P_adapt = waterfill(M, 1 ./ gamma(1: M-1));
CSIT(n, snrdb) = mean(log2(1 + gamma(1:M-1) .* P_adapt));

% Constant Power 
CSIR_only(n, snrdb) = mean(log2(1 + gamma));

% Truncated Inversion
Pout = zeros(1, 500);
sigma_0 = zeros(1, 500);
capacity = zeros(1, 500);
for k = 1: 500
    thr = 10 ^ (snrdb/10 - 2 + k * 0.004 ); 
    
    Pout(k) = numel(gamma(gamma > thr))/ 10000;
    sigma_0(k) = 1 ./ mean(1 ./ gamma(gamma > thr));
    capacity(k)= Pout(k) * log2(1+ sigma_0(k)); 
end
truncated_inv(n, snrdb) = max(capacity);



sigma = 1 ./ mean(1 ./ gamma);
zero_outage(n, snrdb) = log2(1+sigma); 
end
end

%% Plot Capacity

capacity_rayleigh = [AWGN(1, :); CSIT(1, :); CSIR_only(1, :); truncated_inv(1, :); zero_outage(1, :)];
capacity_nakagami = [AWGN(2, :); CSIT(2, :); CSIR_only(2, :); truncated_inv(2, :); zero_outage(2, :)];
capacity_lognormal = [AWGN(3, :); CSIT(3, :); CSIR_only(3, :); truncated_inv(3, :); zero_outage(3, :)];


markers = {'none', 'o', '*', 'diamond', 'pentagram'};
name = {'AWGN', 'Waterfiling', 'CSIR only', 'Truncated Inversion', 'Channel Inversion'};

plt = plot(SNRdB, capacity_rayleigh);
set(plt, {'Marker'}, markers(:))
title('Rayleigh Fading')
xlabel('SNR(dB)')
ylabel('C / B (bps)')
legend(name, 'Location', 'northwest')

figure
plt2 = plot(SNRdB, capacity_nakagami);
set(plt2, {'Marker'}, markers(:))
title('Nakagami Fading (m=2)')
xlabel('SNR(dB)')
ylabel('C / B (bps)')
legend(name, 'Location', 'northwest')

figure
plt3 = plot(SNRdB, capacity_lognormal);
set(plt3, {'Marker'}, markers(:))
title('LogNormal Distribution')
xlabel('SNR(dB)')
ylabel('C / B (bps)')
legend(name, 'Location', 'northwest')