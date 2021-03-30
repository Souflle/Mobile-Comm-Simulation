%% Parameter Setteing

clear; clc;
T = 2000;
fs = 1000;              % Sampling Frequency
t = 1/fs: 1/fs: 2;
              
f_D = [1, 10, 100];     % Max Doppler frequency 
N = 34;
sigma = sqrt(1 / 2);
%% Jake Fading Model
% Received signal : $r\left(t\right)=\textrm{Re}\left\lbrace u\left(t\right)\cdot 
% e^{\;\textrm{j2}\pi f_c t} \cdot \left\lbrack \sum_{n=0}^N \alpha_{n\;} \left(t\right)e^{-j\phi_n 
% \left(t\right)} \right\rbrack \right\rbrace =\textrm{Re}\left\lbrace T\left(t\right)e^{\textrm{j2}\pi 
% f_{c\;} t} \right\rbrace$
% 
% Phase of each path : $\phi {\;}_n \left(t\right)=-2\pi f_D t\cdot \cos \left(\alpha_n 
% \right)-\beta \;{\;}_n$
%% 
% * Assume Arrival Angle  =  Uniformly Distributed : $\alpha {\;}_n =\frac{2\pi 
% \;n}{\;N}$
%% 
% Fading Channel: $T\left(t\right)=\frac{E_0 }{\;\sqrt{\;N}}\left\lbrace \sqrt{\;2}\;\sum_{n=1}^{N_0 
% } \left\lbrack e^{j\left(2\pi \;f_D t\cdot \cos \left(\alpha {\;}_n \right)+\phi_n 
% \right)} +e^{-j\left(2\pi \;f_D t\cdot \cos \left(\alpha {\;}_n \right)+\phi_{-n} 
% \right)} \right\rbrack +e^{j\left(2\pi \;f_D t+\phi_N \right)} {+e}^{-j\left(2\pi 
% \;f_D t+{\phi -}_N \right)} \right\rbrace$ 
%% 
% * $N_0 \;=\;\frac{1}{2}\left(\frac{N}{2}-1\right)$
% * Uniform Distribution of Initial Phase : $\beta {\;}_n =\frac{\pi n\;}{N_0 
% }$
% * Doppler term: $e^{j\phi_n } =e^{-j\phi_{-n} } \Rightarrow {\frac{\sqrt{\;2}}{2}e}^{j\beta 
% {\;}_n } =\frac{\sqrt{\;2}}{2}\left\lbrack \cos \left(\beta_n \right)+\textrm{jsin}\left(\beta 
% {\;}_n \right)\right\rbrack$
% * Max Doppler term:  $e^{j\phi_N } =e^{-j\phi_{-N} } \Rightarrow {\frac{\sqrt{\;2}}{2}e}^{j\alpha 
% \;} =\frac{\sqrt{\;2}}{2}\left\lbrack \cos \left(\alpha \;\right)+\textrm{jsin}\left(\alpha 
% \;\right)\right\rbrack$
%% 
% In-phase Component : $h_c \left(t\right)=2\sum_{n=1}^{N_0 } \cos \beta {\;}_n 
% \cos \left(2\pi \left(f_D \cos {\;\alpha }_n \right)\cdot t\right)+\sqrt{\;2}\cos 
% \left(\alpha \right)\cos \left(2\pi f_D t\right)\;$
% 
% Quadratic Component: $h_s \left(t\right)=2\sum_{n=1}^{N_0 } \sin \beta {\;}_n 
% \cos \left(2\pi \left(f_D \cos {\;\alpha }_n \right)\cdot t\right)+\sqrt{\;2}\sin 
% \left(\alpha \right)\cos \left(2\pi f_D t\right)\;$

N_0 = 1 / 2 * (N / 2 - 1);
h_c = zeros(3, T);
h_s = zeros(3, T);

for e = 1 : 3
for n = 1 : N_0
    beta_n = pi * n / N_0;
    alpha_n = 2 * pi * n / N;
    % alpha = pi / 4
    h_c(e, :) = h_c(e, :) + 2 * cos(beta_n) * cos(2 * pi * f_D(e) * cos(alpha_n) * t) + ...
                cos(2*pi*f_D(e) * t);
    h_s(e, :) = h_s(e, :) + 2 * sin(beta_n) * cos(2 * pi * f_D(e) * cos(alpha_n) * t) + ...
                cos(2*pi*f_D(e) * t);
end
end

h_c = h_c * sigma / sqrt(N);
h_s = h_s * sigma / sqrt(N);

%% Plot

h_env = sqrt(h_c.^2 + h_s.^2);
h_norm = h_env ./ mean(h_env, 2);
h_norm_dB = 10 * log10(h_norm) ;

subplot(3, 1, 1)
plot(t, h_norm_dB(1, :))
ylabel('Envelope(dB)')
title('f_D = 1 Hz')

subplot(3, 1, 2)
plot(t, h_norm_dB(2, :))
ylabel('Envelope(dB)')
title('f_D = 10 Hz')

subplot(3, 1, 3)
plot(t, h_norm_dB(3, :))
ylabel('Envelope(dB)')
title('f_D = 100 Hz')
%% Auto-Correlation

auto_corr_1 = xcorr(h_c(1, :), 'coeff');
auto_corr_2 = xcorr(h_c(2, :), 'coeff');
auto_corr_3 = xcorr(h_c(3, :), 'coeff');

r_1 = [auto_corr_1(2000:3999); besselj(0, 2*pi*f_D(1) * t)];
r_2 = [auto_corr_2(2000:3999); besselj(0, 2*pi*f_D(2) * t)];
r_3 = [auto_corr_3(2000:3999); besselj(0, 2*pi*f_D(3) * t)];

figure
subplot(3, 1, 1)
plot(t, r_1)
title('f_D = 1 Hz')
subplot(3, 1, 2)
plot(t, r_2)
title('f_D = 10 Hz')
subplot(3, 1, 3)
plot(t, r_3)
title('f_D = 100 Hz')