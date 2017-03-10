%% Set signal length L = 100 and num_samples sample points.
num_samples = 1024;
N = num_samples / 2;
L = 100; % Lenth of signal
Ls = L/num_samples; % Time step
Fs = 1/L; % Sampling freq.

%% Plot sinc
sig = zeros(1, num_samples);
for j = 1:num_samples
    sig(j) = sinc((-N+j-1)*Ls);
end
Sig = fftshift(fft(sig));

v = Fs * (-N:N-1);
x = Ls * (-N:N-1);
figure
subplot(1, 2, 1)
plot(x, sig)
title('Time-domain sinc function')
xlabel('Time (sec)')

subplot(1, 2, 2)
plot(v, abs(Sig))
title('Spectrum of sinc function')
xlabel('Frequency (Hz)')