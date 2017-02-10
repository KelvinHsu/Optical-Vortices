%% Set signal length T = 100 and 2N + 1 sample points.
N = 500;
T = 100; % Lenth of signal
Ts = T/(2*N + 1); % Time step
Fs = 1/Ts; % Sampling freq.

%% Plot sinc
sig = zeros(1, 2*N + 1);
for j = 1:2*N+1
    sig(j) = sinc((-N+j-1)*Ts);
end
Sig = fftshift(fft(sig));

f = Fs * (-N:N)/N;

figure
subplot(1, 2, 1)
plot((-N:N)*Ts, sig)
title('Time-domain sinc function')
xlabel('Time (sec)')

subplot(1, 2, 2)
plot(f, abs(Sig))
title('Spectrum of sinc function')
xlabel('Frequency (Hz)')

%% Plot the 2D function exp(i \phi) and its fft2
f = @(x, y) exp(1i * atan(y/x));

sig2 = zeros(2*N + 1, 2*N + 1);
for j = 1:2*N+1
    for k = 1:2*N+1
        if j ~= N + 1
            sig2(j, k) = f((j - N - 1), (k - N - 1));
        else
            sig2(j, k) = exp(1i * pi/2 * sign(k - N - 1));
        end
    end
end
Sig2 = fftshift(fft2(sig2));

t = 1:2*N+1;
s = 1:2*N+1;

figure
subplot(1, 2, 1)
plot3((t - N - 1) * Ts, (s - N - 1) * Ts, abs(sig2(t, s)))
title('Amplitude (Time)')
xlabel('Time (sec)')
ylabel('Time (sec)')

subplot(1, 2, 2)
plot3((t - N - 1) * Ts, (s - N - 1) * Ts, angle(sig2(t, s)))
title('Phase angle (Time)')
xlabel('Time (sec)')
ylabel('Time (sec)')

figure
subplot(1, 2, 1)
plot3((t - N - 1)/N * Fs, (s - N - 1)/N * Fs, abs(Sig2(t, s)))
title('Amplitude (Frequency)')
xlabel('Frequency (Hz)')
ylabel('Frequency (Hz)')

subplot(1, 2, 2)
plot3((t - N - 1)/N * Fs, (s - N - 1)/N * Fs, angle(Sig2(t, s)))
title('Phase angle (Frequency)')
xlabel('Frequency (Hz)')
ylabel('Frequency (Hz)')