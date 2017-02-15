clear

%% Set parameters.
num_samples = 1024; % Number of samples
N = num_samples / 2;
L = 100; % Lenth of signal
Ls = L/num_samples; % Time step
Fs = 1/L; % Freq. step
v = Fs * (-N:N-1); % Freq. axis
x = Ls * (-N:N-1); % Time axis
lambda = 1; % Wavelength. Treated as unit
k = 2 * pi / lambda; % Wavenumber
z = 0:0.01:0.5 * lambda; % Propagation along z axis
log_z = logspace(0, 2, 3); % z axis in log-scale

%% Plot the 2D function exp(i \phi) and its fft2.
f = @(x, y) exp(1i * atan(y/x));
% Produce the wanted signal. Since the range of arctan is [-pi/2, pi/2], one
% needs to extend the function so that exp(i \phi) can be properly sampled.
field = zeros(num_samples, num_samples);
for i = 1:num_samples
    for j = 1:num_samples
        if i > N + 1
            field(i, j) = f((i - N - 1), (j - N - 1));
        elseif i < N + 1 && j > N + 1
            field(i, j) = f((i - N - 1), (j - N - 1)) * exp(1i * pi);
        elseif i < N + 1 && j < N + 1
            field(i, j) = f((i - N - 1), (j - N - 1)) * exp(-1i * pi);
        else
            field(i, j) = exp(1i * pi/2 * sign(j - N - 1));
        end
    end
end
Field = fftshift(fft2(field));

%% Calculate the field expression after propagation along z.
g = @(vx, vy) sqrt(k^2 - 4 * pi^2 * (vx.^2 + vy.^2));

mu = zeros(num_samples, num_samples);
for i = 1:num_samples
    for j = 1:num_samples
        mu(i, j) = g(v(i), v(j));
    end
end

% Field function corresponds to different z
field_z = zeros(num_samples, num_samples, length(z));
for i = 1:length(z)
    field_z(:, :, i) =  ifft2(ifftshift(exp(1i * mu * z(i)) .* Field));
end

% Version of log_z
log_field_z = zeros(num_samples, num_samples, length(log_z));
for i = 1:length(log_z)
    log_field_z(:, :, i) =  ifft2(ifftshift(exp(1i * mu * log_z(i)) .* Field));
end

%% Analyse the phase
% These are x, y values corresponding to r = 0.1, 0.3, 0.5, 1 lambda.
r_x = N + 1 + [1 3 5 10];
r_y = N + 1 + [0 0 1 2];

Phi = zeros(length(r_x), length(z));
dip = zeros(length(r_x), length(z));
for i = 1:length(r_x)
    for j = 1:length(z)
    I = angle(field_z(:, :, j) ./ field);
    Phi(i, j) = I(r_x(i), r_y(i));
    dip(i, j) = k * z(j) - Phi(i, j);
    % I is exactly the I in the paper, which depends on r and z.
    % To extract the data at r^2 = x^2 + y^2 = 1, for example,
    % one needs to set n^2 + m^2 = Fs^2, which corresponds approximately to the sample
    % point (N + 11, N + 1) = (N + 10, N + 2) = ... = (N + 1, N + 11).
    end
end

%% Archive
clear f g i j mu
save('data.mat')

