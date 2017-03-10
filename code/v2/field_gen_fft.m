clear
close all

%% Set parameters.
% Wave parameters
lambda = 633e-9; % Wavelength of He-Ne Laser (633 nm)
k = 2 * pi / lambda; % Wavenumber
W0 = 3.88e-4; % Beam waist 0.388 mm
z0 = pi * W0^2 / lambda;

% Sampling parameters
num_samples = 1024; % Number of samples
N = num_samples / 2;
L = 1000e-6; % Width of laser beam, unit in meter
Ls = L/num_samples; % Sample step
Fs = 1/L; % Spectrum resolution (freq. step)
n = num_samples * lambda / L; % Samples per wavelength

% Coordinate parameters
x = Ls * (-N:N-1); % Space axis
v = Fs * (-N:N-1); % Spectrum axis
%z = (0:0.01:0.49) * lambda; % Propagation along z axis
z = logspace(0, 4, 3) * lambda; % z axis in log-scale
[X, Y] = meshgrid(x);
R = X.^2 + Y.^2;
PHI = atan2(Y, X);
[VX, VY] = meshgrid(v);
idx_eff = find(abs(x) <= L/2); % Indices of elements at which the field will not be truncated.

%% Initial field generation
% Aperture field expression
f_init = @(r, phi) 2 .* r ./ (sqrt(pi) * W0^2) .* exp(- r.^2 ./ W0^2 + 1i .* phi);

% LG(1, 0) field
f_ref = @(r, phi, z) 2 .* r ./ (sqrt(pi) * W0^2 .* (1 + (z/z0)^2)) ...
    .* exp(- r^2 ./ (W0^2 .* (1 + (z/z0)^2))) .*  exp(1i .* phi) ...
    .* exp(1i .* (2.*atan(z/z0) + 0.5 .* k .* r^2 .* z ./ (z^2 + z0^2)) );

% Spectrum transport
h = @(vx, vy, z) exp(-1i .* 2 .* pi.^2 .* z .* (vx.^2 + vy.^2) ./ k );
 
% Generate initial field
U0 = zeros(num_samples, num_samples);
temp = f_init(R, PHI);
U0(idx_eff, idx_eff) = temp(idx_eff, idx_eff);

%% Fresnel diffraction
A = zeros(num_samples, num_samples, length(z)); % spectrum at z
U = zeros(num_samples, num_samples, length(z)); % field at z

A0 = fftshift(fft2(U0));
for i = 1:length(z)
   % Spectrum at z
   A(:, :, i) = A0 .* h(VX, VY, z(i));
   % Field at z
   U(:, :, i) = ifft2(ifftshift(A(:, :, i)));
end

% Reference field LG(1, 0)
U_ref = zeros(num_samples, num_samples, length(z));
for i = 1:length(z)
   U_ref(:, :, i) = f_ref(R, PHI, z(i)); 
end

%% Plot figures
figure('units', 'pixels', 'position', [10, 10, 1200, 300])
title('Phase evolution')
for i = 1:length(z)
    subplot(1, length(z), i)
    pcolor(X ./ lambda, Y ./ lambda, angle(U(:, :, 1)));
    shading interp
    xlabel('$x$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
    ylabel('$y$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
    colorbar
end

% figure('units', 'pixels', 'position', [10, 10, 1200, 300])
% title('Phase evolution (reference)')
% for i = 1:length(z)
%     subplot(1, length(z), i)
%     pcolor(X ./ lambda, Y ./ lambda, angle(U_ref(:, :, 1)));
%     shading interp
%     xlabel('$x$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
%     ylabel('$y$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
%     colorbar
% end

figure('units', 'pixels', 'position', [10, 10, 1200, 300])
title('Amplitude evolution')
for i = 1:length(z)
    subplot(1, length(z), i)
    surf(X ./ lambda, Y ./ lambda, abs(U(:, :, 1)));
    shading interp
    xlabel('$x$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
    ylabel('$y$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
    colorbar
end

% figure('units', 'pixels', 'position', [10, 10, 1200, 300])
% title('Amplitude evolution (reference)')
% for i = 1:length(z)
%     subplot(1, length(z), i)
%     surf(X ./ lambda, Y ./ lambda, abs(U_ref(:, :, 1)));
%     shading interp
%     xlabel('$x$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
%     ylabel('$y$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
%     colorbar
% end




