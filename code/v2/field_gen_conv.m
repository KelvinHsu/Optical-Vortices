clear

%% Set parameters.
num_samples = 1024; % Number of samples
N = num_samples / 2;
L = 1e-5; % Sampling width, unit in meter
Ls = L/num_samples; % Sample step
x = Ls * (-N:N-1); % Spectrum resolution (freq. step)
lambda = 633e-9; % Wavelength of He-Ne Laser (633 nm)
k = 2 * pi / lambda; % Wavenumber
%z = (0:0.01:0.49) * lambda; % Propagation along z axis
z = logspace(0, 2, 3) * lambda; % z axis in log-scale

[X, Y] = meshgrid(x);
R = X.^2 + Y.^2;
PHI = atan2(Y, X);

W0 = 3.88e-4; % Beam waist 0.388 mm
z0 = pi * W0^2 / lambda;

%% Initial field generation
% Aperture field expression
f_init = @(r, phi) 2 .* r ./ (sqrt(pi) * W0^2) .* exp(- r.^2 ./ W0^2 + 1i .* phi);

% LG(1, 0) field
f_ref = @(r, phi, z) 2 .* r ./ (sqrt(pi) * W0^2 .* (1 + (z/z0)^2)) ...
    .* exp(- r^2 ./ (W0^2 .* (1 + (z/z0)^2))) .*  exp(1i .* phi) ...
    .* exp(1i .* (2.*atan(z/z0) + 0.5 .* k .* r^2 .* z ./ (z^2 + z0^2)) );

U0 = f_init(R, PHI);

h = @(x, y, z) exp(-1i .* k .* (x.^2 + y.^2) ./ (2*z)) / (1i .* lambda .* z);
fresnel = zeros(num_samples, num_samples, length(z));
for i = 1:length(z)
    fresnel(:, :, i) = h(X, Y, z(i));
end

U = zeros(num_samples, num_samples, length(z));
for i = 1:length(z)
    U = conv2(U0, fresnel(:, :, i), 'same');
end

%% Plot figures
figure('units', 'pixels', 'position', [10, 10, 1200, 300])
for i = 1:length(z)
    subplot(1, length(z), i)
    surf(X, Y, abs(U(:, :, 1)));
    shading interp
    xlabel('$x$ (m)', 'interpreter', 'LaTeX')
    ylabel('$y$ (m)', 'interpreter', 'LaTeX')
    colorbar
end

figure('units', 'pixels', 'position', [10, 10, 1200, 300])
for i = 1:length(z)
    subplot(1, length(z), i)
    surf(X, Y, angle(U(:, :, 1)));
    shading interp
    xlabel('$x$ (m)', 'interpreter', 'LaTeX')
    ylabel('$y$ (m)', 'interpreter', 'LaTeX')
    colorbar
end

