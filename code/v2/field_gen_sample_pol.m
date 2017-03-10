clear

%% Set parameters.
num_samples = 1024; % Number of samples
N = num_samples / 2;
L = 1e-5; % Radius of laser beam, unit in meter
Ls = L/num_samples; % Time step
Ds = 2 * pi /num_samples; % Freq. step
lambda = 633e-9; % Wavelength of He-Ne Laser (633 nm)
k = 2 * pi / lambda; % Wavenumber
%z = (0:0.01:0.49) * lambda; % Propagation along z axis
z = logspace(0, 2, 3) * lambda; % z axis in log-scale

r = Ls * (0:num_samples - 1); % Freq. axis
phi = Ds * (-N:N-1); % Time axis
[R, PHI] = meshgrid(r, phi);
X = R .* cos(PHI);
Y = R .* sin(PHI);

W0 = 3.88e-4; % Beam waist 0.388 mm
z0 = pi * W0^2 / lambda;

%% Initial field generation
% Aperture field
f_init = @(r, phi) 2 .* r ./ (sqrt(pi) * W0^2) .* exp(- r.^2 ./ W0^2 + 1i .* phi);

% LG(1, 0) field
f_ref = @(r, phi, z) 2 .* r ./ (sqrt(pi) * W0^2 .* (1 + (z/z0)^2)) ...
    .* exp(- r^2 ./ (W0^2 .* (1 + (z/z0)^2))) .*  exp(1i .* phi) ...
    .* exp(1i .* (2.*atan(z/z0) + 0.5 .* k .* r^2 .* z ./ (z^2 + z0^2)) );

field = f_init(R, PHI);

figure
surf(X, Y, abs(field));
shading interp

figure
surf(X, Y, angle(field));
shading interp

%% Convolution
%h = @(x, y, z) exp(1i .* k .* (x^2 + y^2) ./ (2*z)) / (1i * lambda * z);
%fresnel = h(X, Y, );

% How to convert field to cart. cor.?
for i=1:length(z)
   %field_z = conv2 
end

%% Reference field
field_ref = zeros(num_samples, num_samples, length(z));
for i = 1:length(z)
   field_ref(:, :, i) = f_ref(R, PHI, z(i)); 
end

figure
subplot(1, 2, 1)
surf(X, Y, abs(field_ref(:, :, 1)));
shading interp

subplot(1, 2, 2)
surf(X, Y, angle(field_ref(:, :, 1)));
shading interp



