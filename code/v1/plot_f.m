lambda = 1; % Wavelength, unit in meters.
k = 2 * pi / lambda; % Wavenumber
z = logspace(0, 2, 3) * lambda; % Propagation distance
rho_max = 2 / lambda; % Integration upper bound
r = 0.1 * lambda;

%% Plot f
f = @(rho) (rho.^(-1) .* exp(1i .* z(2) .* sqrt(k^2 - 4 .* pi^2 .* rho.^2)) ...
        .* besselj(1, 2 .* pi .* rho .* r));
fplot(f, [0 rho_max])
xlabel('$\rho$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')