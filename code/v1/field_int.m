%%%%%%% Simulation by numerical integral %%%%%%%
%% Parameters
lambda = 672e-9; % Wavelength, unit in meters.
k = 2 * pi / lambda; % Wavenumber
r = linspace(0, 16, 100) * lambda; % r-coordinates
theta = linspace(0, 2*pi, 100); % theta-coordinate
z = logspace(0, 2, 3) * lambda; % Propagation distance
[R, THETA] = meshgrid(r, theta);
rho_max = 10 / lambda; % Integration upper bound

%% Generate the plot
figure
for i = 1:length(z)
    f = @(rho) rho.^(-1) * exp(1i .* z(i) .* sqrt(k^2 - 4 * pi^2 .* rho.^2)) ...
        .* besselj(1, 2 * pi .* rho * r);
    I = integral(f, 0, rho_max, 'ArrayValued', true);
    [I, ~] = meshgrid(I);
    subplot(3, 2, -1 + 2*i)
    surf(R .* cos(THETA) / lambda, R .* sin(THETA) / lambda, angle(I))
    shading interp
    colorbar
    xlabel('$x$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
    ylabel('$y$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
    zlabel('$\Phi$ (rad)', 'interpreter', 'LaTeX')
    
    subplot(3, 2, 2*i)
    pcolor(R .* cos(THETA) / lambda, R .* sin(THETA) / lambda, ...
        atan2(R .* sin(THETA), R .* cos(THETA)) + angle(I))
    shading interp
    colorbar
    xlabel('$x$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
    ylabel('$y$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
end