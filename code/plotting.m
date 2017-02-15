load('data.mat')

%% Plot space function and spectrum
figure
subplot(1, 2, 1)
surf(x, x, abs(field))
title('Amplitude (space)')
xlabel('$x$ (m)', 'interpreter', 'LaTeX')
ylabel('$y$ (m)', 'interpreter', 'LaTeX')
shading interp
colorbar

subplot(1, 2, 2)
surf(x, x, angle(field))
title('Phase angle (space)')
xlabel('$x$ (m)', 'interpreter', 'LaTeX')
ylabel('$y$ (m)', 'interpreter', 'LaTeX')
zlabel('$\phi$ (rad)', 'interpreter', 'LaTeX')
view(225, 45)
shading interp
colorbar

figure
subplot(1, 2, 1)
surf(v, v, abs(Field))
title('Amplitude (frequency)')
xlabel('$\nu_x$ (1/m)', 'interpreter', 'LaTeX')
ylabel('$\nu_y$ (1/m)', 'interpreter', 'LaTeX')
shading interp
colorbar

subplot(1, 2, 2)
surf(v, v, angle(Field))
title('Phase angle (frequency)')
xlabel('$\nu_x$ (1/m)', 'interpreter', 'LaTeX')
ylabel('$\nu_y$ (1/m)', 'interpreter', 'LaTeX')
zlabel('$\phi$ (rad)', 'interpreter', 'LaTeX')
view(135, 45)
shading interp
colorbar

%% Plot (Phi, z) curves for different r's.

figure
hold on
for i = 1:length(r_x)
    plot(z, Phi(i, :))
end
plot(z, k * z)
hold off
xlabel('Propagation distance $z$ (normalized by $\lambda$)', 'interpreter', 'LaTeX')
ylabel('$\Phi$ (rad)', 'interpreter', 'LaTeX')
legend({'$r = 0.1 \lambda$', '$r = 0.3 \lambda$', ...
    '$r = 0.5 \lambda$', '$r = 1.0 \lambda$', 'plane wave'}, 'interpreter', 'LaTeX')

%% Plot (kz - Phi, z) curves for different r's.
figure
hold on
for i = 1:length(r_x)
    plot(z, dip(i, :))
end
hold off
xlabel('Propagation distance $z$ (normalized by lambda)', 'interpreter', 'LaTeX')
ylabel('$kz - \Phi$ (rad)', 'interpreter', 'LaTeX')
legend({'$r = 0.1 \lambda$', '$r = 0.3 \lambda$', ...
    '$r = 0.5 \lambda$', '$r = 1.0 \lambda$'}, 'interpreter', 'LaTeX')

%% Plot evolution of total phase
figure
for i = 1:length(log_z)
    subplot(3, 1, i)
    str = strcat('$z = ', num2str(log_z(i)), ' \lambda$');
    title(str, 'interpreter', 'LaTeX')
    pcolor(x, x, angle(log_field_z(:, :, i)))
    xlabel('$x$ (m)', 'interpreter', 'LaTeX')
    ylabel('$y$ (m)', 'interpreter', 'LaTeX')
    zlabel('$\phi + \Phi$ (rad)', 'interpreter', 'LaTeX')
    shading interp
    colorbar
end

%% Plot evplution of retarded phase
figure
for i = 1:length(log_z)
    subplot(3, 1, i)
    str = strcat('$z = ', num2str(log_z(i)), ' \lambda$');
    title(str, 'interpreter', 'LaTeX')
    surf(x, x, angle(log_field_z(:, :, i)) - angle(field))
    xlabel('$x$ (m)', 'interpreter', 'LaTeX')
    ylabel('$y$ (m)', 'interpreter', 'LaTeX')
    zlabel('$\Phi$ (rad)', 'interpreter', 'LaTeX')
    shading interp
    colorbar
end


