function dy = raketODE(t, y, k, g, theta, F_t)
    % Parametrar
    % t - aktuell tid
    % y - nuvarande tillstånd i form av [x; y; vx; vy; m]
    % k - motståndskoefficient
    % g - gravitationsacceleration
    % theta - avfyrningsvinkel
    % F_t - skjutkraft

    vx = y(3);
    vy = y(4);
    m = y(5);

    % Skjutkraftens riktning
    Fx = F_t * cos(theta);
    Fy = F_t * sin(theta);

    % Luftmotstånd
    D = k * sqrt(vx^2 + vy^2);
    Dx = D * vx;
    Dy = D * vy;

    % Tillståndsderivat
    dvx = (Fx - Dx) / m;
    dvy = (Fy - Dy) / m - g;
    dm = -k;

    % Om massan blir mindre än noll sätts derivatan av massan till noll
    if m < 0
        dm = 0;
    end

    dy = [vx; vy; dvx; dvy; dm];
end

function [T, Y] = runge_kutta_raket(derivFunc, tspan, y0, dt, k, g, theta, F)
    T = tspan(1):dt:tspan(2);
    Y = zeros(length(y0), length(T));
    Y(:, 1) = y0;

    for i = 1:(length(T) - 1)
        t = T(i);
        y = Y(:, i);

        % Skjutkraft beräknas inuti derivatfunktionen baserat på tid och massförbrukning
        F_t = (y(5) > k*t) * F; % Skjutkraft endast när det finns bränsle

        k1 = derivFunc(t, y, k, g, theta, F_t);
        k2 = derivFunc(t + dt/2, y + k1*(dt/2), k, g, theta, F_t);
        k3 = derivFunc(t + dt/2, y + k2*(dt/2), k, g, theta, F_t);
        k4 = derivFunc(t + dt, y + k3*dt, k, g, theta, F_t);

        Y(:, i+1) = y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end

% Startvillkor
y0 = [0; 0; 20*cos(80*pi/180); 20*sin(80*pi/180); 0.050]; % [x, y, vx, vy, m]
tspan = [0, 10]; % Tidsintervall
dt = 0.01; % Tidssteg
g = 9.82; % m/s^2, gravitation
k = 0.08; % s^-1, motståndskoefficient
theta = 80 * (pi/180); % avfyrningsvinkel i radianer
F = 1; % N, skjutkraft

% Anrop till Runge-Kutta integratorn
[T, Y] = runge_kutta_raket(@raketODE, tspan, y0, dt, k, g, theta, F);

% Plotta resultatet
plot(Y(1,:), Y(2,:));
title('Raketens bana');
xlabel('Horisontellt avstånd (m)');
ylabel('Vertikal höjd (m)');
grid on;