function rocket_simulation
    % Initial conditions
    angle = deg2rad(80);
    e_x_0 = cos(angle) * 20;
    e_y_0 = sin(angle) * 20;

    % Time settings and discretization
    t_tot = 5;
    h = 0.0001;
    N = t_tot / h;

    % Calculate velocities using the original Runge-Kutta method for velocities
    [X, Y] = RungeKutta(@eDeriv, e_x_0, e_y_0, h, N);

    % Calculate positions using the new Runge-Kutta method for positions
    [pos_X, pos_Y] = RungeKuttaPosition(X, Y, 0, 0, h, N);

    % Plotting the results
    figure;
    subplot(2, 1, 1);
    plot(pos_X, pos_Y);
    xlabel('X position');
    ylabel('Y Position');
    title('Rocket Trajectory');
    grid on;

    subplot(2, 1, 2);
    plot(linspace(0, t_tot, N+1), sqrt(X.^2 + Y.^2));
    xlabel('Time (s)');
    ylabel('Speed (m/s)');
    title('Rocket Speed vs Time');
    grid on;
end

function [X, Y] = RungeKutta(f, vx0, vy0, h, N)
    X = zeros(1, N+1);
    Y = zeros(1, N+1);
    X(1) = vx0;
    Y(1) = vy0;
    t = 0;
    for i = 1:N
        [K1x, K1y] = f(X(i), Y(i), t);
        [K2x, K2y] = f(X(i) + h/2 * K1x, Y(i) + h/2 * K1y, t + h/2);
        [K3x, K3y] = f(X(i) + h/2 * K2x, Y(i) + h/2 * K2y, t + h/2);
        [K4x, K4y] = f(X(i) + h * K3x, Y(i) + h * K3y, t + h);
        X(i+1) = X(i) + (h/6) * (K1x + 2*K2x + 2*K3x + K4x);
        Y(i+1) = Y(i) + (h/6) * (K1y + 2*K2y + 2*K3y + K4y);
        t = t + h;
    end
end

function [posX, posY] = RungeKuttaPosition(vx, vy, pos_x0, pos_y0, h, N)
    posX = zeros(1, N+1);
    posY = zeros(1, N+1);
    posX(1) = pos_x0; % Initial position x
    posY(1) = pos_y0; % Initial position y
    for i = 1:N
        posX(i+1) = posX(i) + h * vx(i);
        posY(i+1) = posY(i) + h * vy(i);
    end
end

function [e_x_prim, e_y_prim] = eDeriv(x_i, y_i, t)
    k_x = 0.001;
    k_y = 0.001;
    g = 9.82;
    m_0 = 0.05;
    k = 0.08;
    
    V = sqrt(x_i^2 + y_i^2);
    phi = atan2(y_i, x_i);
    
    if t <= 0.08
        m = m_0 - (k * t);
        F = 1;
    else
        F = 0;
        m = m_0 - (k * 0.08);
    end

    e_x_prim = (F * cos(phi) - k_x * x_i * V) / m;
    e_y_prim = ((F * sin(phi) - k_y * y_i * V) / m) - g;
end
