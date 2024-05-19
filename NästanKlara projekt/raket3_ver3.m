function rocket_simulation()
    % Initial conditions
    angle = deg2rad(80);
    e_x_0 = cos(angle) * 20;
    e_y_0 = sin(angle) * 20;
    pos_x0 = 0; % Initial position x-component
    pos_y0 = 0; % Initial position y-component

    % Time settings and discretization
    t_tot = 5;
    h = 0.0001;
    N = round(t_tot / h); % Total steps

    % Target landing position
    target_x = 8.5;

    % Newton-Raphson Shooting Method for finding `F`
    F_initial = 0.0; % Initial guess for `F`
    F_tolerance = 0.01; % Tolerance for landing position
    F = F_initial; % Starting value for iteration

    while true
        % Run simulation with current `F`
        [posX, posY, vx, vy] = simulate_trajectory(F, e_x_0, e_y_0, pos_x0, pos_y0, h, N);

        % Landing error
        final_x = posX(end);
        error = final_x - target_x;

        if abs(error) <= F_tolerance
            break; % Close enough, exit loop
        end

        % Numerical derivative approximation
        delta_F = 0.01; % Small increment for finite difference
        [posX_delta, ~, ~, ~] = simulate_trajectory(F + delta_F, e_x_0, e_y_0, pos_x0, pos_y0, h, N);
        error_delta = posX_delta(end) - target_x;
        derivative = (error_delta - error) / delta_F;

        % Newton-Raphson update for `F`
        F = F - error / derivative;

        fprintf("Trying F=%f, Error=%f\n", F, error);
    end

    fprintf("Final F=%f yields landing x=%f\n", F, final_x);

    % Plot final trajectory
    figure;
    subplot(2, 1, 1);
    plot(posX, posY);
    xlabel('X position');
    ylabel('Y Position');
    title('Rocket Trajectory');
    grid on;

    subplot(2, 1, 2);
    plot(linspace(0, t_tot, N+1), sqrt(vx.^2 + vy.^2));
    xlabel('Time (s)');
    ylabel('Speed (m/s)');
    title('Rocket Speed vs Time');
    grid on;
end

function [posX, posY, vx, vy] = simulate_trajectory(F, vx0, vy0, pos_x0, pos_y0, h, N)
    % Initialize arrays for positions and velocities
    posX = zeros(1, N+1);
    posY = zeros(1, N+1);
    vx = zeros(1, N+1);
    vy = zeros(1, N+1);
    posX(1) = pos_x0;
    posY(1) = pos_y0;
    vx(1) = vx0;
    vy(1) = vy0;

    % Run RK4 integration
    t = 0;
    for i = 1:N
        [K1vx, K1vy] = accel_func(vx(i), vy(i), F, t);
        [K2vx, K2vy] = accel_func(vx(i) + h/2 * K1vx, vy(i) + h/2 * K1vy, F, t + h/2);
        [K3vx, K3vy] = accel_func(vx(i) + h/2 * K2vx, vy(i) + h/2 * K2vy, F, t + h/2);
        [K4vx, K4vy] = accel_func(vx(i) + h * K3vx, vy(i) + h * K3vy, F, t + h);

        vx(i+1) = vx(i) + (h/6) * (K1vx + 2*K2vx + 2*K3vx + K4vx);
        vy(i+1) = vy(i) + (h/6) * (K1vy + 2*K2vy + 2*K3vy + K4vy);

        posX(i+1) = posX(i) + h * vx(i);
        posY(i+1) = posY(i) + h * vy(i);

        t = t + h;
    end
end

function [ax, ay] = accel_func(vx, vy, F, t)
    k_x = 0.001;
    k_y = 0.001;
    g = 9.82;

    V = sqrt(vx^2 + vy^2);

    if t <= 0.08
        ax = (F * cos(atan2(vy, vx)) - k_x * vx * V);
        ay = ((F * sin(atan2(vy, vx)) - k_y * vy * V) - g);
    else
        ax = -k_x * vx * V;
        ay = -k_y * vy * V - g;
    end
end
