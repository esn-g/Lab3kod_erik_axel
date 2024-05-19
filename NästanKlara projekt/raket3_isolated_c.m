function rocket_simulation()
    % Initial conditions
    angle = deg2rad(80);
    vx0 = cos(angle) * 20; % Initial velocity x-component
    vy0 = sin(angle) * 20; % Initial velocity y-component
    pos_x0 = 0; % Initial position x-component
    pos_y0 = 0; % Initial position y-component

    % Time settings and discretization
    t_tot = 5;
    h = 0.0001;
    N = round(t_tot / h); % Number of steps

    % Target landing position
    target_x = 8.5;

    % Newton-Raphson Shooting Method for finding `F`
    F_initial = 0.5; % Initial guess for `F`
    F_tolerance = 0.01; % Tolerance for landing position
    F = F_initial; % Starting value for iteration

    while true
        % Run simulation with current `F`
        [posX, ~, ~, ~] = RK4_positions_second_order(@accel_func, vx0, vy0, pos_x0, pos_y0, h, N, F);

        % Landing error
        final_x = posX(end);
        error = final_x - target_x;

        if abs(error) <= F_tolerance
            break; % Close enough, exit loop
        end

        % Numerical derivative approximation
        delta_F = 0.01; % Small increment for finite difference
        [~, ~, ~, ~] = RK4_positions_second_order(@accel_func, vx0, vy0, pos_x0, pos_y0, h, N, F + delta_F);
        error_delta = posX(end) - target_x;
        derivative = (error_delta - error) / delta_F;

        % Newton-Raphson update for `F`
        F = F - error / derivative;

        fprintf("Trying F=%f, Error=%f\n", F, error);
    end

    fprintf("Final F=%f yields landing x=%f\n", F, final_x);

    % Plot final trajectory
    [posX, posY, velX, velY] = RK4_positions_second_order(@accel_func, vx0, vy0, pos_x0, pos_y0, h, N, F);

    figure;
    subplot(2, 1, 1);
    plot(posX, posY);
    xlabel('X position');
    ylabel('Y Position');
    title('Rocket Trajectory');
    grid on;

    subplot(2, 1, 2);
    plot(linspace(0, t_tot, N+1), sqrt(velX.^2 + velY.^2));
    xlabel('Time (s)');
    ylabel('Speed (m/s)');
    title('Rocket Speed vs Time');
    grid on;
end

function [posX, posY, velX, velY] = RK4_positions_second_order(accel_func, vx0, vy0, pos_x0, pos_y0, h, N, F)
    % Initialize arrays for positions and velocities
    posX = zeros(1, N+1);
    posY = zeros(1, N+1);
    velX = zeros(1, N+1);
    velY = zeros(1, N+1);
    posX(1) = pos_x0;
    posY(1) = pos_y0;
    velX(1) = vx0;
    velY(1) = vy0;

    for i = 1:N
        t = (i - 1) * h;

        % k1: initial accelerations and increments
        [ax1, ay1] = accel_func(velX(i), velY(i), posX(i), posY(i), t, F);
        k1vx = h * ax1;
        k1vy = h * ay1;
        k1x = h * velX(i);
        k1y = h * velY(i);

        % k2: mid-step
        [ax2, ay2] = accel_func(velX(i) + k1vx/2, velY(i) + k1vy/2, posX(i) + k1x/2, posY(i) + k1y/2, t + h/2, F);
        k2vx = h * ax2;
        k2vy = h * ay2;
        k2x = h * (velX(i) + k1vx/2);
        k2y = h * (velY(i) + k1vy/2);

        % k3: another mid-step
        [ax3, ay3] = accel_func(velX(i) + k2vx/2, velY(i) + k2vy/2, posX(i) + k2x/2, posY(i) + k2y/2, t + h/2, F);
        k3vx = h * ax3;
        k3vy = h * ay3;
        k3x = h * (velX(i) + k2vx/2);
        k3y = h * (velY(i) + k2vy/2);

        % k4: end-step
        [ax4, ay4] = accel_func(velX(i) + k3vx, velY(i) + k3vy, posX(i) + k3x, posY(i) + k3y, t + h, F);
        k4vx = h * ax4;
        k4vy = h * ay4;
        k4x = h * (velX(i) + k3vx);
        k4y = h * (velY(i) + k3vy);

        % Update velocities and positions
        velX(i+1) = velX(i) + (1/6) * (k1vx + 2*k2vx + 2*k3vx + k4vx);
        velY(i+1) = velY(i) + (1/6) * (k1vy + 2*k2vy + 2*k3vy + k4vy);
        posX(i+1) = posX(i) + (1/6) * (k1x + 2*k2x + 2*k3x + k4x);
        posY(i+1) = posY(i) + (1/6) * (k1y + 2*k2y + 2*k3y + k4y);
    end
end

function [ax, ay] = accel_func(vx, vy, posX, posY, t, F)
    % Define accelerations based on velocities and positions
    k_x = 0.001;
    k_y = 0.001;
    g = 9.82;

    V = sqrt(vx^2 + vy^2);

    ax = (F * vx / V - k_x * vx * V);
    ay = (F * vy / V - k_y * vy * V) - g;
end
