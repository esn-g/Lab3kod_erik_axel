format longE

% Initial parameters and tolerances
T = 0.0001;
x_target = 7.5;
y_target = 15;

F = 0.9; 
F_last = 1;
angle_pre = deg2rad(80);
angle = angle_pre;
angle_last = 1.7;
E_x_last = 1;
E_y_last = 1;

for i = 1:1
    % Simulation parameters
    t_tot = 12;
    h = 0.00001;
    N = ceil(t_tot / h);

    % Initial velocities
    e_x_0 = cos(angle) * 20;
    e_y_0 = sin(angle) * 20;

    % Runge-Kutta integration
    [X, Y] = RungeKutta(@eDeriv, e_x_0, e_y_0, h, N, F);
    [pos_X, pos_Y] = Integrate(X, Y, 0, 0, h, N);
    
    % Detect landing position
    idx = find(pos_Y(1:end-1) .* pos_Y(2:end) < 0, 1);
    pos_X_landing = pos_X(idx);
    pos_Y_max = max(pos_Y);
    
    % Calculate errors
    E_x = pos_X_landing - x_target;
    E_y = pos_Y_max - y_target;
    E = [E_x; E_y];
    
    % Jacobian matrix
    J = [dG1dF, dG1dphi; dG2dF, dG2dphi];
    
    % Correction term
    t = J \ (-E);
    
    % Update parameters
    F_last = F;
    E_x_last = E_x;
    E_y_last = E_y;
    angle_last = angle;

    F = F + t(1);
    angle = angle + t(2);
    
    % Print debug information
    disp(['Iteration: ', num2str(i)]);
    disp(['F: ', num2str(F)]);
    disp(['angle: ', num2str(rad2deg(angle))]);
    disp(['E_x: ', num2str(E_x)]);
    disp(['E_y: ', num2str(E_y)]);
    disp(['pos_X_landing: ', num2str(pos_X_landing)]);
    disp(['pos_Y_max: ', num2str(pos_Y_max)]);
    disp('----------------------------');
    
    % Check for convergence
    if norm(E) < T
        disp('Converged');
        disp(['pos_X_landing: ', num2str(pos_X_landing)]);
        disp(['pos_Y_max: ', num2str(pos_Y_max)]);
        disp(['F: ', num2str(F)]);
        disp(['angle: ', num2str(rad2deg(angle))]);
        break
    end
end

% Runge-Kutta method
function [X, Y] = RungeKutta(f, vx0, vy0, h, N, F)
    X = zeros(1, N+1);
    Y = zeros(1, N+1);
    X(1) = vx0;
    Y(1) = vy0;
    t = 0;
    for i = 1:N
        [K1x, K1y] = f(X(i), Y(i), t, F);
        [K2x, K2y] = f(X(i) + h/2 * K1x, Y(i) + h/2 * K1y, t + h/2, F);
        [K3x, K3y] = f(X(i) + h/2 * K2x, Y(i) + h/2 * K2y, t + h/2, F);
        [K4x, K4y] = f(X(i) + h * K3x, Y(i) + h * K3y, t + h, F);
        X(i+1) = X(i) + (h/6) * (K1x + 2*K2x + 2*K3x + K4x);
        Y(i+1) = Y(i) + (h/6) * (K1y + 2*K2y + 2*K3y + K4y);
        t = t + h;
    end
end

% Integrate positions from velocities
function [posX, posY] = Integrate(vx, vy, pos_x0, pos_y0, h, N)
    posX = zeros(1, N+1);
    posY = zeros(1, N+1);
    posX(1) = pos_x0;
    posY(1) = pos_y0;
    for i = 1:N
        posX(i+1) = posX(i) + h * vx(i);
        posY(i+1) = posY(i) + h * vy(i);
    end
end

% Derivative function
function [e_x_prim, e_y_prim] = eDeriv(x_i, y_i, t, F)
    k_x = 0.001;
    k_y = 0.001;
    g = 9.82;
    m_0 = 0.05;
    k = 0.08;
    V = sqrt(x_i^2 + y_i^2);
    phi = atan2(y_i, x_i);
    if t <= 0.08
        m = m_0 - (k * t);
    else
        F = 0;
        m = m_0 - (k * 0.08);
    end
    e_x_prim = (F * cos(phi) - k_x * x_i * V) / m;
    e_y_prim = ((F * sin(phi) - k_y * y_i * V) / m) - g;
end