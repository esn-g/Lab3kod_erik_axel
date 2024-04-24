% initialvärden
angle = deg2rad(80);
e_x_0 = cos(angle)*20;
e_y_0 = sin(angle)*20;

% tidsinställningar och diskretisering
t_tot = 5;
h = 0.0001;
N = t_tot/h;

[X, Y, pos_X, pos_Y] = RungeKutta(@eDeriv, e_x_0, e_y_0, h, N);

Y_pos = trapQueen(Y, h);
X_pos = trapQueen(X, h);

max(Y_pos)
max(X_pos)

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


function trap_out = trapQueen(vector, h)
    trap_out = zeros(size(vector)); % allokera 
    cumulative_sum = 0; % Initialize cumulative sum
    for i = 1:length(vector)
        cumulative_sum = cumulative_sum + vector(i); % Add current value to cumulative sum
        trap_out(i) = cumulative_sum * h; % Multiply cumulative sum by step size and assign to integrated vector
    end
end

function [e_x_prim, e_y_prim] = eDeriv(x_i, y_i, t)
    k_x = 0.001;    % ?? ?
    k_y = 0.001;    % konstant
    g = 9.82;       % graav
    m_0 = 0.05;     % startmassa
    k = 0.08;       % någon viktkonstant???
    
    % hastighetsmagnitud
    V = sqrt((x_i^2) + (y_i^2));
    % vinkel
    phi = atan2(y_i,x_i);
    
    % medan bränslet varar, KRAFT
    if t <= 0.08 
        m = m_0 - (k*t);
        F = 1;
    else
        F = 0;
        m = m_0-k*0.08;
    end

    
    % beräkning av derivator
    e_x_prim = (F * cos(phi) - k_x * x_i * V) / m;
    e_y_prim = (((F * sin(phi) - k_y * y_i * V) / m) -g);

end

function [X, Y, posX, posY] = RungeKutta(f, vx0, vy0, h, N)
    X = zeros(1, N+1);
    Y = zeros(1, N+1);
    posX = zeros(1, N+1);  % Position in x
    posY = zeros(1, N+1);  % Position in y
    X(1) = vx0;
    Y(1) = vy0;
    posX(1) = 0; % Initial position x
    posY(1) = 0; % Initial position y
    t = 0;
    for i = 1:N
        [K1x, K1y] = f(X(i), Y(i), t);
        [K2x, K2y] = f(X(i) + h/2 * K1x, Y(i) + h/2 * K1y, t + h/2);
        [K3x, K3y] = f(X(i) + h/2 * K2x, Y(i) + h/2 * K2y, t + h/2);
        [K4x, K4y] = f(X(i) + h * K3x, Y(i) + h * K3y, t + h);
        X(i+1) = X(i) + (h/6) * (K1x + 2*K2x + 2*K3x + K4x);
        Y(i+1) = Y(i) + (h/6) * (K1y + 2*K2y + 2*K3y + K4y);
        posX(i+1) = posX(i) + h * X(i);  % Integrate velocity to get position
        posY(i+1) = posY(i) + h * Y(i);  % Integrate velocity to get position
        t = t + h;
    end
end

