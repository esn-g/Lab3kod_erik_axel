% initialvärden
angle = deg2rad(80);
e_x_0 = cos(angle)*20;
e_y_0 = sin(angle)*20;

% tidsinställningar och diskretisering
t_tot = 10;
h = 0.001;
N = t_tot/h;

[Y, X] = RungeKutta(@eDeriv, e_x_0, e_y_0, h, N);

Y_pos = trapQueen(Y, h);
X_pos = trapQueen(X, h);

max(Y_pos)
max(X_pos)

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
        F = 1;
    else
        F = 0;
    end
    
    % minskning av massa
    m = m_0 - (k*t);
    
    % beräkning av derivator
    e_x_prim = (F * cos(phi) - k_x * x_i * V) / m;
    e_y_prim = (F * sin(phi) - k_y * y_i * V - g) / m;

end

function [X, Y] = RungeKutta(f, x0, y0, h, N)
    X = zeros(1, N+1);
    Y = zeros(1, N+1);
    X(1) = x0;
    Y(1) = y0;
    t = 0;
    for i = 1:N
        K1 = f(X(i), Y(i),t);
        K2 = f(X(i) + h/2, Y(i) + h/2 * K1, t);
        K3 = f(X(i) + h/2, Y(i) + h/2 * K2, t);
        K4 = f(X(i) + h, Y(i) + h * K3, t);
        y_next = Y(i) + h/6 * (K1 + 2*K2 + 2*K3 + K4);
        x_next = X(i) + h;
        X(i+1) = x_next;
        Y(i+1) = y_next;
        t = i*h;
    end
end
