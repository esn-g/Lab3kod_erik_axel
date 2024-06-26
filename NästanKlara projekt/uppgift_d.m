format longE
%e_x_0 = cos(angle) * 20;
%e_y_0 = sin(angle) * 20;


% start
pos_x0 = 0;
pos_y0 =0;

% inskjut initial, T TOlerans
T = 0.01;  
x_target = 7.5;
y_target = 15;

% Gissad F
F = 0; 
F_last = 0;
angle_pre = deg2rad(79);
angle = angle_pre; %gissad vinkel
angle_last = deg2rad(89);
E_x_last = 3;
E_y_last = 4;

for i = 1:1
    % Tid and diskert
    t_tot = 18;
    h = 0.00001;
    N = ceil(t_tot / h);

    % vinkelberäkning beror
    e_x_0 = cos(angle) * 20;
    e_y_0 = sin(angle) * 20;
   

    % Runge kutta
    [X, Y] = RungeKutta(@eDeriv, e_x_0, e_y_0, h, N, F);
    
    % poisitoner
    [pos_X, pos_Y] = Integrate(X, Y, 0, 0, h, N);
    
    % Detect indices where sign changes KOLLA ÖVER
    idx = find(pos_Y(1:end-1) .* pos_Y(2:end) < 0);

    pos_X_landing = pos_X(idx)
    pos_Y_max = max(pos_Y)
    
    % E blir target funktionen f(x) = pos(F) - x_target = 0 är problemet
    % skapa f(x) = 0 funktioner för newtons metod
    E_x = pos_X_landing-x_target;
    E_y = pos_Y_max-y_target;
    
    % funktionsvärdesvektor
    E = [E_x, E_y]';
    
    % skapa jacobiankomponenter
    dG1dF = (E_x_last - E_x)/(F_last - F);
    dG1dphi = (E_x_last - E_x)/(angle_last - angle);
    dG2dF = (E_y_last - E_y)/(F_last - F);
    dG2dphi = (E_y_last - E_y)/(angle_last - angle);
    
    % jakobianen
    J = [dG1dF,dG1dphi; dG2dF,dG2dphi];
    
    cond(J)
    % korrektionsterm enligt Sauer
    korr = J\E;
    
    % spara nuvarande värden innan uppdatering
    F_last = F;
    E_x_last = E_x;
    E_y_last = E_y;
    angle_last = angle;

    % uppdatera F och angle för att hitta bättre lösnig
    F = F - korr(1)
    angle = angle - korr(2)
    
    disp(['Iteration: ', num2str(i)]);
    disp(['F: ', num2str(F)]);
    disp(['angle: ', num2str(rad2deg(angle))]);
    disp(['E_x: ', num2str(E_x)]);
    disp(['E_y: ', num2str(E_y)]);
    disp(['pos_X_landing: ', num2str(pos_X_landing)]);
    disp(['pos_Y_max: ', num2str(pos_Y_max)]);

    if abs(E) < abs(T)
        disp(pos_X_landing)
        disp(pos_Y_max)
        disp(F)
        disp(angle)
        break
    end
end


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

function [posX, posY] = Integrate(vx, vy, pos_x0, pos_y0, h, N)
    posX = zeros(1, N+1);
    posY = zeros(1, N+1);
    posX(1) = pos_x0; % start position x
    posY(1) = pos_y0; % start position y
    for i = 1:N
        posX(i+1) = posX(i) + h * vx(i);
        posY(i+1) = posY(i) + h * vy(i);
    end
end

function [e_x_prim, e_y_prim] = eDeriv(x_i, y_i, t, F)
    % F = F;
    k_x = 0.001;
    k_y = 0.001;
    g = 9.82;
    m_0 = 0.05;
    k = 0.08;

    % size(x_i)
    V = sqrt(x_i^2 + y_i^2);
    phi = atan2(y_i, x_i);
    if t <= 0.08
        m = m_0 - (k * t);
        F = F;
    else
        F = 0;
        m = m_0 - (k * 0.08);
    end

    e_x_prim = (F * cos(phi) - k_x * x_i * V) / m;
    e_y_prim = ((F * sin(phi) - k_y * y_i * V) / m) - g;
end