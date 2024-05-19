format longE

angle = deg2rad(80);
v_0 = 20;


params = {'none','angle', 'g', 'k_x', 'k_y', 'm_0', 'bryttid', 'v_0', 'k'};
base_values = [angle, 9.82, 0.001, 0.001, 0.05, 0.08, v_0, 1, 0.08];
perbutation = 0.01; % 1% perturbation

for param_index = 1:length(params)+1
    perturbed_values = base_values;
    perturbed_values(param_index) = base_values(param_index) * (1 + perbutation);
    angle = perturbed_values(2);
    g = perturbed_values(3);
    k_x = perturbed_values(4);
    k_y = perturbed_values(5);
    m_0 = perturbed_values(6);
    bryttid = perturbed_values(7);
    v_0 = perturbed_values(8);
    k = perturbed_values(9);
    

    e_x_0 = cos(angle) * v_0;
    e_y_0 = sin(angle) * v_0;
    % start
    pos_x0 = 0;
    pos_y0 =0;
    
    % inskjut initial, T TOlerans
    T = 0.000001;  
    x_target = 8.5;
    
    % Gissad F
    F = 0.9; 
    E_last = 2;
    F_last = 1;
    
    for i = 1:100 
        % Tid and diskert
        t_tot = 5;
        h = 0.00001;
        N = ceil(t_tot / h);
       
    
        % Runge kutta
        [X, Y] = RungeKutta(@(x_i, y_i, t) eDeriv(x_i, y_i, t, g, k_x, k_y, m_0, bryttid, F, k), e_x_0, e_y_0, h, N, F);
        
        % poisitoner
        [pos_X, pos_Y] = Integrate(X, Y, 0, 0, h, N);
        
        % Detect indices where sign changes KOLLA ÖVER 
        for j = 1:length(pos_Y)-1
            if pos_Y(j) * pos_Y(j+1) < 0
                idx = j; % Append the index to idx
            end
        end
    
        x0 = pos_X(idx);
        y0 = pos_Y(idx);
        x1 = pos_X(idx+1);
        y1 = pos_Y(idx+1);
    
        t_factor = -y0 / (y1 - y0); % Linjär interpolation
        pos_X_landing = x0 + (x1 - x0) * t_factor;
        
        % E blir target funktionen f(x) = pos(F) - x_target = 0 är problemet
        % som
        E = pos_X_landing-x_target;
        
        dFdx = (E_last - E)/(F_last - F);
        
        F_last = F;
    
        F = F - E/dFdx;
        
        E_last = E;
    
        if abs(E) < abs(T)
            disp(['Iteration: ', num2str(i)]);
            disp(['F: ', num2str(F)]);
            disp(['pos_X_landing: ', num2str(pos_X_landing)]);
            disp('----------------------------');
            break
        end
    end
    if param_index == 1
        F_ostord = F;
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

function [e_x_prim, e_y_prim] = eDeriv(x_i, y_i, t, g, k_x, k_y, m_0, bryttid, F_max, k)
    V = sqrt(x_i^2 + y_i^2);
    phi = atan2(y_i, x_i);

    if t <= bryttid
        m = m_0 - (k * t);
        F = F_max;
    else
        F = 0;
        m = m_0 - (k * bryttid);
    end

    e_x_prim = (F * cos(phi) - k_x * x_i * V) / m;
    e_y_prim = ((F * sin(phi) - k_y * y_i * V) / m) - g;
end
