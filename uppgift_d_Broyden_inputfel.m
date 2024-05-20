format longE

% base values
v_0 = 20;
g = 9.82;
k_x = 0.001;
k_y = 0.001;
m_0 = 0.05;
bryttid = 0.08;
k = 0.08;
params = {'g', 'k_x', 'k_y', 'm_0', 'bryttid', 'v_0', 'k'};
base_values = [g, k_x, k_y, m_0, bryttid, v_0, k];
perturbation = 0.01; % 1% perturbation

% Store the original values
original_F = 0;

F_pert_list = zeros(1,length(params));
input_error = 0;

% Loop over each parameter to perturb and analyze the effect
for param_index = 0:length(params)
    perturbed_values = base_values;
    
    if param_index > 0
        perturbed_values(param_index) = base_values(param_index) * (1 + perturbation);
    end
    
    % Extract perturbed parameters
    g = perturbed_values(1);
    k_x = perturbed_values(2);
    k_y = perturbed_values(3);
    m_0 = perturbed_values(4);
    bryttid = perturbed_values(5);
    v_0 = perturbed_values(6);
    k = perturbed_values(7);

    % start
    pos_x0 = 0;
    pos_y0 =0;

    % inskjut initial, T TOlerans
    T = 1e-13; %lägsta
    x_target = 7.5;
    y_target = 15;
    
    % Gissad F
    F = 0; 
    F_last = 0;
    angle = deg2rad(81); %gissad vinkel
    angle_last = deg2rad(80);
    E_x_last = 3;
    E_y_last = 4;
    E_last = [E_x_last;E_y_last];
    z = [F;angle];
    z_last = [F_last;angle_last];
    A_last  = eye(2,2);
    
    for i = 1:100
        disp(i)
        % Tid and diskert
        t_tot = 4;
        h = 0.000001; %lägsta
        N = ceil(t_tot / h);
    
        % vinkelberäkning beror
        e_x_0 = cos(z(2)) * v_0;
        e_y_0 = sin(z(2)) * v_0;
       
    
        % Runge kutta
        %[X, Y] = RungeKutta(@(x_i, y_i, t) eDeriv(x_i, y_i, t, g, k_x, k_y, m_0, bryttid, z(1), k), e_x_0, e_y_0, h, N);

        %Euler
        [X, Y] = ForwardEuler(@(x_i, y_i, t) eDeriv(x_i, y_i, t, g, k_x, k_y, m_0, bryttid, z(1), k), e_x_0, e_y_0, h, N);
        
        % poisitoner
        [pos_X, pos_Y] = Integrate(X, Y, 0, 0, h, N);
        
        % Detect indices where sign changes
        for j = 1:length(pos_Y)-1
            if pos_Y(j) * pos_Y(j+1) < 0
                idx = j; % Append the index to idx
            end
        end
        
        x0 = pos_X(idx);
        y0 = pos_Y(idx);
        x1 = pos_X(idx + 1);
        y1 = pos_Y(idx + 1);

        t_factor = -y0 / (y1 - y0); % Linear interpolation
        pos_X_landing = x0 + (x1 - x0) * t_factor;
        

        [pos_Y_max, maxIndex] = max(pos_Y);
        %----------------------------interpolera-----------------------------
        % välj ett fönster att interpolera i 
        half_window_size = 20;
        start_index = max(1, maxIndex - half_window_size);
        end_index = min(length(pos_X), maxIndex + half_window_size - 1);
        
        % välj punkter
        x_points = pos_X(start_index:end_index);
        y_points = pos_Y(start_index:end_index);
        
        % kvadratisk ansats
        % y = ax^2 + bx + c
        A = [x_points'.^2, x_points', ones(length(x_points), 1)];
        b = y_points';
        
        % mkm lösning abc
        coeffs = A \ b;
        
        a = coeffs(1);
        b = coeffs(2);
        c = coeffs(3);
        
        % hitta x för max pos analytiskt
        x_vertex = -b / (2 * a);
        
        % beräkna y_max
        y_max_interpolated = a * x_vertex^2 + b * x_vertex + c;
        
        % skriv över
        pos_Y_max = y_max_interpolated;
        %-------------------------------------------------------------------
        
        % E blir target funktionen f(x) = pos(F) - x_target = 0 är problemet
        % skapa f(x) = 0 funktioner för newtons metod
        E_x = pos_X_landing-x_target;
        E_y = pos_Y_max-y_target;
        
        % funktionsvärdesvektor E = G
        E = [E_x, E_y]';
        
        %Broyden beräkningar
        % stora delta
        DELTA = E - E_last;
    
        % lilla delta
        delta = z-z_last;
    
        % Broyden
        A = A_last+(((DELTA-A_last*delta)*(delta'))/((delta')*delta));
        
        
        % korrektionsterm enligt Sauer
        korr = inv(A)*E;
        
        % spara nuvarande värden innan uppdatering
        F_last = F;
        E_x_last = E_x;
        E_y_last = E_y;
        angle_last = angle;
        z_last = [F_last;angle_last];
        E_last = [E_x_last;E_y_last];
        A_last = A;
    
        % uppdatera F och angle för att hitta bättre lösnig
        F = z_last(1) - korr(1);
        angle = z_last(2) - korr(2);
    
        z = [F;angle];
        
        %isp(['Iteration: ', num2str(i)]);
        % disp(['F: ', num2str(z(1))]);
        % disp(['angle: ', num2str(rad2deg(z(2)))]);
        % disp(['E_x: ', num2str(E_x)]);
        % disp(['E_y: ', num2str(E_y)]);
        % disp(['pos_X_landing: ', num2str(pos_X_landing)]);
        % disp(['pos_Y_max: ', num2str(pos_Y_max)]);
        % disp('---------------------------------------------')
    
        if abs(E) < abs(T)
            if param_index == 0
                original_F = z(1);
                original_angle = z(2)
                disp('Original (no perturbation):');
            else
                disp(['Perturbed parameter: ', params{param_index}]);
            end
            disp(['Tolerance acuireq'])
            disp(['Iteration: ', num2str(i)]);
            disp(['F: ', num2str(z(1),10)]);
            disp(['angle: ', num2str(rad2deg(z(2)),10)]);
            disp(['E_x: ', num2str(E_x,10)]);
            disp(['E_y: ', num2str(E_y,10)]);
            disp(['pos_X_landing: ', num2str(pos_X_landing,15)]);
            disp(['pos_Y_max: ', num2str(pos_Y_max,15)]);
            error_F = abs(sum(z(1) - original_F));
            error_angle = abs(sum(z(2) - original_angle));
            disp(['Inputfel F: ', num2str(error_F,15)])
            disp(['Inputfel vinkel: ', num2str(error_angle,15)])
            disp('---------------------------------------------')
            break
        end
    end

    input_error_F = input_error + abs(sum(z(1) - original_F));
    input_error_ang = input_error + abs(sum(z(2) - original_angle));
end
fprintf('input störningsanalys:\n');
disp(['F: ', num2str(input_error_F,15)])
disp(['Vinkel: ', num2str(input_error_ang,15)])

function [X, Y] = RungeKutta(f, vx0, vy0, h, N)
    X = zeros(1, N + 1);
    Y = zeros(1, N + 1);
    X(1) = vx0;
    Y(1) = vy0;
    t = 0;
    for i = 1:N
        [K1x, K1y] = f(X(i), Y(i), t);
        [K2x, K2y] = f(X(i) + h / 2 * K1x, Y(i) + h / 2 * K1y, t + h / 2);
        [K3x, K3y] = f(X(i) + h / 2 * K2x, Y(i) + h / 2 * K2y, t + h / 2);
        [K4x, K4y] = f(X(i) + h * K3x, Y(i) + h * K3y, t + h);
        X(i + 1) = X(i) + (h / 6) * (K1x + 2 * K2x + 2 * K3x + K4x);
        Y(i + 1) = Y(i) + (h / 6) * (K1y + 2 * K2y + 2 * K3y + K4y);
        t = t + h;
    end
end

% 
% function [X, Y] = RungeKutta(f, vx0, vy0, h, N, F)
%     X = zeros(1, N+1);
%     Y = zeros(1, N+1);
%     X(1) = vx0;
%     Y(1) = vy0;
%     t = 0;
%     for i = 1:N
%         [K1x, K1y] = f(X(i), Y(i), t, F);
%         [K2x, K2y] = f(X(i) + h/2 * K1x, Y(i) + h/2 * K1y, t + h/2, F);
%         [K3x, K3y] = f(X(i) + h/2 * K2x, Y(i) + h/2 * K2y, t + h/2, F);
%         [K4x, K4y] = f(X(i) + h * K3x, Y(i) + h * K3y, t + h, F);
%         X(i+1) = X(i) + (h/6) * (K1x + 2*K2x + 2*K3x + K4x);
%         Y(i+1) = Y(i) + (h/6) * (K1y + 2*K2y + 2*K3y + K4y);
%         t = t + h;
%     end
% end

function [X, Y] = ForwardEuler(f, vx0, vy0, h, N)
    X = zeros(1, N + 1);
    Y = zeros(1, N + 1);
    X(1) = vx0;
    Y(1) = vy0;
    t = 0;
    for i = 1:N
        [Kx, Ky] = f(X(i), Y(i), t);
        X(i+1) = X(i) + h * Kx;
        Y(i+1) = Y(i) + h * Ky;
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

function [e_x_prim, e_y_prim] = eDeriv(x_i, y_i, t, g, k_x, k_y, m_0, bryttid, F, k)
    V = sqrt(x_i^2 + y_i^2);
    phi = atan2(y_i, x_i);

    if t <= bryttid
        m = m_0 - (k * t);
    else
        F = 0;
        m = m_0 - (k * bryttid);
    end

    e_x_prim = (F * cos(phi) - k_x * x_i * V) / m;
    e_y_prim = ((F * sin(phi) - k_y * y_i * V) / m) - g;
end


% function [e_x_prim, e_y_prim] = eDeriv(x_i, y_i, t, F)
%     % F = F;
%     k_x = 0.001;
%     k_y = 0.001;
%     g = 9.82;
%     m_0 = 0.05;
%     k = 0.08;
% 
%     % size(x_i)
%     V = sqrt(x_i^2 + y_i^2);
%     phi = atan2(y_i, x_i);
%     if t <= 0.08
%         m = m_0 - (k * t);
%         F = F;
%     else
%         F = 0;
%         m = m_0 - (k * 0.08);
%     end
% 
%     e_x_prim = (F * cos(phi) - k_x * x_i * V) / m;
%     e_y_prim = ((F * sin(phi) - k_y * y_i * V) / m) - g;
% end