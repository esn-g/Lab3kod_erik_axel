format longE

h_values = [0.000001, 0.0000005]; %lägsta med euler
error_num = zeros(2);

for h_index = 1:length(h_values)

    h = h_values(h_index);
    %e_x_0 = cos(angle) * 20;
    %e_y_0 = sin(angle) * 20;
    
    
    % start
    pos_x0 = 0;
    pos_y0 =0;
    
    % inskjut initial, T TOlerans
    T = 1e-13; %lägs
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
        N = ceil(t_tot / h);
    
        % vinkelberäkning beror
        e_x_0 = cos(z(2)) * 20;
        e_y_0 = sin(z(2)) * 20;
       
    
        % Runge kutta
        %[X, Y] = RungeKutta(@eDeriv, e_x_0, e_y_0, h, N, z(1));

        %Euler
        [X, Y] = ForwardEuler(@eDeriv, e_x_0, e_y_0, h, N, z(1));
        
        % poisitoner
        [pos_X, pos_Y] = Integrate(X, Y, 0, 0, h, N);
        
        % Detect indices where sign changes KOLLA ÖVER
        % !!!!!!!!!!!!!!!!!!!!!!!!!
        % idx = find(pos_Y(1:end-1) .* pos_Y(2:end) < 0);
        
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
    
        [pos_Y_max, maxIndex] = max(pos_Y);
        %----------------------------interpolera-----------------------------
        % välj ett fönster att interpolera i 
        half_window_size = 1;
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

        error_num(h_index,:) = z; % Anta att den första intersektionen är landningspunkten
  
        if abs(E) < abs(T)
            %Plotta resultat
            figure;
            subplot(2, 1, 1);
            plot(pos_X, pos_Y);
            xlabel('X position');
            ylabel('Y Position');
            title('Raketbana');
            grid on;

            subplot(2, 1, 2);
            plot(linspace(0, t_tot, N+1), sqrt(X.^2 + Y.^2));
            xlabel('Tid (s)');
            ylabel('Hastighet (m/s)');
            title('Rakhastighet vs Tid');
            grid on;
            disp(['Tolerance aquired'])
            disp(['Iteration: ', num2str(i)]);
            disp(['F: ', num2str(z(1),10)]);
            disp(['angle: ', num2str(rad2deg(z(2)),10)]);
            disp(['E_x: ', num2str(E_x,10)]);
            disp(['E_y: ', num2str(E_y,10)]);
            disp(['pos_X_landing: ', num2str(pos_X_landing,15)]);
            disp(['pos_Y_max: ', num2str(pos_Y_max,15)]);
            disp('---------------------------------------------')
            break
        end
     
    end

end
error = abs(error_num(2,:)- error_num(1,:));
disp(['Numeriskt fel för F: ', num2str(error(1),15)])
disp(['Numeriskt fel för vinkel: ', num2str(error(2),15)])

% Beräkna felordningen för varje metod
p_Y = abs(log2(errors_num(1, 1)) / log2(errors_num(1, 2)));
p_X = abs(log2(errors_num(2, 1)) / log2(errors_num(2, 2)));

% Presentera resultaten
disp(['Felordning för Max Y-position: ', num2str(p_Y)]);
disp(['Felordning för X-intersektion: ', num2str(p_X)]);


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

function [X, Y] = ForwardEuler(f, vx0, vy0, h, N, F)
    X = zeros(1, N+1);
    Y = zeros(1, N+1);
    X(1) = vx0;
    Y(1) = vy0;
    t = 0;
    for i = 1:N
        [Kx, Ky] = f(X(i), Y(i), t, F);
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