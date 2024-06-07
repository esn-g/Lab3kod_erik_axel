
% Initialparametrar och initiala tillstånd
angle = deg2rad(80);
v_0 = 20;
e_x_0 = cos(angle) * v_0;
e_y_0 = sin(angle) * v_0;

% Startpunkter
pos_x0 = 0;
pos_y0 = 0;

% Total tid och diskretisering
t_tot = 4; %5
h = 0.000001; % lägsta möjliga 
N = round(t_tot / h,0);

% Beräkna hastigheter med RK4
%[X, Y] = RungeKutta(@(x_i, y_i, t) eDeriv(x_i, y_i, t, 9.82, 0.001, 0.001, 0.05, 0.08, 1, 0.08), e_x_0, e_y_0, h, N);

% Euler
[X, Y] = ForwardEuler(@(x_i, y_i, t) eDeriv(x_i, y_i, t, 9.82, 0.001, 0.001, 0.05, 0.08, 1, 0.08), e_x_0, e_y_0, h, N);

% Beräkna positioner med Integrator
[pos_X, pos_Y] = Integrator(X, Y, pos_x0, pos_y0, h, N);

% Plotta resulatat
    figure;
    subplot(2, 1, 1);
    plot(pos_X, pos_Y);
    xlabel('X position');
    ylabel('Y Position');
    title('Raketbana');
    grid on;

    subplot(2, 1, 2);
    plot(linspace(0, t_tot, N+1), sqrt(X.^2 + Y.^2));
    xlabel('Time (s)');
    ylabel('Speed (m/s)');
    title('Rockethastighet vs tid');
    grid on;


%% a) Max Y position
[max_posY, maxIndex] = max(pos_Y);
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
max_posY = y_max_interpolated;
fprintf('a) Max Y interpolerat:\n')
disp(y_max_interpolated)
%-------------------------------------------------------------------


%% b) X-intersektioner när Y = 0
zero_crossings = [];
for i = 1:N
    if pos_Y(i) * pos_Y(i+1) <= 0
        zero_crossings = [zero_crossings, i];
    end
end

x_intersections = [];
for i = 1:length(zero_crossings)
    idx = zero_crossings(i);
    x0 = pos_X(idx);
    y0 = pos_Y(idx);
    x1 = pos_X(idx+1);
    y1 = pos_Y(idx+1);

    t_factor = -y0 / (y1 - y0); % Linjär interpolation
    x_cross = x0 + (x1 - x0) * t_factor;
    x_intersections = [x_intersections, x_cross];
end

fprintf('b) X-intersektioner när Y = 0:\n');
disp(x_intersections);

%% Felanalys
% Steg 1: Numeriska metodens fel
h_values = [0.000001, 0.0000005]; % lägsta möjliga
errors_num = zeros(2, length(h_values));

for h_index = 1:length(h_values)
    h = h_values(h_index);
    N = round(t_tot / h,0);

    % Beräkna hastigheter med RK4
    %[X, Y] = RungeKutta(@(x_i, y_i, t) eDeriv(x_i, y_i, t, 9.82, 0.001, 0.001, 0.05, 0.08, 1, 0.08), e_x_0, e_y_0, h, N);
    
    %Euler
    [X, Y] = ForwardEuler(@(x_i, y_i, t) eDeriv(x_i, y_i, t, 9.82, 0.001, 0.001, 0.05, 0.08, 1, 0.08), e_x_0, e_y_0, h, N);

    % Beräkna positioner med Integrator
    [pos_X, pos_Y] = Integrator(X, Y, pos_x0, pos_y0, h, N);

    % Beräkna max Y-position och X-intersektioner
    %-------------------------------------------------------------------
    [max_posYnum, maxIndex] = max(pos_Y);

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
    max_posYnum = y_max_interpolated;
    %------------------------------------------------------------------

    zero_crossings = [];
    for i = 1:N
        if pos_Y(i) * pos_Y(i+1) <= 0
            zero_crossings = [zero_crossings, i];
        end
    end

    x_intersectionsNum = [];
    for i = 1:length(zero_crossings)
        idx = zero_crossings(i);
        x0 = pos_X(idx);
        y0 = pos_Y(idx);
        x1 = pos_X(idx+1);
        y1 = pos_Y(idx+1);

        t_factor = -y0 / (y1 - y0); % Linjär interpolation
        x_cross = x0 + (x1 - x0) * t_factor;
        x_intersectionsNum = [x_intersectionsNum, x_cross];
    end

    errors_num(1, h_index) = max_posYnum;
    errors_num(2, h_index) = x_intersectionsNum(2); % Anta att den första intersektionen är landningspunkten
end
% Skriv ut num_error_bound
num_error_bound = max(abs(diff(errors_num, 1, 2)), [], 2);
disp('Numeriska felgränser:');
disp(num_error_bound);


% Beräkna felordningen för varje metod
p_Y = abs(log2(errors_num(1, 1)) / log2(errors_num(1, 2)));
p_X = abs(log2(errors_num(2, 1)) / log2(errors_num(2, 2)));

% Presentera resultaten
disp(['Felordning för Max Y-position: ', num2str(p_Y, 20)]);
disp(['Felordning för X-intersektion: ', num2str(p_X, 20)]);


%%
% Steg 2: Osäkerhet i indata
params = {'angle', 'g', 'k_x', 'k_y', 'm_0', 'bryttid', 'v_0', 'F_max', 'k'};
base_values = [angle, 9.82, 0.001, 0.001, 0.05, 0.08, v_0, 1, 0.08];
perbutation = 0.01; % 1% perturbation

input_error_bound = zeros(2, length(params));

for param_index = 1:length(params)
    perturbed_values = base_values;
    perturbed_values(param_index) = base_values(param_index) * (1 + perbutation);
    angle_pert = perturbed_values(1);
    g_pert = perturbed_values(2);
    k_x_pert = perturbed_values(3);
    k_y_pert = perturbed_values(4);
    m_0_pert = perturbed_values(5);
    bryttid_pert = perturbed_values(6);
    v_0_pert = perturbed_values(7);
    F_max_pert = perturbed_values(8);
    k_pert = perturbed_values(9);

    e_x_0_pert = cos(angle_pert) * v_0_pert;
    e_y_0_pert = sin(angle_pert) * v_0_pert;

    % Beräkna hastigheter med RK4
    %[X, Y] = RungeKutta(@(x_i, y_i, t) eDeriv(x_i, y_i, t, g_pert, k_x_pert, k_y_pert, m_0_pert, bryttid_pert, F_max_pert, k_pert), e_x_0_pert, e_y_0_pert, h, N);

    %Euler
    [X, Y] = ForwardEuler(@(x_i, y_i, t) eDeriv(x_i, y_i, t, g_pert, k_x_pert, k_y_pert, m_0_pert, bryttid_pert, F_max_pert, k_pert), e_x_0_pert, e_y_0_pert, h, N);
    % Beräkna positioner med Integrator
    [pos_X, pos_Y] = Integrator(X, Y, pos_x0, pos_y0, h, N);

    % Beräkna max Y-position och X-intersektioner
        %-------------------------------------------------------------------
    [max_posY_perturbed, maxIndex] = max(pos_Y);

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
    max_posY_perturbed = y_max_interpolated;
    %------------------------------------------------------------------

    zero_crossings = [];
    for i = 1:N
        if pos_Y(i) * pos_Y(i+1) <= 0
            zero_crossings = [zero_crossings, i];
        end
    end

    x_intersections_perturbed = [];
    for i = 1:length(zero_crossings)
        idx = zero_crossings(i);
        x0 = pos_X(idx);
        y0 = pos_Y(idx);
        x1 = pos_X(idx+1);
        y1 = pos_Y(idx+1);

        t_factor = -y0 / (y1 - y0); % Linjär interpolation
        x_cross = x0 + (x1 - x0) * t_factor;
        x_intersections_perturbed = [x_intersections_perturbed, x_cross];
    end

    input_error_bound(:, param_index) = abs([max_posY_perturbed - max_posY; x_intersections_perturbed(2) - x_intersections(2)]);
end

total_input_error_bound= sum(input_error_bound,2);
total_input_error_bound_max_y = total_input_error_bound(1);
size(total_input_error_bound_max_y)
total_input_error_bound_intersection = total_input_error_bound(2);
size(total_input_error_bound_intersection)

% Total felgräns

fprintf('Metodens felgränser:\n');
fprintf('max Y:\n')
disp(num_error_bound(1));
fprintf('Landningsposition:\n')
disp(num_error_bound(2));
disp('------------------------')
fprintf('input störningsanalys:\n');
fprintf('Max Y:\n')
disp(total_input_error_bound_max_y);
fprintf('Landningsposition:\n')
disp(total_input_error_bound_intersection);
disp('------------------------')

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

function [X, Y] = ForwardEuler(f, vx0, vy0, h, N)
    X = zeros(1, N+1);
    Y = zeros(1, N+1);
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

function [posX, posY] = Integrator(vx, vy, pos_x0, pos_y0, h, N)
    posX = zeros(1, N+1);
    posY = zeros(1, N+1);
    posX(1) = pos_x0; % Initial position x
    posY(1) = pos_y0; % Initial position y
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
