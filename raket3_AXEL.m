function [X, Y] = RungeKutta(f, x0, y0, h, n)
    X = zeros(1, n+1);
    Y = zeros(1, n+1);
    X(1) = x0;
    Y(1) = y0;
    for i = 1:n
        K1 = f(X(i), Y(i));
        K2 = f(X(i) + h/2, Y(i) + h/2 * K1);
        K3 = f(X(i) + h/2, Y(i) + h/2 * K2);
        K4 = f(X(i) + h, Y(i) + h * K3);
        y_next = Y(i) + h/6 * (K1 + 2*K2 + 2*K3 + K4);
        x_next = X(i) + h;
        X(i+1) = x_next;
        Y(i+1) = y_next;
    end
end

function dydt = f(t, y)
    % Definiera konstanter
    m = 100; % massan av raketten (kg)
    g = 9.81; % tyngdaccelerationen (m/s^2)
    k_x = 0.1; % friktionskoefficient i x-riktning
    k_y = 0.2; % friktionskoefficient i y-riktning
    theta = pi/4; % vinkeln hos raketens rörelse (45 grader)
    F_0 = 1; %start krut
    
    % Extrahera variablerna från vektorn y
    x = y(1); % position i x-riktning
    y = y(2); % position i y-riktning
    x_prime = y(3); % hastighet i x-riktning
    y_prime = y(4); % hastighet i y-riktning
    
    % Beräkna V
    V = sqrt(x_prime^2 + y_prime^2);
    
    % Beräkna krafterna
    if t<0.08 
        F=F_0;
    else
        F=0;
    end

    F_x = F * cos(theta);
    F_y = F * sin(theta);
    
    % Beräkna accelerationerna
    x_double_prime = (F_x - k_x * x_prime * V) / m;
    y_double_prime = (F_y - k_y * y_prime * V - m * g) / m;
    
    % Returnera derivatorna
    dydt = [x_prime; y_prime; x_double_prime; y_double_prime];
end

% Ange initialvärden och parametrar för Runge-Kutta-metoden
angle = deg2rad(80);
x0 = cos(angle)*20;
y0 = sin(angle)*20;
h = 0.001;
n = 1/h;

% Anropa Runge-Kutta-funktionen för att lösa differentialekvationen
[X, Y] = RungeKutta(f(t,y), x0, y0, h, n);

% Skriv ut eller använd X och Y för att analysera resultatet
disp(X);
disp(Y);

