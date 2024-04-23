

% Konstanter

m0 = 0.050;    
k = 0.08;      
F = 1;         
g = 9.82;      
k_x = 0.001;
k_y = 0.001;


% Omskrivning till ett system av första ordningens diff eq.

% x_prim = v_x;
% y_prim = v_y;

% hjälp funktioner 
V = sqrt(v_x^2 + v_y^2);
phi = atan2(v_y, v_x);

% andra derivatorna
v_x_prim = (F * cos(phi) - k_x * v_x * V) / m;
v_y_prim = (F * sin(phi) - k_y * v_y * V - g) / m;

state_vector = [x_prim, y_prim, v_x_prim, v_y_prim];

if t <= 0.08 % medan bränslet varar
    F_t = 1;
else
    F_t = 0;
end


function [t, y] = rungeKuttaSolver(f, t0, y0, tEnd, h)
    % f - funk tions hand tag för f(t,y)
    % t0 - start tid
    % y0 - initi alt värde av y
    % tEnd - sluttid
    % h - steglängd

    % Minnesallokering
    t = t0:h:tEnd;
    y = zeros(1, length(t));
    y(1) = y0;

    % Runge-Kutta-metoden
    for i = 1:(length(t) - 1)
        % k1 är den initiala lutningen, vilket är lutningen i början av intervallet
        k1 = h * f(t(i), y(i));
        
        % k2 är lutningen vid intervallets mitt, beräknad med k1 för att uppskatta y vid t + h/2
        k2 = h * f(t(i) + 0.5 * h, y(i) + 0.5 * k1);
        
        % k3 är en annan lutningsuppskattning vid mittintervallet, men nu använder vi k2 för att uppskatta y vid t + h/2
        k3 = h * f(t(i) + 0.5 * h, y(i) + 0.5 * k2);
        
        % k4 är lutningen vid intervallets slut, beräknad med k3 för att uppskatta y vid t + h
        k4 = h * f(t(i) + h, y(i) + k3);
        
        % Kombinera de viktade medelvärdena av k1, k2, k3 och k4 för att beräkna nästa värde av y
        y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end

    % Returnera arrayerna av t och y värden
end
