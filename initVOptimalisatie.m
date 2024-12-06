% Metingen
metingen = [0.123, 0.111, 0.129];
gem_metingen = sum(metingen)/length(metingen);

% CONSTANTEN
g = 9.81; % N/kg
rho = 1.23; % kg/m^3 (lucht)

A = 0.0009073; % m^2, oppervlakte
m = 23 * 10^-3; % kg, massa van de cilinder
%C_d = 4.612742; % N/m, wrijvingsconstante, berekend in script: KKO
C_d = 0.7;
%v_max = 3.6; % m / s, maximale snelheid (dit is een placeholder waarde, werd echt berekenend in experiment 1)

theta = pi/4; % radiaal, de hoek waarop de cilinder gelanceerd wordt

% stap-gerelateerde constanten (zie Euler-methode)
dt = 1/5000;
te_bekijken_aftand = 5;
aantal_stappen = te_bekijken_aftand / dt;

v_init = lsqnonlin(@(params) schuineWorpVerschilMetMetingen(params(1), C_d, rho, gem_metingen, A, m, g, theta, dt, aantal_stappen), 5);
fprintf("V_init geoptimaliseerd: %f\n", v_init);
afstand = schuineWorpVerschilMetMetingen(v_init, C_d, rho, gem_metingen, A, m, g, theta, dt, aantal_stappen) + gem_metingen;
fprintf("afstand: %f\n", afstand);

function [verschil] = schuineWorpVerschilMetMetingen(v_max, C_d, rho, gem_metingen, A, m, g, theta, dt, aantal_stappen)
    % snelheid bij t = 0
    v_ix = v_max * cos(theta);
    v_iy = v_max * sin(theta);
    
    % positie bij t = 0
    x_i = 0;
    y_i = 0;

    for i = 1:aantal_stappen
        % nieuwe wrijving, resulterende versnelling, snelheid & punt berekenen
        F_wx = -1/2 * C_d * A * rho * abs(v_ix) * v_ix;
        F_wy = -1/2 * C_d * A * rho * abs(v_iy) * v_iy;
        a_x = F_wx / m;
        a_y = F_wy / m - g;
        v_ix = v_ix + a_x * dt;
        v_iy = v_iy + a_y * dt;
        x_i = x_i + v_ix * dt;
        y_i = y_i + v_iy * dt;
    
        % stoppen zodra de cilinder op de grond is geland, en het verschil
        % tussen de metingen, en het model berekenen (= uitvoer)
        if y_i < 0
            verschil = x_i - gem_metingen;
            break;
        end
    end
end