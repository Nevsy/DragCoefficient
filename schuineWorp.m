% CONSTANTEN
g = 9.81; % N/kg
rho = 1.23; % kg/m^3 (lucht)

A = 0.0009073; % m^2, oppervlakte
m = 23 * 10^-3; % kg, massa van de cilinder
C_d = 4.612742; % N/m, wrijvingsconstante, berekend in script: KKO
v_max = 5; % m / s, maximale snelheid (dit is een placeholder waarde, werd echt berekenend in experiment 1)

theta = pi/4; % radiaal, de hoek waarop de cilinder gelanceerd wordt

v_x = v_max * cos(theta);
v_y = v_max * sin(theta);

% Oorspronkelijke wrijvingskracht
F_wx = -1/2 * C_d * A * rho * abs(v_x) * v_x;
F_wy = -1/2 * C_d * A * rho * abs(v_y) * v_y;

% versnelling
a_x = F_wx / m;
a_y = F_wy / m - g;
a_x_noFw = 0;
a_y_noFw = -g;

% stap-gerelateerde constanten (zie Euler-methode)
stapgrootte = .005;
te_bekijken_aftand = 2;
aantal_stappen = te_bekijken_aftand / stapgrootte;

% snelheid bij t = 0
v_ix = v_x;
v_iy = v_y;
v_ix_noFw = v_x;
v_iy_noFw = v_y;

% positie bij t = 0
x_i = 0;
y_i = 0;
x_i_noFw = 0;
y_i_noFw = 0;

% initialisatie van de matrix met de te plotten punten
punten = zeros(2, aantal_stappen);
punten_noFw = zeros(2, aantal_stappen);

for i = 1:aantal_stappen
    % proceduraal plotten (animatie, vertraagt de berekening heel wat)
    plot(punten(1, 1:i), punten(2, 1:i), 'b');
    hold on;
    plot(punten_noFw(1, 1:i), punten_noFw(2, 1:i), 'r');
    legend("Met wrijving", "Zonder wrijving");
    xlabel("Afstand x (m)")
    ylabel("Hoogte h (m)")
    title("Simulatie parabool, met, en zonder wrijving")

    % nieuwe wrijving, resulterende versnelling, snelheid & punt berekenen
    F_wx = -1/2 * C_d * A * rho * abs(v_ix) * v_ix;
    F_wy = -1/2 * C_d * A * rho * abs(v_iy) * v_iy;
    a_x = F_wx / m;
    a_y = F_wy / m - g;
    v_ix = v_ix + a_x * stapgrootte;
    v_iy = v_iy + a_y * stapgrootte;
    x_i = x_i + v_ix * stapgrootte;
    y_i = y_i + v_iy * stapgrootte;

    v_iy_noFw = v_iy_noFw + a_y_noFw * stapgrootte; % snelheid van x veranderd niet zonder wrijving
    x_i_noFw = x_i_noFw + v_ix_noFw * stapgrootte;
    y_i_noFw = y_i_noFw + v_iy_noFw * stapgrootte;
    
    % toevoegen van het berekend punt aan de puntenmatrix
    punten(1:2, i + 1) = [x_i ; y_i];
    punten_noFw(1:2, i + 1) = [x_i_noFw; y_i_noFw];

    % tekenen van het plot
    drawnow;

    % stoppen met het plotten zodra de cilinder op de grond is geland, maar
    % doorgaan met het plotten van de andere curve, nl. dat y negatief wordt
    % (er zou een functie gemaakt moeten worden om herhaling van code te vermijden)
    if y_i < 0
        x_i
        for j = i:aantal_stappen
            plot(punten(1, 1:i), punten(2, 1:i), 'b');
            hold on;
            plot(punten_noFw(1, 1:j), punten_noFw(2, 1:j), 'r');
            legend("Met wrijving", "Zonder wrijving");

            v_iy_noFw = v_iy_noFw + a_y_noFw * stapgrootte; % snelheid van x veranderd niet zonder wrijving
            x_i_noFw = x_i_noFw + v_ix_noFw * stapgrootte;
            y_i_noFw = y_i_noFw + v_iy_noFw * stapgrootte;
            
            punten_noFw(1:2, j + 1) = [x_i_noFw; y_i_noFw];
            drawnow;
            if y_i_noFw < 0
                x_i_noFw
                break;
            end
        end
        break;
    end
end