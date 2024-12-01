metingen_echt = [ 0.027 0.033333333 0.066666667 0.1 0.133333333 0.166666667 0.2 0.233333333 0.266666667 0.3 0.333333333 0.366666667 0.4 0.433333333 0.466666667 0.5 0.533333333 0.566666667 0.6 0.633333333 0.666666667 0.7;
             0 0.0042 0.0196 0.049 0.091 0.147 0.217 0.301 0.399 0.5075 0.6265 0.7525 0.8855 1.0255 1.1725 1.3265 1.4875 1.6555 1.8305 2.0125 2.2015 2.3975
]; % rij 1: t (s), rij 2: afstand (m)

g = 9.81; % N/kg
A = 0.0009073; % m^2, oppervlakte
m = 23 * 10^-3; % kg, massa van de cilinder
rho = 1.23;

% Parameters simulatie
dt = 1/3000;
aantal_stappen = 7 / dt;

% Echte data berkening cd
bestFit = lsqnonlin(@(params) vrijeValVerschilMetMetingen(params(1), rho, metingen_echt, A, m, g, dt), 0.5);
fprintf("Real data best fit:\n")
fprintf("\tCd = %f\n", bestFit)

% Plotten resultaat (curve)
[punten, eind_index] = vrijeValPunten(bestFit, rho, A, m, g, dt, metingen_echt(2, end), aantal_stappen);
plot((1:eind_index)*dt ,-punten(1:eind_index), 'r');
hold on

% Verschuiven data (.0313 komt van script: Tijd verschuiving optimisatie)
metingen_aangepast = metingen_echt;
metingen_aangepast(1, :) = metingen_aangepast(1, :) + 0.0313;
% Verschoven data berekening cd
bestFit = lsqnonlin(@(params) vrijeValVerschilMetMetingen(params(1), rho, metingen_aangepast, A, m, g, dt), 0.5);
fprintf("Real data + .0313 best fit:\n")
fprintf("\tCd = %f\n", bestFit)

% Plotten resultaat (curve)
[punten, eind_index] = vrijeValPunten(bestFit, rho, A, m, g, dt, metingen_echt(2, end), aantal_stappen);
plot((1:eind_index)*dt, -punten(1:eind_index), 'b');
hold on

% Plotten datapunten, veschoven en initieel
plot(metingen_aangepast(1, :), -metingen_aangepast(2, :), 'bo');
hold on
plot(metingen_echt(1, :), -metingen_echt(2, :), 'ro');
legend(["best fit: verschoven data: cd = 4.6", "best fit: initiële data: cd = -1.6", "verschoven data", "initiële"])
xlabel("Tijd t (s)");
ylabel("Hoogte h (m)")
title('Vergelijking ideale verschuiving metingen, met de gefitte curve')

% Functie gebruikt om de lsqnonlin functie te gebruiken (mooiere versie
% hiervan gebruikt in twee andere scripts)
% Output is een vector met het verschil tussen de datapunten en het
% dichtsbijzijnde berekende punt
function [verschil_lijst] = vrijeValVerschilMetMetingen(cd, rho, metingen, A, m, g, dt)
    % positie & snelheid bij t = 0
    vyi = 0;
    yi = 0;

    aantal_stappen = round(metingen(1, end)/dt);
    verschil_lijst = zeros(1, size(metingen, 2));
    punten_lijst = zeros(1, aantal_stappen + 1); % afgelegde afstand

    for i = 1:aantal_stappen
        % nieuwe wrijving, resulterende versnelling, snelheid & punt berekenen
        F_w = -1/2 * cd * A * rho * vyi^2;
        ay = F_w / m + g; % F_w negatief tov de zin van versnelling
        vyi = vyi + ay * dt;
        yi = yi + vyi * dt;
        
        % toevoegen van het berekend punt aan de puntenmatrix
        punten_lijst(i + 1) = yi;
    end

    % voor elke meting het verschil tussen de verwachte waarde en de
    % meting berekenen en dit in een vector plaatsen
    for k = 1:size(metingen, 2)
        % t = dt * i (i = #iteraties) => i = t/dt om de metingen op
        % dezelfde schaal te plaatsen als onze berekeningen
        % we vinden de tijd waar onze meting zich bevindt terug in de berekende waarden
        %l = min(round(metingen(1, k) / dt) + 1, aantal_stappen + 1);
        sim_at_measurements = interp1((0:aantal_stappen) * dt, punten_lijst, metingen(1, :), 'linear', 'extrap');
        %verschil_lijst(k) = (punten_lijst(l) - metingen(2, k));
        verschil_lijst = (sim_at_measurements - metingen(2, :));
    end
end

% enkel voor het plotten van functies, simulatie van een val met bepaalde
% constanten tot een bepaalde hoogte
% Output is een matrix met de berkende punten (en overeenkomende tijd) en
% de index (~ tijd) dat het voorwerp een bepaalde verticale afstand heeft
% afgelegd
function [punten_lijst, geindigd_bij] = vrijeValPunten(cd, rho, A, m, g, dt, floor, max_stappen)
    % positie & snelheid bij t = 0
    vyi = 0;
    yi = 0;

    punten_lijst = zeros(1, max_stappen + 1); % afgelegde afstand

    for i = 1:max_stappen
        % nieuwe wrijving, resulterende versnelling, snelheid & punt berekenen
        F_w = -1/2 * cd * A * rho * vyi^2;
        ay = F_w / m + g; % F_w negatief tov de zin van versnelling
        vyi = vyi + ay * dt;
        yi = yi + vyi * dt;
        
        % toevoegen van het berekend punt aan de puntenmatrix
        punten_lijst(i + 1) = yi;

        if yi > floor
            geindigd_bij = i;
            break
        end
    end
    % als de code hier komt, zonder geindigd bij te hebben gecreerd, is er een fout met de 'max stappen', zal dus
    % een fout maken (ideaal gezien zouden we hier 'geindigd_bij = -1'
    % zetten en later telkens checken of we out of bounds zijn)
end