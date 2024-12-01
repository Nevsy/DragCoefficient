metingen_echt = [ 0.027 0.033333333 0.066666667 0.1 0.133333333 0.166666667 0.2 0.233333333 0.266666667 0.3 0.333333333 0.366666667 0.4 0.433333333 0.466666667 0.5 0.533333333 0.566666667 0.6 0.633333333 0.666666667 0.7;
             0 0.0042 0.0196 0.049 0.091 0.147 0.217 0.301 0.399 0.5075 0.6265 0.7525 0.8855 1.0255 1.1725 1.3265 1.4875 1.6555 1.8305 2.0125 2.2015 2.3975
];

g = 9.81; % N/kg

A = 0.0009073; % m^2, oppervlakte
m = 23 * 10^-3; % kg, massa van de cilinder

dt = 1/3000;
aantal_stappen = 500 / dt;
rho = 1.23;

val_afstand = 2.4;

function [best_cd, best_shift, lowest_std] = optimize_time_shift(metingen_echt, rho, A, m, g, dt, val_afstand)
    % verschuivingen
    shift_range = linspace(-1/30, 1/30, 200); % Adjust range as needed
    
    % initialisatie variabelen
    std_afw = zeros(size(shift_range));
    cd_values = zeros(size(shift_range));
    
    % Optimisatie
    for i = 1:length(shift_range)
        % Verschuiven metingen over tijd
        metingen_shifted = metingen_echt;
        metingen_shifted(1, :) = metingen_echt(1, :) + shift_range(i);
        
        % de beste cd vinden (de uitvoer doen verdwijnen)
        options = optimoptions('lsqnonlin', 'Display', 'off');
        cd_best = lsqnonlin(@(cd) afwijkingMetingModel(cd, metingen_shifted, rho, A, m, g, dt), 0.5, [], [], options);
        
        % Model toepassen met gevonden cd waarde
        max_tijd = metingen_shifted(1, end);
        [punten, ~] = vrijeVal_afstand(cd_best, rho, A, m, g, dt, val_afstand, round(max_tijd / dt));
        
        % Interpoleren tussen de metingen en het model om correcte punten
        % te vergelijken
        sim_at_measurements = interp1((0:size(punten,2)-1) * dt, punten, metingen_shifted(1, :), 'linear', 'extrap');
        
        % Std afwijking berekenen
        std_afw(i) = std(sim_at_measurements - metingen_shifted(2, :));
        cd_values(i) = cd_best;
    end
    
    % Beste verschuiving bijhouden
    [lowest_std, best_idx] = min(std_afw);
    best_shift = shift_range(best_idx);
    best_cd = cd_values(best_idx);
end

% Functie hierboven toepassen
[best_cd, best_shift, lowest_std] = optimize_time_shift(metingen_echt, rho, A, m, g, dt, val_afstand);

% Beste verschuiving toepassen
metingen_shifted = metingen_echt;
metingen_shifted(1, :) = metingen_echt(1, :) + best_shift;

% Resultaat plotten
[beste_punten, eind_index_beste] = vrijeVal_afstand(best_cd, rho, A, m, g, dt, val_afstand, round(metingen_shifted(1, end) / dt));
plot(1:eind_index_beste , -beste_punten(1:eind_index_beste), 'r');
hold on
plot(metingen_shifted(1, :) / dt, -metingen_shifted(2, :), 'ro');
legend(["best fit", "data"])

% functie die een model voor de vrije val van een voorwerp simuleert, tot
% een bepaalde hoogte
function [punten_lijst, geindigd_bij] = vrijeVal_afstand(cd, rho, A, m, g, dt, floor, max_stappen)
    % positie & snelheid bij t = 0
    vyi = 0;
    yi = 0;

    punten_lijst = zeros(1, max_stappen + 1); % afgelegde afstand
    geindigd_bij = 0;

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
end

% functie die een vrije val van t seconden simuleert
function [punten_lijst] = vrijeVal_tijd(cd, rho, A, m, g, dt, t)
    % positie & snelheid bij t = 0
    vyi = 0;
    yi = 0;

    aantal_stappen = round(t/dt); 
    punten_lijst = zeros(1, aantal_stappen + 1); % afgelegde afstand

    for i = 1:aantal_stappen
        % nieuwe wrijving, resulterende versnelling, snelheid & punt berekenen
        F_w = -1/2 * cd * A * rho * vyi^2;
        ay = F_w / m + g; % F_w negatief tov de zin van versnelling
        vyi = vyi + ay * dt;
        yi = yi + vyi * dt;
        
        % toevoegen van het berekend punt aan de puntenmatrix, met wat ruis
        punten_lijst(i + 1) = yi;
    end
end

% functie die het verchil returnt, tussen de metingen, en overeenkomstige
% punten op een gesimuleerde curve
function [verschillen_lijst] = afwijkingMetingModel(cd, metingen, rho, A, m, g, dt)
    t = metingen(1, end);
    punten = vrijeVal_tijd(cd, rho, A, m, g, dt, t);
    aantal_stappen = round(t/dt);
    % zelfde interpolatie dan bij 'optimize_time_shift' functie
    sim_at_measurements = interp1((0:aantal_stappen) * dt, punten, metingen(1, :), 'linear', 'extrap');
    verschillen_lijst = (sim_at_measurements - metingen(2, :));
end
