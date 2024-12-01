% Monte Carlo Methode voor het Bepalen van de Distributie van de Cd Waarden
clear
% CONSTANTS
g = 9.81; % N/kg
A = 0.0009073; % m^2, oppervlakte
m = 23 * 10^-3; % kg, massa van de cilinder
rho = 1.23; % kg/m^3
cd = 4.612742; % bepaald in script: KKO

% Verschuiven metingen
metingen = load("metingen.mat").metingen;
metingen(1, :) = metingen(1, :) + 0.0313;

% seed voor random number generator
rng('default')

% Simulatie parameters
dt = 1/500;
val_afstand = 2.4;
max_tijd = metingen(1, end);
aantal_stappen = round(max_tijd/dt);

% Std afwijking bepalen oorsprongkelijke metingen (hiervoor werd
% geoptimiseerd in het script: tijd_verschuiving optimisatie, is nu zo klein mogelijk)
[punten, ~] = vrijeVal_afstand(cd, rho, A, m, g, dt, val_afstand, round(max_tijd / dt));
sim_at_measurements = interp1((0:aantal_stappen) * dt, punten, metingen(1, :), 'linear', 'extrap');
sigma = std(sim_at_measurements - metingen(2, :));

% Monte Carlo simulatie parameters
datapoints = 300;
best_cd = zeros(1, datapoints);
t = metingen(1, end);

% Monte Carlo simulatie
for i = 1:datapoints
    % We simuleren metingen met een std af
    punten = vrijeVal_tijd(cd, rho, A, m, g, dt, t, sigma);
    
    % Create simulated measurement structure
    gesimuleerde_meting = [dt*(1:length(punten)) ; punten];
    
    % Find best-fit drag coefficient for noisy trajectory
    best_cd(i) = lsqnonlin(@(params) afwijkingMetingModel(params(1), gesimuleerde_meting, rho, A, m, g, dt, 0), 0.5);
end

% Statistische analyse
sigma_cd = std(best_cd);
mu_cd = sum(best_cd) / datapoints;

% Uitvoer
fprintf("Std distrib. measurements: %f\n", sigma);
fprintf("Distribution of the Cd with measurement's std distrib.:\n")
fprintf("\tStd distrib.: %f\n", sigma_cd)
fprintf("\tAverage: %f (should equal: %f)\n", mu_cd, cd);

% Plotten in een histogram
figure;
histogram(best_cd);
title('Verdeling wrijvingscoëficienten (Cd)');
xlabel('Wrijvingscoëficient Cd');
ylabel('Frequentie');

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
function [punten_lijst] = vrijeVal_tijd(cd, rho, A, m, g, dt, t, sigma)
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
        punten_lijst(i + 1) = yi + normrnd(0, sigma);
    end
end

% functie die het verchil returnt, tussen de metingen, en overeenkomstige
% punten op een gesimuleerde curve
function [verschillen_lijst] = afwijkingMetingModel(cd, metingen, rho, A, m, g, dt, sigma)
    t = metingen(1, end);
    punten = vrijeVal_tijd(cd, rho, A, m, g, dt, t, sigma);
    sim_at_measurements = interp1((0:length(punten)-1) * dt, punten, metingen(1, :), 'linear', 'extrap');
    verschillen_lijst = (sim_at_measurements - metingen(2, :));
end