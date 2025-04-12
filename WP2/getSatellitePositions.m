function satellitePositions = getSatellitePositions(sc, satArray, time)
    % Funzione per calcolare la posizione di ciascun satellite in tempo reale
    % Inputs:
    %   sc - satelliteScenario oggetto
    %   satArray - array di satelliti creati nel satelliteScenario
    %   time - tempo di interesse in formato datetime
    % Outputs:
    %   satellitePositions - matrice con le posizioni (latitudine, longitudine, altitudine) di ciascun satellite
    
    % Impostiamo il fuso orario di 'time' per farlo corrispondere a quello di 'sc.StartTime'
    time.TimeZone = sc.StartTime.TimeZone;  % Imposta il fuso orario di 'time' per farlo corrispondere a 'sc.StartTime'
    
    % Calcoliamo il tempo trascorso rispetto a StartTime in secondi
    elapsedTime = seconds(time - sc.StartTime);
    
    satellitePositions = zeros(length(satArray), 3);  % Per memorizzare lat, lon, alt
    
    % Calcoliamo la posizione per ciascun satellite
    for idx = 1:length(satArray)
        sat = satArray(idx);
        
        % Otteniamo la posizione del satellite al tempo specificato utilizzando satellitePosition
        [posECEF, ~] = satellitePosition(sc, sat, elapsedTime);
        
        % Converte la posizione ECEF in latitudine, longitudine e altitudine
        [lat, lon, alt] = ecef2geodetic(posECEF(1), posECEF(2), posECEF(3));
        
        % Salviamo la posizione
        satellitePositions(idx, :) = [lat, lon, alt];
    end
end


