function [indQRS,sogliaINT,sogliaFIR] = rilevamentoComplessiQRS(derINTecg,INTecg,FIRecg,fs,campioneIniziale,campioneFinale)
% Rilevamento dei complessi QRS

% ----------------------------- PARAMETRI ---------------------------------
% I parametri utilizzati sono tutti citati nella relazione.

TEMPO_DI_REFRAZIONE = 0.25; % 250 ms
% Intervallo di tempo in cui vengono ignorati i successivi picchi rilevati.

PERCENTUALE_SOGLIA_INIZIALE_SEGNALE = 0.3; % 30%
PERCENTUALE_SOGLIA_INIZIALE_RUMORE = 0.5; % 50%
% Percentuali utilizzate per l'inizializzazione della soglia

% Inizializzazine soglia del segnale integrato:
livelloSegnaleINT = PERCENTUALE_SOGLIA_INIZIALE_SEGNALE*max(INTecg(campioneIniziale:campioneFinale));
livelloRumoreINT = PERCENTUALE_SOGLIA_INIZIALE_RUMORE*mean(INTecg(campioneIniziale:campioneFinale));
sogliaINT = livelloRumoreINT + 0.25*(livelloSegnaleINT - livelloRumoreINT);

% Inizializzazione soglia del segnale filtrato dagli artefatti:
livelloSegnaleFIR = PERCENTUALE_SOGLIA_INIZIALE_SEGNALE*max(FIRecg(campioneIniziale:campioneFinale));
livelloRumoreFIR = PERCENTUALE_SOGLIA_INIZIALE_RUMORE*mean(FIRecg(campioneIniziale:campioneFinale));
sogliaFIR = livelloRumoreFIR + 0.25*(livelloSegnaleFIR - livelloRumoreFIR);

% Dichiarazione vettore dei complessi QRS (picchi R) rilevati:
picchi = 1;
% Questo picco verrà rimosso alla fine della funzione poiché non
% rappresenta un complesso QRS.
% L'inizializzazione a 1 causa inoltre un problema del controllo del TEMPO 
% DI REFRAZIONE. Quando controllo il primo picco con la condizione: 
% ((i - picchi(end))/fs) > TEMPO_DI_REFRAZIONE
% picchi(end) vale 1 e di conseguenza ignora il primo picco reale rilevato
% se in un intervallo di tempo minore o uguale al TEMPO_DI_RIFRAZIONE.

% Calcolo del numero di campioni
numeroCampioni = length(FIRecg);

% Ricerca dei picchi:
for i = 2 : numeroCampioni - 1 
   % Controllo sia un massimo locale nella derivata (punto max pendenza)
   if derINTecg(i) > derINTecg(i-1) && derINTecg(i) > derINTecg(i+1)
       % Controllo sia passato il tempo di refrazione
       if ((i - picchi(end))/fs) > TEMPO_DI_REFRAZIONE
           % Controllo sia sopra la soglia INT 
           if INTecg(i) > sogliaINT  
                % Controllo sia sopra la soglia FIR
                if FIRecg(i) > sogliaFIR  
                    % Essendo un complesso QRS, aggiorno i valori dei livelli dei segnali 
                    livelloSegnaleINT = 0.125*INTecg(i) + 0.875*livelloSegnaleINT; 
                    livelloSegnaleFIR = 0.125*FIRecg(i) + 0.875*livelloSegnaleFIR;
                    % Salvo il numero del campione relativo al picco:
                    picchi = [picchi; i];  %#ok<AGROW>
                    % Warning soppresso: dato che il vettore aumenta le 
                    % dimensioni solo in determinate condizioni, la 
                    % preallocazione non è fattibile.
                else
                    % Aggiorno il livello del rumore FIR
                    livelloRumoreFIR = 0.125*FIRecg(i) + 0.875*livelloRumoreFIR;
                end
            else  
                % Aggiorno il livello del rumore INT
                livelloRumoreINT = 0.125*INTecg(i) + 0.875*livelloRumoreINT; 
            end
            % Aggiorno la soglie ogni ciclo sulla base dei valori dei livelli
            sogliaINT = livelloRumoreINT + 0.25*(livelloSegnaleINT - livelloRumoreINT); % aggiorno sogliaINT
            sogliaFIR = livelloRumoreFIR + 0.25*(livelloSegnaleFIR - livelloRumoreFIR); % aggiorno sogliaBP  
       end 
   end
end

% Scartiamo il primo picco
indQRS = picchi(2:end);

end

