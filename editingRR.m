function [intervalliRR_interpolati,t_RR] = editingRR(intervalliRR,istantiRR,flagStampa)
% Correzione dei valori anomali dovuti potenzialmente a errori di misura: 
% Queste modifiche possono compromettere la correttezza delle informazioni,
% garantiscono risultati migliori nell'utilizzo di metodi statistici e lo
% studio in frequenza.

% Sarebbe sempre meglio eseguire un editing manuale del verrore degli
% intervalli RR.

% ----------------------------- PARAMETRI ---------------------------------
% I parametri utilizzati sono tutti citati nella relazione.

SOGLIA_INTERVALLO_CONFIDENZA = 0.5; % 50%
% La perdita di qualche battito porta a valori anomali di intervalli RR.
% Rimuovendo gli intervalli fuori da un determinato intervallo di confidenza
% ci si avvicina ad una distribuzione normale.

FREQUENZA_DI_CAMPIONAMENTO = 4; % Hz
% Per un corretto lo studio in frequenza di un segnale a tempo discreto 
% con passo di campionamento variabile non è adeguato. 
% Scelta di interpolare il segnale ottenendo un passo di campionamento 
% di fs = 4 Hz (Riportata della relazione)

% ------------------------- Editing del segnale ---------------------------

for k = 2:length(intervalliRR)-1
    % Se un intervallo RR si trova fuori dall'intervallo di confidenza
    if intervalliRR(k) < (1 - SOGLIA_INTERVALLO_CONFIDENZA)*mean(intervalliRR) ...
            || intervalliRR(k) > (1 + SOGLIA_INTERVALLO_CONFIDENZA)*mean(intervalliRR)
        % Correzione presunto picco anomalo
        intervalliRR(k) = (intervalliRR(k-1)+intervalliRR(k+1))/2;
        % Se un intervalloRR viene ritenuto anomalo viene sostituito con la
        % media degli intervalli adiacenti.
    end
end

%---------------------- Interpolazione del segnale ------------------------

% Frequenza di campionamento
fs = FREQUENZA_DI_CAMPIONAMENTO; % [Hz]

% Asse temporale campionato a fs
t_RR = istantiRR(1):1/fs:istantiRR(length(istantiRR));
% Diverso dall'asse t su cui è fenitio il segnale ECG, si hanno meno
% campioni poiché il segnale si conclude con l'ultimo picco R rilevato.

t_RR = t_RR'; % Per comodità lo traspongo in vettore colonna.

% % Interpolazione lineare
% intervalliRR_interpolati = interp1(istantiRR,intervalliRR,t_RR','linear')'; 

% Interpolazione cubica
intervalliRR_interpolati = interp1(istantiRR,intervalliRR,t_RR','spline')'; 

% % Visualizzazione segnale post interpolazione:
% figure()
% plot(t_RR,intervalliRR_interpolati,'--')
% figure()
% histogram(intervalliRR_interpolati)

% ---------------------------- STAMPA A VIDEO? ----------------------------

if flagStampa
    
    % Stampa intervalli RR post pulizia
    figure('Name','Distribuzione RR','NumberTitle','Off')
    histogram(intervalliRR)
    xlabel('Tempo [s]')
    ylabel('Numero di ricorrenze')
    title('Distribuzione degli intervalli RR editati')

    % Stampa tacogramma post pulizia
    numeroBattiti = find(intervalliRR); 
    %numeroBattiti ha la stessa cardinalità di indQRS(1:end-1)/fs ma andiamo
    % a rappresentare gli intervalli RR in funzione del numero di battiti
    % cardiaci (ogni campione è un battito cardiaco).

    figure('Name','Tacogramma','NumberTitle','Off')
    plot(numeroBattiti,intervalliRR)
    xlabel('Numero campioni')
    ylabel('intervalli RR [s]')
    title('Tacogramma editato')
    clear numeroBattiti
    
end
end

