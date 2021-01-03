function [intervalliRR,istantiRR] = calcoloRR(indQRS,t,flagStampa)
% Semplice calcolo del vettore degli intervalli RR

% Calcolo della frequenza di campionamento:
fs = 1/mean(diff(t));

% Calcolo degli intervalli RR
intervalliRR = diff(indQRS.*(1/fs));
istantiRR = indQRS(2:end)/fs;
% Escludo il primo valore di indQRS perché il primo valore del vettore
% intervalliRR è l'intervalloRR tra il primo e il secondo battito


% ---------------------------- STAMPA A VIDEO? ----------------------------

if flagStampa
    
    % Stampa della distribuzione degli intervalli RR 
    Istogramma = figure('Name','Distribuzione RR','NumberTitle','Off');
    histogram(intervalliRR)
    xlabel('Tempo [s]')
    ylabel('Numero di ricorrenze')
    title('Distribuzione degli intervalli RR')

    print(Istogramma,'Grafici\Distribuzione intervalli RR','-dpng')
    
    % Stampa tacogramma
    numeroBattiti = find(intervalliRR); 
    %numeroBattiti ha la stessa cardinalità di indQRS(1:end-1)/fs ma andiamo
    % a rappresentare gli intervalli RR in funzione del numero di battiti
    % cardiaci (ogni campione è un battito cardiaco).

    Tacogramma = figure('Name','Tacogramma','NumberTitle','Off');
    plot(numeroBattiti,intervalliRR)
    xlabel('Numero campioni')
    ylabel('intervalli RR [s]')
    title('Tacogramma')
    clear numeroBattiti
    
    print(Tacogramma,'Grafici\Tacogramma','-dpng')

end

end


