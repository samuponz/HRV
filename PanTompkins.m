function [FIRecg,INTecg,indQRS,sogliaINT,sogliaFIR] = PanTompkins(RAWecg,t,flagStampa)
% Implementazione dell'algoritmo di Pan-Tompkins

% ----------------------------- PARAMETRI ---------------------------------
% I parametri utilizzati sono tutti citati nella relazione.

DURATA_MASSIMA_COMPLESSO_QRS = 0.15; %150 ms
% Paramentro utilizzato per la costruzione del filtro integratore.

BP = [5 15]; % [5 - 15] Hz.
% Banda del filtro per la rimozione degli artefatti

DURATAMEDIABATTITO = 1; % 1 s
TOLLERANZA = 3*DURATAMEDIABATTITO; % 3 s 
% Tolleranza per il controllo di intervalliRR anomali. Se viene rilevato un
% intervallo RR con un picco anomalo viene effettuata nuovamente la
% registrazione.

% Calcolo della frequenza di campionamento:
fs = 1/mean(diff(t));

% --------------------- RIMOZIONE DEGLI ARTEFATTI -------------------------

% Calcolo dei coefficienti del filtro: FIR - Window-based - (Kaiser window)
ordine = 80;
a = 1;
b = fir1(ordine, BP./(fs/2), kaiser(ordine+1,8)); %ordine, dimensione finestra, tipo di finestra

% Costruzione del filtro:
% 1-D LCCDE filter: i coefficienti del numeratore e del denominatore (a e b)
% definiscono la funzione di trasferimento del filtro
% a = 1 => FIR
FIRecg = filter(b,a,RAWecg);

% Risposta in frequenza del filtro:
% n è il numero di punti di valutazione (intero positivo maggiore uguale a 2).
% Se n non è specificato è fissato a 512. 
% n deve essere maggiore dell'ordine del filtro.
% Senza parametri di uscita stampa direttamente a video il modulo della risposta in
% frequenza
% freqz(b, 1, 512, fs);

%-----------------------------DERIVAZIONE----------------------------------

% Calcolo dei coefficienti:
a = 8;
b = [2 1 0 -1 -2];

% Costruzione del filtro:
DERecg = filter(b, a, FIRecg);

% Risposta in frequenza del filtro:
% freqz(b, a, 512, fs);

%-----------------------ELEVEAMENTO AL QUADRATO----------------------------

SQRecg= DERecg.^2;

%----------------------------INTEGRAZIONE----------------------------------

% DurataFinestra per la media mobile
DurataFinestra = DURATA_MASSIMA_COMPLESSO_QRS; % ms

% Calcolo del numero di campioni che identificano la finestra:
N = round(DurataFinestra * fs);

% Calcolo dei coefficienti:
a = N;
b = ones(1,N);

% Costruzione del filtro:
INTecg = filter(b, a, SQRecg);

% Risposta in frequenza del filtro:
% freqz(b, a, 512, fs);

%--------------------RILEVAMENTO DEI COMPLESSI QRS-------------------------

% Realizzazione del filtro derivatore per la ricerca dei massimi locali
% dei punti di pendenza massima del segnale integrato.

% Calcolo dei coefficienti:
a = 1;
b = [fs -fs]; % prendo solo due campioni

% Costruzione del filtro derivatore:
derINTecg = filter(b,a,INTecg);


% Finestra di 2 secondi per l'inizializzazione delle soglie usate per la
% corretta rilevazione dei complessi QRS.
campioneIniziale = 1;
campioneFinale = round(campioneIniziale + round(2*fs));

[indQRS,sogliaINT,sogliaFIR] = rilevamentoComplessiQRS(derINTecg,INTecg,FIRecg,fs,campioneIniziale,campioneFinale);
% Il seguente algoritmo rileva i picchi QRS. Se viene rivelata una distanza 
% temporale tra due complessi QRS (picchi R) consecutivi maggiore alla
% tolleranza specificata, viene effettuato nuovamente l'intero rilevamento dei
% complessi QRS basandosi su una nuova finestra iniziale (sempre di 2 
% secondi, semplicemente traslata nel tempo) per l'inizializzazione delle 
% soglie usate per la corretta rilevazione dei complessi QRS.
% Questo controllo viene fatto perché se nella finestra di inizializzazione
% compaiono valori anomali del segnale ECG le soglie vengono inizializzate
% a valori molto alti che portano ad ignorare i successivi complessi QRS,
% prima che la soglia si riporti a valori plausibili.

while any(diff(indQRS).*(1/fs) > TOLLERANZA) % Possibilità di ERRORE se presenti tanti valori anomali del segnale ECG
 
    campioneIniziale = campioneFinale;
    campioneFinale = round(campioneIniziale + round(2*fs));
    [indQRS,sogliaINT,sogliaFIR] = rilevamentoComplessiQRS(derINTecg,INTecg,FIRecg,fs,campioneIniziale,campioneFinale);

end

% ---------------------------- STAMPA A VIDEO? ----------------------------

if flagStampa

    titles = ["ECG grezzo","ECG filtrato dagli artefatti con il passa banda", ...
        "ECG Derivato","ECG Elevato al quadrato","ECG Integrato"];

    segnali = [RAWecg,FIRecg,DERecg,SQRecg,INTecg];

    PanTompkins = figure('Units','Normalized','Position',[0, 0, 1, 1],'Name','Elaborazione ECG', ...
        'NumberTitle','Off'); 

    % Guardo solo una finestra di 5 secondi del segnale per apprezzarne la forma
    finestra = (t > 60) & (t < 65);
    
    for i = 1:length(titles)
        subplot(5,1,i);
        plot(t(finestra), segnali((finestra),i));
        xlabel('Tempo [s]');
        ylabel('ECG [V]');
        title(titles(i));
    end
    
    % Ricavo un grafico con gli istanti dei picchi rilevati
    complessiQRS = zeros(1,length(t));
    complessiQRS(indQRS) = 5*mean(INTecg);
    
    % Stampo il grafico sovrapposto all'ultimo del ciclo for, ossia quello
    % dell'ECG integrato.
    hold on
    h = stem(t(finestra),complessiQRS(finestra));
    set(h, 'Marker', '.','Markersize', 3)
    
    print(PanTompkins,'Grafici/PanTompkins','-dpng')
    
end

end

