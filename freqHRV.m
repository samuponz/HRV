function [infoSNA, varianza] = freqHRV(intervalliRR_interpolati,t_RR,flagStampa)
% Calcolo delle componenti frequenziali del segnale degli intervalli RR, in
% modo da osservare l'attività del sistema nervoso autonomo.

% L'asse dei tempi t_RR serve solo per la stampa dei grafici.

% ----------------------------- PARAMETRI ---------------------------------
% I parametri utilizzati sono tutti citati nella relazione.

% Componenti dello spettro del segnale
VLF = [0.003, 0.04];
LF = [0.04, 0.15];
HF = [0.15, 0.4];

% Frequenza di campionamento del segnale 
fs = 1/(t_RR(2)-t_RR(1));

%------------------- Rimozione media / trend polinomiali ------------------

intervalliRR_interpolati = intervalliRR_interpolati - mean(intervalliRR_interpolati);
% intervalliRR_interpolati = detrend(intervalliRR_interpolati);

%---------------------------- Stima della PSD -----------------------------

[PSD,f] = myWelch(intervalliRR_interpolati,t_RR,fs,flagStampa);
% Implementazione del metodo di Welch, mi ha fatto sudare.

%--------------------------- Analisi della PSD ----------------------------

[infoSNA, varianza] = calcoloPotenze(PSD,f,VLF,LF,HF,flagStampa);
% Analisi della PSD ricavata con il metodo di Welch

end

function [PSD, f] = myWelch(x,t,fs,flagStampa)
% Possiamo stimare la PSD calcolando il modulo quadro della FFT del segnale.
% Per avere la vera PSD dovremmo avere un numero di campioni infinito, in
% questo modo la media della stima del valore atteso della PSD coincide 
% con la media PSD reale (non c'è BIAS). Nel caso della varianza invece
% anche con un numero di campioni infinto abbiamo una varianza consistente.
% Metodo di Welch: di ridurre la varianza passando alla popolazione delle
% medie, mediando più periodogrammi. Quando passiamo alla popolazione delle 
% medie, se mediamo su segmenti che sono statisticamente indipendenti
% otteniamo una riduzione della varianza. Eseguendo la media di k variabili 
% casuali statisticamente indipendenti otteniamo una riduzione della 
% varianza di un fattore k.

% Dividiamo il segnale di lunghezza Q in segmenti di lunghezza L.
% Eseguiamo il periodogramma di ogni segmento, che ha la forma:
% x[n + mR]*w[n] con m = 0,1,...,K-1 e n = 0,1,...,L-1.
% Se R (hop size) è minore di L (dimensione segmento) abbaimo un overlap.
% La finestra w[n] non è strettamente necessaria, se introdotta però è da
% considerare anche nella normalizzazione della PSD.

% ----------------------------- PARAMETRI ---------------------------------
% I parametri utilizzati sono tutti citati nella relazione.

% Lughezza del singolo segmento del segnale, coincide con la dimenione della finestra.
L = 256;

SOVRAPPOSIZIONE_PERCENTUALE = 0.5; % Sovrapposizione del 50%
% Determina la sovrapposizione tra i segmenti del segnale.

% ------------------- ---- IMPLEMENTAZIONE METODO -------------------------

% x[m] con m = 0,1,...,Q-1 %intervalli RR [s]
Q = length(x); % Per comodità

% TIPOLOGIA FINESTRA:
w = hanning(L); 
% w = hamming(L);
% Altre tipologie di fineatre: Blackman, Blackman-Harris, Kaiser-Bessel.

% w = ones(1,L); 
% finestra rettangolare. Se selezionata, è neceaario aggiustare la
% normalizzazione della PSD. 

% Hop size
R = ceil(L*(1 - SOVRAPPOSIZIONE_PERCENTUALE));

% Costruzione segmenti (periodogrammi):
% Di fatto è un generico finestramento; se non si moltiplica il segmento 
% per una funzione finestra, sarà come aver finestrato con un rect.

% Numero di finestre totali:
K = floor((Q - L)/R) + 1;
% Se Q non è divisibile per K le finestre non coprono tutto il segnale e
% vendono tagliati fuori dei valori in fondo al segnale, al massimo 
% R(Hop size) - 1 campioni. 

segmenti = zeros(L,K); % preallocazione della matrice che conterrà
% ogni segmento. Ho K segmenti di lunghezza L.

% Dato che dovremo effettuare una moltiplicazione campione per campione tra
% x e w questi dovranno avere la stessa diemnsione. Traspongo se
% necessario.

if isrow(x) 
    x = x';
    % voglio lavorare con vettori colonna. la funzione fft() eseguita su
    % matrice riconosce come segnali i vettori colonna.
end

% Creazione figura per la stampa dei segmenti finestrati
if flagStampa
    Segmenti = figure('Name','Segmenti','NumberTitle','Off');
    xlabel('Tempo [s]')
    ylabel('intervalli RR [s] (a media nulla e campionaeti uniformemente)')
    title('Segmenti del segnale degli intervalli RR ricampionati uniformemente')
    hold on
end

% Creazione dei segmenti
for r = 0:K-1  
    
    campioneIniziale = r*R + 1;
    campioneFinale = campioneIniziale + (L - 1) ;
    
    tr = t(campioneIniziale:campioneFinale);
    xr = x(campioneIniziale:campioneFinale);
    
    % Segmento del segnale finestrato
    v = xr .* w; % entrambi vettori colonna
    
    % Stampa dei segmenti finestrati
    if flagStampa 
        plot(tr,v)
    end
    
    % Raccolgo i segmenti in una matrice
    segmenti(:, r + 1) = v; 
    
end

if flagStampa
    print(Segmenti,'Grafici\SegmentiFinestrati','-dpng')
end

nfft = 2^nextpow2(L);
    
% FFT dei segmenti (si passa dal tempo discreto alla frequenza discreta) 
V = fft(segmenti,nfft);  
% Viene effettuata la DFT su ogni segmento del segnale, di lunghezza L. 

% % Normalizzazione FFT
V = V/L; 
% La normalizzazione viene effettuata dividendo la fft per la
% lunghezza del segnale a tempo discreto originale, nel nostro caso il
% segmento finestrato v (NON si tiene conto dello zero padding).

%definizione asse delle frequenze: 
f = (0:nfft - 1)*(fs/nfft); % [Hz]
% La DFT su nfft punti ha un passo di campionamento di 1/nfft.
% Passiamo all'asse delle frequenze moltiplicando ogni campione per la 
% frequenza di campionamento. 

% Stampa FFT normalizzata in funzione delle frequenze
if flagStampa 
    figure('Name','Singole FFT','NumberTitle','Off')
    plot(f,abs(V))
    xlabel('Frequenza [Hz]')
    ylabel('FFT normalizzata [s]')
    title('FFT normalizzate dei segmenti del segnale')
end

% Potenza del segmento in funzione della frequenza (periodogramma)
P = abs(V).^2; 

% Normalizzazione Potenza della FFT
P = (1/((w'*w)/L))*(1/(fs/L))*P; %(????? CONROLLA ?????)
% Si normalizza per il fattore sum(w.^2)/L per compensare la perdita di
% energia dovuta all'utilizzo della finestra.

% P = (1/L)*P;
% Normalizzazione nel caso in cui w = ones(L,1) (finestra rettangolare)


% Stampa della PSD normalizzata di ogni segmento
if flagStampa
    figure('Name','Singole PSD','NumberTitle','Off')
    plot(f,P)
    xlabel('Frequenza [Hz]')
    ylabel('PSD [s^2/Hz]')
    title('PSD normalizzate dei segmenti del segnale')
end

% Stima della PSD come media delle PSD dei segmenti
PSD = (1/K)*sum(P,2);

% PSD unilatera
puntiSingoli = round(nfft/2) + 1;
PSD = PSD(1:puntiSingoli);
PSD(2:end-1) = 2* PSD(2:end-1);
% Essendo il segnale originale un segnale reale la sua PSD è pari, è 
% possibile valutare la PSD unilatera.
% Per la conservazione della potenza è necessario raddoppiare l'area della 
% PSD unilatera.

f = (0:puntiSingoli - 1)*(fs/nfft);
% Asse delle frequenze da 0 a fs/2 - 1.

end

function [infoSNA, varianza] = calcoloPotenze(PSD,f,VLF,LF,HF,flagStampa)
% Analisi della PSD stiamata con il metodo di Welch

% Bande delle varie componenti in frequenza
intervalloVLF = (f >= VLF(1)) & (f <= VLF(2));
intervalloLF = (f >= LF(1)) & (f <= LF(2));
intervalloHF = (f >= HF(1)) & (f <= HF(2));

% Calcolo delle aree (potenze assolute) delle varie componenti in frequenza
areaVLF = trapz(f(intervalloVLF),PSD(intervalloVLF));
areaLF = trapz(f(intervalloLF),PSD(intervalloLF));
areaHF = trapz(f(intervalloHF),PSD(intervalloHF));
areaTotale = areaVLF + areaLF + areaHF;

varianza = areaTotale;

% Calcolo dei massimi valori di potenza delle componenti frequenziali

% Osservazione dell'attività del sistema nervoso autonomo
infoSNA.simpatico = areaLF/(areaLF+areaHF);
% Mostra l'attività simpatica e vagale del sistema nervoso autonomo.
infoSNA.parasimpatico = areaHF/(areaLF+areaHF);
% Mostra l'attività vagale del sistema nervoso autonomo.

% Rapporto LF/HF
infoSNA.rapportoLFHF = areaLF/areaHF;
% Mostra la predominanza dell'attività simpatica su quella vagale.

if flagStampa
    
    % Rilevamento dei punti del grafico in cui stampare i testi
    [~, indice] = min(abs(f - (VLF(2) - VLF(1))/2));
    puntoVLF = [double(f(indice)),double(max(PSD(intervalloVLF)))];
    [~, indice] = min(abs(f - (LF(1) + (LF(2) - LF(1))/2)));
    puntoLF = [double(f(indice)),double(max(PSD(intervalloLF)))];
    [~, indice] = min(abs(f - (HF(1) + (HF(2) - HF(1))/2))); 
    puntoHF = [double(f(indice)),double(max(PSD(intervalloHF)))];

    % Stampa della stima della PSD
    StimaPSD = figure('Name','PSD','NumberTitle','Off');
    hold on
    xline(VLF(1),':','LineWidth',1) % Per questo serve MATLAB R2018b o una versione successiva
    xline(VLF(2),':','LineWidth',1) % Per questo serve MATLAB R2018b o una versione successiva
    xline(LF(2),':','LineWidth',1) % Per questo serve MATLAB R2018b o una versione successiva
    xline(HF(2),':','LineWidth',1) % Per questo serve MATLAB R2018b o una versione successiva
    xlim([0 0.6])
    plot(f,PSD)
    text(puntoVLF(1),puntoVLF(2),'VLF')
    text(puntoLF(1),puntoLF(2),'LF')
    text(puntoHF(1),puntoHF(2),'HF')
    xlabel('Frequenza [Hz]')
    ylabel('PSD [s^2/Hz]')
    title('PSD degli intervalli RR')
    
    print(StimaPSD,'Grafici\Stima PSD','-dpng')

end

end

