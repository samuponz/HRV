% SAMUELE PONZIN - 719596 - INGEGNERIA ELETTRONICA E DELLE TELECOMUNICAZIONI 
% DII - DIPARTIMENTO DI INGEGNERIA DELL'INFORMAZIONE
% UNIVERSITA' DEGLI STUDI DI BRESCIA


%  LA VARIABILITA' DELL'HEART RATE: COME IL SISTEMA NERVOSO AGISCE SULLA
%                          FREQUENZA CARDIACA


% --------- SETUP: Estrazione dei dati selezionati dal database -----------

% WaveForm DataBase (WFDB) Toolbox

% Non è necessario installare il toolbox, tutti i segnali sono presenti
% nella cartella del progetto in formato *.mat
% I segnali sono stati recuperati dal database del sito Physionet.org, che
% mette a disposizione un particolare tool per scaricare direttamente i
% segnali trovati nel database e convertirli in *.mat per poterli caricare
% facilmente nel proprio workspace.

% Comandi utili:

% wfdb:      - apre la documentazione del tool
% pbsearch:  - apre la pagina del sito per cercare i segnali secondo determinati parametri
% rdsamp:    - importa il segnale in formato standard *.dat nel workspace
% wfdb2mat:  - converte un segnale dal formato standard *.dat al formato *.mat. 
% rdmat:     - importa il segnale in formato *.mat (genereato dal comando WFDB2MAT)

% Sono stati analizzati i segnali ECG provenienti dal database FANTASIA,
% presente su PhysioNet.org. Tutte le informazioni sono presenti al seguente 
% link: https://physionet.org/content/fantasia/1.0.0/

% Prendiamo ad esempio tre segnali ECG:

% ECG1) sesso: F; età: 32 anni.

% wfdb2mat('f1y09') 
% % Per il comando precedente devono essere presenti della cartella del 
% % progetto i file "f1y09.dat" e "f1y09.hea".
% [tm, signal, ~, ~] = rdmat('f1y09m'); 
% % Gli output ignorati sono rispettivamente la frequenza di campionamento 
% % e una struttura contentente altre informazioni sui segnali.
% t = tm'; 
% RAWecg = signal(:,2); 
% % Estrazione dei dati di interesse, nella prima colonna è presente
% % l'attività respiratoria del soggetto durante la registrazione.
% save('ECG1','t','RAWecg'); 
% % Salvataggio dei  dati di interesse in formato *.mat.

% ECG2) sesso: M; età: 76 anni.

% wfdb2mat('f1o05')
% [tm,signal, ~, ~]=rdmat('f1o05m');
% t = tm'; 
% RAWecg = signal(:,2);
% save('ECG2','t','RAWecg'); 

% ECG3) sesso: F; età: 21 anni.

% wfdb2mat('f1y10')
% [tm,signal, ~, ~]=rdmat('f1y10m');
% t = tm'; 
% RAWecg = signal(:,2);
% save('ECG3','t','RAWecg');

% -------------------------- ANALISI DELL'HRV -----------------------------

clear
close all

ECG1 = load('ECG1');
ECG2 = load('ECG2');
ECG3 = load('ECG3'); 

% Selezionare quale segnale ECG analizzare:
list = {'ECG1 (sesso: F; età: 32 anni)','ECG2 (sesso: M; età: 76 anni)','ECG3 (sesso: F; età: 21 anni)'};
[indx,tf] = listdlg('PromptString',{'Seleziona un segnale ECG;',...
    'Puoi selezionare un solo ECG alla volta.',''},'SelectionMode','single','ListString',list);

switch indx
    case 1
        fprintf('\nECG1 (sesso: F; età: 32 anni):\n');
        studioHRV(ECG1); 
    case 2
        fprintf('\nECG2 (sesso: M; età: 76 anni):\n');
        studioHRV(ECG2);
    case 3
        fprintf('\nECG3 (sesso: F; età: 21 anni):\n');
        studioHRV(ECG3);
    otherwise
        fprintf('Scelta non valida, si assume come scelta ECG1.');
        fprintf('\nECG1 (sesso: F; età: 32 anni):\n');
        studioHRV(ECG1);
end

clear indx list tf

function studioHRV(ECG)

% ------------------------- Estrazione dei dati ---------------------------

RAWecg = ECG.RAWecg;
t = ECG.t;

fs = 1/mean(diff(t));

% ------------------------------- STANDARD --------------------------------
% Decidere se effettuare lo studio sull'intero segnale o su una finestra di
% 5 minuti. Necessario valutare registrazioni della stessa durata per poter 
% confrontare i risultati dello studio della HRV di diversi segnali ECG.

standard = false;

prompt = 'Analizzare solo i primi 5 minuti del segnale? s/n: ';
str = input(prompt,'s');
if str == 's'
    standard = true;
end

clear prompt

istanteIniziale = 1; % [s]

% Si dà per scontato che le registrazioni siano minimo di 5 minuti.
if standard
    durataStandard = 300; % [s]
    campioneIniziale = istanteIniziale*fs;
    numeroCampioniStandard = durataStandard*fs;
    t = t(campioneIniziale:(campioneIniziale + numeroCampioniStandard - 1));
    RAWecg = RAWecg(campioneIniziale:(campioneIniziale + numeroCampioniStandard - 1));
end

% Durata segnale analizzato
durataSegnaleIntero = (length(RAWecg)*(1/fs)/60); % [minuti]
fprintf("Registrazione di segnale ECG di %.2f minuti.  \n", durataSegnaleIntero);


% ---------------------------- STAMPA A VIDEO? ----------------------------

flagStampa = false;

prompt = 'Visualizzare tutti i grafici? s/n: ';
str = input(prompt,'s');
if str == 's'
    flagStampa = true;
end

clear prompt

% ------------------ ELABORAZIONE DEL SEGNALE DIGITALE --------------------

% Determinazione dei complessi QRS: algoritmo di Pan-Tompkins
[FIRecg, INTecg, indQRS, sogliaINT, sogliaFIR] = PanTompkins(RAWecg, t, flagStampa);

% -------------------- ANALISI DEL SEGNALE DIGITALE -----------------------

[intervalliRR,istantiRR] = calcoloRR(indQRS,t,flagStampa);

% Calcolo HR
HR = round(60/mean(intervalliRR)); %frequenza cardiaca nell'unità di battiti al minuto

% Editing intervalli RR
[intervalliRR_interpolati, t_RR] = editingRR(intervalliRR, istantiRR, flagStampa);
% Questo passo è necessario per garantire le migliori prestazioni dei
% metodi statitici e nel dominio della
% frequenza.

% Calcolo SDNN
SDNN = std(intervalliRR_interpolati*1000); % [ms]

% Analisi delle componenti frequenziali della HRV
[infoSNA, varianza] = freqHRV(intervalliRR_interpolati, t_RR, flagStampa);

% ----------------- Stampa informazioni principali ------------------------

% Stampa delle informazioni principali ottenute dal segnale ECG, per 
% visualizzare tutte le informazioni ottenute dal segnle è necessario
% abilitare il flag di stampa.

fprintf("Risultati:\n");

% Stampa HR
fprintf("Heart rate: %d battiti al minuto. \n", HR);

% Stampa SDNN
fprintf("SDNN: %.3f ms \n", SDNN);

% Stampa varianza dalla PSD
fprintf("Varianza del segnale ricavata dalla PSD degli intervalli RR: %.3f ms \n", varianza*1000);

% Stampa deviazione standard ricavata dalla PSD
fprintf("Deviazione standard ricavata dalla PSD del\nsegnale degli intervalli RR: %.3f ms \n", sqrt(varianza*1000));

% Stampa rapporto LF/HF
fprintf("Rapporto LF/HF: %.2f \n", infoSNA.rapportoLFHF);

% ------------------ Confronto segnali nel tempo --------------------------

% Grafico interattivo per confrontare direttamente alcuni dei segnali di
% interesse.

pause(1) % Per assicurarsi che i grafici vengano stampati in ordine

% Per comodità è possibile visualizzare una finestra temporale del segnale.
inizioFinestra = istanteIniziale + 20; % [s]
intervalloTempoVisualizzato = 3; % [s]
% Valori arbitrari, è solo per permettere una buona visualizzazione dei
% segnali.

Confronto = figure('Units','Normalized','Position',[0.6, 0.1, 0.8, 0.8],'Name','Confronto','NumberTitle','Off');
xlim([inizioFinestra inizioFinestra+intervalloTempoVisualizzato])
hold on
plotbrowser('on')

yline(mean(intervalliRR),'--','LineWidth',2) % Per questo serve MATLAB R2018b o una versione successiva
yline(sogliaINT,':','LineWidth',2) % Per questo serve MATLAB R2018b o una versione successiva
yline(sogliaFIR,'-.','LineWidth',2) % Per questo serve MATLAB R2018b o una versione successiva

h = stem(istantiRR,intervalliRR); % Visualizzati con stem solo per comodità
set(h, 'Marker', '*','Markersize', 3)
plot(t, FIRecg)
plot(t, INTecg)

xlabel('Tempo [s]')
ylabel('Segnale [mV] & Intervalli RR [s]')
title(['Confronto nel tempo dei segnali in una finestra di ', num2str(intervalloTempoVisualizzato), ' secondi'])
legend('Media RR','sogliaINT','sogliaFIR','intervalliRR','FIRecg','INTecg','Location','northeast','NumColumns',3)

print(Confronto,'Grafici\ConfrontoSegnaliTempo','-dpng')

clear sogliaINT sogliaFIR h

% ------------- Stampa attività sistema nervoso autonomo ------------------

pause(1) % Per assicurarsi che i grafici vengano stampati in ordine

dati = [double(infoSNA.simpatico),double(infoSNA.parasimpatico)];
SNA = figure('Units','Normalized','Position',[0.3, 0.3, 0.4, 0.4],'Name','SNA','NumberTitle','Off');
explode = [0 1];
pie(dati, explode);
legend({'Attività simpatica','Attività vagale'},'Location','southeast')
title('Blinciamento delle attività dei rami del sistema nervoso autonomo')

print(SNA,'Grafici\BilanciamentoSNA','-dpng')

end





