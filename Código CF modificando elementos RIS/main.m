% Vaciar espacio de trabajo y cerrar figuras
close all;
clear;

%% Setup de simulación
nbrOfSetups = 100;   % Número de escenarios
nbrOfRealizations = 100;    % Número de realizaciones

L = 100;         % Número de APs
N_AP = 1;        % Antenas por AP
K = 10;          % Número de UEs
S = 20;          % Número de RIS
tau_c = 200;     % Longitud del bloque de coherencia
tau_p = 10;      % Longitud del piloto
p = 100;         % Potencia de transmisión (mW)
fc = 3.5;         % Frecuencia (GHz)
LoS = 2;         % Linea de visión directa
% Desviación estándar angular en el modelo de dispersión local (en radianes)
ASD_varphi = deg2rad(15);  % angulo de azimut 

% Arreglos 3D para guardar resultados por tipo de canal 
SE_PMMSE_DCC = zeros(K, nbrOfSetups, 3);  
%SE_MR_DIST   = zeros(K, nbrOfSetups, 3);

%% Numero de elementos por RIS
N_RIS_values = [50, 100, 200];

for n_ris = 1:length(N_RIS_values)
    N_RIS = N_RIS_values(n_ris);
    for n = 1:nbrOfSetups
        disp(['Setup ' num2str(n) '/' num2str(nbrOfSetups) ' asistido por ' num2str(N_RIS) ' elementos por RIS']);
    
        % Generar escenario
        [R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,pilotIndex,D,HMean_AP_UE, HMean_AP_RIS, HMean_RIS_UE, probLoS_AP_UE, probLoS_RIS_UE] = setup(L,K,N_AP,N_RIS,tau_p,n,ASD_varphi,LoS,fc,S);
        
        % Asignacion de RIS
        risAssignment = assignRIS(probLoS_AP_UE, probLoS_RIS_UE);
    
        % Estimar canales
        [Hhat,H_eq,R_eq,B,C] = channelEstimates(R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,nbrOfRealizations,L,K,S,N_AP,N_RIS,tau_p,pilotIndex,p,HMean_AP_UE,HMean_AP_RIS, HMean_RIS_UE,risAssignment);
    
        % Calcular SE
        [SE_P_MMSE, SE_MR_dist] = SE_uplink(Hhat,H_eq,D,B,C,tau_c,tau_p,nbrOfRealizations,N_AP,K,L,p,R_eq,pilotIndex);
    
        % Guardar resultados en la dimensión
        SE_PMMSE_DCC(:,n,n_ris) = SE_P_MMSE;
        %SE_MR_DIST(:,n,n_ris)  = SE_MR_dist;
    
        clear Hhat H_eq B C R_eq;
    end
end

%% Graficar resultados
figure; hold on; box on;
set(gca,'fontsize',16);

% P-MMSE 
%plot(sort(SE_PMMSE_DCC(:)), linspace(0,1,K*nbrOfSetups), 'k-', 'LineWidth', 2);

% P-MMSE 
aux1 = SE_PMMSE_DCC(:,:,1); % 0 RIS
aux2 = SE_PMMSE_DCC(:,:,2); % 5 RIS
aux3 = SE_PMMSE_DCC(:,:,3); % 10 RIS

plot(sort(aux1(:)), linspace(0,1,K*nbrOfSetups), 'k-', 'LineWidth', 2);
plot(sort(aux2(:)), linspace(0,1,K*nbrOfSetups), 'r-',  'LineWidth', 2);
plot(sort(aux3(:)), linspace(0,1,K*nbrOfSetups), 'g-', 'LineWidth', 2);

% MR Distribuido 
%plot(sort(SE_MR_DIST(:)), linspace(0,1,K*nbrOfSetups), 'r-', 'LineWidth', 2);

% Ejes y leyenda
xlabel('Spectral efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
ylabel('CDF', 'Interpreter', 'Latex');
%legend({'P-MMSE ', 'MR dist'}, 'Interpreter', 'Latex', 'Location', 'SouthEast');
legend({'P-MMSE 50 N-RIS', 'P-MMSE 100 N-RIS', 'P-MMSE 200 N-RIS'}, 'Interpreter', 'Latex', 'Location', 'SouthEast');
xlim([0 25]);
