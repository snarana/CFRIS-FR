% Vaciar espacio de trabajo y cerrar figuras
close all;
clear;

%% Setup de simulación
nbrOfSetups = 100;   % Número de escenarios 
nbrOfRealizations = 1000;    % Número de realizaciones

L = 100;         % Número de APs
N = 1;           % Antenas por AP
K = 10;          % Número de UEs
tau_c = 200;     % Longitud del bloque de coherencia
tau_p = 10;      % Longitud del piloto
p = 100;         % Potencia de transmisión (mW)
fc = 8;         % Frecuencia (GHz)

% Desviación estándar angular en el modelo de dispersión local (en radianes)
ASD_varphi = deg2rad(15);  % angulo de azimut 
ASD_theta = deg2rad(15);   % angulo de elevación

% Arreglos 3D para guardar resultados por tipo de canal (LoS=0,1,2)
SE_PMMSE_DCC = zeros(K, nbrOfSetups, 3);  
SE_MR_DIST   = zeros(K, nbrOfSetups, 3);

%% Bucle sobre los tipos de canal
for LoS = 0:2
    for n = 1:nbrOfSetups
        disp(['Setup ' num2str(n) '/' num2str(nbrOfSetups) ' for LoS = ' num2str(LoS)]);

        % Generar escenario
        [R,pilotIndex,D,HMeanWithoutPhase] = setup(L,K,N,tau_p,n,ASD_varphi,ASD_theta,LoS,fc);

        % Estimar canales
        [Hhat,H,B,C] = channelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p,HMeanWithoutPhase);

        % Calcular SE
        [SE_P_MMSE, SE_MR_dist] = SE_uplink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex);

        % Guardar resultados en la dimensión LoS+1
        SE_PMMSE_DCC(:,n,LoS+1) = SE_P_MMSE;
        SE_MR_DIST(:,n,LoS+1)  = SE_MR_dist;

        clear Hhat H B C R;
    end
end

%% Graficar resultados
figure; hold on; box on;
set(gca,'fontsize',16);

% P-MMSE (negro)
aux1 = SE_PMMSE_DCC(:,:,1); % NLOS
aux2 = SE_PMMSE_DCC(:,:,2); % LOS
aux3 = SE_PMMSE_DCC(:,:,3); % NLOS/LOS

plot(sort(aux1(:)), linspace(0,1,K*nbrOfSetups), 'k--', 'LineWidth', 2);
plot(sort(aux2(:)), linspace(0,1,K*nbrOfSetups), 'k-',  'LineWidth', 2);
plot(sort(aux3(:)), linspace(0,1,K*nbrOfSetups), 'k-.', 'LineWidth', 2);

% MR Distribuido (rojo)
aux1 = SE_MR_DIST(:,:,1);
aux2 = SE_MR_DIST(:,:,2);
aux3 = SE_MR_DIST(:,:,3);

plot(sort(aux1(:)), linspace(0,1,K*nbrOfSetups), 'r--', 'LineWidth', 2);
plot(sort(aux2(:)), linspace(0,1,K*nbrOfSetups), 'r-',  'LineWidth', 2);
plot(sort(aux3(:)), linspace(0,1,K*nbrOfSetups), 'r-.', 'LineWidth', 2);

% Ejes y leyenda
xlabel('Spectral efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
ylabel('CDF', 'Interpreter', 'Latex');
legend({ ...
    'P-MMSE NLOS', 'P-MMSE LOS', 'P-MMSE NLOS/LOS', ...
    'MR dist NLOS', 'MR dist LOS', 'MR dist NLOS/LOS' ...
    }, 'Interpreter', 'Latex', 'Location', 'SouthEast');
xlim([0 20]);
