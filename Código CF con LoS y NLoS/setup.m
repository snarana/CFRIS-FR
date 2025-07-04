function [R,pilotIndex,D,HMean] = setup(L,K,N,tau_p,seed,ASD_varphi,ASD_theta,LoS,fc)
%% Definir escenario
HMean = zeros(N*L,K);

%Establecer el número de semilla si se especifica distinto de cero
if (nargin>4)&&(seed>0)
    rng(seed)
end

squareLength = 1000;    %Tamaño del área de cobertura
B = 20e6;   %Ancho de banda (Hz)
noiseFigure = 7;    % dB
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;    % Potencia de ruido calculada (dBm)
sigma_sf = 4;   %Desviación estándar del desvanecimiento de shadow
decorr = 9;    % Distancia de decorrelación del desvanecimiento de shadow

% Definir alturas (metros)
h_BS = 10;  % Estación base de antena
h_E = 1;    % Entorno efectivo de antena
h_UT = 1.5;  % Antena usuario

% Definir alturas efectivas
h_BS_eff = h_BS - h_E;
h_UT_eff = h_UT - h_E;

% Diferencia de altura entre AP y UE (metros)
distanceVertical = h_BS - h_UT;

% Definir el espaciado de la antena (en número de longitudes de onda)
antennaSpacing = 1/2;

c = 3e8;

%Prepararar para guardar resultados
gainOverNoisedB = zeros(L,K);   % SNR en dB por enlace
R = zeros(N,N,L,K);             % Matriz de correlación
distances = zeros(L,K);         % Distancias AP-UE
maxdistances = 300;             
probLoS = zeros(L,K);           % Probabilidad de LoS
pilotIndex = zeros(K);
D = zeros(L,K);                 % Atenuación del canal
betaLoS = zeros(L,K);
masterAPs = zeros(K,1);         % AP maestro por UE

%% Setups
    
    % Posiciones aleatorias de los APs con distribución uniforme
    APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;
    
    % Inicializar posiciones de los usuarios (UEs)
    UEpositions = zeros(K,1);
        
    % Calcular ubicaciones alternativas de los APs con wraparound
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
    
    % Matriz de correlación para shadowing
    shadowCorrMatrix = sigma_sf^2*ones(K,K);
    shadowAPrealizations = zeros(K,L);

    % Inicializar probabilidad de LoS
    if LoS==1
        probLoS=ones(L,K);
    elseif LoS == 0
        probLoS = zeros(L,K);
    end

    % Por cada usuario
    for k = 1:K
    
        % Posición aleatoria del usuario
        UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;
    
        % Calcular distancias 3D con altura de APs
        [distanceAPstoUE, whichpos] = min(abs(APpositionsWrapped - repmat(UEposition, size(APpositionsWrapped))), [], 2);
        distances(:,k) = sqrt(distanceVertical^2 + distanceAPstoUE.^2);
    
        % Modelo de path-loss UMi
        d_BP_eff = 4*h_BS_eff*h_UT_eff*((fc*1e9)./c);

        for l = 1:L
            d = distances(l,k);
            if d >= 10 && d <= d_BP_eff
                betaLoS(l) = 32.4 + 21 * log10(d) + 20 * log10(fc);
            elseif d > d_BP_eff && d <= 5000
                betaLoS(l) = 32.4 + 40 * log10(d) + 20 * log10(fc) - 9.5 * log10(d_BP_eff^2 + (h_BS - h_UT)^2);
            end
        end

        % Path-loss sin LoS
        PathNLoss_aux = 35.3 * log10(distances(:,k)) + 22.4 + 21.3 * log10(fc) - 0.3 * (h_UT - 1.5);
        betaNLoS = max(betaLoS,PathNLoss_aux);

        % Probabilidad aleatoria de LoS si LoS==2
        if LoS==2
            probLoS(:,k)=(rand<((maxdistances-distances(:,k))./maxdistances));
        end
    
        % Generación de shadowing
        if k-1 > 0
            shortestDistances = zeros(k-1,1);
            for i = 1:k-1
                shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
            end
            
            % Media y varianza condicional (teorema de estimación de Kay)
            newcolumn = sigma_sf^2 * 2.^(-shortestDistances / decorr);
            term1 = newcolumn' / shadowCorrMatrix(1:k-1,1:k-1);
            meanvalues = term1 * shadowAPrealizations(1:k-1,:);
            stdvalue = sqrt(sigma_sf^2 - term1 * newcolumn);

        else 
            % Primer usuario
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];
        end 

        % Realización del shadowing
        shadowing = meanvalues + stdvalue * randn(1, L);

        % Inicializar ganancias del canal
        channelGaindB = zeros(L,1);
        gain_LoS = zeros(L,1);
        gain_NLoS = zeros(L,1);
    
        % Asignar ganancias según si hay LoS
        channelGaindB(probLoS(:,k)==1)=-betaLoS(probLoS(:,k)==1);
        channelGaindB(probLoS(:,k)==0)=-betaNLoS(probLoS(:,k)==0);
        channelGain = db2pow(channelGaindB) .* db2pow(shadowing') ./ db2pow(noiseVariancedBm); % Convert to linear scale and include shadowing and noise
    
        % Cálculo del factor de Rician
        ricianFactor = 10.^(1.3 - 0.003 * distances(:,k));
        gain_LoS(probLoS(:,k)==1) = sqrt(ricianFactor(probLoS(:,k)==1) ./ (ricianFactor(probLoS(:,k)==1) + 1)) .* channelGain(probLoS(:,k)==1);
        gain_NLoS(probLoS(:,k)==1) = sqrt(1 ./ (ricianFactor(probLoS(:,k)==1) + 1)) .* (channelGain(probLoS(:,k)==1));
        gain_NLoS(probLoS(:,k)==0)=channelGain(probLoS(:,k)==0); %note that probLoS is always one in the manuscript
    
        % Almacenar ganancia sobre ruido en dB
        gainOverNoisedB(:,k) = pow2db(channelGain);
    
        % Actualizar matrices de shadowing
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
        shadowAPrealizations(k,:) = shadowing;
    
        UEpositions(k) = UEposition;    % Posición del usuario

        % Asignar AP maestro (con mejor canal)
        [~,master] = max(gainOverNoisedB(:,k));
        D(master,k) = 1;
        masterAPs(k) = master;
    
        % Asignar pilotos ortogonales
        if k <= tau_p
            pilotIndex(k) = k;
        else 
            % Elegir el piloto con menor interferencia
            pilotinterference = zeros(tau_p,1);

            for t = 1:tau_p
                pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1) == t)));
            end

            [~,bestpilot] = min(pilotinterference);
            pilotIndex(k) = bestpilot;
        end
    

        % Para cada AP
        for l = 1:L

            % Ángulos horizontal y vertical al UE
            angletoUE_varphi = angle(UEpositions(k) - APpositionsWrapped(l, whichpos(l)));
            angletoUE_theta = asin(distanceVertical / distances(l,k)); 
            
            % Matriz de correlación espacial (modelo de dispersión local)
            if nargin > 5
                R(:,:,l,k) = gain_NLoS(l) * Rlocalscattering(N, angletoUE_varphi, angletoUE_theta, ASD_varphi, ASD_theta, antennaSpacing);
            else
                R(:,:,l,k) = gain_NLoS(l) * eye(N);
            end
    
            % Componente determinista LoS
            arrayResp = exp(1i * pi * (0:N-1).' * sin(angletoUE_varphi)) / sqrt(N);
            HMean(l,k) = sqrt(gain_LoS(l)) * arrayResp;
        end
    end

 % Cada AP sirve al UE con mejor canal para cada piloto
    for l = 1:L
        
        for t = 1:tau_p
            
            pilotUEs = find(t==pilotIndex(:));
            [~,UEindex] = max(gainOverNoisedB(l,pilotUEs));
            D(l,pilotUEs(UEindex)) = 1;
           
        end
        
    end