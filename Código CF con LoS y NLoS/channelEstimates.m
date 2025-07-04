function [Hhat,H,B,C] = channelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p,HMeanWithoutPhase)
%% Generar realizaciones del canal

% Generar canal Rician
M = L*N;
H = zeros(M, nbrOfRealizations, K); % Canal resultante
W = (randn(M, nbrOfRealizations, K) + 1i * randn(M, nbrOfRealizations, K));   % Ruido gaussiano complejo

% Media del canal
HMean=zeros(M,nbrOfRealizations,K); 
HMeanx=reshape(repmat(HMeanWithoutPhase,nbrOfRealizations,1),M,nbrOfRealizations,K); 

% Fase aleatoria para componente LoS
angles= -pi + 2*pi*rand(M,nbrOfRealizations,K);
phaseMatrix=exp(1i*angles);
  
% Canal Rician con correlación espacial 
for l = 1:L
    for k = 1:K
        
        HMean(:,:,k)= phaseMatrix(:,:,k).*HMeanx(:,:,k);  % Aplicar fase aleatoria
        Rsqrt = sqrtm(R(:,:,l,k));
        H((l-1)*N+1:l*N,:,k) = sqrt(0.5)*Rsqrt*W((l-1)*N+1:l*N,:,k) + HMean((l-1)*N+1:l*N,:,k);
       
    end
end
    
%% Estimación del canal

% Matriz identidad de tamaño NxN
eyeN = eye(N);

% Ruido normalizado para la estimación de piloto
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));

% Preparar para almacenar los resultados 
Hhat = zeros(L*N,nbrOfRealizations,K);

% Reservar matriz de correlación de estimación si se requiere
if nargout>2
    B = zeros(size(R));
end

% Reservar matriz de error de estimación si se requiere
if nargout>3
    C = zeros(size(R));
end

% Para cada AP
for l = 1:L
    
    % Para cada piloto
    for t = 1:tau_p
        
        % Señal recibida procesada
        yp = sqrt(p)*tau_p*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,l,t);
        yMean = zeros(N, nbrOfRealizations);  % Inicializar la media esperada de y
        
        % Matriz a invertir en el estimador MMSE
        PsiInv = (p*tau_p*sum(R(:,:,l,t==pilotIndex),4) + eyeN);
        
        % Para cada usuario que usa el piloto t
        for k = find(t==pilotIndex)'
            % Calcular estimación del canal
            RPsi = R(:,:,l,k) / PsiInv;
            yMean = yMean + sqrt(p)*tau_p*HMean((l-1)*N+1:l*N,:,k);
            Hhat((l-1)*N+1:l*N,:,k) = sqrt(p)*RPsi*(yp - yMean) + HMean((l-1)*N+1:l*N,:,k);
            
            % Correlación de estimación
            if nargout>2
                B(:,:,l,k) = p*tau_p*RPsi*R(:,:,l,k);
            end
            
            % Correlación del error de estimación
            if nargout>3
                C(:,:,l,k) = R(:,:,l,k) - B(:,:,l,k);
            end
            
        end
        
    end
    
end
