function R = Rlocalscattering(N,varphi,theta,ASD_varphi,ASD_theta,antennaSpacing)
%Generar la matriz de correlación espacial para el modelo de dispersión local, para la distribución angular gaussiana y ULA
% Establecer el espaciado de antena si no se especifica en la entrada
if  nargin < 6    
    antennaSpacing = 1/2;    % Mitad de distancia de longitud de onda   
end

%Preparar para calcular la primera fila de la matriz
firstRow  = zeros(N,1);
firstRow(1) = 1;

%% Recorrer todas las columnas de la primera fila
for column = 2:N  
    distance = antennaSpacing*(column-1);    % Distancia desde la primera antena
    
    %Considerar diferentes casos dependiendo de la extensión del ángulo
    if (ASD_theta>0) && (ASD_varphi>0)

        F = @(delta,epsilon)exp(1i*2*pi*distance*sin(varphi+delta).*cos(theta+epsilon)).*...
            exp(-delta.^2/(2*ASD_varphi^2))/(sqrt(2*pi)*ASD_varphi).*...
            exp(-epsilon.^2/(2*ASD_theta^2))/(sqrt(2*pi)*ASD_theta);

        firstRow(column) = integral2(F,-20*ASD_varphi,20*ASD_varphi,-20*ASD_theta,20*ASD_theta);
        
    elseif ASD_varphi>0

        F = @(delta)exp(1i*2*pi*distance*sin(varphi+delta).*cos(theta)).*...
            exp(-delta.^2/(2*ASD_varphi^2))/(sqrt(2*pi)*ASD_varphi);

        firstRow(column) = integral(F,-20*ASD_varphi,20*ASD_varphi);
        
    elseif ASD_theta>0

        F = @(epsilon)exp(1i*2*pi*distance*sin(varphi).*cos(theta+epsilon)).*...
            exp(-epsilon.^2/(2*ASD_theta^2))/(sqrt(2*pi)*ASD_theta);

        firstRow(column) = integral(F,-20*ASD_theta,20*ASD_theta);
        
    else
        
        firstRow(column) = exp(1i*2*pi*distance*sin(varphi).*cos(theta));
        
    end
    
end

%Calcular la matriz de correlación espacial utilizando la estructura de
%Toeplitz y escala
R = toeplitz(firstRow);
R = R*(N/trace(R));
