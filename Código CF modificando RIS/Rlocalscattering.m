function R = Rlocalscattering(M,angle,ASDdeg,antennaSpacing)
% Genera la matriz de correlación espacial para un canal con dispersión local.
%
% Notas:
%   - La matriz de correlación se calcula utilizando un modelo de dispersión
%     local gaussiano.
%   - El ángulo de llegada y la desviación estándar angular determinan la
%     correlación entre las señales recibidas en diferentes antenas.
%   - El espaciado entre antenas afecta la fase de la correlación.

ASD = ASDdeg;
diff_angles = repmat(0:(M-1), M, 1) - repmat((0:(M-1))', 1, M);
R = exp(-((pi * antennaSpacing * diff_angles * sin(angle)).^2) / (2 * ASD^2));
% Rotate the correlation matrix according to the nominal angle of departure
varphi = (0:M-1)'*2*pi*antennaSpacing*sin(angle);
rotationVector = exp(1i*varphi);
R = rotationVector*rotationVector'.*R;
end