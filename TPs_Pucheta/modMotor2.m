function [X] = modMotor2(t_etapa, xant, accion, Tl)
% modMotor (funci�n obtenida de los apuntes de clase)
%   Esta funci�n modela de forma lineal el funcionamiento de un motor de CC
%   a partir de las ecuaciones
%   1) ia_p   = -Ra/Laa*ia - Km/Laa*wr + 1/Laa*Va
%   2) wr_p   = Ki/J*ia - Bm/J*wr - 1/J*TL
%   3) tita_p = wr
%
%   Args:
%   - t_etapa: duraci�n de la etapa actual
%   - xant   : valores anteriores de las salidas de inter�s
%   - accion : entrada de control sobre el motor
%
%   Output:
%   - X: valores de las salidas de inter�s en el tiempo actual

% Par�metros del motor a modelar
Ki = 16.52;
J  = 4.13e-06;
B  = 7.465e-16;
Laa = 2.061e-9/J;
Ra = (8.25e-5 - Laa*B)/J;
%Ra = 28.13            % valor calculado a partir de los datos de corriente
                       % en excel
Km = (1 - Ra*B)/Ki;

Va = accion;          % tensi�n de armadura
h  = 1e-7;            % paso de integraci�n

TL = Tl;

% Asignaci�n de los valores actuales de las variables de inter�s para
% iniciar el proceso de c�lculo
wr   = xant(1);
wrp  = xant(2); 
ia   = xant(3);
tita = xant(4);

for ii=1:t_etapa/h
    % C�lculo de las derivadas
    wrpp = (-wrp*(Ra*J + Laa*B) - wr*(Ra*B + Ki*Km) + Va*Ki)/(J*Laa);
    iap  = (-Ra*ia - Km*wr + Va)/Laa;

    % Aplicaci�n del m�todo de Euler para la estimaci�n de los valores de
    % las variables de inter�s en el tiempo posterior
    wrp  = wrp + h*wrpp;
    wrp  = wrp - ((1/J)*TL);          % reducci�n por acci�n del torque
    ia   = ia + iap*h;                % c�lculo de la corriente de armadura
    wr   = wr + h*wrp;                % c�lculo de la velocidad angular
    tita = tita + h*wr;               % c�lculo de la posici�n 
end

% Vector de salida
X = [wr, wrp, ia, tita];

end
