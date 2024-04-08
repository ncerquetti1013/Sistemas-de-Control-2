function [X] = modMotor2(t_etapa, xant, accion, Tl)
% modMotor (función obtenida de los apuntes de clase)
%   Esta función modela de forma lineal el funcionamiento de un motor de CC
%   a partir de las ecuaciones
%   1) ia_p   = -Ra/Laa*ia - Km/Laa*wr + 1/Laa*Va
%   2) wr_p   = Ki/J*ia - Bm/J*wr - 1/J*TL
%   3) tita_p = wr
%
%   Args:
%   - t_etapa: duración de la etapa actual
%   - xant   : valores anteriores de las salidas de interés
%   - accion : entrada de control sobre el motor
%
%   Output:
%   - X: valores de las salidas de interés en el tiempo actual

% Parámetros del motor a modelar
Ki = 16.52;
J  = 4.13e-06;
B  = 7.465e-16;
Laa = 2.061e-9/J;
Ra = (8.25e-5 - Laa*B)/J;
%Ra = 28.13            % valor calculado a partir de los datos de corriente
                       % en excel
Km = (1 - Ra*B)/Ki;

Va = accion;          % tensión de armadura
h  = 1e-7;            % paso de integración

TL = Tl;

% Asignación de los valores actuales de las variables de interés para
% iniciar el proceso de cálculo
wr   = xant(1);
wrp  = xant(2); 
ia   = xant(3);
tita = xant(4);

for ii=1:t_etapa/h
    % Cálculo de las derivadas
    wrpp = (-wrp*(Ra*J + Laa*B) - wr*(Ra*B + Ki*Km) + Va*Ki)/(J*Laa);
    iap  = (-Ra*ia - Km*wr + Va)/Laa;

    % Aplicación del método de Euler para la estimación de los valores de
    % las variables de interés en el tiempo posterior
    wrp  = wrp + h*wrpp;
    wrp  = wrp - ((1/J)*TL);          % reducción por acción del torque
    ia   = ia + iap*h;                % cálculo de la corriente de armadura
    wr   = wr + h*wrp;                % cálculo de la velocidad angular
    tita = tita + h*wr;               % cálculo de la posición 
end

% Vector de salida
X = [wr, wrp, ia, tita];

end
