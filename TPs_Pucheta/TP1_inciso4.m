% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr. Ing. Pucheta, Julian
% Alumno: Cerquetti, Narella
% Actividad Práctica N1 - Caso de estudio 2
%   Inciso 4  
%   Obtener el torque máximo que puede soportar el motor modelado mediante las Ecs. (5)
%   (6) y (7) cuando se lo alimenta con 12V, graficando para 5 segundos de tiempo la
%   velocidad angular y corriente ia para establecer su valor máximo como para dimensionar 
%   dispositivos electrónicos.
% 

%%
clear all;
%%
% PARÁMETROS DEL MOTOR
La  = 366e-6; 
J   = 5e-9; 
Ra  = 55.6; 
Bm  = 0; 
Ki  = 6.49e-3; 
Km  = 6.53e-3;     
E   = 12;                     % tensión de actuación sobre el motor

%%
% MODELADO DEL MOTOR EN EL ESPACIO DE ESTADOS
% Este modelado permite obtener fácilmente las funciones de transferencia
% para las salidas de interés respecto a la tensión aplicada y el torque,
% considerando a éste último como una entrada de perturbación

% Variables de estado.
% x1 = ia, x2 = wr, x3 = theta
% Entradas
% u1 = E, u2 = TL
% Salidas.
% y1 = wr
A = [-Ra/La -Km/La 0 ; Ki/J -Bm/J 0 ; 0 1 0];    % matriz de estados
B = [1/La 0 ; 0 -1/J ; 0 0];                     % matriz de entrada 
C = [0 1 0];                                     % matriz de salida                                                       
D = [0 0];                                       % matriz de transmisión directa   
%%
% SIMULACIÓN

% Condiciones iniciales
% x = [ia; wr; theta]
Xop      = [0 0 0]';                   % punto de operación
wr(1)    = 0;
ia(1)    = 0;
theta(1) = 0;

x = [ia(1) wr(1) theta(1)]';
%%
% VECTORES DE SIMULACIÓN
tF   = 5;                              % tiempo de finalización de la simulación en secs
tint = 10e-6;                          % paso de integración
t    = 0:tint:tF;                      % vector de tiempo para gráfica 
um   = E*ones(1, length(t));           % vector de acción sobre el motor

% se genera el vector de torque
TLmax    = 0.0014;
TL       = zeros(1, length(t));
for i=1:tF/tint
    tt = t(i);
    if (tt>=0.7 && tt<=1.5)
        TL(i) = TLmax;
    else
        TL(i) = 0;
    end
end

plot(t, TL); grid; hold on;
%%
% SIMULACIÓN DEL SISTEMA
for ii=1:(length(t)-1)
    % Cálculo de las derivadas
    xp = A*(x - Xop)+B*[um(ii) TL(ii)]';
    x = x + xp*tint;
    wr(ii+1) = C*x;
    ia(ii+1) = x(1);
    theta(ii+1) = x(3);
end

%%
% Gráficas

% Velocidad angular
subplot(5,1,1);hold on;
plot(t,wr,'r');title('\omega, [rad/seg]');
% Posición
subplot(5,1,2);hold on;
plot(t,theta,'r');title('\theta, [rad]');
% Corriente, ia
subplot(5,1,3);hold on;
plot(t,ia,'r');title('i_a, [A]');
% Tensión, Va
subplot(5,1,4);hold on;
plot(t,um,'r');title('Accion de control, [V]');
%Torque
subplot(5,1,5);hold on;
plot(t,TL,'k');title('Torque, [N.m]');

xlabel('Tiempo [Seg.]');





