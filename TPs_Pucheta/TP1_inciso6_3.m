% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr. Ing. Pucheta, Julian
% Alumno: Cerquetti, Narella
% Actividad Práctica N1 - Caso de estudio 2
%   Inciso 6
%   Implementar un PID en tiempo discreto para que el ángulo del motor permanezca en una 
%   referencia de 1radian sometido al torque dado por la hoja de excel (Tip: partir de KP=0,1; Ki=0,01; KD=5).
%   El PID se implementará sobre el motor cuyos parámetros se obtuvieron en
%   el inciso 5
%%
clc; clear all; close all;
%%
% PARÁMETROS DEL MOTOR
% Se importan del archivo generado en el inciso anterios
filename = 'motorParam.txt';
delimiter = ',';
header = 1;

motParams = importdata(filename,delimiter,header);

Ki = motParams.data(1,:);
Jm = motParams.data(2,:);
Bm = motParams.data(3,:);
La = motParams.data(4,:);
Ra = motParams.data(5,:);
Km = motParams.data(6,:);
Ei = 12;                    % tensión aplicada

%%
% MODELADO DEL MOTOR EN EL ESPACIO DE ESTADOS

% Variables de estado
% x1 = ia, x2 = wr, x3 = theta
% Salidas
% y1 = wr
% Entradas
% u1 = Ei, u2 = TL
As = [-Ra/La -Km/La 0 ; Ki/Jm -Bm/Jm 0 ; 0 1 0]; % matriz de estados      
Bs = [1/La 0 ; 0 -1/Jm ; 0 0];                   % matriz de entrada
Cs = [0 1 0];                                    % matriz de salida
Ds = [0 0];                                      % matriz de transmisión directa

% Condiciones iniciales y punto de operación
ia(1)    = 0;
wr(1)    = 0;
theta(1) = 0;
x        = [ia(1) wr(1) theta(1)]';             % vector variables de estado
xop      = [0 0 0]';                            % punto de operación
wr_ref   = 1;                                   % consigna

% Funciones de transferencia del sistema
% Wr(s)/Ei(s)
[num1 den1] = ss2tf(As,Bs,Cs,Ds,1);
G_w_e       = tf(num1,den1);
% Wr(s)/TL(s)
[num2 den2] = ss2tf(As,Bs,Cs,Ds,2);
G_w_tl      = tf(num2,den2);

%%
% PARÁMETROS DE SIMULACIÓN
poles  = pole(G_w_e);                           % se calculan los polos para obtener la dinámica rápida del sistema
tR     = log(0.95)/real(poles(2));              % dinámica rápida
Ts     = tR/2;                                  % se calcula el tiempo de muestreo << tR
tF     = 0.6;                                   % tiempo final de simulación
t      = 0:Ts:tF;                               % vector de tiempo para las gráficas
tL     = zeros(1, length(t));                   % vector de torque
TL     = 1000e-3;                               % acción de torque

for i=1:tF/Ts
    t_ = t(i);
 
    if (((0.1869<=t_) && (t_<0.3372)) || t_>=0.4872)
        tL(i) = TL;
    else
        tL(i) = 0;
    end
end

e = zeros(1,length(t));         % vector de error
u = zeros(1,length(t));         % vector acción de control PID

%%
% PID discreto
% u(i+1) = u(i) + A*e(i) + B*e(i-1) + C*e(i-2)

% Parámetros PID
Kp = 10; Ki = 4; Kd = 0.1;
A  = ((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts); 
B  = (-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts); 
C  = Kd/Ts;

for i=1 : tF/Ts
    e(i) = wr_ref - theta(i);
    
    % Aplicación del PID discreto 
    if i==1
        u = u + A*e(i);
    end
    if i==2
        u = u + A*e(i) + B*e(i-1);
    end
    if i>2
        u = u + A*e(i) + B*e(i-1) + C*e(i-2);
    end
    
    xp = As*x + Bs*[u(i) tL(i)]';
    x = x + xp*Ts;
    ia(i+1) = x(1);
    wr(i+1) = Cs*x;
    theta(i+1) = x(3);
end

%%
% GRÁFICAS

figure(1)
subplot(4,1,1)
plot(t, wr, 'k'); grid; hold on;
title('\omega_r [Rad/s]'); xlabel('Tiempo [s]'); ylabel('\omega_r [Rad/s]'); ylim([-10 10]);
subplot(4,1,2)
plot(t, ia, 'k'); grid; hold on;
title('Corriente [A]'); xlabel('Tiempo [s]'); ylabel('Corriente [A]'); ylim([0 0.1]); grid;
subplot(4,1,3)
plot(t, theta, 'k'); grid; hold on;
title('\Theta [Rad]');xlabel('Tiempo [s]'); ylabel('\Theta [Rad]'); ylim([0 2]); grid;
subplot(4,1,4)
plot(t, tL, 'k'); grid; hold on;
title('Entradas'); xlabel('Tiempo [s]'); ylabel('[N.m][V]'); legend('Torque')





