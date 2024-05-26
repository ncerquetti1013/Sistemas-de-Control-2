% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Cerquetti, Narella
% Tp N° 2 - Caso de estudio 1 - 
%   Item 1 - Inciso 1
%   Objetivo: evitar que la tensión de la acción de control supere los
%   24[V]
% 
%%
close all; clear all; clc;

%%
% LECTURA DE DATOS
values = xlsread('Curvas_Medidas_Motor_2024');

tt     = values(1:end,1);       % tiempo
W      = values(1:end,2);       % velocidad angular
Ia     = values(1:end,3);       % corriente de armadura
Vin    = values(1:end,4);       % tensión de entrada aplicada al motor 
TL_    = values(1:end,5);       % torque de carga

%%
% Gráficas de datos medidos
%figure(1)
% Velocidad angular 
%subplot(4,1,1);hold on
%plot(tt, W, 'b'); title('Velocidad angular , \omega[rad/seg]'); grid on;hold on;
% Tensión de excitación de entrada
%subplot(4,1,2)
%plot(tt, Vin, 'b'); title('Tensión,[V]'); grid on; hold on;
% Corriente de armadura
%subplot(4,1,3)
%plot(tt, Ia, 'b'); title('Corriente, I_a[A]'); grid on; hold on;
% Torque
%subplot(4,1,4)
%plot(tt, TL_, 'b'); title('Torque [N.m]'); grid on; hold on;

%%
% PASO I: DEFINICIÓN DEL MODELO DE ESTADOS PARA EL SISTEMA CONTINUO.

%-------------------------------------------------------------------------%
%------------------------ Parámetros del motor CC ------------------------%
%-------------------------------------------------------------------------%
% Estos parámetros son los resultantes de adecuar los valores obtenidos por
% el método de Chen a un proceso con un torque de carga.
Ki = 0.01162;
Jm = 2.0628e-9;
Bm = 0;
La = 7.0274e-4;
Ra = 28.13;
Km = 0.0605;

%-------------------------------------------------------------------------%
%-------------------- Modelado en espacio de estados ---------------------%
%-------------------------------------------------------------------------%

% Variables de estado
% x1 = ia, x2 = wr, x3 = theta
% Salidas
% y1 = wr
% Entradas
% u1 = Ei, u2 = TL

% Matriz de estados
Ac = [-Ra/La -Km/La 0 ; 
       Ki/Jm -Bm/Jm 0 ; 
       0      1     0]; 
   
% Matriz de salida (considerando torque)
Bc = [1/La  0 ; 
      0  -1/Jm; 
      0     0];   
  
% Matriz de salida - velocidad y ángulo
% y1 = wr
% y2 = theta
% Se toman estas salidas ya que son las que, según la consigna, se pueden 
% medir
Cc = [0 1 0; 
      0 0 1];   
  
% Matriz de transmisión directa.  
Dc = [0 0; 
      0 0];                                       
%-------------------------------------------------------------------------%
%---------------------- Sistema continuo en VE a LA ----------------------%
%-------------------------------------------------------------------------%
sysCont = ss(Ac,Bc,Cc,Dc);

%%
% PASO II: CONVERSIÓN A TIEMPO DISCRETO.

%-------------------------------------------------------------------------%
%------------------- Cálculo de dinámicas del sistema --------------------%
%-------------------------------------------------------------------------%
val = real(eig(Ac));

% Dinámica rápida
tR   = log(0.95)/val(2); 
tInt = tR/4;   % tiempo de integración máximo (puede ser más pequeño)

% Dinámica lenta
tL   = log(0.05)/val(3); 
tsim = tL*5;   % tiempo de simulación mínimo

%-------------------------------------------------------------------------%
%------------------- Definición de tiempo de muestreo --------------------%
%-------------------------------------------------------------------------%
Ts      = tInt;              % Tiempo de muestreo

%-------------------------------------------------------------------------%
%--------------------- Obtención del sistema discreto --------------------%
%-------------------------------------------------------------------------%
sysDisc = c2d(sysCont, Ts, 'zoh');
Ad     = sysDisc.a;     % matriz de estados
Bd     = sysDisc.b;     % matriz de entrada considerando control sobre el torque
Bd_aux = Bd(:,1);       % matriz de entrada sin considerar control sobre el torque
Cd     = sysDisc.c;     % matriz de salida
Dd     = sysDisc.d;     % matriz de acoplamiento directo

%%
% PASO III: ANÁLISIS DE ALCANZABILIDAD Y CONTROLABILIDAD.

%-------------------------------------------------------------------------%
%------- Verificación de controlabilidad y alcanzabilidad del SD ---------%
%-------------------------------------------------------------------------%
% Controlabilidad
Mc = [Bd_aux Ad*Bd_aux Ad^2*Bd_aux Ad^3*Bd_aux]; 
controlabilidad = rank(Mc); 

% Alcanzabilidad
Ma = [Bd_aux Ad*Bd_aux Ad^2*Bd_aux]; 
alcanzabilidad = rank(Ma);

if(alcanzabilidad == controlabilidad)
    fprintf('\nSistema controlable')
else
    fprintf('\nSistema NO controlable')
end

%%
% PASO IV: DEFINICIÓN DE LA LEY DE CONTROL QUE SE DESEA APLICAR.

%-------------------------------------------------------------------------%
%--------------------- Ley de control con integrador ---------------------%
%-------------------------------------------------------------------------%

% Se define Ley de control con integrador a fin de anular el error en
% estado estable del sistema para una referencia no nula.
%  ___________________ 
% | u = -K*x + Ki*psi |
%  -------------------

%%
% PASO V: DEFINICIÓN DEL SISTEMA AMPLIADO.

%-------------------------------------------------------------------------%
%--------------------------- Sistema ampliado ----------------------------%
%-------------------------------------------------------------------------%
Aa = [Ad , zeros(3,1) ; -Cd(2,:)*Ad, eye(1)];
Ba = [Bd_aux; -Cd(2,:)*Bd_aux];

%-------------------------------------------------------------------------%
%------- Verificación de controlabilidad y alcanzabilidad del SDa --------%
%-------------------------------------------------------------------------%
Maa = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba]; 
Mca = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4]; 

if(rank(Maa) == rank(Mca))
    fprintf('\nSistema ampliado controlable')
else
    fprintf('\nSistema ampliado NO controlable')
end

%%
% PASO VI: IMPLEMENTACIÓN LQR.

%-------------------------------------------------------------------------%
%------------------- Determinación de matrices Q y R ---------------------%
%-------------------------------------------------------------------------%

% CONSIDERACIONES PARA LA DEFINICIÓN DE LAS MATRICES Q Y R.
% 1) Impacto matriz Q.
%       Penaliza los estados del sistema; aumentar los elementos de la
%       diagonal que se corresponden a estados específicos hace que el
%       controlador trate de mantenerlos cerca de 0.
% 2) Impacto matriz R.
%       Penaliza las entradas de control; aumentar los valores de los
%       elementos de R hace que el controlador use menos energía de
%       control.

% PAUTAS RECOMENDADAS.
% 1) Recomendable inicializar Q y R con matrices identidad.
% 2) Calcular el controlador y simular el sistema a LC.
% 3) Evaluar aspectos de importancia del proceso simulado:
%       - Velocidad de la respuesta.
%       - Amplidud máxima de la acción de control.
%       - Oscilaciones.
% 4) Ajustar Q y R en base a los resultados de simulación:
%       - Respuesta lenta --> Aumentar valores de Q para penalizar más los
%       estados.
%       - Acción de control muy grande --> Aumentar valores de R para
%       penalizar más la magnitud de la acción de control.
% 5) Repetir procedimiento hasta obtener respuesat deseada.

% Definición de matrices Q y R.
d = [100 6.25e-6 0.40528 0.1];       % elementos de la diagonal de Q (ia, wr, theta, psita)
Q = diag(d);
R = 1e4;

% Cálculo del controlador
[Klqr, ~, ~] = dlqr(Aa, Ba, Q, R); 
K  = Klqr(1:3);     % Ganancia para los estados originales
Ki = -Klqr(4);      % Ganancia para el integrador

%%
% PASO VII: VERIFICACIÓN.

%-------------------------------------------------------------------------%
%------------------------- Definición de tiempos -------------------------%
%-------------------------------------------------------------------------%
Tf = 0.2;                 % Tiempo total de simulación
dT = tInt;              % Tiempo de integración
p_max = floor(Tf/dT);   % Cantidad de pasos de simulación
deadZone = 1;           % Alinealidad del actuador

%-------------------------------------------------------------------------%
%------------------------ Definición de vectores -------------------------%
%-------------------------------------------------------------------------%

% Tiempo
t = 0:dT:p_max*dT;

% Variables de estado
ia    = zeros(1, p_max+1);
wr    = zeros(1, p_max+1);
theta = zeros(1, p_max+1);
ei    = zeros(1, p_max+1);

% Torque
TL_v = zeros(1, p_max+1);

% Referencia de entrada
ref = zeros(1, p_max+1);

% Acción de control
u = zeros(1, p_max+1);

% Salida
y_out = zeros(1, p_max+1);
%%
%-------------------------------------------------------------------------%
%---------------------------- Inicialización -----------------------------%
%-------------------------------------------------------------------------%

% Inicialización señal de referencia
thetaRef = pi/2; 
ref(1) = thetaRef;
iCounter = 0;
tSwitch  = 5;         % tiempo en el que cambia la referencia en [s]
for i=1:p_max
    iCounter = iCounter + dT;
    if(iCounter > tSwitch)
        thetaRef = -1*thetaRef;
        iCounter = 0;
    end
    ref(i) = thetaRef;
end

%plot(t,ref,'g','LineWidth',1.5);

% Inicialización señal de torque
TLref = 1e-3; 
TL_v(1) = TLref;
iCounter = 0;
tSwitch  = 0.3;         % tiempo en el que cambia la referencia en [s]
for i=1:p_max
    iCounter = iCounter + dT;
    if(iCounter > tSwitch)
        TLref = -1*TLref;
        iCounter = 0;
    end
    TL_v(i) = TLref;
end
TL_v(TL_v<0) = 0;

%plot(t,TL_v,'g','LineWidth',1.5);

% Inicialización de estados
x = [ia(1) wr(1) theta(1)]';        % vector de estados
ei = 0;                             % error de integración

%%
%-------------------------------------------------------------------------%
%------------------------------ Simulación -------------------------------%
%-------------------------------------------------------------------------%
for i=1:p_max
    
    eip = ref(i) - Cd(2,:)*x;       % error 
    ei = ei + eip;                  % salida: en TD no se integra, es un retardo
    
    % Acción de control
    u(i) = -K*x + Ki*ei;             % sin observador
    
    % Alinealidad del actuador
    if(abs(u) < deadZone)
        u(i) = 0;
    else
        u(i) = sign(u(i))*(abs(u(i)) - deadZone );
    end
    
    iap = -(Ra/La)*ia(i) -(Km/La)*wr(i) + (1/La)*u(i);
    ia(i+1) = ia(i) + iap*dT;
    wrp = (Ki/Jm)*ia(i) - (Bm/Jm)*wr(i) - (1/Jm)*TL_v(i);
    wr(i+1) = wr(i) + wrp*dT;
    theta(i+1) = theta(i) + wr(i)*dT;
    
    x = [ia(i+1) wr(i+1) theta(i+1)]';
end

%%
%-------------------------------------------------------------------------%
%------------------------------- Gráficas --------------------------------%
%-------------------------------------------------------------------------%
color='b';
figure(2)
subplot(3,2,1);hold on;
plot(t,ref,'g','LineWidth',1.5);title('Angulo  \phi [rad]'); grid on;hold on;
plot(t,theta, color ,'LineWidth',1.5)
subplot(3,2,2);hold on;
plot(t,wr, color ,'LineWidth',1.5);title('Velocidad Angular  \omega [rad/s]');grid on;hold on; 
subplot(3,2,3);hold on;
plot(t,ia, color ,'LineWidth',1.5);title('Corriente  Ia [A]');grid on;hold on; 
subplot(3,2,4);hold on;
plot(t,u, color ,'LineWidth',1.5);title('Acción de control  V [V]');grid on;hold on;
subplot(3,1,3);hold on;
plot(t,TL_v, color ,'LineWidth',1.5);title('Torque  TL [N/m]');grid on;hold on;

figure(3)
plot(theta, wr, color ,'LineWidth',1.5);title('Plano de Fases'); grid on;hold on;



