% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumna: Cerquetti, Narella
% Tp N� 2 - Caso de estudio 2 -

%   Para el caso del avi�n, emplear un tiempo de integraci�n por Euler 
%   adecuado y un tiempo de simulaci�n de 70seg. 
%   Los par�metros son:
%       - a=0.07; 
%       - omega=9; 
%       - b=5; 
%       - c=150
%   Hallar un controlador para que los polos de lazo cerrado se ubican en 
%       - mu_i=-15+/-15j 
%       - mu_i=-0.5+/-0.5j, 
%   para referencias de 100 y -100 metros en altura, ambas con alturas 
%   iniciales de -500 y 500.
%
%   Item 2 - Inciso 1.
%   Objetivo: Proponer un controlador en tiempo discreto en variables de estado 
%   para que el proceso evolucione en los rangos de validez del modelo, 
%   es decir donde los �ngulos y el valor de la acci�n de control en valor 
%   absoluto son menores a la unidad.
%
%   Item 2 - Inciso 2.
%   Objetivo: Asumiendo que no puede medirse el �ngulo alfa, pero s� el 
%   �ngulo phi y la altura, proponer un esquema que permita lograr el objetivo 
%   de control.
%
%   Item 2 - Inciso 3.
%   Objetivo: Establecer el valor del tiempo de muestreo m�s adecuado para 
%   implementar el dise�o en un sistema micro controlado.
%
%   Item 2 - Inciso 4.
%   Objetivo: Determinar el efecto de la nolinealidad en la acci�n de control, 
%   descripta en la Fig. 4, y verificar cu�l es el m�ximo valor admisible
%   de la nolinealidad.

%%
%clear all; clc; close all;
%%
% PASO II: DEFINICI�N DEL MODELO DE ESTADOS PARA EL SISTEMA CONTINUO.

%-------------------------------------------------------------------------%
%------------------------- Par�metros del avi�n --------------------------%
%-------------------------------------------------------------------------%
a = 0.07;
b = 5;
c = 150;    % velocidad de vuelo [m/s]
w = 9;      % frecuencia natural

%-------------------------------------------------------------------------%
%-------------------- Modelado en espacio de estados ---------------------%
%-------------------------------------------------------------------------%
% Variables de estado - condiciones iniciales
alpha(1) = 0;
phi(1)   = 0;
phi_p(1) = 0;
h(1)     = 500;
%h(1)     = -500;
x = [alpha(1) phi(1) phi_p(1) h(1)]';

% Matriz de estados
Ac = [-a    a   0 0;    % alpha
       0    0   1 0;    % phi
      w^2 -w^2  0 0;    % phi_p
        c   0   0 0];   % h

% Matriz de entrada
Bc = [0 0 b*w^2 0]';

% Matriz de salida
Cc = [0 0 0 1;
      0 1 0 0];             %salida altura y phi

% Matriz de acoplamiento directo
Dc = 0;

%-------------------------------------------------------------------------%
%---------------------- Sistema continuo en VE a LA ----------------------%
%-------------------------------------------------------------------------%
sysCont  = ss(Ac,Bc,Cc,Dc);     


%%
% PASO III: AN�LISIS CONTROLABILIDAD.

%-------------------------------------------------------------------------%
%---------------- Verificaci�n de controlabilidad del SC -----------------%
%-------------------------------------------------------------------------%

% PARA
% Controlabilidad
Mc = [Bc Ac*Bc Ac^2*Bc Ac^3*Bc]; 
controlabilidad = rank(Mc); 

numVE = size(x);

if(numVE(1) == controlabilidad)
    fprintf('\nSistema controlable')
else
    fprintf('\nSistema NO controlable')
end

%%
% PASO IV: DEFINICI�N DE LA LEY DE CONTROL QUE SE DESEA APLICAR.

%-------------------------------------------------------------------------%
%--------------------- Ley de control con integrador ---------------------%
%-------------------------------------------------------------------------%

% Se define Ley de control con una prealimentaci�n para la referencia
%  ___________________ 
% | u = -K*x + G*ref |
%  -------------------
% Al implementar una prealimentaci�n de referencia, el controlador se
% dise�a como un regulador a 0 y luego se a�ade la prealimentaci�n a la
% acci�n de control

%%
% PASO V: DEFINICI�N DE POLOS DESEADOS A LAZO CERRADO.

%-------------------------------------------------------------------------%
%-------------------- Ubicaci�n polos de lazo cerrado --------------------%
%-------------------------------------------------------------------------%

% Polos del sistema original
p1 = -15+15i;
p2 = -15-15i;
p3 = -0.5+0.5i;
p4 = -0.5-0.5i;    

poles = [p1 p2 p3 p4];

%%
% PASO VI: IMPLEMENTACI�N DEL CONTROLADOR.

%-------------------------------------------------------------------------%
%------------------ C�lculo de la matriz del controlador -----------------%
%-------------------------------------------------------------------------%

% Transformaci�n a FCC
auto_val = eig(Ac);
c_ai = poly(auto_val);

W = [c_ai(4) c_ai(3) c_ai(2)  1;
     c_ai(3) c_ai(2)    1     0;
     c_ai(2)    1       0     0;
        1       0       0     0];

T = Mc * W;

A_controlable = inv(T)*Ac*T;

alpha_i = poly(poles);

K_fcc = fliplr([alpha_i(2:end) - c_ai(2:end)])*inv(T);
G_fcc = -inv(Cc(1,:)*inv(Ac-Bc*K_fcc)*Bc);

% Ackerman
K_ack = acker(Ac, Bc, poles);
G_ack = -inv(Cc(1,:)*inv(Ac-Bc*K_ack)*Bc);

% LQR
q = [1 1000000 1 1];
Q = diag(q);
R = 1000000;

Klqr   = lqr(Ac, Bc, Q, R);
G_lqr  = -inv(Cc(1,:)*inv(Ac-Bc*Klqr)*Bc);

eig(Ac - Bc*Klqr)
%%
% PASO VII: SIMULACI�N.

%-------------------------------------------------------------------------%
%---------------------------- Inicializaciones ---------------------------%
%-------------------------------------------------------------------------%

% Referencias
ref = -100;
%ref = 100;

% Tiempos
Tsim  = 70;      % tiempo de simulaci�n
dt    = 1e-4;    % tiempo de integraci�n
steps = Tsim/dt; % pasos de simulaci�n

% Vectores
t = 0:dt:(Tsim-dt);
u1(1)  = 0;        % acci�n de control efectiva

% Punto de operaci�n del sistema
xop =[0 0 0 ref]';

% Integrador
psi(1) = 0;

% Salida
%y_out(1) = [0; 0];

%%
%-------------------------------------------------------------------------%
%------------------------------- Simulaci�n ------------------------------%
%-------------------------------------------------------------------------%

for i=1:steps
   
   % Estados para el instante actual
   x = [alpha(i); phi(i); phi_p(i); h(i)];
   
   % Salida para el instante actual
   y_out = Cc*x;
   
 
   % Acci�n de control
   %u(i) = -K_fcc*x+G_fcc*ref;       % K por transformaci�n can�nica
   %u(i) = -K_ack*x+G_ack*ref;      % K por Ackerman
   u(i) = -Klqr*x+G_lqr*ref;       % K por lqr
   
   % SISTEMA LINEAL
   % Actualizaci�n de variables
   xp = Ac*x + Bc*u(i);      % ecuaci�n de estados
   x  = x + xp*dt;           % actualizaci�n de los estados para el tiempo siguiente
   
   if(i+1<=steps)
       alpha(i+1) = x(1);
       phi(i+1)   = x(2);
       phi_p(i+1) = x(3);
       h(i+1)     = x(4);    
   end
  
    
end    

%%
%-------------------------------------------------------------------------%
%-------------------------------- Gr�ficas -------------------------------%
%-------------------------------------------------------------------------%

color = 'k';

figure(1);

subplot(3,2,1); grid on; hold on;
plot(t, alpha, color);
title('�ngulo con la horizontal, \alpha');
ylabel('\alpha [rad]');
xlabel('Tiempo [s]');

subplot(3,2,2); grid on; hold on;
plot(t, phi, color);
title('�ngulo de cabeceo, \phi');
ylabel('\phi [rad]');
xlabel('Tiempo [s]');

subplot(3,2,3); grid on; hold on;
plot(t, phi_p, color);
title('Velocidad de cambio del �ngulo de cabeceo, \phi_p');
ylabel('\phi_p [rad/s]');
xlabel('Tiempo [s]');

subplot(3,2,4); grid on; hold on;
plot(t, h, color);
title('Altura de vuelo, h');
ylabel('h [m]');
xlabel('Tiempo [s]');

subplot(3,1,3); grid on; hold on;
plot(t, u, color);
title('Acci�n de control, u_t');
ylabel('u_t [V]');
xlabel('Tiempo [s]');