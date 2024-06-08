% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumna: Cerquetti, Narella
% Tp N° 2 - Caso de estudio 2 -

%   Para el caso del avión, emplear un tiempo de integración por Euler 
%   adecuado y un tiempo de simulación de 70seg. 
%   Los parámetros son:
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
%   es decir donde los ángulos y el valor de la acción de control en valor 
%   absoluto son menores a la unidad.
%
%   Item 2 - Inciso 2.
%   Objetivo: Asumiendo que no puede medirse el ángulo alfa, pero sí el 
%   ángulo phi y la altura, proponer un esquema que permita lograr el objetivo 
%   de control.
%
%   Item 2 - Inciso 3.
%   Objetivo: Establecer el valor del tiempo de muestreo más adecuado para 
%   implementar el diseño en un sistema micro controlado.
%
%   Item 2 - Inciso 4.
%   Objetivo: Determinar el efecto de la nolinealidad en la acción de control, 
%   descripta en la Fig. 4, y verificar cuál es el máximo valor admisible
%   de la nolinealidad.

%%
%clear all; clc; close all;
%%
% PASO II: DEFINICIÓN DEL MODELO DE ESTADOS PARA EL SISTEMA CONTINUO.

%-------------------------------------------------------------------------%
%------------------------- Parámetros del avión --------------------------%
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
% PASO III: ANÁLISIS CONTROLABILIDAD.

%-------------------------------------------------------------------------%
%---------------- Verificación de controlabilidad del SC -----------------%
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
% PASO IV: DEFINICIÓN DE LA LEY DE CONTROL QUE SE DESEA APLICAR.

%-------------------------------------------------------------------------%
%--------------------- Ley de control con integrador ---------------------%
%-------------------------------------------------------------------------%

% Se define Ley de control con integrador a fin de anular el error en
% estado estable del sistema para una referencia no nula.
%  ___________________ 
% | u = -K*x + Ki*psi |
%  -------------------

%-------------------------------------------------------------------------%
%--------------------------- Sistema ampliado ----------------------------%
%-------------------------------------------------------------------------%
% Aa =  A   0
%      -C   0
%
% Ba = B
%      0

Aa = [Ac  zeros(4,1) ; -Cc(1,:) 0];
Ba = [Bc; 0];

%-------------------------------------------------------------------------%
%------- Verificación de controlabilidad para el sistema ampliado --------%
%-------------------------------------------------------------------------%
Mca = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba]; 

numVEa = numVE(1)+1;

if(numVEa == rank(Mca))
    fprintf('\nSistema ampliado controlable')
else
    fprintf('\nSistema ampliado NO controlable')
end

%%
% PASO V: DEFINICIÓN DE POLOS DESEADOS A LAZO CERRADO.

%-------------------------------------------------------------------------%
%-------------------- Ubicación polos de lazo cerrado --------------------%
%-------------------------------------------------------------------------%

% Polos del sistema original
p1 = -15+15i;
p2 = -15-15i;
p3 = -0.5+0.5i;
p4 = -0.5-0.5i;
% Polo del integrador
p5 = -.000001;     

poles = [p1 p2 p3 p4 p5];

%%
% PASO VI: IMPLEMENTACIÓN DEL CONTROLADOR.

%-------------------------------------------------------------------------%
%------------------ Cálculo de la matriz del controlador -----------------%
%-------------------------------------------------------------------------%

% Transformación a FCC
auto_val = eig(Aa);
c_ai = poly(auto_val);

Wa = [c_ai(5) c_ai(4) c_ai(3) c_ai(2)  1;
      c_ai(4) c_ai(3) c_ai(2)    1     0;
      c_ai(3) c_ai(2)    1       0     0;
      c_ai(2)    1       0       0     0;
      1          0       0       0     0];

T = Mca * Wa;

A_controlable = inv(T)*Aa*T;

alpha_i = poly(poles);

K_fcc = fliplr([alpha_i(2:end) - c_ai(2:end)])*inv(T)

Kcc = K_fcc(1:4);
KIcc = -K_fcc(end);

% Ackerman
K_ack = acker(Aa, Ba, poles);
Kack = K_ack(1:4);
KIack = -K_ack(end);

% LQR
q = [1 1000000 1 1 .1];
Q = diag(q);
R = 1000000;

Klqr   = lqr(Aa, Ba, Q, R);
K_lqr  = Klqr(1:4);
KI_lqr = -Klqr(5);


eig(Aa - Ba*Klqr)
%%
% PASO VII: SIMULACIÓN.

%-------------------------------------------------------------------------%
%---------------------------- Inicializaciones ---------------------------%
%-------------------------------------------------------------------------%

% Referencias
ref = -100;
%ref = 100;

% Tiempos
Tsim  = 70;      % tiempo de simulación
dt    = 1e-4;    % tiempo de integración
steps = Tsim/dt; % pasos de simulación

% Vectores
t = 0:dt:(Tsim-dt);
u(1)  = 0;        % acción de control efectiva
uu(1) = 0;        % acción de control afectada por alinealidad

% Punto de operación del sistema
xop =[0 0 0 ref]';

% Integrador
psi(1) = 0;

% Salida
%y_out(1) = [0; 0];

%%
%-------------------------------------------------------------------------%
%------------------------------- Simulación ------------------------------%
%-------------------------------------------------------------------------%

for i=1:steps
   
   % Estados para el instante actual
   x = [alpha(i); phi(i); phi_p(i); h(i)];
   
   % Salida para el instante actual
   y_out = Cc*x;
   
   % Error de control
   psi_p = ref - y_out(1); %Cc(1,:)*x;
   % Integral del error de control
   psi(i+1) = psi(i) + psi_p*dt;
   
   % Acción de control
   %u(i) = -K_lqr*x + KI_lqr*psi(i+1);          % LQR
   %u(i) = -Kack*x + KIack*psi(i+1);             % Ackerman
   u(i) = -Kcc*x + KIcc*psi(i+1);               % FCC
   
   % SISTEMA LINEAL
   % Actualización de variables
   xp = Ac*x + Bc*u(i);      % ecuación de estados
   x  = x + xp*dt;          % actualización de los estados para el tiempo siguiente
   
   if(i+1<=steps)
       alpha(i+1) = x(1);
       phi(i+1)   = x(2);
       phi_p(i+1) = x(3);
       h(i+1)     = x(4);    
   end
  
    
end    

%%
%-------------------------------------------------------------------------%
%-------------------------------- Gráficas -------------------------------%
%-------------------------------------------------------------------------%

color = 'm';

figure(1);

subplot(3,2,1); grid on; hold on;
plot(t, alpha, color);
title('Ángulo con la horizontal, \alpha');
ylabel('\alpha [rad]');
xlabel('Tiempo [s]');

subplot(3,2,2); grid on; hold on;
plot(t, phi, color);
title('Ángulo de cabeceo, \phi');
ylabel('\phi [rad]');
xlabel('Tiempo [s]');

subplot(3,2,3); grid on; hold on;
plot(t, phi_p, color);
title('Velocidad de cambio del ángulo de cabeceo, \phi_p');
ylabel('\phi_p [rad/s]');
xlabel('Tiempo [s]');

subplot(3,2,4); grid on; hold on;
plot(t, h, color);
title('Altura de vuelo, h');
ylabel('h [m]');
xlabel('Tiempo [s]');

subplot(3,1,3); grid on; hold on;
plot(t, u, color);
title('Acción de control, u_t');
ylabel('u_t [V]');
xlabel('Tiempo [s]');




