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
clear all; clc; %close all;
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
%h(1)     = 500;
h(1)     = -500;
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
%p1 = -15+15i;
%p2 = -15-15i;
%p3 = -0.5+0.5i;
%p4 = -0.5-0.5i;
%p5 = -0.00001;   

% Polos obtenidos del LQR
p1 = -15.7714+18.1582i;
p2 = -15.7714-18.1582i;
p3 = -0.0427+0.0492i;
p4 = -0.0427-0.0492i;
p5 = -0.0775+0.0000i;

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
q = [1000 10000000 1 1 .01];
Q = diag(q);
R = 10000000/2;

Klqr   = lqr(Aa, Ba, Q, R);
K_lqr  = Klqr(1:4);
KI_lqr = -Klqr(5);

eig(Aa - Ba*Klqr)

%%
% PASO VIII: DISEÑO DEL OBSERVADOR.

%-------------------------------------------------------------------------%
%---------------------------- Sistema dual -------------------------------%
%-------------------------------------------------------------------------%
Ao = Ac';
Bo = Cc';       
Co = Bc';

%-------------------------------------------------------------------------%
%------------------ Determinación de matrices Qo y Ro --------------------%
%-------------------------------------------------------------------------%

% Definición de las matrices Qo y Ro
do = [1 100 1 100];
Qo = 100*diag(do);
Ro = diag([100000000 100000000]);

% Observador
Ko = dlqr(Ao, Bo, Qo, Ro)';
%%
% PASO VIII: SIMULACIÓN.

%-------------------------------------------------------------------------%
%---------------------------- Inicializaciones ---------------------------%
%-------------------------------------------------------------------------%

% Referencias
%ref = -100;
ref = 100;

% Tiempos
Tsim  = 150;      % tiempo de simulación
dt    = 1e-4;    % tiempo de integración
steps = Tsim/dt; % pasos de simulación

% Vectores
t = 0:dt:(Tsim-dt);
u(1)  = 0;        % acción de control efectiva

% Punto de operación del sistema
xop =[0 0 0 ref]';

% Estados observados
x_obs = [0 0 0 h(1)]';
%alpha_obs(1) = 0;
%phi_obs(1) = 0;
%phip_obs(1) = 0;
%h_obs(1) = 0;

% Integrador
psi = 0;

% Zona muerta
deadZone = .1;


%%
%-------------------------------------------------------------------------%
%------------------------------- Simulación ------------------------------%
%-------------------------------------------------------------------------%

for i=1:steps
    
    y_obs = Cc*(x_obs);
    y_out = Cc*x;
    
    % Error de control
    psi_p = ref - y_out(1); 
    % Integral del error de control
    psi = psi + psi_p*dt;
    
    % Acción de control
    %u(i)= -Kcc*(x_obs)+KIcc*psi;     % Con observador
    u(i) = -Kcc*x + KIcc*psi;        % Sin observador
    
    if abs(u(i))<deadZone
        u(i)=0;
    else
        u(i)=sign(u(i))*(abs(u(i))-deadZone);
    end
    
    % Sistema no lineal
    if(i+1<=steps)
        alpha_p = a*(phi(i) - alpha(i));
        phi_pp  = -w^2*(phi(i) - alpha(i) -b*u(i));
        h_p     = c*alpha(i);
        alpha(i+1) = alpha(i) + dt*alpha_p;
        phi_p(i+1) = phi_p(i) + dt*phi_pp;
        phi(i+1)   = phi(i) + dt*phi_p(i);
        h(i+1)     = h(i) + dt*h_p;  
        
        % actualizacuón de estados
        x = [alpha(i+1) phi(i+1) phi_p(i+1) h(i+1)]';
   end
    
    
    x_obsp = Ac*(x_obs) + Bc*u(i) +  Ko*(y_out - y_obs);
    x_obs = x_obs + dt*x_obsp;
  
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




