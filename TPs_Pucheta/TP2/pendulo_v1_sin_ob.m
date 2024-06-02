% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Cerquetti, Narella
% Tp N° 2 - Caso de estudio 3 -

%   Calcular sistema controlador que haga evolucionar al péndulo en el equilibrio 
%   estable.

%   Objetivo de control
%   Partiendo de una condición inicial nula en el desplazamiento y el ángulo 
%   en pi, hacer que el carro se desplace a 10 metros evitando las oscilaciones 
%   de la masa m, considerando que es una grúa. Una vez que delta=10 modificar 
%   a m a un valor 10 veces mayor y volver al origen evitando oscilaciones.

%   Item 3 - Inciso 1
%   Objetivo: Considerar que sólo puede medirse el desplazamiento y el ángulo.
%
%   Item 3 - Inciso 2
%   Objetivo: Especificar el rango posible para el tiempo de muestreo para 
%   implementar el sistema en un microcontrolador.
%
%   Item 3 - Inciso 3
%   Objetivo: Determinar el efecto de la nolinealidad en la acción de control, descripta en la Fig. 4, y verificar cuál
%   es el máximo valor admisible de ésa no linealidad.

%%
close all; clear all; clc;

%%
% PASO I: DEFINICIÓN DEL MODELO DE ESTADOS PARA EL SISTEMA CONTINUO.

%-------------------------------------------------------------------------%
%------------------------ Parámetros del péndulo -------------------------%
%-------------------------------------------------------------------------%

m  = 0.1;
m_ = m*10;
F  = 0.1;
l  = 1.6;
g  = 9.8;
M  = 1.5;

%-------------------------------------------------------------------------%
%-------------------- Modelado en espacio de estados ---------------------%
%-------------------------------------------------------------------------%

% Variables de estado
% x1 = (delta - desp), x2 = delta_p, x3 = (phi - áng), x4 = phi_p
% Salidas
% y1 = delta, y2 = phi
% Entradas
% u1 = u


% Matriz de estados
% RECORRIDO DE 0 A 10: considera m
Ac = [0       1               0        0;  
      0    -F/M            -m*g/M      0;  
      0       0               0        1;  
      0   -F/(l*M)     -g*(m+M)/(l*M)  0];

% RECORRIDO DE 10 A 0: considera m_
Ac_ = [0       1             0         0;           
       0    -F/M         -m_*g/M       0;              
       0       0             0         1;               
       0   -F/(l*M)  -g*(m_+M)/(l*M)   0];
   
% Matriz de entrada 
Bc = [   0   ;
        1/M  ;
         0   ;
      1/(l*M)]; 
  
% Matriz de salida - desplazamiento y ángulo
% y1 = wr
% y2 = theta
% Se toman estas salidas ya que son las variables de estado que se pueden 
% medir
Cc = [1 0 0 0;
      0 0 1 0];  
  
% Matriz de transmisión directa.  
Dc = [0; 
      0];            
  
%-------------------------------------------------------------------------%
%---------------------- Sistema continuo en VE a LA ----------------------%
%-------------------------------------------------------------------------%
sysCont  = ss(Ac,Bc,Cc,Dc);         % m
sysCont_ = ss(Ac_,Bc,Cc,Dc);        % m_
%%
% PASO II: CONVERSIÓN A TIEMPO DISCRETO.

%-------------------------------------------------------------------------%
%------------------- Cálculo de dinámicas del sistema --------------------%
%-------------------------------------------------------------------------%
val = real(eig(Ac));

% Dinámica rápida
tR   = log(0.95)/val(2); 
tInt = tR/4;            % tiempo de integración máximo 

% Dinámica lenta
tL   = log(0.05)/val(3); 
tsim = tL*5;            % tiempo de simulación mínimo

%-------------------------------------------------------------------------%
%------------------- Definición de tiempo de muestreo --------------------%
%-------------------------------------------------------------------------%
Ts      = 1e-2;         % Tiempo de muestreo

%-------------------------------------------------------------------------%
%--------------------- Obtención del sistema discreto --------------------%
%-------------------------------------------------------------------------%
sysDisc = c2d(sysCont, Ts, 'zoh');
Ad     = sysDisc.a;     % matriz de estados
Bd     = sysDisc.b;     % matriz de entrada 
Cd     = sysDisc.c;     % matriz de salida
Dd     = sysDisc.d;     % matriz de acoplamiento directo

sysDisc_ = c2d(sysCont_, Ts, 'zoh');
Ad_     = sysDisc_.a;     % matriz de estados
Bd_     = sysDisc_.b;     % matriz de entrada 
Cd_     = sysDisc_.c;     % matriz de salida
Dd_     = sysDisc_.d;     % matriz de acoplamiento directo

%%
% PASO III: ANÁLISIS DE ALCANZABILIDAD Y CONTROLABILIDAD.

%-------------------------------------------------------------------------%
%------- Verificación de controlabilidad y alcanzabilidad del SD ---------%
%-------------------------------------------------------------------------%

% PARA m
% Controlabilidad
Mc = [Bd Ad*Bd Ad^2*Bd Ad^3*Bd Ad^4*Bd]; 
controlabilidad = rank(Mc); 

% Alcanzabilidad
Ma = [Bd Ad*Bd Ad^2*Bd Ad^3*Bd]; 
alcanzabilidad = rank(Ma);

if(alcanzabilidad == controlabilidad)
    fprintf('\nSistema (m) controlable')
else
    fprintf('\nSistema (m) NO controlable')
end

% PARA m_
% Controlabilidad
Mc_ = [Bd_ Ad_*Bd_ Ad_^2*Bd_ Ad_^3*Bd_ Ad_^4*Bd_]; 
controlabilidad_ = rank(Mc_); 

% Alcanzabilidad
Ma_ = [Bd_ Ad_*Bd_ Ad_^2*Bd_ Ad_^3*Bd_]; 
alcanzabilidad_ = rank(Ma_);

if(alcanzabilidad_ == controlabilidad_)
    fprintf('\nSistema (m_) controlable')
else
    fprintf('\nSistema (m_) NO controlable')
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
% Para el sistema ampliado se toma la parte de la matriz Cd de interés, es
% decir la que se va a comparar con la referencia. como en este caso la
% referencia es el desplazaminiento, se toma la primera fila de Cd
Aa = [Ad , zeros(4,1) ; -Cd(1,:)*Ad, eye(1)];
Ba = [Bd; -Cd(1,:)*Bd];

Aa_ = [Ad_ , zeros(4,1) ; -Cd_(1,:)*Ad_, eye(1)];
Ba_ = [Bd_; -Cd_(1,:)*Bd_];

%-------------------------------------------------------------------------%
%------- Verificación de controlabilidad y alcanzabilidad del SDa --------%
%-------------------------------------------------------------------------%
% PARA m
%Maa = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba Aa^5*Ba]; 
%Mca = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4*Ba Aa^5*Ba Aa^6]; 

%if(rank(Maa) == rank(Mca))
%    fprintf('\nSistema ampliado (m) controlable')
%else
%    fprintf('\nSistema ampliado (m) NO controlable')
%end

% PARA m_
%Maa_ = [Ba_ Aa_*Ba_ Aa_^2*Ba_ Aa_^3*Ba_ Aa_^4*Ba_ Aa_^5*Ba_]; 
%Mca_ = [Ba_ Aa_*Ba_ Aa_^2*Ba_ Aa_^3*Ba_ Aa_^4*Ba_ Aa_^5*Ba_ Aa_^6]; 

%if(rank(Maa_) == rank(Mca_))
%    fprintf('\nSistema ampliado (m_) controlable')
%else
%    fprintf('\nSistema ampliado (m_) NO controlable')
%end

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

% CÁLCULO LQR1 - Trayectoria de 0 a 10 m
% Definición de matrices Q y R.
d1 = [1 50 500 .1 .0001];         % elementos de la diagonal: (delta, delta_p, phi, phi_p, ei)       
Q1 = diag(d1);
R1 = .5; %0.1; 

% Cálculo del controlador
[Klqr, ~, ~] = dlqr(Aa, Ba, Q1, R1); 
K  = Klqr(1:4);     % Ganancia para los estados originales
Ki = -Klqr(5);      % Ganancia para el integrador

% CÁLCULO LQR2 - Trayectoria de 10 a 0
% Definición de matrices Q y R.
d2 = [10 50 90 .01 .0001];         % elementos de la diagonal: (delta, delta_p, phi, phi_p, ei)       
Q2 = diag(d2);
R2 = .5; %0.01; 

% Cálculo del controlador
[Klqr_, ~, ~] = dlqr(Aa_, Ba_, Q2, R2); 
K_  = Klqr_(1:4);     % Ganancia para los estados originales
Ki_ = -Klqr_(5);      % Ganancia para el integrador


%%
% PASO VII: DISEÑO DEL OBSERVADOR.

%-------------------------------------------------------------------------%
%---------------------------- Sistema dual -------------------------------%
%-------------------------------------------------------------------------%
% PARA m
Ao = Ad';
Bo = Cd';       
Co = Bd';

% PARA m_
Ao_ = Ad_';
Bo_ = Cd_';       
Co_ = Bd_';

%-------------------------------------------------------------------------%
%------------------ Determinación de matrices Qo y Ro --------------------%
%-------------------------------------------------------------------------%

% CONSIDERACIONES PARA LA DEFINICIÓN DE LAS MATRICES Qo Y Ro.
% 1) Matriz Qo.
%       Matriz de covarianza del ruido del proceso. Representa la 
%       incertidumbre en el modelo del sistema. Si los estados son conocidos 
%       con precisión, los valores en Qo deben ser pequeños. Si hay 
%       incertidumbre en los estados, los valores deben ser mayores.
% 2) Matriz Ro.
%       Matriz de covarianza del ruido de la medición. Representa la 
%       incertidumbre en las mediciones de salida. Si las mediciones son 
%       precisas, los valores en Ro deben ser pequeños. Si las mediciones 
%       tienen ruido significativo, los valores deben ser mayores.

% Definición de las matrices Qo y Ro

% PARA m
%do = [.01 .01 .01 .0001]; 
%Qo = diag(do); 
%Ro = diag([10000 100000]);
   
% Cálculo del controlador del sistema dual (m)
%Ko = (dlqr(Ao,Bo,Qo,Ro))';

% PARA m_
% Las ponderaciones de Q y R para el sistema corespondiente a m_ en
% principio debieran ser las mismas dado que no hay información que
% implique que la precisión con las que se conocen los estados ni la
% incertidumbre de las mediciones de la salida cambie.
%do_ = [.01 .01 .01 .0001]; 
%Qo_ = diag(do_); 
%Ro_ = diag([10000 100000]);
   
% Cálculo del controlador del sistema dual (m_)
%Ko_ = (dlqr(Ao_,Bo_,Qo_,Ro_))';
   
%----------------------------------------|   
% REVISAR Y AJUSTAR VALORES DE OBSERVADOR|
%--------------- ------------------------|  
%%
% PASO VIII: VERIFICACIÓN.

%-------------------------------------------------------------------------%
%------------------------- Definición de tiempos -------------------------%
%-------------------------------------------------------------------------%
dT    = 1e-4;            % tiempo de integración de Euler
Tsim  = 20;              % tiempo de simulación
p_max = floor(Tsim/Ts);  % cantidad de puntos a simular.

%-------------------------------------------------------------------------%
%---------------------------- Inicialización -----------------------------%
%-------------------------------------------------------------------------%
% Condiciones iniciales
phiIn  = pi;            % ángulo inicial del péndulo
posRef = 10;            % referencia de posición a donde se quiere desplazar el carro           

% Cond iniciales VE
d(1)     = 0;
d_p(1)   = 0;
phi(1)   = phiIn;
phi_p(1) = 0;
phi_pp(1) = 0;

x   = [d(1) d_p(1) phi(1) phi_p(1)]';     % estado
xop = [0 0 pi 0]';                        % pto de operación

h = Ts/Tsim;
u = [];                                   % vector acc de control (auxiliar)

ref   = 10;                               % posición de referencia inicial
flag  = 0;

ei(1) = 0;
x_hat = [0 0 pi 0]';

% Definición de zona muerta
deathZone = 0.5;

%-------------------------------------------------------------------------%
%------------------------ Definición de vectores -------------------------%
%-------------------------------------------------------------------------%
t = 0:h:Tsim;

%%
%-------------------------------------------------------------------------%
%------------------------------ Simulación -------------------------------%
%-------------------------------------------------------------------------%
K_c  = K;
Ki_c = Ki;
%Ko_c = Ko;
m_c  = m;
A    = Ad;
B    = Bd;

i    = 1;
for ki=1:p_max
    
    % Salida de dos componentes
    y_out   = Cd*x;             % salida del sistema
    %y_out_o = Cd*(x_hat-xop);   % salida del observador
   
    ei(ki+1)= ei(ki)+ref-y_out(1);
    
    %Ley de control
    u1(ki) = -K_c*(x - xop) + Ki_c*ei(ki+1);          % sin observador
    %u1(ki)  = -K_c*(x_hat - xop) + Ki_c*ei(ki+1);     % con observador
    
    % Zona Muerta
    if(abs(u1(ki)) < deathZone)
        u1(ki) = 0;               
    else
        u1(ki) = sign(u1(ki))*(abs(u1(ki)) - deathZone);
    end
    %-----------------------------------------------------
    
    for kii=1:Ts/h
        
        u(i) = u1(ki);
        
        % Cálculo por sistema no lineal
        d_pp       = (1/(M+m_c))*(u(i) - m_c*l*phi_pp*cos(phi(i)) + m_c*l*phi_p(i)^2*sin(phi(i)) - F*d_p(i));
        phi_pp     = (1/l)*(g*sin(phi(i)) - d_pp*cos(phi(i)));
        d_p(i+1)   = d_p(i) + h*d_pp;
        d(i+1)     = d(i) + h*d_p(i);
        phi_p(i+1) = phi_p(i) + h*phi_pp;
        phi(i+1)   = phi(i) + h*phi_p(i);
        
        % Cambia de sistema una vez que alcanza los 10[m]
        if(d(i) >= 9.99)
            if(flag == 0)
                ref  = 0;
                m_c  = m_;
                flag = 1;
                K_c  = K_;
                Ki_c = Ki_;
                %Ko_c = Ko_;
                A    = Ad_;
                B    = Bd_;
            end
        end
        i=i+1;
    end
    
    % Actualización de los estados
    % Estados del sistema
    x     = [d(i-1) d_p(i-1) phi(i-1) phi_p(i-1)]';   
    % Estados estimados por el observador
    %x_hat = A*x_hat + B*u1(ki) + Ko_c*(y_out - y_out_o) + xop;
end

u(i) = u1(ki);

%%
%-------------------------------------------------------------------------%
%------------------------------- Gráficas --------------------------------%
%-------------------------------------------------------------------------%
color = 'b';

figure(1);
subplot(3,2,1); grid on; hold on;
plot(t,phi_p,color,'LineWidth',1.5);grid on; title('Velocidad angular \phi_p');

subplot(3,2,2); grid on; hold on;
plot(t,phi,color,'LineWidth',1.5); title('Ángulo \phi');xlabel('Tiempo');

subplot(3,2,3); grid on; hold on;
plot(t,d,color,'LineWidth',1.5);title('Posición grúa \delta');xlabel('Tiempo');

subplot(3,2,4); grid on; hold on;
plot(t,d_p,color,'LineWidth',1.5);title('Velocidad de grúa \delta_p');

subplot(3,1,3); grid on; hold on;
plot(t,u,color,'LineWidth',1.5);title('Acción de control u');xlabel('Tiempo en Seg.');
 
figure(2);
subplot(2,1,1);grid on; hold on;
plot(phi,phi_p,color,'LineWidth',1.5);
title('Ángulo vs Velocidad angular');
xlabel('Ángulo');ylabel('Velocidad angular');
 
subplot(2,1,2);grid on; hold on;
plot(d,d_p,color,'LineWidth',1.5);
title('Distancia vs velocidad');
xlabel('Distancia');ylabel('Velocidad');
   
   

