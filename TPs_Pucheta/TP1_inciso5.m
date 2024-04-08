% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr. Ing. Pucheta, Julian
% Alumno: Cerquetti, Narella
% Actividad Práctica N1 - Caso de estudio 2
%   Inciso 5
%   A partir de las curvas de mediciones de las variables graficadas en la Fig. 1-3, se requiere 
%   obtener el modelo del sistema considerando como entrada un escalón de 12V, como salida 
%   a la velocidad angular y como perturbación al toque TL aplicado. En el archivo 
%   Curvas_Medidas_Motor.xls están las mediciones, en la primer hoja
%   los valores y en la segunda los nombres. Se requiere obtener el modelo dinámico, para 
%   establecer las constantes de la corriente.
%%
clear all;
%% 
% LECTURA DE DATOS
values = xlsread('Curvas_Medidas_Motor_2024');

tt     = values(1:end,1);       % tiempo
W      = values(1:end,2);       % velocidad angular
Ia     = values(1:end,3);       % corriente de armadura
Vin    = values(1:end,4);       % tensión de entrada aplicada al motor 
TL_    = values(1:end,5);       % torque de carga

% Gráficas de señales obtenidas de la plantilla
figure(1)
% Velocidad angular 
subplot(4,1,1);hold on
plot(tt, W, 'b'); title('Velocidad angular , \omega[rad/seg]'); grid on;hold on;
% Tensión de excitación de entrada
subplot(4,1,2)
plot(tt, Vin, 'b'); title('Tensión,[V]'); grid on; hold on;
% Corriente de armadura
subplot(4,1,3)
plot(tt, Ia, 'b'); title('Corriente, I_a[A]'); grid on; hold on;
% Torque
subplot(4,1,4)
plot(tt, TL_, 'b'); title('Torque [N.m]'); grid on; hold on;

%%
% PARÁMETROS DE SIMULACIÓN
tint = 1e-7; 
tF   = 0.6;
Ts   = tint;
E    = 12;
t    = 0:Ts:tF;                         % vector de tiempo para gráfica 
um   = zeros(1, length(t));             % vector de acción sobre el motor

% se genera el vector de torque
TLmax    = 1500e-3;
TL       = zeros(1, length(t));
for i=1:tF/Ts
    t_ = t(i);
    if(t_>=0.0351)
        um(i) = E;
    end
    if (((0.1869<=t_) && (t_<0.3372)) || t_>=0.4872)
        TL(i) = TLmax;
    else
        TL(i) = 0;
    end
end

subplot(2,1,1);
plot(t, um); grid; hold on;
subplot(2,1,2);
plot(t, TL); grid; hold on;
%%
% MÉTODO DE CHEN
% G(s) = wr(s)/U(s)
 
% 1) Elección de valores sobre el plot de Wr
t1_W = values(703,1); y1_W = (values(703,2)); 
t2_W = values(704,1); y2_W = (values(704,2));
t3_W = values(705,1); y3_W = (values(705,2));

%figure(3)
%plot([t1_W t2_W t3_W],[y1_W, y2_W, y3_W],'+');hold on;

% 2) Cálculo de ganancia en estado estacionario
K_W = values(16683,2)/12;     % ganancia en ess

% 3) Normalización de valores de la señal
% Se dividen los valores obtenidos del gráfico por el valor de la entrada
% dado que con el método se busca la respuesta del sistema a una entrada
% escalón unitario
y1_W = y1_W/12; y2_W = y2_W/12; y3_W = y3_W/12;

% 4) Aplicación del método
% 4.1) Se define ki = y(ti)/K - 1 para cada ecuación resultante de los 3
% puntos tomados en el inciso 1 de esta sección
k1_W = (y1_W/K_W) - 1;   
k2_W = (y2_W/K_W) - 1;   
k3_W = (y3_W/K_W) - 1;   

% 4.2) Se despejan los valores de alfa1, alfa2 y beta
% alf1 = [k1*k2 + k3 - sqrt(b)]/[2*(k1^2 + k2)]               -- Chen ec.21
% alf2 = [k1*k2 + k3 + sqrt(b)]/[2*(k1^2 + k2)]               -- Chen ec.22
% beta = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3 -- Chen ec.23
b_W     = 4*(k1_W^3)*k3_W - 3*(k1_W^2)*(k2_W^2) - 4*(k2_W^3) + (k3_W^2) + 6*k1_W*k2_W*k3_W; 
alfa1_W = (k1_W*k2_W + k3_W - sqrt(b_W))/(2*(k1_W^2 + k2_W)); 
alfa2_W = (k1_W*k2_W + k3_W + sqrt(b_W))/(2*(k1_W^2 + k2_W)); 
beta_W  = (k1_W+alfa2_W)/(alfa1_W-alfa2_W);                 % -- Chen ec.20 

% 4.3) Cálculo de las constantes de tiempo de la FT
% T1 = -t1/ln(alfa1)                                          -- Chen ec.25
% T2 = -t1/ln(alfa2)                                          -- Chen ec.25
% T3 = beta*(T1 - T2) + T1                                    -- Chen ec.25
T1_W = -(t1_W - 0.0351)/log(alfa1_W);           
T2_W = -(t1_W - 0.0351)/log(alfa2_W);           
T3_W =(beta_W*(T1_W - T2_W)) + T1_W;      

% 4.4) Definición de la función de transferencia obtenida con el método
G_W = tf(K_W*[T3_W 1],conv([T1_W 1],[T2_W 1]));

% 5) Gráfica de los resultado obtenidos
% Se compara la señal dada por los datos de la hoja de cálculo con la
% salida del sistema al hacer uso de la FT obtenida con el método de Chen
[y_G_W,t_G_W] = lsim(G_W,um,t);
figure(4)
plot(tt, W, 'k' ); grid on; hold on;
plot(t_G_W, y_G_W,'r'); title('Respuesta al sistema aplicando M. de Chen vs. Señal graficada de la tabla de datos');
legend({'W_t de excel','W_t método de Chen'},'Location','southeast')

%% 
% MÉTODO DE CHEN
% G(s) = Ia(s)/U(s)
 
% 1) Elección de valores sobre el plot de Wr
t1_ia = values(703,1); y1_ia = (values(703,3)); 
t2_ia = values(704,1); y2_ia = (values(704,3));
t3_ia = values(705,1); y3_ia = (values(705,3));

%figure(5)
%plot([t1_W t2_W t3_W],[y1_W, y2_W, y3_W],'+');hold on;

% 2) Cálculo de ganancia en estado estacionario
K_ia = values(16572,3)/12;     % ganancia en ess

% 3) Normalización de valores de la señal
% Se dividen los valores obtenidos del gráfico por el valor de la entrada
% dado que con el método se busca la respuesta del sistema a una entrada
% escalón unitario
y1_ia = y1_ia/12; y2_ia = y2_ia/12; y3_ia = y3_ia/12;

% 4) Aplicación del método
% 4.1) Se define ki = y(ti)/K - 1 para cada ecuación resultante de los 3
% puntos tomados en el inciso 1 de esta sección
k1_ia = (y1_ia/K_ia) - 1;   
k2_ia = (y2_ia/K_ia) - 1;   
k3_ia = (y3_ia/K_ia) - 1;   

% 4.2) Se despejan los valores de alfa1, alfa2 y beta
% alf1 = [k1*k2 + k3 - sqrt(b)]/[2*(k1^2 + k2)]               -- Chen ec.21
% alf2 = [k1*k2 + k3 + sqrt(b)]/[2*(k1^2 + k2)]               -- Chen ec.22
% beta = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3 -- Chen ec.23
b_ia     = 4*(k1_ia^3)*k3_ia - 3*(k1_ia^2)*(k2_ia^2) - 4*(k2_ia^3) + (k3_ia^2) + 6*k1_ia*k2_ia*k3_ia; 
alfa1_ia = (k1_ia*k2_ia + k3_ia - sqrt(b_ia))/(2*(k1_ia^2 + k2_ia)); 
alfa2_ia = (k1_ia*k2_ia + k3_ia + sqrt(b_ia))/(2*(k1_ia^2 + k2_ia)); 
beta_ia  = (k1_ia + alfa2_ia)/(alfa1_ia - alfa2_ia);        % -- Chen ec.20 

% 4.3) Cálculo de las constantes de tiempo de la FT
% T1 = -t1/ln(alfa1)                                          -- Chen ec.25
% T2 = -t1/ln(alfa2)                                          -- Chen ec.25
% T3 = beta*(T1 - T2) + T1                                    -- Chen ec.25
T1_ia = -(t1_ia - 0.0351)/log(alfa1_ia);           
T2_ia = -(t1_ia - 0.0351)/log(alfa2_ia);           
T3_ia =(beta_ia*(T1_ia - T2_ia)) + T1_ia;      

% 4.4) Definición de la función de transferencia obtenida con el método
G_ia = tf(K_ia*[T3_ia 1],conv([T1_ia 1],[T2_ia 1]));

% 5) Gráfica de los resultado obtenidos
% Se compara la señal dada por los datos de la hoja de cálculo con la
% salida del sistema al hacer uso de la FT obtenida con el método de Chen
[y_G_ia, t_G_ia] = lsim(G_ia, um, t);
figure(6)
plot(tt, Ia, 'k'); grid on; hold on;
plot(t_G_ia, y_G_ia,'r'); title('Respuesta al sistema aplicando M. de Chen vs. Señal graficada de la tabla de datos');
legend({'Ia_t de excel','IA_t método de Chen'},'Location','southeast')

%%
% CÁLCULO DE LOS PARÁMETROS DEL MOTOR A PARTIR DE G_ia
Ki = G_W.num{1}(3);
J  = G_ia.num{1}(2);
Bm = G_ia.num{1}(3);                         % se lo puede aproximar a 0
La = G_ia.den{1}(1)/J;
Ra = (G_ia.den{1}(2)-La*Bm)/J;
Km = (1 - Ra*Bm)/Ki;
%%
% Archivo .txt para guardar los parámetros del motor
%e0 = {'Ki';'J';'Bm';'La';'Ra';'Km'};
%e1 = [Ki;J;Bm;La;Ra;Km];
%tb = table(e0,e1);
%writetable(tb, 'motorParam.txt');
%type motorParam.txt

%%
% SIMULACIÓN DEL FUNCIONAMIENTO DEL MOTOR A PARTIR DE LOS PARÁMETROS 
% OBTENIDOS

% Modelado en el espacio de estados
A = [-Ra/La -Km/La 0 ; Ki/J -Bm/J 0 ; 0 1 0];
B = [1/La 0 ; 0 -1/J ; 0 0];
C = [0 1 0];                               
D = [0 0];

% Condiciones iniciales
Xop      = [0 0 0]';                    % punto de operación
wr(1)    = 0;
ia(1)    = 0;
theta(1) = 0;
    
x = [ia(1) wr(1) theta(1)]';            % variables de estado, x = [ia; wr; theta]

%%
% SIMULACIÓN
for ii=1:(length(t)-1)
    % Cálculo de las derivadas
    xp = A*(x - Xop)+B*[um(ii) TL(ii)]';
    x  = x + xp*Ts;
    wr(ii+1)    = C*x;
    ia(ii+1)    = x(1);
    theta(ii+1) = x(3);
end


%%
% GRÁFICAS
%t=0:t_etapa:tF;

figure(7)
% Velocidad angular
subplot(4,1,1);
plot(tt, W, 'g' );title('Velocidad angular , \omega[rad/seg]'); grid on;hold on; 
plot(t,wr,'k');hold on;
legend({'w de excel','w aproximada'},'Location','southeast')
% Tensión de entrada
subplot(4,1,2);
plot(t,um,'r');title('Tension de Entrada'); grid on; hold on;
xlabel('Tiempo [Seg.]');hold on;
% Corriente
subplot(4,1,3);hold on;
plot(tt, Ia, 'g'); title('Corriente, I_a[A]'); grid on; hold on;
plot(t,ia,'k');
xlabel('Tiempo [Seg.]');hold on;
legend({'ia de excel','ia aproximada'},'Location','southeast')
% Posición
subplot(4,1,4);
plot(tt, TL_, 'g'); title('Torque, [N.m/V]'); grid on; hold on;
plot(t, TL, 'k'); grid on; hold on;

xlabel('Tiempo [Seg.]');hold on;

