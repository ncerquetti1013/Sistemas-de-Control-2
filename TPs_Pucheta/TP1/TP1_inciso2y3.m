% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr. Ing. Pucheta, Julian
% Alumno: Cerquetti, Narella
% Actividad Práctica N1 - Caso de estudio 1
%   Inciso 3  
%   En el archivo Curvas_Medidas_RLC.xls (datos en la hoja 1 y etiquetas en la hoja 2)
%   encontrarán las series de datos que deberían emplear para deducir los valores de R, L y C 
%   del circuito. Emplear el método de la respuesta al escalón, tomando como salida la tensión 
%   en el capacitor.
%
%   Inciso 4
%   Una vez determinados los parámetros R, L y C, emplear la serie de corriente desde 
%   0.05seg en adelante para validar el resultado.
% 

%%
clear all;
%%
% Lectura de datos 
values   = xlsread('Curvas_Medidas_RLC_2024.xls'); 
time     = values(1:end,1);
iSignal  = values(1:end,2);
vcSignal = values(1:end,3);
uSignal  = values(1:end,4);

% Gráficas de datos importados
% 1) Tensión en el capacitor
figure(1)
%subplot(3,1,1);
%plot(time,vcSignal, 'b');title('Tension Capacitor , V_t');grid on;hold on; 
% 2) Corriente en el circuito
%subplot(3,1,2);
plot(time, iSignal, 'b');title('Corriente , I_t');grid on;hold on;
% 3) Tensión de entrada
%subplot(3,1,3);
%plot(time, uSignal, 'b');title('Tensión de entrada , U_t');grid on;hold on;

%%
% Definición de señal para verificación del método
t   = linspace(0,0.1,1000);
u   = linspace(0,0,1000);
vin = 12;
ii  = 0;

for i=1:1000-1
 ii=ii+1;
  if ii<100
      u(i)=0;
 elseif ii>=100 && ii<=500
      u(i)=vin;
 else
      u(i)=vin*-1;
  end
end

figure(4)
plot(t,u, 'm');title('Tensión de Entrada, u_t');grid on
%%
% MÉTODO DE CHEN
% G(s) = Vc(s)/U(s)

% Elección de las tres muestras para aplicar el método, teniendo en cuenta
% que las mismas deben tomarse sobre el transitorio de la señal y siguiendo
% las condiciones del método (ver artículo)

% 1) Elección de valores sobre el plot de Vc
t1_v = values(121,1); y1_v = values(121,3); 
t2_v = values(141,1); y2_v = values(141,3);
t3_v = values(161,1); y3_v = values(161,3);

figure(5)
plot([t1_v t2_v t3_v],[y1_v,y2_v,y3_v],'+');hold on;

% 2) Cálculo de ganancia en estado estacionario
% Se obtienen desde la gráfica ploteada con los datos del excel una vez que
% la respuesta del sistema se estabiliza
K_v = values(470,3)/vin;     

% 3) Normalización de valores
% Se dividen los valores obtenidos del gráfico por el valor de la entrada
% dado que con el método se busca la respuesta del sistema a una entrada
% escalón unitario
y1_v = y1_v/vin; y2_v = y2_v/vin; y3_v = y3_v/vin;

% 4) Aplicación del método
% 4.1) Se define ki = y(ti)/K - 1 para cada ecuación resultante de los 3
% puntos tomados en el inciso 1
k1_v = (y1_v/K_v) - 1;   
k2_v = (y2_v/K_v) - 1;   
k3_v = (y3_v/K_v) - 1;   

% 4.2) Se despejan los valores de alfa1, alfa2 y beta
% alf1 = [k1*k2 + k3 - sqrt(b)]/[2*(k1^2 + k2)]               -- Chen ec.21
% alf2 = [k1*k2 + k3 + sqrt(b)]/[2*(k1^2 + k2)]               -- Chen ec.22
% beta = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3 -- Chen ec.23
b_v     = 4*(k1_v^3)*k3_v - 3*(k1_v^2)*(k2_v^2) - 4*(k2_v^3) + (k3_v^2) + 6*k1_v*k2_v*k3_v; 
alfa1_v = (k1_v*k2_v + k3_v - sqrt(b_v))/(2*(k1_v^2 + k2_v)); 
alfa2_v = (k1_v*k2_v + k3_v + sqrt(b_v))/(2*(k1_v^2 + k2_v)); 
beta_v  = (k1_v+alfa2_v)/(alfa1_v-alfa2_v);                 % -- Chen ec.20 

% 4.3) Cálculo de las constantes de tiempo de la FT
% T1 = -t1/ln(alfa1)                                          -- Chen ec.25
% T2 = -t1/ln(alfa2)                                          -- Chen ec.25
% T3 = beta*(T1 - T2) + T1                                    -- Chen ec.25
T1_v = -(t1_v - 0.01)/log(alfa1_v);           
T2_v = -(t1_v - 0.01)/log(alfa2_v);           
T3_v =(beta_v*(T1_v - T2_v)) + T1_v;          

% 4.4) Definición de la función de transferencia obtenida con el método
G_v = tf(K_v*[T3_v 1],conv([T1_v 1],[T2_v 1]));

% 5) Gráfica de los resultado obtenidos
% Se compara la señal dada por los datos de la hoja de cálculo con la
% salida del sistema al hacer uso de la FT obtenida con el método de Chen
[y_G_v,t_G_v] = lsim(G_v,u,t);
figure(7)
plot(time,vcSignal, 'k' ); grid on; hold on;
plot(t_G_v,y_G_v,'r'); title('Respuesta al sistema aplicando M. de Chen vs. Señal graficada de la tabal de datos');
legend({'Vc_t de excel','Vc_t método de Chen'},'Location','southeast')

%%
% MÉTODO DE CHEN
% G(s) = I(s)/U(s)
% Se repite el método de Chen pero en esta instancia se calcula FT respecto
% a la corriente del sistema, dado que del desarrollo anterior no se puede
% despejar directamente los valores de R, C y L sin fijar o asumir uno de
% ellos.
% El procedimiento a seguir es el mismo que el del caso anterior

% 1) Elección de valores sobre el plot de Vc
t1_i = values(103,1); y1_i = values(103,2); 
t2_i = values(105,1); y2_i = values(105,2);
t3_i = values(107,1); y3_i = values(107,2);

figure(8)
plot([t1_i t2_i t3_i],[y1_i,y2_i,y3_i],'+');hold on;

% 2) Cálculo de ganancia en estado estacionario
% Se obtienen desde la gráfica ploteada con los datos del excel una vez que
% la respuesta del sistema se estabiliza
K_i = values(495,2)/vin;     

% 3) Normalización de valores
% Se dividen los valores obtenidos del gráfico por el valor de la entrada
% dado que con el método se busca la respuesta del sistema a una entrada
% escalón unitario
y1_i = y1_i/vin; y2_i = y2_i/vin; y3_i = y3_i/vin;

% 4) Aplicación del método
% 4.1) Se define ki = y(ti)/K - 1 para cada ecuación resultante de los 3
% puntos tomados en el inciso 1
k1_i= (y1_i/K_i) - 1;
k2_i= (y2_i/K_i) - 1;
k3_i= (y3_i/K_i) - 1;

% 4.2) Se despejan los valores de alfa1, alfa2 y beta
% alf1 = [k1*k2 + k3 - sqrt(b)]/[2*(k1^2 + k2)]               -- Chen ec.21
% alf2 = [k1*k2 + k3 + sqrt(b)]/[2*(k1^2 + k2)]               -- Chen ec.22
% beta = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3 -- Chen ec.23
b_i     = (4*k1_i^3*k3_i) - (3*k1_i^2*k2_i^2)-(4*k2_i^3) + (k3_i^2 )+ (6*k1_i*k2_i*k3_i);
alf1_i  = ((k1_i*k2_i) + k3_i - sqrt(b_i))/(2*(k1_i^2 + k2_i)); 
alf2_i  = ((k1_i*k2_i) + k3_i + sqrt(b_i))/(2*(k1_i^2+k2_i)); 
beta_i  = (k1_i + alf2_i)/(alf1_i - alf2_i);                % -- Chen ec.20

% 4.3) Cálculo de las constantes de tiempo de la FT
% T1 = -t1/ln(alfa1)                                          -- Chen ec.25
% T2 = -t1/ln(alfa2)                                          -- Chen ec.25
% T3 = beta*(T1 - T2) + T1                                    -- Chen ec.25
T1_i = -(t1_i - 0.01)/log(alf1_i);
T2_i = -(t1_i - 0.01)/log(alf2_i);
T3_i = (beta_i*(T1_i - T2_i)) + T1_i;         

% 4.4) Definición de la función de transferencia obtenida con el método
G_i = tf(K_i*[T3_i 1],conv([T1_i 1],[T2_i 1]));

% 5) Gráfica de los resultado obtenidos
% Se compara la señal dada por los datos de la hoja de cálculo con la
% salida del sistema al hacer uso de la FT obtenida con el método de Chen
[y_G_i,t_G_i] = lsim(G_i,u,t);
figure(9)
plot(time, iSignal, 'k'); grid on; hold on;
plot(t_G_i,y_G_i,'m'); title('Respuesta al sistema aplicando M. de Chen vs. Señal graficada de la tabal de datos');
legend({'I_t de excel','I_t método de Chen'},'Location','southeast')

%%
% CÁLCULO DE LOS VALORES DE CAPACIDAD, INDUCTANCIA Y RESISTENCIA
% G_i(s) = sCap / [LCap*s^2 + RCap*s + 1]
Cap =  G_i.num{1}(2);
L   = (G_i.den{1}(1))/Cap;
R   = (G_i.den{1}(2))/Cap;

%%
% Verifiación de los resultados

% Matrices del espacio de estados
A = [-R/L -1/L; 1/Cap 0];
B = [1/L; 0];
C = [1 0];
D = 0;

%Definicion de la ecuación de estado y de salida (salida de corriente)
G1=ss(A,B,C,D);
[yout,yt]=lsim(G1,(u),t);
figure(6)
plot(yt,yout,'b');grid on; hold on;
plot(values(:,1),values(:,2),'r'); title('Comparación de corriente');
legend({'i(t) aproximada con valores RLC calculados','i(t) de excel'},'Location','southeast')