% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr. Ing. Pucheta, Julian
% Alumno: Cerquetti, Narella
% Actividad Pr�ctica N1 - Caso de estudio 1 
%   Inciso 1  
%   Asignar valores a R=47ohm, L=1uHy, y C=100nF. Obtener simulaciones que permitan 
%   estudiar la din�mica del sistema, con una entrada de tensi�n escal�n de 12V, que cada 
%   1ms cambia de signo.
%   
%%
clear all;
%%
% Representaci�n del sistema en variables de estado
% [xp] = A*[x] + B*[u]
% [y]  = C*[x] + D*[u] 

% Par�metros de dise�o
R = 47; L = 1e-6; Cap = 100e-9; vin = 12;   

% Punto de operaci�n
% Es necesario definir el punto de operaci�n para encontrar las ecuaciones
% de estado
Il(1) = 0;
Vcl(1)= 0;
y(1)  = 0;
Xop   = [0 0]';          % punto de operaci�n
x     = [Il(1) Vcl(1)]'; % variables de estado

% Espacio de estados
A=[-R/L -1/L; 1/Cap 0];   % matriz de estados
B=[1/L; 0];               % matriz de entrada
C=[R 0];                  % matriz de salida
D=[0];                    % matriz de transmision directa
%%
% Obtenci�n de FT a partir del espacio de estados para el estudio del
% comportamiento (din�mica del sistema)

% C�lculo del numerador y denominador de FT a partir de las ecs. de estado
[numG, denG] = ss2tf(A,B,C,D);
% Funci�n de transferencia del sistema
G = tf(numG, denG);
% Polos de la FT
poles = roots(denG);
% Tama�o de paso y tiempo de simulaci�n
tR   = log(0.95)/poles(1);    % din�mica r�pida
tint = tR/10;                 % tiempo de integraci�n (10 veces menor)
tsim = 2e-3                   % no se escoje en funci�n de la din�mica 
                              % lenta del sistema dado que el t_est << que 
                              % el tiempo de cambio de la entrada
step = tsim/tint;
t    = linspace(0, tsim, step);
u    = linspace(0, 0, step);    % vector inicial de entrada (todos 0s)

%%
% flag para controlar switcheo de la entrada 
swTime = 0;     
z      = size(t);
z      = z(2)/2;

for i=1:step-1
 swTime = swTime+1;
 if (swTime >= z)
     swTime = 0;
     vin = vin*(-1);
 end
 
 u(i) = vin;
 
 %Variables de sistema lineal 
 xp = A*(x-Xop)+B*u(i); % calcula derivadas de las VE
 x = x + xp*tint;       % actualiza el valor de las VE para el pr�ximo
                        % paso por integraci�n por Euler
 Y = C*x;               % calcula salida
 % Actualizaci�n de las variables de inter�s para los nuevos instantes de
 % tiempo
 y(i+1)=Y(1);           % actualiza valor de la salida
 Il(i+1)=x(1);          % actualiza valor de la corriente en el circuito
 Vcl(i+1)=x(2);         % actualiza valor de la tensi�n en C
end

%%
% Gr�ficas

% Gr�fico de Il, Vc y Vin
figure(1)
%subplot(4,1,1);
%plot(t,Il, 'b' );title('Corriente , i'); grid on; 
%subplot(4,1,2);
%plot(t,Vcl, 'r' );title('Tensi�n Capacitor , v_c');grid on
%subplot(4,1,3); 
%plot(t,u, 'm' );title('Tensi�n de Entrada, u');grid on
%subplot(4,1,4);
plot(t,y, 'm' );title('Tensi�n de Salida, v_r');grid on
