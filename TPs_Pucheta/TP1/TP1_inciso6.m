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
clc;clear;close all;
%%
X        = -[0; 0; 0; 0];
t_etapa  = 1e-7;
tF       = 1.0; %0.6
thetaRef = 1;
color_='r';

%%
% PID
% Ley del PID para tiempo discreto
% u(k+1) = u(k) + A*e(k) + B*e(k-1) + C*e(k-2) 

% Constantes del PID
Kp = 10; Ki = 20; Kd = 0; % prueba basada en valores usuales implementados para este tipo de PIDs
%Kp=10;Ki=30;Kd=2;

% Definición de coeficientes para el cálculo de la acción de control
Ts = t_etapa;
A1 = ((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1 = (-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1 = Kd/Ts;
e  = zeros(tF/t_etapa,1); 
u=0;        % salida del PID

ii = 0;
for t=0:t_etapa:tF
    ii = ii+1; k = ii+2;
    TL = 0;
    
    if((0.4<=t)&& (t<0.7))
        TL = 6.0e-3;
    end

    X = modMotor2(t_etapa, X, u, TL);

    e(k) = thetaRef - X(4); %ERROR
    u    = u + A1*e(k) + B1*e(k-1) + C1*e(k-2); %PID

    if u>12         %limito accion de control a +-12
        u=12;
    end
    if u<-12
        u=-12;
    end
    
    x1(ii)=X(1);    %wr
    x2(ii)=X(2);    %wrp
    x3(ii)=X(3);    %ia
    x4(ii)=X(4);    %tita
    acc(ii)=u;      %accion de control
    torque(ii)=TL;
end

%%
% Gráficos

t=0:t_etapa:tF;

figure(1)
subplot(5,1,1);hold on;
plot(t,x4,color_);title('\Theta_t [rad/sec]');hold on;

subplot(5,1,2);hold on;
plot(t,acc,color_);title('Accion de control');hold on;

subplot(5,1,3);hold on;
plot(t,x3,color_);title('Ia [A]');hold on;

subplot(5,1,4);hold on;
plot(t,x1,color_);title('\omega_t');hold on;

subplot(5,1,5);hold on;
plot(t,torque,color_);title('T_L');hold on;