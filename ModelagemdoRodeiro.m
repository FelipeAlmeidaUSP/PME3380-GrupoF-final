clear;
clc;

#parametros
m = 1751;
b0 = 0.7175;
e0 = 0.61;
lambda = 0.1; #conicidade
r0 = 0.4572;
theta = 800; #momento de inercia
kx = 4.5687e+6;
ky = 100;
cx = 20;
cy = 27.6;


#matrizes
[M] = [m, 0; 0, theta];
[C] = [2*cy, 0; 0, 2*cx*(b0^2)];
[K] = [ky, -2*cy; 2*cx*b0*lambda/r0, kx*(exp(2))];


#preparacao pra plotar
i = 1;
n = 20;
parte_real = zeros(1, 11);
parte_imag = zeros(1, 11);
parte_real2 = zeros(1, 11);
parte_imag2 = zeros(1, 11);
parte_real3 = zeros(1, 11);
parte_imag3 = zeros(1, 11);
parte_real4 = zeros(1, 11);
parte_imag4 = zeros(1, 11);

#calculo dos autovalores para diferentes valores de V0
while i<=n
  V0(i) = 2 +4*i;  #velocidade do trem

  #matrizes da equacao geral
  [A] = [-(inv([M])*[C]*(1/V0(i))), -(inv([M])*[K]); eye(size(M)), zeros(size(M))];
  [B] = [inv([M]); zeros(size(M))];

  #Autovalores de A
  [autovetores, autovalores] = eig([A]);
  autoval = [autovalores(1, 1), autovalores(2, 2), autovalores(3, 3), autovalores(4, 4)];

  parte_real1(i) = real(autoval(1))*1000;
  parte_imag1(i) = imag(autoval(1))*1000;

  parte_real2(i) = real(autoval(2))*1000;
  parte_imag2(i) = imag(autoval(2))*1000;

  parte_real3(i) = real(autoval(3))*1000;
  parte_imag3(i) = imag(autoval(3));

  parte_real4(i) = real(autoval(4))*1000;
  parte_imag4(i) = imag(autoval(4));


  i =i+1;
endwhile



figure()
subplot(2, 2, 1);
scatter(parte_real1, parte_imag1, 50, V0, 'x');
legend(arrayfun(@(val) ['Autovalor'], unique(V0), 'UniformOutput', false));
xlabel('parte real');
ylabel('parte imaginária');
title('Autovalores do primeiro modo de vibrar');
legend('off');
grid on;

subplot(2, 2, 3);
scatter(parte_real2, parte_imag2, 50, V0, '+');
legend(arrayfun(@(val) ['Autovalor'], unique(V0), 'UniformOutput', false));
xlabel('parte real');
ylabel('parte imaginária');
legend('off');
grid on;

subplot(2, 2, 2);
scatter(parte_real3, parte_imag3, 50, V0, 'x');
legend(arrayfun(@(val) ['Autovalor'], unique(V0), 'UniformOutput', false));
xlabel('parte real');
ylabel('parte imaginária');
title('Autovalores do segundo modo de vibrar');
legend('off');
axis([-2.5 0 205 206]);
grid on;

subplot(2, 2, 4);
scatter(parte_real4, parte_imag4, 50, V0, '+');
legend(arrayfun(@(val) ['Autovalor'], unique(V0), 'UniformOutput', false));
xlabel('parte real');
ylabel('parte imaginária');
legend('off');
axis([-2.5 0 -206 -205]);
grid on;


