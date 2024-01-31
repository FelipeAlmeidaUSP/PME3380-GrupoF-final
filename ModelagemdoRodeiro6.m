clear;
clc;

#Parâmetros
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

#Matrizes
[M] = [m, 0; 0, theta];
[C] = [2*cy, 0; 0, 2*cx*(b0^2)];
[K] = [ky, -2*cy; 2*cx*b0*lambda/r0, kx*(exp(2))];

#Velocidade do trem e raio de curvatura
V0 = 40; #m/s
Raio0 = 10e+40; #m
Raio1 = -1200; #m
Raio2 = 300; #m
Raio3 = 800; #m
Raio4 = -400; #m

#Matrizes da equação geral
[A] = [-(inv(M)*C*(1/V0)), -(inv(M)*K); eye(size(M)), zeros(size(M))];
[B] = [inv(M); zeros(size(M))];

#Autovalores e Autovetores de A
[autovetores, autovalores] = eig(A);
R = autovetores;

#condição inicial x0
x0 = [0; 0; 0; 0];

#discretização do tempo
tsim = 50; #tempo de simulação
n = 200; #num de partições
Tn = tsim/n; #periodo de cada divisao

#preparando para o loop
Tp = zeros(1, n+1);
Tp(1, 1) = 1;
x = zeros(length(x0), n+1);
x(:, 1) = x0;
i = 1;

while i<=20 #primeiro percurso reto
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio0); Mpsi + (2*kx*(b0^2)/Raio0)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile


while i<=50 #primeira curva
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio1); Mpsi + (2*kx*(b0^2)/Raio1)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile


while i<=70 #segundo percurso reto
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio0); Mpsi + (2*kx*(b0^2)/Raio0)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile


while i<=90 #segunda curva
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio2); Mpsi + (2*kx*(b0^2)/Raio2)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile


while i<=120 #terceiro percurso reto
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio0); Mpsi + (2*kx*(b0^2)/Raio0)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile


while i<=140 #terceira curva
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio3); Mpsi + (2*kx*(b0^2)/Raio3)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile


while i<=150 #quarto percurso reto
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio0); Mpsi + (2*kx*(b0^2)/Raio0)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile


while i<=170 #quarta curva
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio4); Mpsi + (2*kx*(b0^2)/Raio4)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile


while i<=200 #quinto percurso reto
  T = Tp(i)
  #função de entrada
  Fy = 0*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio0); Mpsi + (2*kx*(b0^2)/Raio0)];
  u = [[F]];

  #matriz da exponencial no intervalo T
  [V] = (autovalores * T);
  [V] = expm([V]);
  #matriz fundamental no intervalo T
  matfund = real([R] * [V] * inv([R]));

  #solução
  x(:, i+1) = (matfund * x(:, i)) + (inv([A]) * (matfund - eye(size(matfund))) * [B] * u);

  Tp(i+1) = Tp(i) + 1;
  i = i + 1;
endwhile

#tempo total
t = Tp*Tn

figure();
subplot(2, 1, 1);
plot(t, x(3, :)/100);
xlabel('Tempo (s)');
ylabel('Deslocamento (cm)');
title('Deslocamento uy');

subplot(2, 1, 2);
plot(t, x(4, :)*100);
xlabel('Tempo (s)');
ylabel('Deslocamento (rad)');
title('Deslocamento angular');

figure();
subplot(2, 1, 1);
plot(t, x(1, :)/100);
xlabel('Tempo (s)');
ylabel('Velocidade (cm/s)');
title('Velocidade em y');

subplot(2, 1, 2);
plot(t, x(2, :));
xlabel('Tempo (s)');
ylabel('Velocidade (rads/s)');
title('Velocidade angular');

#velocidade de cada roda
w = V0/r0; #velocidade angular
r1 = r0 + (x(3, :)/1000)*lambda; #raio da roda externa
r2 = r0 - (x(3, :)/1000)*lambda; #raio da roda interna
v1 = w*r1;
v2 = w*r2;

figure();
plot(t,v1, 'DisplayName', 'Roda Direita');
hold on
plot(t,v2, 'DisplayName', 'Roda Esquerda');
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');
title('Velocidade de cada roda');
legend('show');

#representação do trajeto percorrido
x1_rode = zeros(1, n+1);
x2_rode = zeros(1, n+1);
theta_rode = zeros(1, n+1);
i = 1;

while i <= n
  T = Tn;

  w = V0/r0;
  r1 = r0 + (x(3, i)/1000)*lambda;
  r2 = r0 - (x(3, i)/1000)*lambda;
  v1 = w * r1;
  v2 = w * r2;

  #deslocamento das rodas para cada T
  delta_x1 = v1 * T;
  delta_x2 = v2 * T;
  delta_theta = (delta_x2 - delta_x1) / (2 * b0);

  #calculando a posicao do rodeiro
  x1_rode(i+1) = x1_rode(i) + delta_x1 * cos(theta_rode(i)) - delta_x2 * sin(theta_rode(i));
  x2_rode(i+1) = x2_rode(i) + delta_x1 * sin(theta_rode(i)) + delta_x2 * cos(theta_rode(i));
  theta_rode(i+1) = theta_rode(i) + delta_theta;

  i = i + 1;
endwhile

#visualizando o trajeto no plano xy
figure();
h = plot(x1_rode, x2_rode, '-o');
set(h, 'MarkerSize', 2);
xlabel('Deslocamento em x (m)');
ylabel('Deslocamento em y (m)');
title('Trajeto no Plano XY do Rodeiro');
grid on;
