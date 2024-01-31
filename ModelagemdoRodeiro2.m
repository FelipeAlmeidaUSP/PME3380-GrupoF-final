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
Raio = 500; #m

#Matrizes da equação geral
[A] = [-(inv(M)*C*(1/V0)), -(inv(M)*K); eye(size(M)), zeros(size(M))];
[B] = [inv(M); zeros(size(M))];

#Autovalores e Autovetores de A
[autovetores, autovalores] = eig(A);
R = autovetores;

#condição inicial x0
x0 = [0; 0; 20; 0];

#discretização do tempo
tsim = 80; #tempo de simulação
n = 200; #num de partições
Tn = tsim/n; #periodo de cada divisao

#preparando para o loop
Tp = zeros(1, n+1);
Tp(1, 1) = 1;
x = zeros(length(x0), n+1);
x(:, 1) = x0;
i = 1;

while i<=n
  T = Tp(i)
  #função de entrada
  Fy = 120*sin(30*T); #força oscilatória de entrada senoidal
  Mpsi = 0; #momento de entrada
  [F] = [Fy + (m*(V0^2)/Raio); Mpsi + (2*kx*(b0^2)/Raio)];
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
plot(t,v1, 'DisplayName', 'Roda Externa');
hold on
plot(t,v2, 'DisplayName', 'Roda Interna');
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');
title('Velocidade de cada roda');
legend('show');
