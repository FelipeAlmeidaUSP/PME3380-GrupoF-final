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


#Velocidade do trem e raio de curvatura
V0 = 40;
Raio = 500; #m

#condição inicial x0
x0 = [0; 0; 0.4; 0];


#discretização do tempo
tsim = 800; #tempo de simulação
n = 200; #num de partições
Tn = tsim/n; #periodo de cada divisao

#preparando para o loop
Tp = zeros(1, n+1);
Tp(1, 1) = 1;
x = zeros(length(x0), n+1);
x(:, 1) = x0;
i = 1;
k=1;
j=8;
resultados_uy = zeros(n+1, j);

while k<=j
  cx = 0.02 * 10*k;
  cy = 0.276 * 10*k;

    #Matrizes
    [M] = [m, 0; 0, theta];
    [C] = [2*cy, 0; 0, 2*cx*(b0^2)];
    [K] = [ky, -2*cy; 2*cx*b0*lambda/r0, kx*(exp(2))];
    #Matrizes da equação geral
    [A] = [-(inv(M)*C*(1/V0)), -(inv(M)*K); eye(size(M)), zeros(size(M))];
    [B] = [inv(M); zeros(size(M))];

    #Autovalores e Autovetores de A
    [autovetores, autovalores] = eig(A);
    R = autovetores;
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
   resultados_uy(:, k) = x(3, :)/100;
   k = k +1;
endwhile

#tempo total
t = Tp*Tn
#para plotar
k = linspace(1, 8, 8);

figure()
surf(k, t, resultados_uy, 'EdgeColor', 'interp');
axis ij
xlabel('Expoente do amortecimento z');
ylabel('Tempo (s)');
zlabel('Deslocamento uy (cm)');
title('Deslocamento uy em função do amortecimento e do tempo');

