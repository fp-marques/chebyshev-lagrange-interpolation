function pp = lagrange_interp(x,y,xp)
% LAGRANGE_INTERP  Polinómio interpolador na forma de Lagrange
% pp=lagrange_interp(x,y,xp)
% x=nós de interpolação (vetor coluna, n+1 valores)
% y=valores da função nos nós, y = f(x) (vetor coluna)
% xp=pontos onde se quer avaliar o polinómio (vetor coluna)
% pp=valores do polinómio P_n nos pontos xp (vetor coluna)

nnodes=size(x,1); %nº de nós.
nnodesp=size(xp,1); %nº de pontos de avaliação.

pp=zeros(nnodesp,1); %inicializa vetor coluna com zeros com tamanho de nnodesp.
ii=1:nnodes; %cria vetor linha com índices dos nós.
for k=1:nnodes %inicializa o ciclo que constrói L_k (de L_1 até L_nnodes).
    L(:,k)=ones(nnodesp,1); %criamos uma matriz coluna a coluna onde cada uma tem tamanho nnodesp e com valor 1 em todas as entradas (para permitir a multiplicação).
    ink=setdiff(ii,k); %cria vetor linha retirando o k do vetor linha ii (os j ≠ k da fórmula de L_k).
    for jj=ink %ciclo interno para cálculo dos L_k.
        L(:,k)=L(:,k).*(xp-x(jj))/(x(k)-x(jj)); %produtório para obter as novas entradas de cada L_k.
    end
    pp=pp+L(:,k)*y(k); %soma a contribuição de cada nó (o seu L_k multiplicado pelo valor da função nesse nó); no fim do ciclo, pp contém o valor do polinómio em todos os pontos xp.
end
end