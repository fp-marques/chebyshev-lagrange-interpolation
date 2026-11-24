function D = derivada_b(x,m)
% Derivada de ordem m de f(x) = exp(x)*ln(1+exp(-x)) 12-10-2026
% sem cálculo simbólico, a partir de duas relações obtidas à mão:
%   f'(x) = f(x) - s(x) , com s(x) = 1/(1+exp(-x))
%   s'(x) = s(x)*(1 - s(x))
% derivando m vezes: f^(m) = f - (s + s' + s'' + ... + s^(m-1))
% cada derivada de s é um polinómio em s: se s^(k) = P(s) então
%   s^(k+1) = P'(s)*s' = P'(s)*(s - s^2)
% os polinómios guardam-se como vetores de coeficientes (polyval, polyder, conv)
% Nota: verificou-se que o resultado é exato até m = 25, exceto perto de x = 37, onde f calculada como no enunciado dá 0
% perde-se precisão por cancelamento subtrativo (erros de arredondamento)
%
% x : pontos onde se calcula a derivada (vetor coluna)
% m : ordem da derivada
% D : valores de f^(m) nos pontos x
%
s=1./(1+exp(-x));
f=exp(x).*log(1+exp(-x));   % f tal como no enunciado (igual ao main_b_exp)
P=[1 0];                          % s^(0) = s  ->  P(s) = s
soma=zeros(size(x));
for k=0:m-1
    soma=soma+polyval(P,s);       % soma s^(k)
    P=conv(polyder(P),[-1 1 0]);  % s^(k+1) = P'(s)*(s - s^2)
end
D=f-soma;
end
