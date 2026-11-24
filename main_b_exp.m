clear all
clc
close all
% f(x) calculada tal como esta no enunciado
% atencao: para x grande exp(-x) fica abaixo de eps/2 e o computador faz 1+exp(-x) = 1,
% por isso perto de x = 37 os valores de f saem errados (vemos isso nos resultados)
f=@(x) exp(x).*log(1+exp(-x));
a=0;
b=37;
nvalores=[2:2:60]';       % graus que vamos testar (cada um usa n+1 nos)
nshow=[4 8 16];           % so estes 3 vao para os graficos (a) e (b), com todos nao se via nada
ncomp=16;                 % grau usado no grafico (c) para comparar com o spline
mmax=25;                  % ate esta ordem a derivada_b da valores certos (depois perde precisao)
ntab=[2 4 8 16 24 40 60]';   % graus que aparecem na tabela (os mesmos do relatorio)

%malha fina 
% pontos onde calculamos o polinomio e medimos o erro (1000, como na alinea a)
xp=[a:(b-a)/999:b]';
nnodesp=size(xp,1);
fp=f(xp);                 % valor "exato" da funcao (com o mesmo problema perto de x = 37)
% malha ainda mais fina so para procurar o maximo da derivada (a derivada muda muito perto de x=0)
xd=[a:(b-a)/36999:b]';

nn=size(nvalores,1);
% um valor por cada n, vao sendo preenchidos dentro do ciclo
Emax=zeros(nn,1);         % erro maximo de Lagrange
Smax=zeros(nn,1);         % erro maximo do spline
Maj=zeros(nn,1);          % majorante (5)
Dmax=zeros(nn,1);         % max|f^(n+1)| em [a,b]
h=zeros(nn,1);            % maior distancia entre dois nos seguidos
for in=1:nn
    n=nvalores(in);       % in=posicao onde guardamos, n = grau do polinomio
    % tabela x,y: nos de Chebyshev e f nesses nos
    x=nos_chebyshev(a,b,n);
    y=f(x);
    % Lagrange e spline na malha toda, os dois com a mesma tabela
    pp=lagrange_interp(x,y,xp);
    sp=interp1(x,y,xp,'spline');
    % erro em cada ponto e depois o maior em valor absoluto
    E=fp-pp;
    Es=fp-sp;
    Emax(in)=max(abs(E));
    Smax(in)=max(abs(Es));
    % majorante (5): aqui o max da derivada de ordem n+1 vem da derivada_b
    % so confiamos nela ate a ordem 25, depois fica NaN (nao aparece no grafico)
    if n+1<=mmax
        Dmax(in)=max(abs(derivada_b(xd,n+1)));
    else
        Dmax(in)=NaN;
    end
    Maj(in)=(b-a)^(n+1)/(2^(2*n+1)*factorial(n+1))*Dmax(in);
    % h = maior espacamento entre nos (como estao ordenados o diff ja da positivo)
    h(in)=max(diff(x));
end

% so as linhas de ntab, que sao as que estao no relatorio quando o majorante nao foi calculado (n>24) aparece um traco
fprintf('   n      h     max|f^(n+1)|  erro max Lagr.  majorante (5)  erro max spline\n')
for in=1:nn
    if ismember(nvalores(in),ntab)
        if isnan(Maj(in))
            fprintf('%4d  %6.3f        -          %9.2e          -           %9.2e\n', ...
                nvalores(in), h(in), Emax(in), Smax(in))
        else
            fprintf('%4d  %6.3f    %9.2e      %9.2e      %9.2e      %9.2e\n', ...
                nvalores(in), h(in), Dmax(in), Emax(in), Maj(in), Smax(in))
        end
    end
end

%o problema do arredondamento a funcao verdadeira cresce de ln2 = 0.693 para 1 e nunca passa de 1 exp(-x) < eps/2 quando x > -log(eps/2), e a partir dai 1+exp(-x) da exatamente 1 e f da 0
fprintf('\n-log(eps/2) = %.2f\n', -log(eps/2))
xs=[30 35 36 36.5 37];
for i=1:5
    fprintf('f(%4.1f) = %.4f\n', xs(i), f(xs(i)))
end

azul=[42 120 214]/255;
laranja=[235 104 52]/255;
verde=[27 175 122]/255;
azuis=[134 182 239; 42 120 214; 16 66 129]/255;

figure('Position',[100 100 1000 750],'DefaultAxesFontSize',11)
% (a) a funcao e os polinomios para os graus de nshow (ve-se o ruido de f perto de x = 37) os polinomios voltam a ser calculados aqui, e mais simples do que guarda-los no ciclo
subplot(2,2,1)
plot(xp,fp,'k-','LineWidth',2); hold on
for i=1:3
    x=nos_chebyshev(a,b,nshow(i));
    pp=lagrange_interp(x,f(x),xp);
    plot(xp,pp,'-','Color',azuis(i,:),'LineWidth',1.5)
end
legend('f(x)',['P_{' num2str(nshow(1)) '}(x)'],['P_{' num2str(nshow(2)) '}(x)'], ...
    ['P_{' num2str(nshow(3)) '}(x)'],'Location','best')
xlabel('x'); ylabel('y'); title('(a) f(x) e P_n(x)'); grid on
% (b) erro em cada ponto, em escala log (os bicos para baixo sao os nos)
subplot(2,2,2)
for i=1:3
    x=nos_chebyshev(a,b,nshow(i));
    pp=lagrange_interp(x,f(x),xp);
    semilogy(xp,abs(fp-pp),'-','Color',azuis(i,:),'LineWidth',1.5); hold on
end
legend(['n = ' num2str(nshow(1))],['n = ' num2str(nshow(2))],['n = ' num2str(nshow(3))],'Location','best')
xlabel('x'); ylabel('|E(x)|'); title('(b) erro em cada ponto'); grid on
% (c) Lagrange contra spline com exatamente os mesmos nos e valores
subplot(2,2,3)
x=nos_chebyshev(a,b,ncomp);
y=f(x);
pp=lagrange_interp(x,y,xp);
sp=interp1(x,y,xp,'spline');
semilogy(xp,abs(fp-pp),'-','Color',azul,'LineWidth',1.5); hold on
semilogy(xp,abs(fp-sp),'-','Color',laranja,'LineWidth',1.5)
legend('Lagrange','spline','Location','best')
xlabel('x'); ylabel('|E(x)|'); title(['(c) Lagrange vs spline, n = ' num2str(ncomp)]); grid on
% (d) erro maximo em funcao de h, log10 nos dois eixos o majorante so tem pontos ate n = 24 (os NaN nao sao desenhados)
subplot(2,2,4)
plot(log10(h),log10(Emax),'o-','Color',azul,'LineWidth',1.5); hold on
plot(log10(h),log10(Smax),'s-','Color',laranja,'LineWidth',1.5)
plot(log10(h),log10(Maj),'--','Color',verde,'LineWidth',1.5)
legend('erro máx. Lagrange','erro máx. spline',['majorante (5), n \leq ' num2str(mmax-1)],'Location','best')
xlabel('log_{10}(h)'); ylabel('log_{10}(erro máximo)'); title('(d) erro máximo vs h'); grid on

% guarda a figura na pasta figuras (cria a pasta se ainda nao existir)
if ~exist('figuras','dir'), mkdir('figuras'); end
exportgraphics(gcf,fullfile('figuras','b_exp.png'),'Resolution',200)
