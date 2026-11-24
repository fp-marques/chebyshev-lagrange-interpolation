clear all
clc
close all
%dados da funcao a) 
f=@(x) cos(x);
a=0;
b=2*pi;
nvalores=[2:2:30]';       % graus que vamos testar (cada um usa n+1 nos)
nshow=[2 4 8];            % so estes 3 vao para o grafico (a), com todos nao se via nada
ncomp=8;                  % grau usado no grafico (c) para comparar com o spline
ntab=[2 4 6 8 18 20 30]'; % graus que aparecem na tabela (os mesmos do relatorio)
% pontos onde calculamos o polinomio e medimos o erro
% e o mesmo xp das aulas mas com 1000 pontos em vez de 100, para apanhar bem o erro maximo
xp=[a:(b-a)/999:b]';
nnodesp=size(xp,1);
fp=f(xp);                 % valor exato do cos nesses pontos, serve de referencia

%ciclo nos graus n 
nn=size(nvalores,1);
% um valor por cada n, vao sendo preenchidos dentro do ciclo
Emax=zeros(nn,1);         % erro maximo de Lagrange
Smax=zeros(nn,1);         % erro maximo do spline
Maj=zeros(nn,1);          % majorante (5)
h=zeros(nn,1);            % maior distancia entre dois nos seguidos
for in=1:nn
    n=nvalores(in);       % in = posicao onde guardamos, n = grau do polinomio
    % tabela x,y: nos de Chebyshev e o cos nesses nos
    x=nos_chebyshev(a,b,n);
    y=f(x);
    % Lagrange e spline na malha toda, os dois com a mesma tabela
    pp=lagrange_interp(x,y,xp);
    sp=interp1(x,y,xp,'spline');
    % erro em cada ponto (pode ser + ou -) e depois o maior em valor absoluto
    E=fp-pp;
    Es=fp-sp;
    Emax(in)=max(abs(E));
    Smax(in)=max(abs(Es));
    % majorante (5): qualquer derivada do cos e +-sin ou +-cos, por isso o max e 1
    Dmax=1;
    Maj(in)=(b-a)^(n+1)/(2^(2*n+1)*factorial(n+1))*Dmax;
    % h = maior espacamento entre nos (como estao ordenados o diff ja da positivo)
    h(in)=max(diff(x));
end
%
%tabela 
% uma linha por grau de ntab, da para ver logo se o erro fica abaixo do majorante (os outros graus so aparecem no grafico (d))
fprintf('   n      h     erro max Lagrange   majorante (5)   erro max spline\n')
for in=1:nn
    if ismember(nvalores(in),ntab)
        fprintf('%4d   %6.3f      %9.2e         %9.2e        %9.2e\n', ...
            nvalores(in), h(in), Emax(in), Maj(in), Smax(in))
    end
end

azul=[42 120 214]/255;
laranja=[235 104 52]/255;
verde=[27 175 122]/255;
azuis=[134 182 239; 42 120 214; 16 66 129]/255;

figure('Position',[100 100 1000 750],'DefaultAxesFontSize',11)
% (a) o cos e os polinomios para os graus de nshow, para ver a olho quao boa e a aproximacao
% os polinomios voltam a ser calculados aqui, e mais simples do que guarda-los no ciclo
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
% (b) erro em cada ponto, em escala log porque os erros dos varios n sao muito diferentes
% os bicos para baixo sao os nos (ai o erro e zero)
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
% (d) erro maximo em funcao de h, log10 nos dois eixos (como no grafico do erro da aula de dif. finitas)
% cada ponto e um valor de n, os da direita tem poucos nos e os da esquerda muitos
subplot(2,2,4)
plot(log10(h),log10(Emax),'o-','Color',azul,'LineWidth',1.5); hold on
plot(log10(h),log10(Smax),'s-','Color',laranja,'LineWidth',1.5)
plot(log10(h),log10(Maj),'--','Color',verde,'LineWidth',1.5)
legend('erro máx. Lagrange','erro máx. spline','majorante teórico (5)','Location','best')
xlabel('log_{10}(h)'); ylabel('log_{10}(erro máximo)'); title('(d) erro máximo vs h'); grid on

% guarda a figura na pasta figuras (cria a pasta se ainda nao existir)
if ~exist('figuras','dir'), mkdir('figuras'); end
exportgraphics(gcf,fullfile('figuras','a_cos.png'),'Resolution',200)
