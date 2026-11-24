clear all
clc
close all
% Interpolação de Lagrange com nós de Chebyshev
% f(x) = 1/(1+25x^2) (função de Runge) , x em [-1, 1]

% Dados da função c
f=@(x) 1./(1+25*x.^2);
a=-1;
b=1;
nvalores=[2:2:60]';   %graus n do polinómio que vamos testar (sequência 2, 4, 6, ..., 60)
nshow=[4 8 16];       %graus mostrados nos gráficos de f, Pn e erro
ncomp=16;             %grau usado na comparação com o spline
nrunge=10;            %grau usado na comparação com nós igualmente espaçados

% Malha fina onde se calcula o erro 
xp = [a:(b-a)/999:b]';      % 1000 pontos igualmente espaçados em [a,b]
% nnodesp = size(xp,1);       % nº de pontos da malha (1000); não é usado neste script, o lagrange_interp calcula o seu
fp = f(xp);                 % valores exatos da função na malha

% Ciclo nos graus n 
nn = size(nvalores,1);      

% Iniciação de vetores com zeros
Emax = zeros(nn,1);         % erro máximo de Lagrange (Chebyshev)
Smax = zeros(nn,1);         % erro máximo do spline
Qmax = zeros(nn,1);         % erro máximo de Lagrange com nós igualmente espaçados
Maj = zeros(nn,1);          % majorante teórico, equação (5) do relatório
h = zeros(nn,1);            % maior espaçamento entre nós

for in=1:nn
    n=nvalores(in);   % Usam-se duas variáveis porque os graus vão de 2 em 2, mas as posições dos vetores vão de de 1 em 1
    
    % tabela x,y com n+1 nós de Chebyshev
    x=nos_chebyshev(a,b,n);
    y=f(x); % O par (x, y) é a tabela que se vai interpolar

    % polinómio de Lagrange e spline nos pontos da malha
    pp=lagrange_interp(x,y,xp);
    sp=interp1(x,y,xp,'spline');

    % erro em cada ponto da malha e erro máximo
    E = fp-pp;
    Es = fp-sp;
    Emax(in) = max(abs(E));
    Smax(in) = max(abs(Es));

    % para comparar: n+1 nós igualmente espaçados (como na aula, A04 efeito Runge)
    xe=[a:(b-a)/n:b]';
    pe=lagrange_interp(xe,f(xe),xp);
    Qmax(in)=max(abs(fp-pe)); % erro máximo com nós igualmente espaçados, para comparar com Emax (fenómeno de Runge)

    % majorante (5): f(x) = Re[1/(1+5ix)] -> f^(m)(x) = Re[m!(-5i)^m/(1+5ix)^(m+1)]
    % logo |f^(m)(x)| <= m!*5^m (igualdade em x=0 para m par)
    Dmax = factorial(n+1)*5^(n+1);
    Maj(in) = (b-a)^(n+1)/(2^(2*n+1)*factorial(n+1))*Dmax;

    % maior espaçamento entre nós consecutivos
    h(in)=max(abs(diff(x)));
end

% Tabela de resultados
fprintf('%4s  %8s  %18s  %14s  %16s  %26s\n', ...
    'n','h','erro max Lagrange','majorante (5)','erro max spline','erro max igualm. espaçados')
for in=1:nn
    fprintf('%4d  %8.4f  %18.3e  %14.3e  %16.3e  %26.3e\n', ...
        nvalores(in), h(in), Emax(in), Maj(in), Smax(in), Qmax(in))
end

% Gráficos 
azul=[42 120 214]/255;                            % Lagrange (Chebyshev)
laranja=[235 104 52]/255;                         % spline
verde=[27 175 122]/255;                           % majorante
violeta=[74 58 167]/255;                          % Lagrange (nós igualmente espaçados)
azuis=[134 182 239; 42 120 214; 16 66 129]/255;   % vários n (claro -> escuro)


figure('Units','normalized','OuterPosition',[0 0.05 1 0.95])
theme(gcf,'light')
set(gcf,'Color','w')        % fundo branco
% (1) função e polinómios interpoladores
subplot(2,2,1)
plot(xp,fp,'r-','LineWidth',2); hold on
for i=1:3
    x=nos_chebyshev(a,b,nshow(i));
    pp=lagrange_interp(x,f(x),xp);
    plot(xp,pp,'-','Color',azuis(i,:),'LineWidth',1.5)
end
legend('f(x)',['P_{' num2str(nshow(1)) '}(x)'],['P_{' num2str(nshow(2)) '}(x)'], ...
    ['P_{' num2str(nshow(3)) '}(x)'],'Location','northeast')
xlabel('x'); ylabel('y'); title('(a) f(x) e polinómio de Lagrange P_n(x)'); grid on
xlim([a b])
% (2) erro em cada ponto
subplot(2,2,2)
for i=1:3
    x=nos_chebyshev(a,b,nshow(i));
    pp=lagrange_interp(x,f(x),xp);
    semilogy(xp,abs(fp-pp),'-','Color',azuis(i,:),'LineWidth',1.5); hold on
end
legend(['n = ' num2str(nshow(1))],['n = ' num2str(nshow(2))],['n = ' num2str(nshow(3))],'Location','best')
xlabel('x'); ylabel('|E(x)|'); title('(b) erro em cada ponto |f(x) - P_n(x)|'); grid on
xlim([a b])
% (3) comparação com o spline (mesmos nós)
subplot(2,2,3)
x=nos_chebyshev(a,b,ncomp);
y=f(x);
pp=lagrange_interp(x,y,xp);
sp=interp1(x,y,xp,'spline');
semilogy(xp,abs(fp-pp),'-','Color',azul,'LineWidth',1.5); hold on
semilogy(xp,abs(fp-sp),'-','Color',laranja,'LineWidth',1.5)
legend('Lagrange','spline (interp1)','Location','best')
xlabel('x'); ylabel('|E(x)|'); title(['(c) Lagrange vs spline, n = ' num2str(ncomp)]); grid on
xlim([a b])
% (4) erro máximo em função do espaçamento dos nós (escala logarítmica)
subplot(2,2,4)
plot(log10(h),log10(Emax),'o-','Color',azul,'LineWidth',1.5); hold on
plot(log10(h),log10(Smax),'s-','Color',laranja,'LineWidth',1.5)
plot(log10(h),log10(Maj),'--','Color',verde,'LineWidth',1.5)
legend('erro máx. Lagrange','erro máx. spline','majorante (5)','Location','best')
xlabel('log_{10}(h)'); ylabel('log_{10}(erro máximo)'); title('(d) erro máximo vs espaçamento h'); grid on
xlim([min(log10(h)) max(log10(h))])

if ~exist('figuras','dir'), mkdir('figuras'); end
exportgraphics(gcf,fullfile('figuras','c_runge.png'),'Resolution',200)

% figura 2: nós de Chebyshev vs nós igualmente espaçados (fenómeno de Runge)
figure('Units','normalized','OuterPosition',[0 0.3 1 0.55])
theme(gcf,'light')
set(gcf,'Color','w')        % fundo branco
subplot(1,2,1)
x=nos_chebyshev(a,b,nrunge);
xe=[a:(b-a)/nrunge:b]';
plot(xp,fp,'r-','LineWidth',2); hold on
plot(xp,lagrange_interp(xe,f(xe),xp),'-','Color',violeta,'LineWidth',1.5)
plot(xp,lagrange_interp(x,f(x),xp),'-','Color',azul,'LineWidth',1.5)
plot(xe,f(xe),'o','Color',violeta,'MarkerSize',7)
plot(x,f(x),'x','Color',azul,'MarkerSize',9,'LineWidth',1.5)
legend('f(x)','nós igualmente espaçados','nós de Chebyshev','Location','north')
xlabel('x'); ylabel('y'); title(['(a) n = ' num2str(nrunge)]); grid on
xlim([a b])
subplot(1,2,2)
semilogy(nvalores,Qmax,'o-','Color',violeta,'LineWidth',1.5); hold on
semilogy(nvalores,Emax,'o-','Color',azul,'LineWidth',1.5)
legend('nós igualmente espaçados','nós de Chebyshev','Location','northwest')
xlabel('n'); ylabel('erro máximo'); title('(b) erro máximo em função de n'); grid on
xlim([min(nvalores) max(nvalores)])
exportgraphics(gcf,fullfile('figuras','c_chebyshev_vs_igualmente_espacados.png'),'Resolution',200)