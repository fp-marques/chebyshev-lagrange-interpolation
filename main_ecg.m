clear all
clc
close all
%leitura dos dados (codigo fornecido no enunciado) 
% o csv tem duas colunas com cabecalho: time_s (tempo em s) e ECG_mV (sinal em mV)
% as amostras estao de 0.0111 em 0.0111 s, de 0 ate 7.99 s
ref=readtable('P1_ECG_25pct.csv');
tref=ref.time_s;
yref=ref.ECG_mV;

mvalores=[25 30 60 100 200 600 616 617 650 720]';  % valores de m a testar (os da tabela do relatorio)
% 25 e 30 porque o QRS entra entre os dois, 616 e 617 porque e ai que o Lagrange comeca a falhar
% (testamos m a m entre 610 e 620: com 616 ainda funciona, com 617 ja aparece o primeiro NaN)
mshow=[30 600 720];      % so estes 3 vao para a figura (entra o pico R, ainda funciona, ja falha)

%ciclo nos valores de m
nm=size(mvalores,1);
% um valor por cada m, vao sendo preenchidos dentro do ciclo
Dmax=zeros(nm,1);         % maior diferenca entre Lagrange e o spline
nfalha=zeros(nm,1);       % quantos pontos deram NaN ou Inf no Lagrange
for im=1:nm
    m=mvalores(im);       % im = posicao onde guardamos, m = numero de nos
    % o intervalo vai da 1a amostra ate a m-esima
    t1=tref(1);
    tm=tref(m);
    % m nos de Chebyshev nesse intervalo (m nos -> grau n = m-1)
    % os nos nao calham em cima das amostras, por isso o y vem do spline
    x=nos_chebyshev(t1,tm,m-1);
    y=interp1(tref,yref,x,'spline');
    % malha fina de 1000 pontos no intervalo, como o xp das outras alineas
    tp=[t1:(tm-t1)/999:tm]';
    pp=lagrange_interp(x,y,tp);
    % nao ha funcao exata, por isso a referencia e o proprio spline das amostras
    sp=interp1(tref,yref,tp,'spline');
    % para m grande o Lagrange da NaN/Inf (numeros maiores que o computador aguenta)
    % contamos esses pontos e so usamos os que deram valor finito para a diferenca
    ok=isfinite(pp);
    nfalha(im)=sum(~ok);
    Dmax(im)=max(abs(pp(ok)-sp(ok)));
end

% uma linha por m: intervalo de tempo, diferenca maxima e pontos que falharam
% se houve pontos NaN a diferenca so dava para calcular nos pontos bons,
% por isso nesses casos nao a mostramos (aparece um traco, como no relatorio)
fprintf('   m    intervalo [s]   max|P - spline| [mV]   pontos NaN/Inf\n')
for im=1:nm
    if nfalha(im)==0
        fprintf('%4d   [0, %5.3f]        %9.2e           %4d\n', ...
            mvalores(im), tref(mvalores(im)), Dmax(im), nfalha(im))
    else
        fprintf('%4d   [0, %5.3f]            -               %4d\n', ...
            mvalores(im), tref(mvalores(im)), nfalha(im))
    end
end

azul=[42 120 214]/255;     % Lagrange, a mesma cor das outras alineas

% um grafico por cada m de mshow, numa so linha (como nas outras alineas) voltamos a calcular o polinomio aqui, como nas outras alineas letra um pouco maior porque a figura vai ficar pequena no relatorio
figure('Position',[50 50 1000 330],'DefaultAxesFontSize',11)
for i=1:3
    m=mshow(i);
    t1=tref(1);
    tm=tref(m);
    x=nos_chebyshev(t1,tm,m-1);
    y=interp1(tref,yref,x,'spline');
    tp=[t1:(tm-t1)/1999:tm]';          % 2000 pontos para a curva sair mais lisa
    pp=lagrange_interp(x,y,tp);
    sp=interp1(tref,yref,tp,'spline');  % spline das amostras nos mesmos pontos (e a referencia, como na tabela)
    subplot(1,3,i)
    % linha preta = spline, linha azul = Lagrange por cima
    % o preto e mais grosso para se ver por baixo do azul quando os dois coincidem
    plot(tp,sp,'k-','LineWidth',2); hold on
    plot(tp,pp,'-','Color',azul,'LineWidth',1.2)
    % eixo y fixo em todos para se poder comparar (e se o polinomio disparar nao estraga a escala)
    xlim([t1 tm]); ylim([-0.4 1.0]); grid on
    xlabel('t [s]'); ylabel('ECG [mV]')
    title(['m = ' num2str(m) '  (grau ' num2str(m-1) ')'])
    if i==1
        legend('spline','Lagrange','Location','northwest')   % legenda so no primeiro
    end
end
% guarda a figura na pasta figuras (cria a pasta se ainda nao existir)
if ~exist('figuras','dir'), mkdir('figuras'); end
exportgraphics(gcf,fullfile('figuras','ecg.png'),'Resolution',200)
