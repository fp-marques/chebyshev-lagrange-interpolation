function x = nos_chebyshev(a,b,n)
%Nós de Chebyshev
%a e b são as pontas do intervalo [a,b] onde queremos os nós;
%n é o grau do polinómio que queremos; (vamos ter n+1 nós)
%x é o que a função devolve: o vetor com os nós.

% constantes que não dependem de k, calculadas uma só vez antes do ciclo
C1=(b-a)/2;              
C2=(a+b)/2;              
C3=pi/(2*(n+1));
x=zeros(n+1,1); % vetor coluna para os n+1 nós
for k=0:n % ciclo k=0 até n
    x(k+1)=C1*cos((2*k+1)*C3)+C2; % índice k+1 porque o MATLAB começa em 1
end
x=sort(x); % ordena-se por ordem crescente

end
