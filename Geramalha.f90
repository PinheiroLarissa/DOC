    subroutine Geramalha

	! Programa Gerador de Malha - 23:07:2016 -Alvim
    ! Versão          : Versão 2
    ! Data da revisão : 20/08/2018
    ! Revisor         : Larissa
    ! Observações da última revisão: Transposição dos dados de entrada, J = radiais = 37 (retas radiais) e I = biconcavas = 15 (circ. concêntricos)
    !                                Retirada do "go to 1" do if de convergência
    !                                As derivadas estão sendo todas calculadas na subrotina: CalculaDerivadas_KsiEta
    !                                As subrotinas CalculaX_ksi, CalculaX_Eta, CalculaY_ksi e CalculaY_Eta foram removidas.
    !                                A malha bicôncava com simetria está sendo rodada.
    !                                A malha foi refinada para 93x200
    !                                A Malha está sendo gerada com 1/4 de simetria.

	USE DadosGerais
    USE double_Precision
    !USE DadosTermicos

    IMPLICIT NONE

	INTEGER      :: k,l,i,j,kmax
    REAL(KIND=DP):: aux2, aux3, auxt, omega
    ALLOCATE(diferencax(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(diferencay(1:biconcavas, 1:radiais), STAT = AllocateStatus)

    !Abertura dos arquivos

	open(unit=4,  file='coeficientes.txt',  status='unknown')
	open(unit=3,  file='derivadas.txt',     status='unknown')
    open(unit=12, file='resultado.txt',     status='unknown')
	open(unit=2,  file='dadosentrada.txt',  status='unknown')
    open(unit=26, file='convergencias.txt', status='unknown')
    OPEN(UNIT=777, FILE = 'ResultadoMalhaMemoria.txt', status = 'unknown') ! Essa unidade salva o x e o y em um txt para não precisar gerar a malha em todas as rodadas
    OPEN(UNIT=888, FILE = 'ResultadoCoeficientesMemoria.txt', status = 'unknown') ! Essa unidade salva o x e o y em um txt para não precisar gerar a malha em todas as rodadas


    !Escreve pontos de entrada lidos
	call EspelhoEntrada

	!Armazena em xold e yold os vetores lidos n entrada 
	 xold = x
	 yold = y
     
     
    WRITE(12,*) 'pontos iniciais xold'
	WRITE(12,*) '                    '
	WRITE(12,2) ((xold(i,j),i=1,biconcavas),j=1,radiais)
	WRITE(12,*) 'pontos iniciais yold'
	WRITE(12,*) '                    '
	WRITE(12,2) ((yold(i,j),i=1,biconcavas),j=1,radiais)

2	FORMAT(32(E15.4),/)

	write(*,*)'Entre com o valor de kmax (tolerancia para a geracao numerica da malha):'

    read(*,*) kmax

	! Começa a iteração

	do k=1,kmax

	!Calcula as derivadas com xold e yold para X_eta, Y_eta, X_ksi e Y_ksi
    call CalculaDerivadas_KsiEta 

    !Calcula coeficientes das equacoes para x e y
    call CalculaCoeficientes

    !Resolve as equacoes diferenciais numericamente para x e y por gauss-Seidl
    call CalculaCoordenadas

    !Armazena as diferencas entre iteracoes consecutivas

     do i=2,biconcavas-1

	 do j=2,radiais-1

	 if((xold(i,j).ne.0.0).and.(yold(i,j).ne.0.0))then

	 diferencax(i,j)=(x(i,j)-xold(i,j))/xold(i,j)
	
 	 diferencay(i,j)=(y(i,j)-yold(i,j))/yold(i,j)

	  endif

	 end do

	 end do

	 !Calcula o maximo valor absoluto das diferencas
	
	 aux2=maxval(abs(diferencax))
	!WRITE(*,*) 'aux2=',aux2
	 aux3=maxval(abs(diferencay))
     !	WRITE(*,*) 'aux3=',aux3
!pause

	 
	 !Critério de parada
	 if(aux2<=1E-5 .and. aux3<=1E-5) then
          
     write(26,*)'As cooordenadas convergiram  em ',k,'iteracoes'
     
     !Escreve as derivadas calculadas
	 call EscreveDerivadas
    
	 !Escreve coeficientes calculados
	 call EscreveCoeficientes
             
	 !Escreve resultados de interesse
	  call EscreveResultado 
      CALL EscreveMalhaEmTXT
      CALL EscreveCoeficientesEmTXT
      
      go to 1
	
      endif

     !Armazena em xold  e yold os novos valores de x e y
      omega=1.0
      xold=omega*x+(1.0-omega)*xold
	  !xold=x
      
      yold=omega*y+(1.0-omega)*yold
	  !yold=y

	   end do
    
      if (k>=kmax) then
  
       endif 

	  !Fechamento dos arquivos

1     close(unit=2)

	  close(unit=3)

	  close(unit=4)

      close(unit=12)

	  close(unit=26)
      
      CLOSE(unit=777)

	end subroutine Geramalha