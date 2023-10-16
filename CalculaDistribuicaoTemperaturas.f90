    subroutine CalculaDistribuicaoTemperaturas
    ! Subrotina principal do cálculo de temperaturas bidimensionais
    ! Malha aumentada para 93x200 e simetria de 1/4
    ! Última revisão: 13/08/2018 - Larissa
   
    

	use dadosGerais
    USE double_Precision
	use DadosTermicos

	implicit none

	integer::i,j,k,KTMAX
	REAL(KIND=DP)::auxT,omega
    
    OPEN(unit=27,file='convergenciaT.txt',status='unknown')
    OPEN(unit=22,file='EscreveDadosTermicos.txt',status='unknown')
    OPEN(unit=666,file='TESTE11122019.txt',status='unknown')
    
    
    CALL LeDadosTermicos
    !calculo das temperaturas
    
    !call LedadosTermicos 
    call EspelhoDadosTermicos
    
	!Armazena em temperaturaOld os vetores lidos na entrada

	TemperaturaOld=Temperatura

	!omega=1.0

	WRITE(*,*) 'Entre com KTMAX, maximo de iteracoes para calculo das temperaturas:'

	READ (*,*) KTMAX
    WRITE(*,*) KTMAX
    
	! Comeca a iteracao
	

	DO k = 1, KTMAX
        WRITE(*,*) 'ESCREVENDO EM 666'
     
        OPEN(unit=666,file='TESTE11122019.txt',status='unknown')
2000    Format(32(F10.3),/)
    write(666,2000) (((temperatura(i,j)),i=1,biconcavas),j=1,radiais)
   
    CLOSE(UNIT=666)
    
    CALL CalculaTemperaturas
        WRITE(*,*) 'Temperatura(1,1)',Temperatura(1,1)
        WRITE(*,*) '  '
        PAUSE
        WRITE(*,*) 'Temperatura(1,50)',Temperatura(1,11)
        WRITE(*,*) '  '
        PAUSE
        WRITE(*,*) 'Temperatura(F1-1,50)',Temperatura(F1-1,11)
        WRITE(*,*) '  '
        PAUSE
        WRITE(*,*) 'Temperatura(F1,50)', Temperatura(F1,11)
        WRITE(*,*) '  '
        PAUSE
        WRITE(*,*) 'Temperatura(F2,50)', Temperatura(F2,11)
        WRITE(*,*) '  '
        PAUSE
        WRITE(*,*) 'Temperatura(F2+1,50)', Temperatura(F2+1,11)
        WRITE(*,*)' '
        PAUSE
        WRITE(*,*) 'Temperatura(biconcavas,50)', Temperatura(biconcavas,11)
        PAUSE
        WRITE(*,*) 'Temperatura(40,1)', Temperatura(8,1)
        PAUSE
         WRITE(*,*) 'Temperatura(40,radiais)', Temperatura(8,radiais)
       PAUSE


	!temperatura=omega*temperatura+(1.0-omega)*temperaturaOld  
    
    !Armazena as diferencas entre iteracoes consecutivas
  
        
    DO i=2,biconcavas-1
	DO j=2,radiais-1

	diferencaT(i,j)=(Temperatura(i,j)-TemperaturaOld(i,j))/TemperaturaOld(i,j)
	
    END DO
    END DO

	 !Calcula o máximo valor absoluto das diferencas
	 auxT=maxval(abs(diferencaT))


	
	 !Critério de parada
	  IF(auxT<=10e-3) THEN

	    WRITE(27,*)'Temperaturas convergiram  em',k,'iteracoes'
      
	   CALL EscreveTemperaturas 
       
       GO TO 10
		END IF

    !Armazena em Temperaturaold os novos valores de temperatura
    Temperaturaold=Temperatura
       
    END DO
    

	IF(k>=KTMAX) THEN

	WRITE(27,*) 'Temperaturas NAO convergiram em' ,k,'iteracoes'
    
       
    CALL EscreveTemperaturas 

	END IF
  
10  CLOSE(unit=27)
    CLOSE(unit=22)
    CLOSE(unit=666)


	end subroutine CalculaDistribuicaoTemperaturas