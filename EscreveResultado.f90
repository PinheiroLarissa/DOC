    SUBROUTINE EscreveResultado

	USE DadosGerais
    USE double_Precision

	IMPLICIT NONE

	INTEGER::i,j

    !Escreve os resultados no resultados.txt, já formatado para ser plotado em Scilab
	WRITE(12,*) 'x=['
    WRITE(12,1) ((x(i,j),i=1,biconcavas),j=1,radiais)
    WRITE(12,*) ']'
    WRITE(12,*) '       '
    
	WRITE(12,*) 'y=['
    WRITE(12,1) ((y(i,j),i=1,biconcavas),j=1,radiais)
	WRITE(12,*) ']'

    WRITE(12,*) 'plot(x,y)'
1   FORMAT(32(E15.4),/)



	END SUBROUTINE EscreveResultado