SUBROUTINE EscreveCoeficientesEmTXT
	USE DadosGerais
    USE double_Precision

	IMPLICIT NONE

	INTEGER::i,j

    !Escreve os coeficientes para quando a malha não for gerada, assim poderá ser lida diretamente do TXT
    WRITE(888,1) ((alfa(i,j),i=1,biconcavas),j=1,radiais)
    WRITE(888,*) '       '
    WRITE(888,1) ((beta(i,j),i=1,biconcavas),j=1,radiais)
    WRITE(888,*) '       '
    WRITE(888,1) ((gama(i,j),i=1,biconcavas),j=1,radiais)
    WRITE(888,*) '       '
    WRITE(888,1) ((Jacobiano(i,j),i=1,biconcavas),j=1,radiais)

1   FORMAT(32(E15.4),/)



	END SUBROUTINE EscreveCoeficientesEmTXT