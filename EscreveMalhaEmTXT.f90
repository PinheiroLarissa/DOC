SUBROUTINE EscreveMalhaEmTXT
	USE DadosGerais
    USE double_Precision

	IMPLICIT NONE

	INTEGER::i,j

    !Escreve os resultados no resultados.txt para quando a malha não for gerada, assim poderá ser lida diretamente do TXT
	WRITE(777,*) biconcavas, radiais, F1, F2
    WRITE(777,1) ((x(i,j),i=1,biconcavas),j=1,radiais)
    WRITE(777,*) '       '

    WRITE(777,1) ((y(i,j),i=1,biconcavas),j=1,radiais)

1   FORMAT(32(E15.4),/)



	END SUBROUTINE EscreveMalhaEmTXT