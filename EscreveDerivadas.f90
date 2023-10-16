    SUBROUTINE EscreveDerivadas

	USE DadosGerais
    USE double_Precision
	IMPLICIT NONE

	INTEGER::i,j

	!escreveas as derivadas em relacao a ksi e eta

	WRITE(3,*) 'Derivadas X_Ksi'
	WRITE(3,1) ((xksi(i,j),i=1,biconcavas),j=1,radiais)

	WRITE(3,*) 'Derivadas X_Eta'
	WRITE(3,1) ((xeta(i,j),i=1,biconcavas),j=1,radiais)

	WRITE(3,*) 'Derivadas Y_Ksi'
	WRITE(3,1) ((yksi(i,j),i=1,biconcavas),j=1,radiais)

	WRITE(3,*) 'Derivadas Y_Eta'
	WRITE(3,1) ((yeta(i,j),i=1,biconcavas),j=1,radiais)

1	FORMAT(32(e18.5),1x,/)
    WRITE(*,*) 'ESCREVEU AS DERIVADAS'

	END SUBROUTINE EscreveDerivadas