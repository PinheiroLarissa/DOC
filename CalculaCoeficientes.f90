    SUBROUTINE CalculaCoeficientes

	USE dadosGerais
    USE double_Precision

	IMPLICIT NONE
    INTEGER::i,j
    
    ALLOCATE(alfa(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(beta(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(gama(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(Jacobiano(1:biconcavas, 1:radiais), STAT = AllocateStatus)

	DO i=1,biconcavas

	    DO j=1,radiais
        
    alfa(i,j)=xeta(i,j)*xeta(i,j)+yeta(i,j)*yeta(i,j)

	beta(i,j)=xksi(i,j)*xeta(i,j)+yksi(i,j)*yeta(i,j)

	gama(i,j)=xksi(i,j)*xksi(i,j)+yksi(i,j)*yksi(i,j)
    
    Jacobiano(i,j)=xksi(i,j)*yeta(i,j)-xeta(i,j)*yksi(i,j)  

        ENDDO
	ENDDO

    end subroutine CalculaCoeficientes
