  SUBROUTINE LeDadosProblema
    USE dadosGerais
    
    IMPLICIT NONE

	INTEGER      :: i,j
    
    biconcavas = 32
    radiais = 63
    F1 = 6
    F2 = 26
    Ifinal = 32 
    OPEN(unit=111,file='dadosproblema.txt',status='old')
     
    ALLOCATE(X(1:biconcavas, 1:2*radiais), STAT = AllocateStatus)
    ALLOCATE(Y(1:biconcavas, 1:2*radiais), STAT = AllocateStatus)
    READ(111,*) ((x(i,j),i=1,biconcavas),j=1,radiais)
    READ(111,*) ((y(i,j),i=1,biconcavas),j=1,radiais)
   

    
 CLOSE (UNIT=111)
        END SUBROUTINE 