SUBROUTINE TesteMatriz
    !========================================================================================!
    !                                                                                        !
    !   Teste para verificar o "COLUNM WISE" do FORTRAN e o WRITE.                           !
    !   Resultado:  ****meu FORTRAN É COLUNM WISE.****                                       !
    !   Resultado2: ****da forma como está o WRITE, a forma colunm wise fica PRESERVADO.**** !
    !                                                                                        !
    !========================================================================================!
    IMPLICIT NONE
    
    
    INTEGER, DIMENSION(:,:), ALLOCATABLE ::AA
    INTEGER:: a,b,i,j,c
    OPEN(UNIT = 1000, FILE = 'TesteMatriz.txt', STATUS = 'unknown') 
    
    a = 6
    b = 4
    c = 0
    
    ALLOCATE(AA(a,b))
    
    DO i=1,a
        DO J=1,b
            c=c+1
            AA(i,j) = c
        ENDDO
    ENDDO
            
100 FORMAT(6(I))    
    WRITE(1000,100) ((AA(i,j),i=1,a),j=1,b)
    
    CLOSE (UNIT=1000)
    END SUBROUTINE TesteMatriz