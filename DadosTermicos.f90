    Module DadosTermicos
    USE double_Precision
    
    IMPLICIT NONE
    
     REAL(KIND=DP) :: Kf, Kc 
     REAL(KIND=DP) :: densidade_de_potencia
     REAL(KIND=DP) :: h_int, h_out, hg
     REAL(KIND=DP) :: T_int, T_out
    
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE:: Temperatura,TemperaturaOld,diferencaT
	REAL(KIND=DP), DIMENSION(:),   ALLOCATABLE:: TemperaturaR
    
    CONTAINS
    !
    SUBROUTINE LedadosTermicos
    
    USE dadosGerais 
    
    IMPLICIT NONE
    INTEGER        :: i,j!,Ifinal
    
    
    ALLOCATE(Temperatura   (1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(TemperaturaOld(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(diferencaT    (1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(TemperaturaR  (0:biconcavas-1), STAT = AllocateStatus)
    OPEN(unit=21,file='dadoscalor.txt',status='old')
    OPEN(unit=211,file='TemperaturaInicial.txt',status='old')


    READ(21,*) Kf, Kc, densidade_de_potencia, h_int, h_out, hg, T_int, T_out
    READ(211,*)((Temperatura(i,j),i=1,biconcavas),j=1,radiais)

        DO j=1,radiais
        Do i=F1,F2
    Temperatura(i,j) =800. ! apenas um chute inicial para atribuir valores à matriz Temperatura
            ENDDO
        ENDDO
                DO j=1,radiais
            Do i=2,F1-1
    Temperatura(i,j) = 350. ! apenas um chute inicial para atribuir valores à matriz Temperatura
            ENDDO
                ENDDO
                
            DO j=1,radiais
            Do i=F2+1,biconcavas
    Temperatura(i,j) =350. ! apenas um chute inicial para atribuir valores à matriz Temperatura
            ENDDO
        ENDDO
        
   ! CLOSE(unit=21)
   ! CLOSE(unit=211)

    
    END SUBROUTINE LedadosTermicos
    

    END MODULE DadosTermicos