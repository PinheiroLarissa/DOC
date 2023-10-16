    subroutine CalculaTemperaturas_SEMOLD
    ! Calcula as temperaturas numericamente por dif. finitas no contorno e nos pontos internos.
    ! As condições do contorno interno e externo são convectivas e das simetrias (reflexivas) é gradiente zero.
    ! As quinas são tratadas como pertencendo às paredes reflexivas, interna e externa. (05/10/2019)
    ! Acerto das condições de contorno após transposição do programa.
    ! Aumento da malha para 45x100 bicôncava com simetria de 1/4 (05/10/2019)
    ! O Balanço de energia das quinas foi refeito (05/10/2019)
    ! Última revisão: 09/10/2019 - Larissa
    
    USE DadosGerais
    USE DadosTermicos

	IMPLICIT NONE
	
	INTEGER::i,j
    REAL::AT,BT,CT
    ALLOCATE(Temperatura   (1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(TemperaturaOld(1:biconcavas, 1:radiais), STAT = AllocateStatus)

!Obs.: eq. na página 21 do caderno
DO j=2,radiais-1
!=================================================================================================================================================================================
!----------------------------------------------------------------   QUINA (1,1)   ------------------------------------------------------------------------------------------------
!=================================================================================================================================================================================
    
	!Cálculo da Temperatura na quina Q1 (1,1) - quina entre PAREDE OESTE e SUL - interna                         ****(CHECADA EM 17/10 - AGORA É REFLEXIVA)****
    Temperatura(1,1) = (1.0/(-1. + (beta(1,1)/alfa(1,1)))) * ((-1./3.)*(4*Temperatura(2,1) - Temperatura(3,1)) + beta(1,1)/(3*alfa(1,1))*(4*Temperatura(1,2) - Temperatura(1,3))  ) 

!=================================================================================================================================================================================
!--------------------------------------------------------  PAREDE SUL (REFLEXIVA INFERIOR)  -------------------------------------------------------------------------------------
!=================================================================================================================================================================================
 	DO i=2,biconcavas-1 !Considera o corte da simetria inferior: essa pareds reflexiva é calculada de uma vez, sem repartir em CLAD, GAP e FUEL.

    ! PAREDE OC2 - Reflexiva "sul" (CHECADO)
	Temperatura(i,1) = (1./3.)*(4*Temperatura(i,2) -TemperaturaOld(i,3)) - (beta(i,1)/(3*gama(i,1))) *(Temperatura(i+1,1) - Temperatura(i-1,1))
    
    ENDDO
!=================================================================================================================================================================================
!--------------------------------------------------------   CLAD1 + GAP1: INTERNA CONVECTIVA   -----------------------------------------------------------------------------------
!=================================================================================================================================================================================
	!DO j=2,radiais-1
    ! PAREDE OC1 - Parede CONVECTIVA interna da pastilha (CHECADO)
    Temperatura(1,j)=(1./( (3*Kc*alfa(1,j)/(2*Jacobiano(1,j)*SQRT(alfa(1,j))*h_int)) + 1 ))&
                      *( (Kc/(2*Jacobiano(1,j)*SQRT(alfa(1,j))*h_int))*(alfa(1,j)*(4*TemperaturaOld(2,j) - Temperatura(3,j))&
                      - beta(1,j)*(Temperatura(1,j+1) - Temperatura(1,j-1))) + T_int) 
            
    !END DO
    !WRITE(*,*)'TemperaturaOld interna NA PAREDE INTERNA CONVECTIVA DO CLAD1:', TemperaturaOld(1,50)
    !PAUSE
    
!==================================================================================================================================================================================
!--------------------------------------------------------------------- CLAD 1 -----------------------------------------------------------------------------------------------------
!==================================================================================================================================================================================
!DO j=2, radiais-1 !Pontos Internos (CHECADO)

!---------------------------------------------------------------------COM GAP------------------------------------------------------------------------------------------------------
    IF (dg .NE. 0) THEN   !Com GAP, o F1-1 é a parede do GAP. 
    DO i=2, F1-2
         
    Temperatura(i,j) = (0.5/(alfa(i,j)+gama(i,j)))*( alfa(i,j)*(Temperatura(i+1,j) + Temperatura(i-1,j))&
                 - 0.5*beta(i,j)*(Temperatura(i+1,j+1)-Temperatura(i-1,j+1)-TemperaturaOld(i+1,j-1)+Temperatura(i-1,j-1))&
                 + gama(i,j)*(Temperatura(i,j+1) + Temperatura(i,j-1)))

    END DO
        
    Temperatura(F1-1,j) =(1./( (-3*Kc*alfa(F1-1,j)*ro/(2*Jacobiano(F1-1,j)*SQRT(alfa(F1-1,j))*hg*rfi)) -1 ))&
                      *( (Kc*ro/(2*Jacobiano(F1-1,j)*SQRT(alfa(F1-1,j))*hg*rfi))*(alfa(F1-1,j)*(-4*Temperatura(F1-2,j) + Temperatura(F1-3,j))&
                      - beta(F1-1,j)*(Temperatura(F1-1,j+1) - Temperatura(F1-1,j-1))) - Temperatura(F1,j))
    ELSE !Sem gap, o último ponto interno é o F1-1 (parede OC3)
!---------------------------------------------------------------------SEM GAP----------------------------------------------------------------------------------------------------
   
    DO i=2, F1-1
        
    Temperatura(i,j) = (0.5/(alfa(i,j)+gama(i,j)))*( alfa(i,j)*(Temperatura(i+1,j) + Temperatura(i-1,j))&
                 - 0.5*beta(i,j)*(Temperatura(i+1,j+1)-Temperatura(i-1,j+1)-Temperatura(i+1,j-1)+Temperatura(i-1,j-1))&
                 + gama(i,j)*(Temperatura(i,j+1) + Temperatura(i,j-1)))

    END DO
    END IF

!END DO 

!===============================================================================================================================================================
!---------------------------------------------------------------------- FUEL -----------------------------------------------------------------------------------
!===============================================================================================================================================================
    !DO j=2,radiais-1 !(CHECADO)
        
!----------------------------------------------------------------------COM GAP-------------------------------------------------------------------------

    IF (dg .NE. 0) THEN !F1 e F2 são pontos de fronteira
    Temperatura(F1,j)=(1./( (3*Kf*alfa(F1,j)*rfi/(2*Jacobiano(F1,j)*SQRT(alfa(F1,j))*hg*ro)) + 1 ))&
                      *( (Kf*rfi/(2*Jacobiano(F1,j)*SQRT(alfa(F1,j))*hg*ro))*(alfa(F1,j)*(4*Temperatura(F1+1,j) - Temperatura(F1+2,j))&
                      - beta(F1,j)*(Temperatura(F1,j+1) - Temperatura(F1,j-1))) + Temperatura(F1-1,j))
        
    DO i=F1+1, F2-1
	
    Temperatura(i,j)=(0.5/(alfa(i,j)+gama(i,j)))*( alfa(i,j)*(Temperatura(i+1,j) + Temperatura(i-1,j))&
                        - 0.5*beta(i,j)*(Temperatura(i+1,j+1)-Temperatura(i-1,j+1)-Temperatura(i+1,j-1)+Temperatura(i-1,j-1))&
                        + gama(i,j)*(Temperatura(i,j+1) + Temperatura(i,j-1))+ ( (Jacobiano(i,j))**2 * densidade_de_potencia/Kf))
    END DO
    
    Temperatura(F2,j)= (1./( (3*Kf*alfa(F2,j)*rfo/(2*Jacobiano(F2,j)*SQRT(alfa(F2,j))*hg*rci)) + 1 ))&
                      *( (-Kf*rfo/(2*Jacobiano(F2,j)*SQRT(alfa(F2,j))*hg*rci))*(alfa(F2,j)*(-4*Temperatura(F2-1,j) + Temperatura(F2-2,j))&
                      - beta(F2,j)*(Temperatura(F2,j+1) - Temperatura(F2,j-1))) + Temperatura(F2+1,j))
    
    ELSE
!----------------------------------------------------------------------SEM GAP------------------------------------------------------------------------- 
    DO i=F1, F2
	
    Temperatura(i,j)=(0.5/(alfa(i,j)+gama(i,j)))*( alfa(i,j)*(Temperatura(i+1,j) + Temperatura(i-1,j))&
                        - 0.5*beta(i,j)*(Temperatura(i+1,j+1)-Temperatura(i-1,j+1)-Temperatura(i+1,j-1)+Temperatura(i-1,j-1))&
                        + gama(i,j)*(Temperatura(i,j+1) + Temperatura(i,j-1))+ ( (Jacobiano(i,j))**2*densidade_de_potencia/Kf))
    END DO
    ENDIF
    
    
!END DO

!===============================================================================================================================================================
!---------------------------------------------------------------------- CLAD 2 ---------------------------------------------------------------------------------
!===============================================================================================================================================================
!DO j=2, radiais-1
    
!-----------------------------------------------------------------------COM GAP---------------------------------------------------------------------------------
   IF (dg .NE. 0) THEN ! Sem gap, o primeiro ponto interno do GAP é o F2+1 (parede LC1), com GAP, o F2+1 é a parede do GAP
        
        ! PAREDE LC3 - Parede entre CLAD2GAP2 ---> Parede do GAP2 SE **HOUVER** GAP (CHECADO)
        Temperatura(F2+1,j) = (1./( (3*Kc*alfa(F2+1,j)*rci/(2*Jacobiano(F2+1,j)*SQRT(alfa(F2+1,j))*hg*rfo)) +1 ))&
                      *( (Kc*rci/(2*Jacobiano(F2+1,j)*SQRT(alfa(F2+1,j))*hg*rfo))*(alfa(F2+1,j)*(4*Temperatura(F2+2,j) - Temperatura(F2+3,j))&
                      - beta(F2+1,j)*(Temperatura(F2+1,j+1) - Temperatura(F2+1,j-1))) + Temperatura(F2,j))
         
        DO i = F2+2, biconcavas - 1
         
        Temperatura(i,j) = (0.5/(alfa(i,j)+gama(i,j)))*( alfa(i,j)*(Temperatura(i+1,j) + Temperatura(i-1,j))&
                        - 0.5*beta(i,j)*(Temperatura(i+1,j+1)-Temperatura(i-1,j+1)-Temperatura(i+1,j-1)+Temperatura(i-1,j-1))&
                        + gama(i,j)*(Temperatura(i,j+1) + Temperatura(i,j-1)))
        
        END DO
        ELSE
!-----------------------------------------------------------------------SEM GAP----------------------------------------------------------------------------------
  
        DO i = F2+1, biconcavas - 1
         
        Temperatura(i,j) = (0.5/(alfa(i,j)+gama(i,j)))*( alfa(i,j)*(Temperatura(i+1,j) + Temperatura(i-1,j))&
                        - 0.5*beta(i,j)*(Temperatura(i+1,j+1)-Temperatura(i-1,j+1)-Temperatura(i+1,j-1)+Temperatura(i-1,j-1))&
                        + gama(i,j)*(Temperatura(i,j+1) + Temperatura(i,j-1)))
        END DO
        END IF
 
!END DO
!===============================================================================================================================================================
!------------------------------------------------------------    QUINA (BINCONCAVAS, 1)    ---------------------------------------------------------------------
!===============================================================================================================================================================
    !Cálculo da Temperatura na quina Q4 (biconcavas,1) - quina entre PAREDE LESTE e NORTE - externa               ****(CHECADA EM 17/10 - AGORA É REFLEXIVA)****
    Temperatura(biconcavas,1) = (1.0/(1. + (beta(biconcavas,1)/alfa(biconcavas,1)))) * ((-1./3.)*(-4*Temperatura(biconcavas-1,1) + Temperatura(biconcavas-2,1)) &
        + beta(biconcavas,1)/(3*alfa(biconcavas,1))*(4*Temperatura(biconcavas,2) - Temperatura(biconcavas,3))  )
!===============================================================================================================================================================
!-------------------------------------------------------    CLAD2 + GAP2: EXTERNA CONVECTIVA    ----------------------------------------------------------------
!===============================================================================================================================================================
!DO j=2,radiais-1
    ! PAREDE LC3 - Parede CONVECTIVA externa da pastilha (CHECADO)
        
        Temperatura(biconcavas,j)=(1./( (3*Kc*alfa(biconcavas,j)/(2*Jacobiano(biconcavas,j)*SQRT(alfa(biconcavas,j))*h_out)) + 1 ))&
                      *( (-Kc/(2*Jacobiano(biconcavas,j)*SQRT(alfa(biconcavas,j))*h_out))*(alfa(biconcavas,j)*(-4*Temperatura(biconcavas-1,j) + Temperatura(biconcavas-2,j))&
                      - beta(biconcavas,j)*(Temperatura(biconcavas,j+1) - Temperatura(biconcavas,j-1))) + T_out) 
!END DO

    
!===============================================================================================================================================================
!------------------------------------------------------------   QUINAS DA PAREDE NORTE   -----------------------------------------------------------------------
!===============================================================================================================================================================

    !Cálculo da Temperatura na quina Q2 (1,radiais) - quina entre PAREDE OESTE e NORTE - interna                  ****(CHECADA EM 17/10 - AGORA É REFLEXIVA)****
    Temperatura(1,radiais) = (1.0/(-1. - (beta(1,radiais)/alfa(1,radiais)))) * ((-1./3.)*(4*Temperatura(2,radiais) - Temperatura(3,radiais))&
        + beta(1,radiais)/(3*alfa(1,radiais))*(-4*Temperatura(1,radiais-1) + Temperatura(1,radiais-2))  )

                        
    !Cálculo da Temperatura na quina Q3 (biconcavas,radiais) - quina entre PAREDE LESTE e SUL - externa           ****(CHECADA EM 17/10 - AGORA É REFLEXIVA)****
    Temperatura(biconcavas,radiais) = (1.0/(1. - (beta(biconcavas,radiais)/alfa(biconcavas,radiais)))) * ((-1./3.)*(-4*Temperatura(biconcavas-1,radiais) + Temperatura(biconcavas-2,radiais)) &
        + beta(biconcavas,radiais)/(3*alfa(biconcavas,radiais))*(-4*Temperatura(biconcavas,radiais-1) + TemperaturaOld(biconcavas,radiais-2))  ) 
!===============================================================================================================================================================
!--------------------------------------------------------  PAREDE NORTE (REFLEXIVA SUPERIOR)  ------------------------------------------------------------------
!===============================================================================================================================================================
DO i=2,biconcavas-1 !Considera o corte da simetria superior e inferior: essa parede reflexiva é calculada de uma vez, sem repartir em CLAD, GAP e FUEL.
        
    ! PAREDE OC4 - Reflexiva "norte" (CHECADO)
    Temperatura(i,radiais) =(-1./3.)*(-4*Temperatura(i,radiais-1) + Temperatura(i,radiais-2)) + (beta(i,radiais)/(3*gama(i,radiais)))*(Temperatura(i+1,radiais) - Temperatura(i-1,radiais))
                        		
END DO  

    ENDDO
    

   END SUBROUTINE CalculaTemperaturas_SEMOLD

