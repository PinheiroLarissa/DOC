    SUBROUTINE  CalculaDerivadas_KsiEta
    ! Calcula as derivadas primeiras xeta, xksi, yeta e yksi com x e y.
    ! Essa derivadas serão usadas para calcular os coeficientes alfa, beta, gama e jacobiana.

	USE dadosGerais
    USE double_Precision
    IMPLICIT NONE 
    INTEGER::i,j 
    
    ALLOCATE(x   (1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(y   (1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(xksi(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(yksi(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(xeta(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(yeta(1:biconcavas, 1:radiais), STAT = AllocateStatus)

    !**************************************************************************************************************
    !--------------------------------------------CONDIÇÕES DE CONTORNO---------------------------------------------
    !**************************************************************************************************************
    
    !=================================================QUINAS=======================================================
    ! Calculo na quina Q1 (- Discretização por diferenças avançadas - conferidas e Ok!

	xeta(1,1) = (1.0/2.0)*( -3.0*x(1,1)+4.0*x(1,2)-x(1,3))
	yeta(1,1) = (1.0/2.0)*( -3.0*y(1,1)+4.0*y(1,2)-y(1,3))
	xksi(1,1) = (1.0/2.0)*( -3.0*x(1,1)+4.0*x(2,1)-x(3,1))
	yksi(1,1) = (1.0/2.0)*( -3.0*y(1,1)+4.0*y(2,1)-y(3,1))
    
    
    ! Calcula da Quina Q2 (1, radiais) - Discretização por dif. avançadas em Ksi e recuadas em Eta - conferidas e OK!
    
    xeta(1,radiais) = (1.0/2.0)*( 3.0*x(1,radiais) -4.0*x(1,radiais-1)+x(1,radiais-2))
	yeta(1,radiais) = (1.0/2.0)*( 3.0*y(1,radiais) -4.0*y(1,radiais-1)+y(1,radiais-2))
	xksi(1,radiais) = (1.0/2.0)*(-3.0*x(1,radiais) +4.0*x(2,radiais)-x(3,radiais))
	yksi(1,radiais) = (1.0/2.0)*(-3.0*y(1,radiais) +4.0*y(2,radiais)-y(3,radiais))
    
    !Cálculo da Quina Q3 (biconcavas, radiais) - Discretização por dif. recuadas - conferidas e OK!

	xeta(biconcavas,radiais) = (1.0/2.0)*( 3.0*x(biconcavas,radiais)-4.0*x(biconcavas,radiais-1)+x(biconcavas,radiais-2))
	yeta(biconcavas,radiais) = (1.0/2.0)*( 3.0*y(biconcavas,radiais)-4.0*y(biconcavas,radiais-1)+y(biconcavas,radiais-2))
	xksi(biconcavas,radiais) = (1.0/2.0)*( 3.0*x(biconcavas,radiais)-4.0*x(biconcavas-1,radiais)+x(biconcavas-2,radiais))
	yksi(biconcavas,radiais) = (1.0/2.0)*( 3.0*y(biconcavas,radiais)-4.0*y(biconcavas-1,radiais)+y(biconcavas-2,radiais))
    
    ! Cálculo da Quina Q4 (biconcavas,1) - Discretização por dif. avançadas em Eta e recuadas em Ksi - conferidas e OK!
	xeta(biconcavas,1) = (1.0/2.0)*( -3.0*x(biconcavas,1)+4.0*x(biconcavas,2)-x(biconcavas,3))
	yeta(biconcavas,1) = (1.0/2.0)*( -3.0*y(biconcavas,1)+4.0*y(biconcavas,2)-y(biconcavas,3))
	xksi(biconcavas,1) = (1.0/2.0)*(  3.0*x(biconcavas,1)-4.0*x(biconcavas-1,1)+x(biconcavas-2,1))
	yksi(biconcavas,1) = (1.0/2.0)*(  3.0*y(biconcavas,1)-4.0*y(biconcavas-1,1)+y(biconcavas-2,1))

    !===============================================PAREDES=======================================================

    ! PAREDE SUL - conferidas e Ok!
	do i=2,biconcavas-1 

	xeta(i,1) = (1.0/2.0)*( -3.0*x(i,1)+4.0*x(i,2)-x(i,3))
	yeta(i,1) = (1.0/2.0)*( -3.0*y(i,1)+4.0*y(i,2)-y(i,3))
	xksi(i,1) = (1.0/2.0)*(x(i+1,1)-x(i-1,1))
	yksi(i,1) = (1.0/2.0)*(y(i+1,1)-y(i-1,1))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!QUINA INTERNA
    
        
    xeta(F1-1,1) = (1.0/2.0)*( -3.0*x(F1-1,1)+4.0*x(F1-1,2)-x(F1-1,3))
	yeta(F1-1,1) = (1.0/2.0)*( -3.0*y(F1-1,1)+4.0*y(F1-1,2)-y(F1-1,3))
	xksi(F1-1,1) = (1.0/2.0)*(  3.0*x(F1-1,1)-4.0*x(F1-2,1)+x(F1-3,1))
	yksi(F1-1,1) = (1.0/2.0)*(  3.0*y(F1-1,1)-4.0*y(F1-2,1)+y(F1-3,1))
	end do
    
    !PAREDE NORTE - conferidas e OK!
	do i=2,biconcavas-1

	xeta(i,radiais) = (1.0/2.0)*( 3.0*x(i,radiais)-4.0*x(i,radiais-1)+x(i,radiais-2))
	yeta(i,radiais) = (1.0/2.0)*( 3.0*y(i,radiais)-4.0*y(i,radiais-1)+y(i,radiais-2))
	xksi(i,radiais) = (1.0/2.0)*(x(i+1,radiais)-x(i-1,radiais))
	yksi(i,radiais) = (1.0/2.0)*(y(i+1,radiais)-y(i-1,radiais))
ENDDO

    ! PAREDE LESTE - conferidas e OK!
	do j=2,radiais-1

    xeta(biconcavas,j) = (1.0/2.0)*(x(biconcavas,j+1)-x(biconcavas,j-1))
	yeta(biconcavas,j) = (1.0/2.0)*(y(biconcavas,j+1)-y(biconcavas,j-1))
	xksi(biconcavas,j) = (1.0/2.0)*( 3.0*x(biconcavas,j)-4.0*x(biconcavas-1,j)+x(biconcavas-2,j))
	yksi(biconcavas,j) = (1.0/2.0)*( 3.0*y(biconcavas,j)-4.0*y(biconcavas-1,j)+y(biconcavas-2,j))
    
	end do

    !PAREDE OESTE - conferidas e OK!
    do j=2,radiais-1
        
    xeta(1,j) = (1.0/2.0)*(x(1,j+1)-x(1,j-1))
	yeta(1,j) = (1.0/2.0)*(y(1,j+1)-y(1,j-1))
	xksi(1,j) = (1.0/2.0)*( -3.0*x(1,j)+4.0*x(2,j)-x(3,j))
	yksi(1,j) = (1.0/2.0)*( -3.0*y(1,j)+4.0*y(2,j)-y(3,j))
	
    end do
    
    !PAREDE CLAD1-GAP1(OC3) 
    do j=2,radiais-1
        
    xeta(F1-1,j) = (1.0/2.0)*(x(F1-1,j+1)-x(F1-1,j-1))
	yeta(F1-1,j) = (1.0/2.0)*(y(F1-1,j+1)-y(F1-1,j-1))
	xksi(F1-1,j) = (1.0/2.0)*( 3.0*x(F1-1,j)-4.0*x(F1-2,j)+x(F1-3,j))
	yksi(F1-1,j) = (1.0/2.0)*( 3.0*y(F1-1,j)-4.0*y(F1-2,j)+y(F1-3,j))
	
    end do
    
    !PAREDE GAP1-FUEL(S1) 
    do j=2,radiais-1
        
    xeta(F1,j) = (1.0/2.0)*(x(F1,j+1)-x(F1,j-1))
	yeta(F1,j) = (1.0/2.0)*(y(F1,j+1)-y(F1,j-1))
	xksi(F1,j) = (1.0/2.0)*( -3.0*x(F1,j)+4.0*x(F1+1,j)-x(F1+2,j))
	yksi(F1,j) = (1.0/2.0)*( -3.0*y(F1,j)+4.0*y(F1+1,j)-y(F1+2,j))
	
    end do
    
    !PAREDE FUEL-GAP2(S3) 
    do j=2,radiais-1
        
    xeta(F2,j) = (1.0/2.0)*(x(F2,j+1)-x(F2,j-1))
	yeta(F2,j) = (1.0/2.0)*(y(F2,j+1)-y(F2,j-1))
	xksi(F2,j) = (1.0/2.0)*( 3.0*x(F2,j)-4.0*x(F2-1,j)+x(F2-2,j))
	yksi(F2,j) = (1.0/2.0)*( 3.0*y(F2,j)-4.0*y(F2-1,j)+y(F2-2,j))
	
    end do
    
    !PAREDE GAP2-CLAD2(LC1) 
    do j=2,radiais-1
        
    xeta(F2+1,j) = (1.0/2.0)*(x(F2+1,j+1)-x(F2+1,j-1))
	yeta(F2+1,j) = (1.0/2.0)*(y(F2+1,j+1)-y(F2+1,j-1))
	xksi(F2+1,j) = (1.0/2.0)*( -3.0*x(F2+1,j)+4.0*x(F2+2,j)-x(F2+3,j))
	yksi(F2+1,j) = (1.0/2.0)*( -3.0*y(F2+1,j)+4.0*y(F2+2,j)-y(F2+3,j))
    
	
    end do
    !**************************************************************************************************************
    !--------------------------------------------PONTOS INTERNOS---------------------------------------------------
    !**************************************************************************************************************
    
    !-------------------------------------------------CLAD1--------------------------------------------------------
    DO i=2, F1-2
    DO j=2, radiais-1
            
    xeta(i,j) = (1.0/2.0)*(x(i,j+1)-x(1,j-1))
	yeta(i,j) = (1.0/2.0)*(y(i,j+1)-y(1,j-1))
    xksi(i,j) = (1.0/2.0)*(x(i+1,j)-x(i-1,j))
	yksi(i,j) = (1.0/2.0)*(y(i+1,j)-y(i-1,j))
    
    END DO
    END DO
    
    !--------------------------------------------------FUEL--------------------------------------------------------

    DO i= F1+1, F2-1
    DO j=2, radiais-1
            
    xeta(i,j) = (1.0/2.0)*(x(i,j+1)-x(1,j-1))
	yeta(i,j) = (1.0/2.0)*(y(i,j+1)-y(1,j-1))
    xksi(i,j) = (1.0/2.0)*(x(i+1,j)-x(i-1,j))
	yksi(i,j) = (1.0/2.0)*(y(i+1,j)-y(i-1,j))
    
    END DO
    END DO
    
    !--------------------------------------------------CLAD2--------------------------------------------------------

    DO i=F2+2, biconcavas-1
    DO j=2, radiais-1
            
    xeta(i,j) = (1.0/2.0)*(x(i,j+1)-x(1,j-1))
	yeta(i,j) = (1.0/2.0)*(y(i,j+1)-y(1,j-1))
    xksi(i,j) = (1.0/2.0)*(x(i+1,j)-x(i-1,j))
	yksi(i,j) = (1.0/2.0)*(y(i+1,j)-y(i-1,j))
    
    END DO
    END DO

	END SUBROUTINE CalculaDerivadas_KsiEta
