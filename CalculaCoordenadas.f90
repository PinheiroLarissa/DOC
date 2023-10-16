    SUBROUTINE CalculaCoordenadas
    !A ortogonalidade não está sendo usada porque não converge
    !Última revisão: 23/07/2018 - Larissa

	USE DadosGerais

	IMPLICIT NONE

	INTEGER::i,j

	!Cálculo das coordenadas x e y

!*********************************************************************************************************************
!--------------------------------------------CONDIÇÕES DE ORTOGONALIDADE----------------------------------------------
!*********************************************************************************************************************
!Obs.: Comente qualquer um dos loops abaixo para ter condição preescrita e não ortogonalidade no contorno da malha
    
    !PAREDE OESTE (CONVECTIVA INTERNA)
	!do j=2,radiais-1
        
    !x(1,j) = ( (-3.0*y(1,j)+4.0*y(2,j)-3.0*y(3,j)) * (y(1,j+1)-y(1,j-1))& 
    !         + (4.0*x(2,j)-x(3,j)) * (x(1,j+1)-x(1,j-1)) )/( 3.0*(x(1,j+1)-x(1,j-1) ))
    
    !y(1,j) = ( (-3.0*x(1,j)+4.0*x(2,j)-3.0*x(3,j)) * (x(1,j+1)-x(1,j-1))& 
    !        + (4.0*y(2,j)-y(3,j)) * (y(1,j+1)-y(1,j-1)) )/( 3.0*(y(1,j+1)-y(1,j-1) ))

    !end do 	
    
    !PAREDE LESTE (CONVECTIVA EXTERNA)
    !do j=2,radiais-1
        
    !x(biconcavas,j) = ( -(3.0*y(biconcavas,j)-4.0*y(biconcavas-1,j)+y(biconcavas-2,j)) * (y(biconcavas,j+1)-y(biconcavas,j-1))& 
    !             - (-4.0*x(biconcavas-1,j)+x(biconcavas-2,j)) * (x(biconcavas,j+1)-x(biconcavas,j-1)) )/( 3.0*(x(biconcavas,j+1)-x(biconcavas,j-1) ))
    
    ! y(biconcavas,j) = ( -(3.0*x(biconcavas,j)-4.0*x(biconcavas-1,j)+x(biconcavas-2,j)) * (x(biconcavas,j+1)-x(biconcavas,j-1))& 
    !             - (-4.0*y(biconcavas-1,j)+y(biconcavas-2,j)) * (y(biconcavas,j+1)-y(biconcavas,j-1)) )/( 3.0*(y(biconcavas,j+1)-y(biconcavas,j-1) ))
    
    !end do
    
    !PAREDE NORTE (REFLEXIVA)
    !do i=2,biconcavas-1
    
    !x(i,radiais) = ( -(3.0*y(i,radiais)-4.0*y(i,radiais-1)+y(i,radiais-2)) * (y(i+1,radiais)-y(i-1,radiais))& 
    !  - (-4.0*x(i,radiais-1)+x(i,radiais-2)) * (x(i+1,radiais)-x(i-1,radiais)) )/( 3.0*(x(i+1,radiais)-x(i-1,radiais) ))
        
    !y(i,radiais) = ( -(3.0*x(i,radiais)-4.0*x(i,radiais-1)+x(i,radiais-2)) * (x(i+1,radiais)-x(i-1,radiais))& 
    ! - (-4.0*y(i,radiais-1)+y(i,radiais-2)) * (y(i+1,radiais)-y(i-1,radiais)) )/( 3.0*(y(i+1,radiais)-y(i-1,radiais) ))
    
    !end do
        
    !PAREDE SUL (REFLEXIVA)
    !do i=2,biconcavas-1
        
    !x(i,1) = ( (-3.0*y(i,1)+4.0*y(i,2)-y(i,3)) * (y(i+1,1)-y(i-1,1))& 
    !        + (4.0*x(i,2)-x(i,3)) * (x(i+1,1)-x(i-1,1)) )/( 3.0*(x(i+1,1)-x(i-1,1) ))
    
    !y(i,1) = ( (-3.0*x(i,1)+4.0*x(i,2)-x(i,3)) * (x(i+1,1)-x(i-1,1))& 
    !      + (4.0*y(i,2)-x(i,3)) * (y(i+1,1)-y(i-1,1)) )/( 3.0*(y(i+1,1)-y(i-1,1) ))
    
    !end do
        
!*********************************************************************************************************************
!-------------------------------------------------PONTOS INTERNOS-----------------------------------------------------
!********************************************************************************************************************* 
!Obs.: equação 28, pg.15 do caderno
    
    !======================================================================================
    !------------------------------------------CLAD1---------------------------------------
    !======================================================================================
    
    DO i=2,F1-2
	DO j=2,radiais-1

	x(i,j)= (0.5/(alfa(i,j)+gama(i,j))*(alfa(i,j)*(xold(i+1,j)+xold(i-1,j))&
	        -0.5*beta(i,j)*(xold(i+1,j+1)-xold(i-1,j+1)-xold(i+1,j-1)+xold(i-1,j-1))&
			+gama(i,j)*(xold(i,j+1)+xold(i,j-1))))



	y(i,j)= (0.5/(alfa(i,j)+gama(i,j))*(alfa(i,j)*(yold(i+1,j)+yold(i-1,j))&
	        -0.5*beta(i,j)*(yold(i+1,j+1)-yold(i-1,j+1)-yold(i+1,j-1)+yold(i-1,j-1))&
			+gama(i,j)*(yold(i,j+1)+yold(i,j-1))))
    ENDDO
    ENDDO
 
    !======================================================================================
    !------------------------------------------FUEL---------------------------------------
    !======================================================================================
    DO i=F1+1,F2-1
	DO j=2,radiais-1

	x(i,j)= (0.5/(alfa(i,j)+gama(i,j))*(alfa(i,j)*(xold(i+1,j)+xold(i-1,j))&
	        -0.5*beta(i,j)*(xold(i+1,j+1)-xold(i-1,j+1)-xold(i+1,j-1)+xold(i-1,j-1))&
			+gama(i,j)*(xold(i,j+1)+xold(i,j-1))))


	y(i,j)= (0.5/(alfa(i,j)+gama(i,j))*(alfa(i,j)*(yold(i+1,j)+yold(i-1,j))&
	        -0.5*beta(i,j)*(yold(i+1,j+1)-yold(i-1,j+1)-yold(i+1,j-1)+yold(i-1,j-1))&
			+gama(i,j)*(yold(i,j+1)+yold(i,j-1))))
    
    ENDDO
    ENDDO
 
    !======================================================================================
    !------------------------------------------CLAD2---------------------------------------
    !======================================================================================
    DO i=F2+2,Ifinal-1
	DO j=2,radiais-1

	x(i,j)= (0.5/(alfa(i,j)+gama(i,j))*(alfa(i,j)*(xold(i+1,j)+xold(i-1,j))&
	        -0.5*beta(i,j)*(xold(i+1,j+1)-xold(i-1,j+1)-xold(i+1,j-1)+xold(i-1,j-1))&
			+gama(i,j)*(xold(i,j+1)+xold(i,j-1))))


	y(i,j)= (0.5/(alfa(i,j)+gama(i,j))*(alfa(i,j)*(yold(i+1,j)+yold(i-1,j))&
	        -0.5*beta(i,j)*(yold(i+1,j+1)-yold(i-1,j+1)-yold(i+1,j-1)+yold(i-1,j-1))&
			+gama(i,j)*(yold(i,j+1)+yold(i,j-1))))
    

    ENDDO
    ENDDO
        
    END SUBROUTINE CalculaCoordenadas