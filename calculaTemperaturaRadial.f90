    subroutine calculaTemperaturaRadial
    ! Versão          : Versão 2
    ! Data da revisão : 13/08/2018
    ! Revisor         : Larissa
    ! Observações da última revisão: Transposição dos dados de entrada, radiais = 200 (200 pontos nas retas radiais)
    !                                Alteração DOs sinais das equações, conforma anotações em caderno
    !                                Alteração da marcha, conforme anotação em caderno.
    !                                Malha bicôncava com simetria de 1/4 e refino de 93x200.
    !                                Malha bicôncava com simetria de 1/4 e refino de 45x100.

    
	USE DadosGerais
    USE double_Precision
	use DadosTermicos

	IMPLICIT NONE

    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE::b,fonte,p,a,h,c
    REAL(KIND=DP) :: deltar,raio_interno,raio_externo
	INTEGER   :: j,Numero_pontos
    
    ! daDOsproblema.txt apresenta x e y com 37 linhas (radiais, retais radiais) e 15 colunas (biconcavas, círculos concêntricos)
    OPEN(unit=1,file='dadosproblema.txt',status='old')
    
    ! Obtém a quantidade dos pontos, principalmente biconcavas que será utilizaDO (quantidade de pontos numa reta radial)
    READ(1,*) radiais,biconcavas
    ALLOCATE(b(0:biconcavas-1))
    ALLOCATE(fonte(0:biconcavas-1))
    ALLOCATE(p(0:biconcavas-1))
    ALLOCATE(a(0:biconcavas-1))
    ALLOCATE(h(0:biconcavas-1))
    ALLOCATE(c(0:biconcavas-1))


    CALL LedadosTermicos 
    
    Numero_pontos = biconcavas-1!porque parte de ZERO são 0 + 14 pontos = 15 pontos

	raio_interno = 1.0
    raio_externo = 1.6
	deltar = (raio_externo - raio_interno)/REAL((biconcavas -1.))
    
    ! Define toDOs os pontos DO vetor TemperaturaR = 1000
	TemperaturaR=1000.0
    
    
    ! Calcula os pontos internos
	DO j=1,Numero_pontos-1! PartinDO DO zero

    a(j)= (1.0 - (1.0/(2.0*(j+(raio_interno/deltar)))))

	c(j)= (1.0 + (1.0/(2.0*(j+(raio_interno/deltar)))))

	b(j)= -2.0

	END DO

	WRITE(*,*)deltar,h_int,Kf

    ! Calcula os pontos DO contorno interno
    b(0) = -2.0*( 1.0+(1.0 - (1.0/(2*raio_interno/deltar)))*(deltar*h_int)/Kf)
	!b(1) = -2.0*(1.0+(1.0-(1.0/(2.0*(1+raio_interno/deltar))*(deltar*h_int)/Kf)))

	c(0) = 2.0
    
    ! Calcula os pontos do contorno externo
	a(Numero_pontos) = 2.0
    
    b(Numero_pontos) = -2.0*( 1.0+(1.0 + (1.0/(2.0*(Numero_Pontos+raio_interno/deltar))))*(deltar*h_out)/Kf)
	!b(Numero_pontos) = -2.0*( 1.0+(1.0-(1.0/(2.0*(Numero_Pontos+raio_interno/deltar))*(deltar*h_out)/Kf)))
    
    ! Define o valor de toDOs os pontos do vetor 'fonte'
	fonte = -(((deltar)**2)*densidade_de_potencia)/Kf
    
    ! Altera o valor da fonte no contorno (0)
	fonte(0) = fonte(0)-2.0*(1.0- (1.0/(2.0*raio_interno/deltar)))*((deltar*h_int*T_int)/Kf)
    
    ! Altera o valor da fonte no contorno (Numero_pontos)
    fonte(Numero_pontos) = fonte(Numero_pontos)-2.0*(1.0 + (1.0/(2.0*(Numero_pontos+raio_interno/deltar))))*((deltar*h_int*T_int)/Kf)
	!fonte(Numero_pontos) = fonte(Numero_pontos)-2.0*(1.0+(1.0/(2.0*(Numero_pontos+raio_interno/deltar)))*((deltar*h_out*T_ext)/Kf))

	! Começa a resolução da matriz tridiagonal (ALgoritmo de Thomas)
	h(0)= c(0)/b(0)

	p(0)=fonte(0)/b(0)

	DO j=1,Numero_pontos-1  ! Pontos internos

	h(j) = (c(j))/(b(j)-a(j)*h(j-1))

	END DO

	DO j =1,Numero_pontos

	p(j) = (fonte(j)-a(j)*p(j-1))/(b(j)-a(j)*h(j-1))

	END DO

	TemperaturaR(Numero_pontos) = p(Numero_pontos)

	DO j = Numero_pontos-1,0,-1

	TemperaturaR(j) = p(j) - h(j)*TemperaturaR(j+1)

	END DO

	open(unit=88,file='TemperaturaRadial.txt',status='unknown')
    
    WRITE(88,*) 'TemperaturaR ='
    DO j = 0, Numero_pontos
    WRITE(88,*)j, ' ', TemperaturaR(j)
    END DO
    
	CLOSE(88)
    
	DO j = 0, Numero_pontos
	WRITE(*,*)TemperaturaR(j)
    END DO

	WRITE(*,*) fonte

	END SUBROUTINE



