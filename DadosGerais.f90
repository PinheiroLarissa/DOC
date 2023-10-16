    MODULE dadosGerais
	! Dados gerais de dimens�o da malha e coeficientes para 'Geramalha'. 
    USE double_Precision
	
   IMPLICIT NONE
    ! ---------------------QUANTIDADE E DEFINI��O DOS PONTOS DA MALHA INICIAL, GERADA POR CASSINI OU ELLIPK (Metodo da Secante)---------------------------    
    INTEGER, PARAMETER                   :: num_bic = 20                              ! N�mero de bic�ncavas conc�ntricas S� do FUEL
    INTEGER                              :: biconcavas, radiais                       ! Total de bic�ncavas da malha de CLAD1+FUEL+CLAD2
    INTEGER                              :: AllocateStatus                            ! STAT do ALLOCATE das matrizes de malha de CLAD e FUEL, e n�mero de pontos das matrizes de CLAD1 e CLAD2
    INTEGER                              :: C1,C2                                     ! C1 = qtdd de pts do CLAD1, C2= qtdd de pts do CLAD2                   
    INTEGER                              :: F1,F2,Ifinal                              ! F1 = Primeiro ponto do FUEL (=C1+1), F2 = �ltimo ponto do FUEL(=F1 + num_bic-1), Ifinal = �ltimo ponto da malha(=F2+1 + C2-1) 
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: X_malha_FUEL, Y_malha_FUEL          ! Matriz final que armazena todos os x e y do FUEL (colunas: bic�ncavas conc�ntricas, linhas: linhas radias).
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: X_Malha_CLAD1,Y_Malha_CLAD1         ! Matrizes que armazenam os valores de x e y do CLAD1 (e to gap1, caso tenha).
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: X_Malha_CLAD2,Y_Malha_CLAD2         ! Matrizes que armazenam os valores de x e y do CLAD2 (e to gap2, caso tenha).    
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
    ! ------------------------------PARAMETRIZ��O E GERA��O DA CURVA MOLDE A PARTIR DA PARAMETRIZA��O E DO METODO DA SECANTE-----------------------------
    REAL(KIND=DP), DIMENSION(0:31)       :: u                                         ! Parte real da parametriza��o
    REAL(KIND=DP), DIMENSION(0:31)       :: v                                         ! Parte imagin�ria da parametriza��o 
    REAL(KIND=DP), DIMENSION(0:SIZE(v)-1):: x_curva_interna                           ! Vetor dos valores da coordenada x da curva interna da bic�ncava da parametriza��o, o molde da geometria (obtida em CASSISNI ou ELLIPK).
    REAL(KIND=DP), DIMENSION(0:SIZE(v)-1):: y_curva_interna                           ! Vetor dos valores da coordenada y da curva interna da bic�ncava da parametriza��o, o molde da geometria (obtida em CASSISNI ou ELLIPK) 
    REAL(KIND=DP), DIMENSION(0:SIZE(v)-1):: x_curva_interna_auxiliar, y_curva_interna_auxiliar ! Curva interna auxiliar para mudar sentido da parametriza��o de anti-hor�rio para hor�rio
    INTEGER                              :: linha_radial, coluna_biconcava            ! Contadores para espelhar as matrizes do 1o QUADRANTE para o 4o QUADRANTE
    REAL(KIND=DP)                        :: step, di                                  ! � o passo entre os pontos da reta radial, 'di': dist�ncia de cada bic6oncava at� a curva interna.
    REAL(KIND=DP)                        :: intervalos                                ! N�mero de intervalos nas retas radiais.
    REAL(KIND=DP)                        :: A_bask, B_bask, C_bask, D_bask, m_incli, Bpp    ! Par�metros da f�rmula de B�skhara para encontrar x_ort. 'm': � a inclina��o da reta secante, 'Bp': coeficiente linear da reta secante.
    REAL(KIND=DP), DIMENSION(0:num_bic-1, 0:SIZE(v)-1):: X_Malha_menor, X_Malha_maior ! Malha auxiliar no processo de c�lculo a partir das retas secante e ortogonal.(est� super dimensionada!)
    REAL(KIND=DP), DIMENSION(0:num_bic-1, 0:SIZE(v)-1):: X_malha1,X_malha2            ! Malhas que armazenam o resultado das duas ra�zes de Bh�skara. (est� super dimensionada para dos CLADS!)
    !------------------------------------------------------------------------------------------------------------------------------------------------------

    ! ------------------------------------------MEDIDAS RELATIVAS �S DIMENS�ES DA PASTILHA---------------------------------------------------------------
    REAL(KIND=DP), PARAMETER             ::dg = 0!0.062e-1                            ! Espessuras do GAP (se dg = 0, ent�o � desconsiderado)
    REAL(KIND=DP), PARAMETER             ::PI =  3.141592653589793238462              ! Valor de pi
    REAL(KIND=DP), PARAMETER             ::dc1_clad1 = 0.826E-1                       ! Espessuras do CLAD1 
    REAL(KIND=DP), PARAMETER             ::dc2_clad2 = 0.826E-1                       ! Espessuras do CLAD2 
    REAL(KIND=DP), PARAMETER             ::df  = 0.2574                            ! Espessura do FUEL, do SpaceClaim e do Yuan
    REAL(KIND=DP), PARAMETER             ::ro  = 3.01E-1                              ! Raio EQUIVALENTE interno do GAP1
    REAL(KIND=DP), PARAMETER             ::rfi = 3.081E-1                             ! Raio EQUIVALENTE interno do FUEL
    REAL(KIND=DP), PARAMETER             ::rfo = 5.988E-1                             ! Raio EQUIVALENTE externo do FUEL
    REAL(KIND=DP), PARAMETER             ::rci = 5.99E-1                              ! Raio EQUIVALENTE externo do GAP2
    REAL(KIND=DP)                        ::dc1                                        ! Espessuras do CLAD1 combinada com DG
    REAL(KIND=DP)                        ::dc2                                        ! Espessuras do CLAD2 combinada com DG    
    INTEGER                              ::razao_fuel_clad1, razao_fuel_clad2                               ! Raz�es entre o FUEL e os CLAD para determinar a quantidade de pontos do CLAD (resultado ap�s arredondamento). 
    REAL(KIND=DP)                        ::razao_fuel_clad1_float,razao_fuel_clad2_float,C1_float,C2_float  ! Raz�es entre a espessura do comburt�vel e do fuel para determinar a qtdd de pts dos CLAD (este resultado � float e por isso � arrendondado). 
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
    ! ------------------------------------------COEFICIENTES E MATRIZES DA GERA��O DA MALHA DE THOMPSON------------------------------------------
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: x,y                                 ! S�o as matrizes x e y da malha completa: CLAD1+FUEL+CLAD2
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: xold,yold,diferencax,diferencay     ! D�o matrizes usadas em Gauss-Seidl como auxiliares para obter a malha convergida de Thompson
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: xksi,yksi,xeta,yeta                 ! Derivadas de x em rela��o a ksi e eta, e y em rela��o a ksi e eta
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: alfa,beta,gama,Jacobiano            ! Coeficientes da malha de Thompson
    !------------------------------------------------------------------------------------------------------------------------------------------------------

   
    CONTAINS 
    !===========================================================================================================================================================================================================================
    !------------------------------------------------------- SUBROUTINE DE LEITURA E ALOCA��O ----------------------------------------------------------------------------------------------------------------------------------
    !===========================================================================================================================================================================================================================

    SUBROUTINE LeDados ! Aloca as matrizes de malha, coeficientes e derivadas do m�todo de Thompson
	
    USE double_Precision
    
    IMPLICIT NONE
             
    !L� os dados do problema: condi��es de contorno e pontos estimados de malha
	integer::i !n�mero de colunas, c�rculos conc�ntricos no dom�nio f�sico
    integer::j !n�mero de linhas, retas radiais que cortam os c�rculos no dom�nio f�sico

	OPEN(unit=1,file='geometria.txt',status='old')

  
    
    ALLOCATE(X(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    ALLOCATE(Y(1:biconcavas, 1:radiais), STAT = AllocateStatus)
    
    ALLOCATE(X_Malha_FUEL (0:num_bic-1, 0:(radiais) -1), STAT = AllocateStatus)
    ALLOCATE(Y_Malha_FUEL (0:num_bic-1, 0:(radiais) -1), STAT = AllocateStatus)
    ALLOCATE(X_Malha_CLAD1(0:C1-1,      0:(radiais) -1), STAT = AllocateStatus)
    ALLOCATE(Y_Malha_CLAD1(0:C1-1,      0:(radiais) -1), STAT = AllocateStatus)
    ALLOCATE(X_Malha_CLAD2(0:C2-1,      0:(radiais) -1), STAT = AllocateStatus)
    ALLOCATE(Y_Malha_CLAD2(0:C2-1,      0:(radiais) -1), STAT = AllocateStatus)

    ALLOCATE(xold         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(yold         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(diferencax   (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(diferencay   (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
   
    ALLOCATE(Xksi         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(yksi         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(xeta         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(yeta         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(alfa         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(beta         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(gama         (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    ALLOCATE(Jacobiano    (1:biconcavas, 1:radiais),   STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Alocacao mal sucedida na SUBROUTINE LeDados do MODULE dadosGerais ***"
    
    ! L� as malhas rec�m calculadas em Cassini_geo ou Ellipk_geo
    !READ(1,*) ((x(i,j),i=1,biconcavas),j=1,radiais)
    !READ(1,*) ((y(i,j),i=1,biconcavas),j=1,radiais)

    CLOSE(unit=1)
    

        END SUBROUTINE Ledados
    END MODULE DadosGerais
    
    