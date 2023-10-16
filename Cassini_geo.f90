   SUBROUTINE Cassini_geo
    
    
    !*********************************************************************************************************************************************************************************************************************   
    ! Programador: Larissa
    ! �ltima revis�o: 14/07/2019
    !
    !Descri��o da subortina Cassini_geo:
    !M�dulo de defini��o geral, contendo: malha em x, y, coeficientes, n�mero de bic�ncavas e retas radiais

    
    USE double_Precision
    USE dadosGerais
    IMPLICIT NONE 
    
    INTEGER      :: i, j, g, indice_bic                         ! Contadores, i -->> linhas, u da parametriza��o, j -->> colunas, v da parametriza��o, g --> contador que corre entre bic�ncavas nas retas radiais.
    INTEGER      :: k, q, p, col                                ! Contadores utilizados nos loops de concatenamento das matrizes de malha
    REAL(KIND=DP):: a, b, ab                                    ! Par�metros de Cassini, 'a', 'b' e 'ab'; 'di': dist�ncia de cada bic6oncava at� a curva interna.
    REAL(KIND=DP):: u0                                          ! u = u0 � o 'u' escolhido para a parametriza��o com base nos a e b dados pelo usu�rio.
    REAL(KIND=DP), DIMENSION(0:SIZE(v)-1):: ro0                 ! Vetor utilizada na obten��o dos da curva interna atrav�s da parametriza��o, com u = u0.

    REAL(KIND=DP), DIMENSION(0:SIZE(v)-1,0:SIZE(u)-1):: rho     ! Matriz utilizada na obten��o dos contornos atrav�s da parametriza��o.
    REAL(KIND=DP), DIMENSION(0:SIZE(v)-1,0:SIZE(u)-1):: x_cas   ! Matriz de coordenadas x dos contornos dados pelas parametriza��o de Cassini.
    REAL(KIND=DP), DIMENSION(0:SIZE(v)-1,0:SIZE(u)-1):: y_cas   ! Matriz de coordenadas y dos contornos dados pelas parametriza��o de Cassini.
     
 
    !***********************************************************************************************************************************************************************************
    !----------------------------------------------------------------- ABERTURA DOS ARQUIVOS: SALVAR A GEOMETRIA E OS CONTORNOS ------------------------------------------------------------
    !***************************************************************************************************************************************************************************************
    
    OPEN(UNIT = 1, FILE = 'geometria.txt', STATUS = 'unknown')           ! Salva a geometra malhada para fornecer como entrada do m�todo de Thompson
    OPEN(UNIT = 200, FILE = 'contornos_cassini.txt', STATUS = 'unknown') ! Salva todos os contornos para que se possa analisar a qualidade da parametriz��o.

    ! Dist�ncia entre os focos (m)
    a  = 0.3151
    
    ! Raz�o a/b que caracteriza as curvas de Cassini (m)
    ab = 0.843
    
    ! A dist�ncia entre um ponto da curva at� os focos � b^2 (m^2)
    b = a/ab
    
    !***********************************************************************************************************************************************************************************
    !--------------------------------------------------------- PARAMETRIZA��O DAS CURVAS DE CASSINI: TODOS OS CONTORNOS ----------------------------------------------------------------
    !***********************************************************************************************************************************************************************************
    
    ! Cria��o do vetor u, que vai de 0 a PI em (SIZE(u) - 1) passos
    DO i = 0, SIZE(u) - 1  
    u(i) =  i*PI/(SIZE(u) - 1)
    ENDDO
    
    ! Cria��o do vetor v, que vai de 0 a PI em (SIZE(v) - 1) passos
    DO j = 0, SIZE(u) - 1 
    v(j) = j*PI/(SIZE(u) - 1)
    ENDDO
    
   ! Cria��o das matrizes x e y de Cassini, para tra��o todos os contornos com o a escolhido, i ->> linhas, u   j->> colunas, v
    DO i = 0, SIZE(u) - 1 
        DO j = 0, SIZE(u) - 1
     

    rho(j,i) = SQRT((EXP(2*u(i)) + 2*EXP(u(i))*COS(v(j)) + 1.0))
    x_cas(j,i) = (a/SQRT(2.0)) * SQRT(ABS((rho(j,i) + (EXP(u(i))*COS(v(j)) + 1.0) )))
    y_cas(j,i) = (a/SQRT(2.0)) * SQRT(ABS((rho(j,i) - (EXP(u(i))*COS(v(j)) + 1.0) )))

        ENDDO
    ENDDO
    
100 FORMAT(51(f14.9),/)
    WRITE(200,*) '======================================================================================'
    WRITE(200,*) 'CURVAS GERADAS UTILIZANDO A PARAMETRIZA��O DE MOON E SPENCER PARA AS CURVAS DE CASSINI'
    WRITE(200,*) '======================================================================================'
    WRITE(200,*) 'a = ',a, 'b = ',b, 'a/b = ', ab,'# colunas: ',j, '# colunas: ',i
    WRITE(200,*) 'Matriz x, linhas--> u, colunas--> v'
    WRITE(200,100) x_cas
    WRITE(200,*) 'Matriz y, linhas--> u, colunas--> v' 
    WRITE(200,100) y_cas
    
!***********************************************************************************************************************************************************************************
!--------------------------------------------------------------------- ESCOLHA DO CONTORNO: U = U0 ---------------------------------------------------------------------------------
!***********************************************************************************************************************************************************************************
    
    
    ! C�lculo de u = u0 para usar na parametriza��o com base no a e no b escolhidos acima
    u0 = log( (a**2 + b**2)/a**2 - 1. )
    
    
    ! Cria��o dos vetores x e y da curva interna (que ser� extrudada com a espessura dada pelo SpaceClaim, conforme m�todo da secante, pg 128 caderno) 
 
        DO j = 0, SIZE(v) - 1
            
    ro0(j) = SQRT(EXP(2*u0) + 2.0*EXP(u0)*COS(v(j)) + 1.0)
    x_curva_interna(j) = a/SQRT(2.0) * SQRT( ABS((ro0(j) + (EXP(u0)*COS(v(j)) + 1.0))) )
    y_curva_interna(j) = a/SQRT(2.0) * SQRT( ABS((ro0(j) - (EXP(u0)*COS(v(j)) + 1.0))) )
        ENDDO
        

!1000   FORMAT(51(f14.9),/)
 !   WRITE(*,*) 'x=['
  !  WRITE(*,1000) x_curva_interna
   !WRITE(*,*) ']' 
   !WRITE(*,*) 'y=['
    !WRITE(*,1000) y_curva_interna
     !WRITE(*,*) ']'
    !WRITE(*,*) 'plot(x,y)'
     !pause
     
         !**ATEN��O: ESCREVENDO A CURVA INTERNA AGORA EM SENTIDO HOR�RIO PARA N�O TER DISCONTINUIDADE QUANDO UTILIZAR A SIMETRIA DE 1/2 DO LADO DIREITO
    x_curva_interna_auxiliar = x_curva_interna
    y_curva_interna_auxiliar = y_curva_interna
    
    k=SIZE(v)
        DO j = 0, SIZE(v) - 1
    k=k-1
    x_curva_interna(j) =  x_curva_interna_auxiliar(k)
    y_curva_interna(j) =  y_curva_interna_auxiliar(k)
        ENDDO

     !   WRITE(*,*) 'x1=['
    !WRITE(*,1000) x_curva_interna
   !WRITE(*,*) ']' 
   !WRITE(*,*) 'y1=['
    !WRITE(*,1000) y_curva_interna
     !WRITE(*,*) ']'
    !WRITE(*,*) 'plot(x1,y1)'
     !pause
!***********************************************************************************************************************************************************************************
!------------------------------------------------------------------ M�TODO DA SECANTE: ENCONTRAM-SE X e Y A CADA RETA ORTOGONAL ----------------------------------------------------
!***********************************************************************************************************************************************************************************
    ! 1. Para verificar os c�lculos, pg. 128 - 129 caderno.
    ! 2. A curva � feita no sentido anti-hor�rio;
    ! 3. A curva interna, gerada pela parametriza��o, � salva na primeira coluna das matrizes X e Y, que ser�o usadas nos c�lculos de gera��o de malha por Thompson;
    ! 4. A curva externa ainda n�o � conhecida, apenas as quinas externas. Essas quinas tamb�m s�o salvas nas matrizes X e Y;
    ! 5. Loop em j: anda ao longo da bic�ncava. Loop em g: anda ao longo da reta ortogonal (ou radial), desenhando a malha da primeira bic�ncava at� a �ltima;
    ! 6. Assim, as colunas das matrizes X e Y s�o bic�ncavas conc�ntricas, e as linhas s�o as retas radiais (desenhadas como retas ortogonais �s retas secantes da curva interna);
    ! 7. Os pontos ao longo de uma reta radial s�o tra�ados acrescentando a espessura "di" � curva interna, sendo di uma fra��o incremental da espessura total "d".
       
        ! Combinando a espessura do CLAD e do GAP para facilitar os c�lculos, e se dg = 0, o primeiro ponto do GAP � bypassado
        dc1 = dc1_clad1 + dg 
        dc2 = dc2_clad2 + dg

!===================================================================================================================================================================================
!--------------------------------------------------------------------------------------------FUEL-----------------------------------------------------------------------------------
!===================================================================================================================================================================================
    ALLOCATE(X_Malha_FUEL(0:num_bic-1, 0:2*(SIZE(v))-2))
    ALLOCATE(Y_Malha_FUEL(0:num_bic-1, 0:2*(SIZE(v))-2))
    
    ! Coordenadas iniciais de X_MALHA_FUEL e Y_MALHA_FUEL  
    X_Malha_FUEL(0,:) = x_curva_interna(:) ! Toda a coluna 0 � prenchida pela bic�ncava definida pela curva interna (primeira bic�ncava conc�ntrica)
    Y_Malha_FUEL(0,:) = y_curva_interna(:) ! Toda a coluna 0 � prenchida pela bic�ncava definida pela curva interna (primeira bic�ncava conc�ntrica)
    
!10000 FORMAT(101(f14.9),/)
 !  WRITE(*,*) 'x2=['
  !  WRITE(*,10000) X_Malha_FUEL(0,:)
   ! WRITE(*,*) ']' 
   ! WRITE(*,*) 'y2=['
    !WRITE(*,10000) Y_Malha_FUEL(0,:)
    ! WRITE(*,*) ']'
     !WRITE(*,*) 'plot(x2,y2)'
    !PAUSE
     
    ! N�o tenho a curva externa ainda, mas tenho as duas "quinas" externas, preciso deles porque uso o ponto da frente e o de tr�s das bic�ncavas interna e externa para o c�culo do coeficiente angular da reta secante.
    X_Malha_FUEL(num_bic-1, 0) = 0                                                  ! O primeiro x da bic. externa � zero. 
    Y_Malha_FUEL(num_bic-1, 0) = y_curva_interna(0) + df                            ! O primeiro y da bic. externa � apenas a soma do primeiro pt. da bic. interna + a espessura.
            
    X_Malha_FUEL(num_bic-1, SIZE(v) - 1) = x_curva_interna(SIZE(v) - 1) + df        ! O �ltimo x da bic. externa � apenas a soma do �ltimo y da bic. interna + a espessura.                        
    Y_Malha_FUEL(num_bic-1, SIZE(v) - 1) = 0                                        ! O �ltimo y da bic. externa � zero.
   
    intervalos = num_bic - 1                ! N�mero total de intervalos entre as bic�ncavas conc�ntricas.
    step    = df/(intervalos)               ! Passo entre os pontos das retas radiais (i. e. bic�ncavas), onde d � a dist�ncia entre a primeira e a �ltima bic�ncava (curva interna e externa).
    ! Loop para encontrar os pontos internos de x_ort e, ap�s isso, y_ort
    DO j = 1, SIZE(v) - 2
        DO g = 1, num_bic - 1                                        ! O 'g' � o �ndice da bic�ncava. A primeira � a curva_interna. A �ltima ser� calculada quando g = intervalos = num_bic-1 = �ltima bic�ncava.
  
    di = step*g                                                      ! Cada ponto da reta radial pertence a uma bic�ncava conc. Para mudar de bic. basta multiplicar g pelo passo 'step'.
              
    X_Malha_FUEL(g,0) = 0
    Y_Malha_FUEL(g,0) = Y_Malha_FUEL(0,0)+di
    
    X_Malha_FUEL(g,SIZE(v)-1)=X_Malha_FUEL(0,SIZE(v)-1)+di
    Y_Malha_FUEL(g,SIZE(v)-1)=0
              
    m_incli  = (Y_Malha_FUEL(0,j+1) - Y_Malha_FUEL(0,j-1)) / ( X_Malha_FUEL(0,j+1) - X_Malha_FUEL(0,j-1))   ! Coeficiente angular da reta secante (expans�o em S�rie de Taylor).

    Bpp  = Y_Malha_FUEL(0,j) + X_Malha_FUEL(0,j)/m_incli                                                    ! C�culo do coeficiente linear da equa��o da curva y_ortogonal(j) = (-1/m)*x_ortogonal(j) + Bpp.
    
    A_bask = 1.0 + 1.0/m_incli**2                                                   ! � o termo "a" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    B_bask = -2*X_Malha_FUEL(0,j) + 2*Y_Malha_FUEL(0,j)/m_incli - 2*Bpp/m_incli                    ! � o termo "b" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    C_bask = X_Malha_FUEL(0,j)**2 + Y_Malha_FUEL(0,j)**2 - 2*Y_Malha_FUEL(0,j)*Bpp + Bpp**2 - di**2    ! � o termo "c" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    D_bask = B_bask**2 - 4*A_bask*C_bask                                            ! � o termo "delta" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    
    X_Malha1(g,j) = (- B_bask + SQRT(ABS(D_bask)) ) /(2*A_bask)                   ! Tenta a ra�z positiva ...                                      
    X_Malha2(g,j) = (- B_bask - SQRT(ABS(D_bask)) ) /(2*A_bask)                   ! ...faz a outra ra�z.                            
   
         IF (X_Malha1(g,j) > X_Malha2(g,j)) THEN
            X_Malha_maior(g,j) = X_Malha1(g,j)
            X_Malha_menor(g,j) = X_Malha2(g,j)
         ELSE
            X_Malha_menor(g,j) = X_Malha1(g,j)
            X_Malha_maior(g,j) = X_Malha2(g,j)
         ENDIF
                     
         IF(m_incli<0) THEN 
            X_Malha_FUEL(g,j) = X_Malha_maior(g,j) 
        ELSE
            X_Malha_FUEL(g,j) = X_Malha_menor(g,j) 
         ENDIF
            Y_Malha_FUEL(g,j)  = (-1.0/m_incli)*(X_Malha_FUEL(g,j)) + Bpp     
    
        ENDDO
    ENDDO
    !=======================================================================================================================================================
    !----------------------------------------------------------------FUEL: 4o QUADRANTE--------------------------------------------------------------
    !=======================================================================================================================================================
    !1) Escreve-se o 4o QUADRANTE com base na simetria do primeiro;
    !2) Para continuar escrevendo a curva no sentido hor�rio, vou espelhar tanto x quanto y;
    !3) O 1o QUADRANTE tem um pt a mais que � o x=0, ou seja, se o primeiro quadrante � SIZE(v) = 51, o total de pontos � 101 e o segundo tem 50 pts.        
   
    DO i = 0, num_bic-1
        linha_radial     = SIZE(v)-2

        DO j = SIZE(v), 2*(SIZE(v))-2
            X_Malha_FUEL(i,j) =  X_Malha_FUEL(i,linha_radial) 
            Y_Malha_FUEL(i,j) = -Y_Malha_FUEL(i,linha_radial)
            linha_radial     = linha_radial - 1 

        ENDDO
    ENDDO
    
    !===================================================================================================================================================================================
    !-----------------------------------------------------------------------------------CLAD1 + GAP1: 1o QUADRANTE--------------------------------------------------------------------
    !===================================================================================================================================================================================    ! Partindo da curva interna do FUEL para obter a malha do clad1 e do gap1 (ou seja, a primeira curva do FUEL e a primeira a ser desenhado do CLAD1 est�o sobrepostas nessa abordagem).
    
    ! Total de pontos do Clad1+Gap1, como fun��o do n�mero de pontos do FUEL
    ! Calcula a raz�o entre a espessura do FUEL e do CLAD+GAP
    razao_fuel_clad1_float = df/dc1_clad1
    
    ! Arredonda a raz�o para baixo, de forma a superestimar a espessura do CLAD+GAP
    razao_fuel_clad1 = FLOOR(razao_fuel_clad1_float)
    
    ! Calcula quantos pontos deve ter o CLAD+GAP de acordo com quantas vezes menor ele � em rela��o ao FUEL
    C1_float = num_bic/razao_fuel_clad1
    
    ! Arrendonda a raz�o para cima, de forma a superestimar o n�mero de pontos do CLAD+GAP
    C1 = CEILING(C1_float) 
    
        
    ! Finalmente alocando as matrizes que v�o armazenar os valores de x e y do CLAD1+GAP1
    ALLOCATE(X_Malha_CLAD1(0:C1-1, 0:2*(SIZE(v))-2))
    ALLOCATE(Y_Malha_CLAD1(0:C1-1, 0:2*(SIZE(v))-2))

    ! N�o tenho a curva interna do CLAD ainda, mas tenho as duas "quinas" internas, preciso delas porque uso o ponto da frente e o de tr�s das bic�ncavas interna e externa para o c�culo do coeficiente angular da reta secante.
    X_Malha_CLAD1(0, 0) = 0                                              ! O primeiro x da bic. interna do CLAD1 � zero.
    Y_Malha_CLAD1(0, 0) = Y_malha_FUEL(0,0) - dc1                        ! O primeiro y da bic. interna do CLAD1 � apenas a subtra��o do primeiro pt. da bic. interna do FUEL - a espessura do CLAD.
            
    X_Malha_CLAD1(0, SIZE(v)-1) = Y_malha_FUEL(0,SIZE(v) - 1) - dc1      ! O �ltimo x da bic. interna do CLAD1 � apenas a subtra��o do �ltimo y da bic. interna do FUEL - a espessura do CLAD.
    Y_Malha_CLAD1(0, SIZE(v)-1) = 0                                      ! O �ltimo y da bic. interna do CLAD1 � zero. 
    
    intervalos = C1 - 1             ! N�mero de intervalos entre as bic�ncavas que comp�em o CLAD1
    step       = dc1/intervalos     ! O passo entre as bic�ncavas do CLAD1
    
    ! Loop para encontrar os pontos internos de x_ort e, ap�s isso, y_ort
    DO j = 1, SIZE(v) - 2

    indice_bic = 0
    DO g = C1 - 1, 0, -1 ! O 'g' � o �ndice da bic�ncava. A primeira calculada � a parede do GAP, se n�o tiver, � a primeira curva mais interna ap�s o FUEL, da direita para a esquerda

        ! Quando o gap est� sendo considerado, dg � diferente de 0 e a �ltima espessura � na verdade a espessura do GAP
        ! Se dg = 0, ent�o pequenas espessuras 'di' s�o diminu�das da curva interna do FUEL para desenhar a malha do CLAD1
        IF(dg .NE. 0) THEN
            IF (g == C1-1) THEN
                 di = dg 
                indice_bic = indice_bic + 1 ! � o �ndice respons�vel por marcar em qual biconcava est� acontecendo a conta

            ELSE 
                di = dg+step*indice_bic        ! Cada ponto da reta radial pertence a uma bic�ncava conc. Para mudar de bic. basta multiplicar g pelo passo 'step'.
                indice_bic = indice_bic + 1 ! � o �ndice respons�vel por marcar em qual biconcava est� acontecendo a conta

            ENDIF
        ELSE                
                indice_bic = indice_bic + 1 ! � o �ndice respons�vel por marcar em qual biconcava est� acontecendo a conta
                di = step*indice_bic        ! Cada ponto da reta radial pertence a uma bic�ncava conc. Para mudar de bic. basta multiplicar g pelo passo 'step'.
        ENDIF   
    
    X_Malha_CLAD1(g,0)         = 0                          !calculando a reta radial que fica sobreposta ao eixo y
    Y_Malha_CLAD1(g,0)         = Y_Malha_FUEL(0,0) - di     !calculando a reta radial que fica sobreposta ao eixo y

    X_Malha_CLAD1(g,SIZE(v)-1) = X_Malha_FUEL(0,SIZE(v)-1)-di           !calculando a reta radial que fica sobreposta ao eixo x
    Y_Malha_CLAD1(g,SIZE(v)-1) =  0                                     !calculando a reta radial que fica sobreposta ao eixo x
                            
    m_incli  = (Y_Malha_FUEL(0,j+1) - Y_Malha_FUEL(0,j-1)) / ( X_Malha_FUEL(0,j+1) - X_Malha_FUEL(0,j-1))    ! Coeficiente angular da reta secante (expans�o em S�rie de Taylor).
    
    Bpp  = Y_Malha_FUEL(0,j) + X_Malha_FUEL(0,j)/m_incli                                            ! C�culo do coeficiente linear da equa��o da curva y_ortogonal(j) = (-1/m)*x_ortogonal(j) + Bpp.
    
    A_bask = 1.0 + 1.0/m_incli**2                                                                   ! � o termo "a" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    B_bask = -2*X_Malha_FUEL(0,j) + 2*Y_Malha_FUEL(0,j)/m_incli - 2*Bpp/m_incli                     ! � o termo "b" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    C_bask = X_Malha_FUEL(0,j)**2 + Y_Malha_FUEL(0,j)**2 - 2*Y_Malha_FUEL(0,j)*Bpp + Bpp**2 - di**2 ! � o termo "c" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    D_bask = B_bask**2 - 4*A_bask*C_bask                                                            ! � o termo "delta" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    
    X_Malha1(g,j) = (- B_bask + SQRT(ABS(D_bask)) ) /(2*A_bask)                       ! Tenta a ra�z positiva ...                                      
    X_Malha2(g,j) = (- B_bask - SQRT(ABS(D_bask)) ) /(2*A_bask)                       ! ...faz a outra ra�z.                            
   
         IF (X_Malha1(g,j) > X_Malha2(g,j)) THEN
            X_Malha_maior(g,j) = X_Malha1(g,j)
            X_Malha_menor(g,j) = X_Malha2(g,j)
         ELSE
            X_Malha_menor(g,j) = X_Malha1(g,j)
            X_Malha_maior(g,j) = X_Malha2(g,j)
         ENDIF
                     
         IF(m_incli<0) THEN 
            X_Malha_CLAD1(g,j) = X_Malha_menor(g,j) 
        ELSE
            X_Malha_CLAD1(g,j) = X_Malha_maior(g,j) 
         ENDIF
            Y_Malha_CLAD1(g,j) = (-1.0/m_incli)*(X_Malha_CLAD1(g,j)) + Bpp     

        ENDDO
    ENDDO
    !=======================================================================================================================================================
    !----------------------------------------------------------------CLAD1+GAP1: 4o QUADRANTE--------------------------------------------------------------
    !=======================================================================================================================================================
    !1) Escreve-se o 4o QUADRANTE com base na simetria do primeiro;
    !2) Para continuar escrevendo a curva no sentido anti-hor�rio, vou espelhar tanto x quanto y;
    !3) O 1o QUADRANTE tem um pt a mais que � o x=0, ou seja, se o primeiro quadrante � SIZE(v) = 51, o total de pontos � 101 e o segundo tem 50 pts.        

    DO i = 0, C1-1
        linha_radial = SIZE(v)-2
        
        DO j = SIZE(v), 2*(SIZE(v))-2
            
            X_Malha_CLAD1(i,j) =  X_Malha_CLAD1(i,linha_radial) 
        
            Y_Malha_CLAD1(i,j) = -Y_Malha_CLAD1(i,linha_radial)
            linha_radial     = linha_radial - 1

        ENDDO
    ENDDO
    !===================================================================================================================================================================================
    !---------------------------------------------------------------------------------------CLAD2 + GAP2-------------------------------------------------------------------------------
    !===================================================================================================================================================================================    ! Partindo da curva externa do FUEL para obter a malha do clad2 e do gap2 (ou seja, a �ltima curva do FUEL e a primeira a ser desenhado do CLAD2 est�o sobrepostas nessa abordagem).
    
    ! Total de pontos do Clad2+Gap2, como fun��o do n�mero de pontos do FUEL
    !Calcula a raz�o entre a espessura do FUEL e do CLAD (a espessura do GAP � muito pequena)
    razao_fuel_clad2_float = df/dc2_clad2
    
    ! Arredonda a raz�o para baixo, de forma a superestimar a espessura do CLAD
    razao_fuel_clad2 = FLOOR(razao_fuel_clad2_float)
    
    ! Calcula quantos pontos deve ter o CLAD de acordo com quantas vezes menor ele � em rela��o ao FUEL
    C2_float = num_bic/razao_fuel_clad2
    
    ! Arrendonda a raz�o para cima, de forma a superestimar o n�mero de pontos do CLAD
    C2 = CEILING(C2_float)
   
    !Finalmente alocando as matrizes que v�o armazenar os valores de x e y do CLAD2
    ALLOCATE(X_Malha_CLAD2(0:C2-1, 0:2*(SIZE(v))-2))
    ALLOCATE(Y_Malha_CLAD2(0:C2-1, 0:2*(SIZE(v))-2))
    
   ! Preenchendo a primeira bic�ncava do CLAD2 com o �ltima curva do FUEL
   !! X_Malha_CLAD2(0,:) = X_Malha_FUEL(num_bic-1,:) ! Toda a coluna 0 � prenchida pela bic�ncava definida pela curva externa do FUEL (primeira bic�ncava conc�ntrica do CLAD2)
    !!Y_Malha_CLAD2(0,:) = Y_Malha_FUEL(num_bic-1,:) ! Toda a coluna 0 � prenchida pela bic�ncava definida pela curva externa do FUEL (primeira bic�ncava conc�ntrica do CLAD2)

    ! N�o tenho a curva externa ainda, mas tenho as duas "quinas" externas, preciso deles porque uso o ponto da frente e o de tr�s das bic�ncavas interna e externa para o c�culo do coeficiente angular da reta secante.
    X_Malha_CLAD2(C2-1, 0) = 0                                                 ! O primeiro x da bic. externa � zero.
    Y_Malha_CLAD2(C2-1, 0) = Y_Malha_FUEL(num_bic-1,0) + dc2                   ! O primeiro y da bic. externa � apenas a soma do primeiro pt. da bic. externa do FUEL + a espessura do CLAD..
            
    X_Malha_CLAD2(C2-1, SIZE(v)-1) = X_Malha_FUEL(num_bic-1,SIZE(v) - 1) + dc2  ! O �ltimo x da bic. externa � apenas a soma do �ltimo y da bic. externa do FUEL + a espessura.
    Y_Malha_CLAD2(C2-1, SIZE(v)-1) = 0                                          ! O �ltimo y da bic. externa � zero.
   
    intervalos = C2 - 1                ! N�mero total de intervalos entre as bic�ncavas conc�ntricas.
    step    = dc2/(intervalos)         ! Passo entre os pontos das retas radiais (i. e. bic�ncavas), onde d � a dist�ncia entre a primeira e a �ltima bic�ncava (curva interna e externa).
    
    ! Loop para encontrar os pontos internos de x_ort e, ap�s isso, y_ort
    DO j = 1, SIZE(v) - 2
        indice_bic=0
        DO g = 0, C2 - 1                                                ! O 'g' � o �ndice da bic�ncava.   
        
        ! Quando o gap est� sendo considerado, dg � diferente de 0 e a primeira espessura � na verdade a espessura do GAP
        ! Se dg = 0, ent�o pequenas espessuras 'di' s�o diminu�das da curva interna para desenhar a malha do CLAD1
        IF(dg .NE. 0.) THEN
            IF (g == 0) THEN
                 di = dg 
                indice_bic = indice_bic + 1 !Para n�o ter "pulo" de espa�o entre parede gap e primeiro ponto do CLAD
            ELSE
                di = dg+step*indice_bic
                indice_bic = indice_bic + 1 !Para n�o ter "pulo" de espa�o entre parede gap e primeiro ponto do CLAD
            ENDIF
        ELSE                     
                indice_bic = indice_bic + 1 ! � o �ndice respons�vel por marcar em qual biconcava est� acontecendo a conta           
                di = step*indice_bic ! Cada ponto da reta radial pertence a uma bic�ncava conc. Para mudar de bic. basta multiplicar indice_bic pelo passo 'step'.
        ENDIF                                                            
       
    X_Malha_CLAD2(g,0)         = 0                                          !calculando a reta radial que fica sobreposta ao eixo y
    Y_Malha_CLAD2(g,0)         = Y_Malha_FUEL(num_bic-1,0) + di             !calculando a reta radial que fica sobreposta ao eixo y

    X_Malha_CLAD2(g,SIZE(v)-1) = X_Malha_FUEL(num_bic-1,SIZE(v)-1)+di       !calculando a reta radial que fica sobreposta ao eixo x    
    Y_Malha_CLAD2(g,SIZE(v)-1) = 0                                          !calculando a reta radial que fica sobreposta ao eixo x
                   
    m_incli  = (Y_Malha_FUEL(num_bic-1,j+1) - Y_Malha_FUEL(num_bic-1,j-1)) / ( X_Malha_FUEL(num_bic-1,j+1) - X_Malha_FUEL(num_bic-1,j-1))   ! Coeficiente angular da reta secante (expans�o em S�rie de Taylor).
    
    Bpp  = Y_Malha_FUEL(num_bic-1,j) + X_Malha_FUEL(num_bic-1,j)/m_incli                                                     ! C�culo do coeficiente linear da equa��o da curva y_ortogonal(j) = (-1/m)*x_ortogonal(j) + Bpp.
    
    A_bask = 1.0 + 1.0/m_incli**2                                                                                            ! � o termo "a" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    B_bask = -2.0*X_Malha_FUEL(num_bic-1,j) + 2*Y_Malha_FUEL(num_bic-1,j)/m_incli - 2*Bpp/m_incli                            ! � o termo "b" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    C_bask = X_Malha_FUEL(num_bic-1,j)**2 + Y_Malha_FUEL(num_bic-1,j)**2 - 2*Y_Malha_FUEL(num_bic-1,j)*Bpp + Bpp**2 - di**2  ! � o termo "c" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    D_bask = B_bask**2 - 4*A_bask*C_bask                                                ! � o termo "delta" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    
    X_Malha1(g,j) = (- B_bask + SQRT(ABS(D_bask)) ) /(2*A_bask)                        ! Tenta a ra�z positiva ...                                      
    X_Malha2(g,j) = (- B_bask - SQRT(ABS(D_bask)) ) /(2*A_bask)                        ! ...faz a outra ra�z.                            
    
         IF (X_Malha1(g,j) > X_Malha2(g,j)) THEN
            X_Malha_maior(g,j) = X_Malha1(g,j)
            X_Malha_menor(g,j) = X_Malha2(g,j)
         ELSE
            X_Malha_menor(g,j) = X_Malha1(g,j)
            X_Malha_maior(g,j) = X_Malha2(g,j)
         ENDIF
                     
         IF(m_incli<0) THEN 
            X_Malha_CLAD2(g,j) = X_Malha_maior(g,j) 
        ELSE
            X_Malha_CLAD2(g,j) = X_Malha_menor(g,j) 
         ENDIF
            Y_Malha_CLAD2(g,j)  = (-1.0/m_incli)*(X_Malha_CLAD2(g,j)) + Bpp     
    
        ENDDO
    ENDDO
    !=======================================================================================================================================================
    !----------------------------------------------------------------CLAD2+GAP2: 4o QUADRANTE--------------------------------------------------------------
    !=======================================================================================================================================================
    !1) Escreve-se o 4o QUADRANTE com base na simetria do primeiro;
    !2) Para continuar escrevendo a curva no sentido anti-hor�rio, vou espelhar tanto x quanto y;
    !3) O 1o QUADRANTE tem um pt a mais que � o x=0, ou seja, se o primeiro quadrante � SIZE(v) = 51, o total de pontos � 101 e o segundo tem 50 pts.        
     
    DO i = 0, C2-1
        linha_radial     = SIZE(v) - 2
        
     DO j = SIZE(v), 2*(SIZE(v))-2
            X_Malha_CLAD2(i,j) =  X_Malha_CLAD2(i,linha_radial) 
            Y_Malha_CLAD2(i,j) = -Y_Malha_CLAD2(i,linha_radial)
            linha_radial     = linha_radial - 1

        ENDDO
    ENDDO    

    !=======================================================================================================================================================
    !----------------------------------------------------------------FORMANDO A MATRIZ FINAL--------------------------------------------------------------
    !=======================================================================================================================================================
    
    ! Definindo o n�mero de bic�ncavas e renomeando para melhor compreens�o nas subrotinas seguintes
    biconcavas = C1 + num_bic + C2 ! Total de bic�ncavas (em n�mero de pontos da parte interna para a externa)
    radiais    = 2*SIZE(v)-1         ! Total de retas radiais (em n�mero pontos longitudinais)
    
    F1      = C1+1                 ! Se C1=15, CLAD1 vai de 0:14 (=0:C1-1) e F1 = C1+1 pts. F1 � LIMITE, C1 � N�MERO DE PONTOS!
    F2      = F1 + num_bic-1       ! F2 � o LIMITE EXTERNO do FUEL, logo FUEL come�a em F1 (=C1+1=ao ponto 16) e anda num_bic pts (F1:num_bic-1), num_bic-1 porque o F1 EST� INCLUSO!
    Ifinal  = F2+1 + C2-1          ! J� o CLAD2 come�a no ponto a frente do final do FUEL e termina C2-1 pts depois, J� QUE O F2+1 EST� INLCUSO tenho que tirar um ponto (como no caso do ZERO...)
 
    ALLOCATE(X(1:biconcavas, 1:radiais))
    ALLOCATE(Y(1:biconcavas, 1:radiais))

    ! Concatenando a matriz final X 
    col = 0
    DO k = 1, C1 
        X(k,:) =X_Malha_CLAD1(col,:)
        col = col + 1
    ENDDO
   
     col = 0
     DO p = F1, F2
        X(p,:) = X_Malha_FUEL(col,:)
        col = col + 1
     ENDDO
     
     col=0
     DO q = F2+1, Ifinal
        X(q,:) =X_Malha_CLAD2(col,:)
        col = col+1
     ENDDO
    
     ! Concatenando a matriz final Y
    col = 0
    DO k = 1, C1
        Y(k,:) =Y_Malha_CLAD1(col,:)
         col = col + 1
    ENDDO
   
     col = 0
     DO p = F1, F2
        Y(p,:) = Y_Malha_FUEL(col,:)
         col = col + 1
     ENDDO
     
     col=0
     DO q = F2+1, Ifinal
        Y(q,:) =Y_Malha_CLAD2(col,:)
         col = col + 1
     ENDDO

     WRITE(*,*) BICONCAVAS
     WRITE(*,*) RADIAIS
     
    !Escreve a geometria malhada em geometria.txt, j� formatado para plot em Scilab
200 FORMAT(32(f14.8),/)
    WRITE(1,*) 'x=['
    WRITE(1,200) ((X(i,j),i=1,biconcavas),j=1,radiais)
    WRITE(1,*) ']'

    WRITE(1,*) '                                                                    '
    WRITE(1,*) '                                                                    '
    WRITE(1,*) 'y=['
    WRITE(1,200) ((Y(i,j),i=1,biconcavas),j=1,radiais)
    WRITE(1,*) ']'
    WRITE(1,*) 'plot(x,y)'
    
    CLOSE(UNIT=1)
    CLOSE(UNIT=200)
    
    END SUBROUTINE