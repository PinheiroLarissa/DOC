SUBROUTINE MetodoSecante
    
    USE double_Precision
    USE dadosGerais
    
    IMPLICIT NONE 
    INTEGER            ::j,g,i,indice_bic, k, q, p, col                          ! Contadores utilizados nos loops de concatenamento das matrizes de malha
    
    OPEN(UNIT = 1,   FILE = 'geometria.txt', STATUS = 'unknown') ! Esse arquivo é usado tanto pela SUBROUTINE Cassini_geo quanto pela SUBROUTINE Eliptico_geo, dependendo de qual seja chamada pelo principal.

    
 !***********************************************************************************************************************************************************************************
    !------------------------------------------------------------------ MÉTODO DA SECANTE: ENCONTRAM-SE X e Y A CADA RETA ORTOGONAL ----------------------------------------------------
    !***********************************************************************************************************************************************************************************
    ! 1. Para verificar os cálculos, pg. 128 - 129 caderno.
    ! 2. A curva é feita no sentido anti-horário;
    ! 3. A curva interna, gerada pela parametrização, é salva na primeira coluna das matrizes X e Y, que serão usadas nos cálculos de geração de malha por Thompson;
    ! 4. A curva externa ainda não é conhecida, apenas as quinas externas. Essas quinas também são salvas nas matrizes X e Y;
    ! 5. Loop em j: anda ao longo da bicôncava. Loop em g: anda ao longo da reta ortogonal (ou radial), desenhando a malha da primeira bicôncava até a última;
    ! 6. Assim, as colunas das matrizes X e Y são bicôncavas concêntricas, e as linhas são as retas radiais (desenhadas como retas ortogonais às retas secantes da curva interna);
    ! 7. Os pontos ao longo de uma reta radial são traçados acrescentando a espessura "di" à curva interna, sendo di uma fração incremental da espessura total "d".
       
        ! Combinando a espessura do CLAD e do GAP para facilitar os cálculos, e se dg = 0, o primeiro ponto do GAP é bypassado
        dc1 = dc1_clad1 + dg 
        dc2 = dc2_clad2 + dg

    !===================================================================================================================================================================================
    !--------------------------------------------------------------------------------------------FUEL-----------------------------------------------------------------------------------
    !===================================================================================================================================================================================
    ALLOCATE(X_Malha_FUEL(0:num_bic-1, 0:2*(SIZE(v) -1)))
    ALLOCATE(Y_Malha_FUEL(0:num_bic-1, 0:2*(SIZE(v) -1)))
    
    ! Coordenadas iniciais de X_MALHA_FUEL e Y_MALHA_FUEL    
    X_Malha_FUEL(0,:) = x_curva_interna(:) ! Toda a coluna 0 é prenchida pela bicôncava definida pela curva interna (primeira bicôncava concêntrica)
    Y_Malha_FUEL(0,:) = y_curva_interna(:) ! Toda a coluna 0 é prenchida pela bicôncava definida pela curva interna (primeira bicôncava concêntrica)
  ! WRITE(*,*) x_curva_interna
   !   WRITE(*,*) Y_curva_interna

   
    
    ! Não tenho a curva externa ainda, mas tenho as duas "quinas" externas, preciso deles porque uso o ponto da frente e o de trás das bicôncavas interna e externa para o cáculo do coeficiente angular da reta secante.
    X_Malha_FUEL(num_bic-1, 0) = x_curva_interna(0) + df                    ! O primeiro x da bic. externa é apenas a soma do primeiro pt. da bic. interna + a espessura.
    Y_Malha_FUEL(num_bic-1, 0) = 0                                          ! O primeiro y da bic. externa é zero.
            
    X_Malha_FUEL(num_bic-1, SIZE(v) - 1) = 0                                ! O último x da bic. externa é zero.
    Y_Malha_FUEL(num_bic-1, SIZE(v) - 1) = y_curva_interna(SIZE(v) - 1) + df ! O último y da bic. externa é apenas a soma do último y da bic. interna + a espessura.
   
    intervalos = num_bic - 1                ! Número total de intervalos entre as bicôncavas concêntricas.
    step    = df/(intervalos)               ! Passo entre os pontos das retas radiais (i. e. bicôncavas), onde d é a distância entre a primeira e a última bicôncava (curva interna e externa).
    ! Loop para encontrar os pontos internos de x_ort e, após isso, y_ort
    DO j = 1, SIZE(v) - 2
        DO g = 1, num_bic - 1                                        ! O 'g' é o índice da bicôncava. A primeira é a curva_interna. A última será calculada quando g = intervalos = num_bic-1 = última bicôncava.
  
    di = step*g                                                      ! Cada ponto da reta radial pertence a uma bicôncava conc. Para mudar de bic. basta multiplicar g pelo passo 'step'.
              
    X_Malha_FUEL(g,0) = X_Malha_FUEL(0,0)+di
    Y_Malha_FUEL(g,SIZE(v)-1)=Y_Malha_FUEL(0,SIZE(v)-1)+di
              
    m_incli  = (Y_Malha_FUEL(0,j+1) - Y_Malha_FUEL(0,j-1)) / ( X_Malha_FUEL(0,j+1) - X_Malha_FUEL(0,j-1))   ! Coeficiente angular da reta secante (expansão em Série de Taylor).
    
    Bpp  = Y_Malha_FUEL(0,j) + X_Malha_FUEL(0,j)/m_incli                                                    ! Cáculo do coeficiente linear da equação da curva y_ortogonal(j) = (-1/m)*x_ortogonal(j) + Bpp.
    
    A_bask = 1.0 + 1.0/m_incli**2                                                   ! É o termo "a" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    B_bask = -2*X_Malha_FUEL(0,j) + 2*Y_Malha_FUEL(0,j)/m_incli - 2*Bpp/m_incli                    ! É o termo "b" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    C_bask = X_Malha_FUEL(0,j)**2 + Y_Malha_FUEL(0,j)**2 - 2*Y_Malha_FUEL(0,j)*Bpp + Bpp**2 - di**2    ! É o termo "c" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    D_bask = B_bask**2 - 4*A_bask*C_bask                                            ! É o termo "delta" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    
    X_Malha1(g,j) = (- B_bask + SQRT(ABS(D_bask)) ) /(2*A_bask)                   ! Tenta a raíz positiva ...                                      
    X_Malha2(g,j) = (- B_bask - SQRT(ABS(D_bask)) ) /(2*A_bask)                   ! ...faz a outra raíz.                            
   
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
    !2) Para continuar escrevendo a curva no sentido anti-horário, vou espelhar tanto x quanto y;
    !3) O 1o QUADRANTE tem um pt a mais que é o x=0, ou seja, se o primeiro quadrante é SIZE(v) = 51, o total de pontos é 101 e o segundo tem 50 pts.        
   
    DO i = 0, num_bic-1
        linha_radial     = SIZE(v)-1

        DO j = SIZE(v), 2*(SIZE(v)-1)
            X_Malha_FUEL(i,j) =  X_Malha_FUEL(i,linha_radial) 
            Y_Malha_FUEL(i,j) = -Y_Malha_FUEL(i,linha_radial)
            linha_radial     = linha_radial - 1 

        ENDDO
    ENDDO
    
    !===================================================================================================================================================================================
    !-----------------------------------------------------------------------------------CLAD1 + GAP1: 1o QUADRANTE--------------------------------------------------------------------
    !===================================================================================================================================================================================    ! Partindo da curva interna do FUEL para obter a malha do clad1 e do gap1 (ou seja, a primeira curva do FUEL e a primeira a ser desenhado do CLAD1 estão sobrepostas nessa abordagem).
    
    ! Total de pontos do Clad1+Gap1, como função do número de pontos do FUEL
    ! Calcula a razão entre a espessura do FUEL e do CLAD+GAP
    razao_fuel_clad1_float = df/dc1_clad1
    
    ! Arredonda a razão para baixo, de forma a superestimar a espessura do CLAD+GAP
    razao_fuel_clad1 = FLOOR(razao_fuel_clad1_float)
    
    ! Calcula quantos pontos deve ter o CLAD+GAP de acordo com quantas vezes menor ele é em relação ao FUEL
    C1_float = num_bic/razao_fuel_clad1
    
    ! Arrendonda a razão para cima, de forma a superestimar o número de pontos do CLAD+GAP
    C1 = CEILING(C1_float) 
    
        
    ! Finalmente alocando as matrizes que vão armazenar os valores de x e y do CLAD1+GAP1
    ALLOCATE(X_Malha_CLAD1(0:C1-1, 0:2*(SIZE(v) -1)))
    ALLOCATE(Y_Malha_CLAD1(0:C1-1, 0:2*(SIZE(v) -1)))

    ! Não tenho a curva interna do CLAD ainda, mas tenho as duas "quinas" internas, preciso delas porque uso o ponto da frente e o de trás das bicôncavas interna e externa para o cáculo do coeficiente angular da reta secante.
    X_Malha_CLAD1(0, 0) = X_malha_FUEL(0,0) - dc1                ! O primeiro x da bic. interna do CLAD é apenas a subtração do primeiro pt. da bic. interna do FUEL - a espessura do CLAD.
    Y_Malha_CLAD1(0, 0) = 0                                          ! O primeiro y da bic. interna do CLAD é zero.
            
    X_Malha_CLAD1(0, SIZE(v)-1) = 0                                  ! O último x da bic. interna do CLAD é zero.
    Y_Malha_CLAD1(0, SIZE(v)-1) = Y_malha_FUEL(0,SIZE(v) - 1) - dc1 ! O último y da bic. interna do CLAD é apenas a subtração do último y da bic. interna do FUEL - a espessura do CLAD.
    
    intervalos = C1 - 1             ! Número de intervalos entre as bicôncavas que compõem o CLAD1
    step       = dc1/intervalos     ! O passo entre as bicôncavas do CLAD1
    
    ! Loop para encontrar os pontos internos de x_ort e, após isso, y_ort
    DO j = 1, SIZE(v) - 2

    indice_bic = 0
    DO g = C1 - 1, 0, -1 ! O 'g' é o índice da bicôncava. A primeira calculada é a parede do GAP, se não tiver, é a primeira curva mais interna após o FUEL, da direita para a esquerda

        ! Quando o gap está sendo considerado, dg é diferente de 0 e a última espessura é na verdade a espessura do GAP
        ! Se dg = 0, então pequenas espessuras 'di' são diminuídas da curva interna do FUEL para desenhar a malha do CLAD1
        IF(dg .NE. 0) THEN
            IF (g == C1-1) THEN
                 di = dg 
                indice_bic = indice_bic + 1 ! é o índice responsável por marcar em qual biconcava está acontecendo a conta

            ELSE 
                di = dg+step*indice_bic        ! Cada ponto da reta radial pertence a uma bicôncava conc. Para mudar de bic. basta multiplicar g pelo passo 'step'.
                indice_bic = indice_bic + 1 ! é o índice responsável por marcar em qual biconcava está acontecendo a conta

            ENDIF
        ELSE                
                indice_bic = indice_bic + 1 ! é o índice responsável por marcar em qual biconcava está acontecendo a conta
                di = step*indice_bic        ! Cada ponto da reta radial pertence a uma bicôncava conc. Para mudar de bic. basta multiplicar g pelo passo 'step'.
        ENDIF   
    
    X_Malha_CLAD1(g,0)                = X_Malha_FUEL(0,0) - di       !calculando a reta radial que fica sobreposta ao eixo x
    Y_Malha_CLAD1(g,0)                 = 0                           !calculando a reta radial que fica sobreposta ao eixo x

    X_Malha_CLAD1(g,SIZE(v)-1)         = 0                           !calculando a reta radial que fica sobreposta ao eixo y
    Y_Malha_CLAD1(g,SIZE(v)-1) = Y_Malha_FUEL(0,SIZE(v)-1)-di        !calculando a reta radial que fica sobreposta ao eixo y
                            
    m_incli  = (Y_Malha_FUEL(0,j+1) - Y_Malha_FUEL(0,j-1)) / ( X_Malha_FUEL(0,j+1) - X_Malha_FUEL(0,j-1))    ! Coeficiente angular da reta secante (expansão em Série de Taylor).
    
    Bpp  = Y_Malha_FUEL(0,j) + X_Malha_FUEL(0,j)/m_incli                                            ! Cáculo do coeficiente linear da equação da curva y_ortogonal(j) = (-1/m)*x_ortogonal(j) + Bpp.
    
    A_bask = 1.0 + 1.0/m_incli**2                                                                   ! É o termo "a" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    B_bask = -2*X_Malha_FUEL(0,j) + 2*Y_Malha_FUEL(0,j)/m_incli - 2*Bpp/m_incli                     ! É o termo "b" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    C_bask = X_Malha_FUEL(0,j)**2 + Y_Malha_FUEL(0,j)**2 - 2*Y_Malha_FUEL(0,j)*Bpp + Bpp**2 - di**2 ! É o termo "c" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    D_bask = B_bask**2 - 4*A_bask*C_bask                                                            ! É o termo "delta" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    
    X_Malha1(g,j) = (- B_bask + SQRT(ABS(D_bask)) ) /(2*A_bask)                       ! Tenta a raíz positiva ...                                      
    X_Malha2(g,j) = (- B_bask - SQRT(ABS(D_bask)) ) /(2*A_bask)                       ! ...faz a outra raíz.                            
   
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
    !2) Para continuar escrevendo a curva no sentido anti-horário, vou espelhar tanto x quanto y;
    !3) O 1o QUADRANTE tem um pt a mais que é o x=0, ou seja, se o primeiro quadrante é SIZE(v) = 51, o total de pontos é 101 e o segundo tem 50 pts.        

    DO i = 0, C1-1
        linha_radial = SIZE(v)-1
        
        DO j = SIZE(v), 2*(SIZE(v)-1)
            
            X_Malha_CLAD1(i,j) =  X_Malha_CLAD1(i,linha_radial) 
        
            Y_Malha_CLAD1(i,j) = -Y_Malha_CLAD1(i,linha_radial)
            linha_radial     = linha_radial - 1

        ENDDO
    ENDDO
    !===================================================================================================================================================================================
    !---------------------------------------------------------------------------------------CLAD2 + GAP2-------------------------------------------------------------------------------
    !===================================================================================================================================================================================    ! Partindo da curva externa do FUEL para obter a malha do clad2 e do gap2 (ou seja, a última curva do FUEL e a primeira a ser desenhado do CLAD2 estão sobrepostas nessa abordagem).
    
    ! Total de pontos do Clad2+Gap2, como função do número de pontos do FUEL
    !Calcula a razão entre a espessura do FUEL e do CLAD (a espessura do GAP é muito pequena)
    razao_fuel_clad2_float = df/dc2_clad2
    
    ! Arredonda a razão para baixo, de forma a superestimar a espessura do CLAD
    razao_fuel_clad2 = FLOOR(razao_fuel_clad2_float)
    
    ! Calcula quantos pontos deve ter o CLAD de acordo com quantas vezes menor ele é em relação ao FUEL
    C2_float = num_bic/razao_fuel_clad2
    
    ! Arrendonda a razão para cima, de forma a superestimar o número de pontos do CLAD
    C2 = CEILING(C2_float)
   
    !Finalmente alocando as matrizes que vão armazenar os valores de x e y do CLAD2
    ALLOCATE(X_Malha_CLAD2(0:C2-1, 0:2*(SIZE(v) -1)))
    ALLOCATE(Y_Malha_CLAD2(0:C2-1, 0:2*(SIZE(v) -1)))
    
   ! Preenchendo a primeira bicôncava do CLAD2 com o última curva do FUEL
   !! X_Malha_CLAD2(0,:) = X_Malha_FUEL(num_bic-1,:) ! Toda a coluna 0 é prenchida pela bicôncava definida pela curva externa do FUEL (primeira bicôncava concêntrica do CLAD2)
    !!Y_Malha_CLAD2(0,:) = Y_Malha_FUEL(num_bic-1,:) ! Toda a coluna 0 é prenchida pela bicôncava definida pela curva externa do FUEL (primeira bicôncava concêntrica do CLAD2)

    ! Não tenho a curva externa ainda, mas tenho as duas "quinas" externas, preciso deles porque uso o ponto da frente e o de trás das bicôncavas interna e externa para o cáculo do coeficiente angular da reta secante.
    X_Malha_CLAD2(C2-1, 0) = X_Malha_FUEL(num_bic-1,0) + dc2               ! O primeiro x da bic. externa é apenas a soma do primeiro pt. da bic. externa do FUEL + a espessura do CLAD.
    Y_Malha_CLAD2(C2-1, 0) = 0                                      ! O primeiro y da bic. externa é zero.
            
    X_Malha_CLAD2(C2-1, SIZE(v)-1) = 0                                ! O último x da bic. externa é zero.
    Y_Malha_CLAD2(C2-1, SIZE(v)-1) = Y_Malha_FUEL(num_bic-1,SIZE(v) - 1) + dc2 ! O último y da bic. externa é apenas a soma do último y da bic. externa do FUEL + a espessura.
   
    intervalos = C2 - 1                ! Número total de intervalos entre as bicôncavas concêntricas.
    step    = dc2/(intervalos)         ! Passo entre os pontos das retas radiais (i. e. bicôncavas), onde d é a distância entre a primeira e a última bicôncava (curva interna e externa).
    
    ! Loop para encontrar os pontos internos de x_ort e, após isso, y_ort
    DO j = 1, SIZE(v) - 2
        indice_bic=0
        DO g = 0, C2 - 1                                                ! O 'g' é o índice da bicôncava.   
        
        ! Quando o gap está sendo considerado, dg é diferente de 0 e a primeira espessura é na verdade a espessura do GAP
        ! Se dg = 0, então pequenas espessuras 'di' são diminuídas da curva interna para desenhar a malha do CLAD1
        IF(dg .NE. 0.) THEN
            IF (g == 0) THEN
                 di = dg 
                indice_bic = indice_bic + 1 !Para não ter "pulo" de espaço entre parede gap e primeiro ponto do CLAD
            ELSE
                di = dg+step*indice_bic
                indice_bic = indice_bic + 1 !Para não ter "pulo" de espaço entre parede gap e primeiro ponto do CLAD
            ENDIF
        ELSE                     
                indice_bic = indice_bic + 1 ! é o índice responsável por marcar em qual biconcava está acontecendo a conta           
                di = step*indice_bic ! Cada ponto da reta radial pertence a uma bicôncava conc. Para mudar de bic. basta multiplicar indice_bic pelo passo 'step'.
        ENDIF                                                            
       
    X_Malha_CLAD2(g,0)         = X_Malha_FUEL(num_bic-1,0) + di          !calculando a reta radial que fica sobreposta ao eixo x
    Y_Malha_CLAD2(g,0)         = 0                               !calculando a reta radial que fica sobreposta ao eixo x

    X_Malha_CLAD2(g,SIZE(v)-1) = 0                               !calculando a reta radial que fica sobreposta ao eixo y
    Y_Malha_CLAD2(g,SIZE(v)-1) = Y_Malha_FUEL(num_bic-1,SIZE(v)-1)+di           !calculando a reta radial que fica sobreposta ao eixo y
                   
    m_incli  = (Y_Malha_FUEL(num_bic-1,j+1) - Y_Malha_FUEL(num_bic-1,j-1)) / ( X_Malha_FUEL(num_bic-1,j+1) - X_Malha_FUEL(num_bic-1,j-1))   ! Coeficiente angular da reta secante (expansão em Série de Taylor).
    
    Bpp  = Y_Malha_FUEL(num_bic-1,j) + X_Malha_FUEL(num_bic-1,j)/m_incli                                                     ! Cáculo do coeficiente linear da equação da curva y_ortogonal(j) = (-1/m)*x_ortogonal(j) + Bpp.
    
    A_bask = 1.0 + 1.0/m_incli**2                                                             ! É o termo "a" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    B_bask = -2.0*X_Malha_FUEL(num_bic-1,j) + 2*Y_Malha_FUEL(num_bic-1,j)/m_incli - 2*Bpp/m_incli                            ! É o termo "b" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    C_bask = X_Malha_FUEL(num_bic-1,j)**2 + Y_Malha_FUEL(num_bic-1,j)**2 - 2*Y_Malha_FUEL(num_bic-1,j)*Bpp + Bpp**2 - di**2    ! É o termo "c" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    D_bask = B_bask**2 - 4*A_bask*C_bask                                                ! É o termo "delta" da equação de Bháskara para achar o x ortogonal (na mesma reta radial).
    
    X_Malha1(g,j) = (- B_bask + SQRT(ABS(D_bask)) ) /(2*A_bask)                        ! Tenta a raíz positiva ...                                      
    X_Malha2(g,j) = (- B_bask - SQRT(ABS(D_bask)) ) /(2*A_bask)                        ! ...faz a outra raíz.                            
    
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
    !2) Para continuar escrevendo a curva no sentido anti-horário, vou espelhar tanto x quanto y;
    !3) O 1o QUADRANTE tem um pt a mais que é o x=0, ou seja, se o primeiro quadrante é SIZE(v) = 51, o total de pontos é 101 e o segundo tem 50 pts.        
     
    DO i = 0, C2-1
        linha_radial     = SIZE(v) - 1
        
     DO j = SIZE(v), 2*(SIZE(v)-1)
            X_Malha_CLAD2(i,j) =  X_Malha_CLAD2(i,linha_radial) 
            Y_Malha_CLAD2(i,j) = -Y_Malha_CLAD2(i,linha_radial)
            linha_radial     = linha_radial - 1

        ENDDO
    ENDDO    

    !=======================================================================================================================================================
    !----------------------------------------------------------------FORMANDO A MATRIZ FINAL--------------------------------------------------------------
    !=======================================================================================================================================================
    
    ! Definindo o número de bicôncavas e renomeando para melhor compreensão nas subrotinas seguintes
    biconcavas = C1 + num_bic + C2 ! Total de bicôncavas (em número de pontos da parte interna para a externa)
    radiais    = 2*SIZE(v)-1       ! Total de retas radiais (em número pontos longitudinais)
    
    F1      = C1+1                 ! Se C1=15, CLAD1 vai de 0:14 (=0:C1-1) e F1 = C1+1 pts. F1 é LIMITE, C1 é NÚMERO DE PONTOS!
    F2      = F1 + num_bic-1       ! F2 é o LIMITE EXTERNO do FUEL, logo FUEL começa em F1 (=C1+1=ao ponto 16) e anda num_bic pts (F1:num_bic-1), num_bic-1 porque o F1 ESTÁ INCLUSO!
    Ifinal  = F2+1 + C2-1          ! Já o CLAD2 começa no ponto a frente do final do FUEL e termina C2-1 pts depois, JÁ QUE O F2+1 ESTÁ INLCUSO tenho que tirar um ponto (como no caso do ZERO...)
 
    ALLOCATE(X(1:biconcavas, 1:2*SIZE(v)-1))
    ALLOCATE(Y(1:biconcavas, 1:2*SIZE(v)-1))

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

    !Escreve a geometria malhada em geometria.txt, já formatado para plot em Scilab
200 FORMAT(75(f14.8),/)
    WRITE(1,*) 'x=['
    WRITE(1,200) X
    WRITE(1,*) ']'

    WRITE(1,*) '                                                                    '
    WRITE(1,*) '                                                                    '
    WRITE(1,*) 'y=['
    WRITE(1,200) Y
    WRITE(1,*) ']'
    WRITE(1,*) 'plot(x,y)'
    
    CLOSE(UNIT=1)
    END SUBROUTINE MetodoSecante