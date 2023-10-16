 SUBROUTINE Ellipk_geo
    !*********************************************************************************************************************************************************************************************************************   
    ! Programador: Larissa
    ! �ltima revis�o: 19/11/2018
    !
    !ATEN��O2: O ALGORITMO USADO FICA INST�VEL QUADO U=K E V=K', POR ISSO ESSE VALORES DEVEM SER TIRADOS DA TABELA EM: https://dlmf.nist.gov/22.20#E1 OU ABRAMOWITZ, PG571
    !
    ! Descri��o da subrotina Ellipk:
    ! 1) Calcula os par�metros da m�dia geom�trica-aritm�tica a e b (ver ABRAMOWITZ, pg. 571);
    ! 2) Calcula as integaris KK (integral el�ptica completa de primeiro tipo) e KKp (integral el�ptica completa de primeiro tipo complementar);
    ! 3) Define os vetores u e v como m�ltiplos de KK e KKp ( de forma a variarem de 0 a KK (ou KKp)) e sn, cn e dn variarem de 0 a 1;
    ! 4) C�culo da amplitude PHI, obtido de forma regressiva de PHI(N) at� PHI(0) (ver ABRAMOWITZ, pg. 571);
    ! 5) C�culo de sn, cn, e dn para o m dado como entrada e os vetores u e v definidos. Para alguns valores de m, u e v o algoritmo fica inst�vel e usam-se valores tabelados;
    !    Quando m = 0 ou m = 1   ---> sn, cn e dn valores tabelados;
    !    Quando u = 0 ou u = KK  ---> sn, cn e dn valores tabelados;
    !    Quando v = 0 ou v = KKp ---> sn, cn e dn valores tabelados;
    ! 6) Obten��o de x e y com a f�rmula da pg. 39 de TU. Esses contornos s�o salvos em 'contornos.txt'. Por�m:
    !    Quando u = KK E v = KKp essa equa��o fornece 0/0 resultando em x e y = NaN (pois ambos os valores s�o tabelados), por isso existe um IF qua transforma esse ponto NaN em 0 (que deveria ser o valor resultante)
    ! 7) Escolhe um contorno fixando o v (colunas de x e y) e defini como a curva interna da pastilha, ap�s o que:
    !    Define o contorno externo como um m�ltiplo do interno (por enquanto);
    !    Loop interno: Desenha as retas radiais de ponto a ponto (bic�ncava a bic�ncava na horizontal) com o passo de X e Y calculados como a "espessura/intervalos" 
    !    Loop externo: muda de "altura" na bic�ncava para tra�ar a pr�xima reta radial, o passo � recalculado para X e Y pois a espessura varia ao longo da geometria
    !    Salva esses vetores X e Y de retais radias como linhas das matrizes finais X_malha, Y_malha no arquivo "geometria.txt"
    !********************************************************************************************************************************************************************************************************************
    
    USE dadosGerais
    USE double_Precision

    IMPLICIT NONE 
    REAL(KIND=DP):: m, k, kc, mc, tol, conv1, conv2                     ! Par�metro m, m�dulo, m�dulo e par�metro complementares, toler�ncias
    REAL(KIND=DP):: KK, KKp                                             ! Integrais el�pticas KK (completa de 1o tipo) e KKp (completa complementar de 1o tipo)
    REAL(KIND=DP), PARAMETER:: a_proporcionalidade = 0.8810933       ! Par�metra a que multiplica a geometria inteira de forma a ter o mesmo di�mtero interno de Cassini e do anular do Yuan (em cm) para v=0.2 e m=0.8
    REAL(KIND=DP), DIMENSION(0:100):: a, b, c, ap, bp, cp               ! Par�metros da m�dia GA (vetores super dimensionados)
    REAL(KIND=DP), DIMENSION(0:SIZE (v)-1, 0:SIZE(u)-1):: x_contorno    ! x(v,u) --> obtido ap�s transforma��o de (x - iy) = const. * cn(u + iv), pg. 39 do Tu
    REAL(KIND=DP), DIMENSION(0:SIZE (v)-1, 0:SIZE(u)-1):: y_contorno    ! y(v,u) --> obtido ap�s transforma��o de (x - iy) = const. * cn(u + iv), pg. 39 do Tu
    
    REAL(KIND=DP), DIMENSION(0:SIZE(u)-1)        :: snu, cnu, dnu       ! Fun��es el�pticas jacobianas
    REAL(KIND=DP), DIMENSION(0:SIZE(v)-1)        :: snv, cnv, dnv       ! Fun��es el�pticas jacobianas 
    REAL(KIND=DP), DIMENSION(0:100, 0:SIZE(u)-1) :: phi_n               ! Matriz PHI: linhas: n+1 valores de PHI (100 - n+1 valores em branco), colunas: valor de PHI para cada u, onde n+1 = numero de itera��es da GA
    REAL(KIND=DP), DIMENSION(0:100, 0:SIZE(v)-1) :: phi_p               ! Matriz PHI: linhas: n+1 valores de PHI (100 - n+1 valores em branco), colunas: valor de PHI para cada v, onde n+1 = numero de itera��es da GA
    INTEGER:: n, p, i, j, w, l, g, t,q,col,indice_bic                   ! Contadores
    INTEGER       :: coluna_escolhida                                   ! coluna de x e y (provenientes do c�lculo com as fun��es el�pticas) escolhida como contorno da malha
     
    
    OPEN(UNIT = 1,   FILE = 'geometria.txt', STATUS = 'unknown')! Esse arquivo � usado tanto pela SUBROUTINE Cassini_geo quanto pela SUBROUTNE Eliptico_geo, dependendo de qual seja chamada pelo principal
    OPEN(UNIT = 100, FILE = 'contornos_eliptico.txt', STATUS = 'unknown')
   
    !WRITE(*,*) "Entre com o parametro m"
    !READ(*,*) m
    m=0.800000000000000000000
    
    ! Defini��o do m�dulo e dos complementares
    mc = 1 - m
    k  = SQRT(m)
    kc = SQRT(1 - k**2)
    

    !************************************************************************************************************************************   
    ! --------------------------- 1) Calcula os valores de a, b e c (Para a integral K (ou KK)) ---------------------------------------
    !************************************************************************************************************************************ 
    
    ! Valores iniciais para a m�dia geom�trica-aritm�tica (GA) para o c�lculo de snu, cnu e dnu e KK
    a(0)  = 1.
    b(0)  = SQRT(mc)
    c(0)  = SQRT(m)
    
    tol = 1E-10
    conv1 = 5*tol  
    n = 0
    DO WHILE (ABS(conv1) .GE. tol)
        n = n + 1
         
        a(n) = 0.5*(a(n-1) + b(n-1))
        b(n) = SQRT(a(n-1)*b(n-1))
        c(n) = 0.5*(a(n-1) - b(n-1))
        
        conv1 = c(n)
        
        a(n-1)  = a(n)
        b(n-1)  = b(n)
        c(n-1)  = c(n)
                                                 
    END DO 
  
    
    
    !************************************************************************************************************************************   
    ! ------------------------ 1.1) Calcula os valores de ap, bp e cp (para a integral complementar K'(KKp)) ----------------------------
    !************************************************************************************************************************************  
        
    ! Valores iniciais para a m�dia geom�trica-aritm�tica (GA) para o c�lculo de snv, cnv e dnv e KKp
    ap(0)  = 1.
    bp(0)  = SQRT(m)
    cp(0)  = SQRT(mc)
        
    
    conv2 = 5*tol
    p = 0
     DO WHILE (ABS(conv2) .GE. tol)        
        p = p + 1 
        
       
        ap(p) = 0.5*(ap(p-1) + bp(p-1))
        bp(p) = SQRT(ap(p-1)*bp(p-1))
        cp(p) = 0.5*(ap(p-1) - bp(p-1)) 
          
        conv2 = cp(p)
        
        ap(p-1)  = ap(p)
        bp(p-1)  = bp(p)
        cp(p-1)  = cp(p)
        
     END DO 
     
    !************************************************************************************************************************************   
    ! ---------------------------------- 2) Integrais el�pticas completas de primeiro tipo ----------------------------------------------
    !************************************************************************************************************************************  
    
    KK = PI/(2.*ap(p))
    KKp = PI/(2.*a(n))   
    
    !************************************************************************************************************************************   
    ! ----------------------------------------- 3) Defini��o dos vetores u e v ----------------------------------------------------------
    !************************************************************************************************************************************ 
    
    DO i = 0, SIZE(u)-1       ! Cria o vetor u com intervalos iguais de 0 at� a dimens�o declarada
        u(i) = i*(1./(SIZE(u)-1))*KK 
        !u(i) = (KK - u(i-1))/3. + u(i-1)
        !u(SIZE(u)-1) = KK
    END DO
    
    
    
    DO j = 0, SIZE(v)-1       ! Cria o vetor v com intervalos iguais de 0 at� a dimens�o declarada
        !v(j) = j*(1./(SIZE(v)-1))*KKp 
        v(j)=0.2*KKp 
        !v(j) = j/10*KKp
        !v(SIZE(v)-1) = KKp
    END DO 
    
    
    !************************************************************************************************************************************   
    ! ---------------------------------------------- 4) C�culo da amplitude PHI ---------------------------------------------------------
    !************************************************************************************************************************************ 
    !Obs: o c�lculo de phi_n-1 est� sempre defasado em rela��o aos vetores a(n) e b(n), pois quando estou calculando phi_n-1 uso o a(n) e 
    !c(n). 
    !Exemplo: se n=4, numero total de elementos em a(n) = n+1=5, phi_n(5,i) = 2**(n+1)*a(n)*u(i), e phi_n(4,1)
    !� calculando no 'ELSE'do loop abaixo utilizando a(4) e c(4) (valores da converg�ncia da m�dia GA) e o phi_n(5,i).
    !O valor de interesse phi_n(0,i) utiliza a(1) e c(1), primeiros dois valores calculados na GA (e n�o a(0) e c(0))
    !************************************************************************************************************************************ 
    
    !----------------------------------------------------------------------------------------------------------------------------------!
    !            u(1)          u(2)          u(3)            u(4)       ...      u(i)                                                  !
    !phi_n =   phi_n(0,1)   phi_n(0,2)     phi_n(0,3)     phi_n(0,4)           phi_n(0,i) ----> pen�ltimo valor de phi (phi_0), quero  !
    !            ...         ...            ...             ...                  ...                                                   !
    !       phi_n(n+1,1)   phi_n(n+1,2)   phi_n(n+1,3)   phi_n(n+1,4)         phi_n(n+1,i) ----> �ltimo valor de phi (phi_N), conhe�o  !
    !----------------------------------------------------------------------------------------------------------------------------------!

    DO i= 0, SIZE(u)-1                      ! Como o �ndice de u parte de 0, e.x. se n = 10, num. elemen. = size (u)= n+1

         DO w = n+1, 0, -1                  ! Soma 1 ao n pois n � o n�mero total de elementos no vetor a(n), que tamb�m parte de 0
             
            IF (w .EQ. n+1) THEN
                 phi_n(n+1,i) = 2**(n+1)*a(n)*u(i)      ! Obtendo o �ltimo phi para cada u (n+1 = n�mero de itera��es para encontrar a GA)
            ELSE   
                 phi_n(w,i) = 0.5*(ASIN(c(w)/a(w)*SIN(phi_n(w+1,i))) + phi_n(w+1,i)) ! Para evitar um ALLOCATE, coloquei phi com apenas dois componentes, j� que s� preciso armazenar phi(0) e phi(1)
            END IF 
            
        END DO
    END DO
    
    !----------------------------------------------------------------------------------------------------------------------------------!
    !            v(1)          v(2)          v(3)            v(4)       ...      v(j)                                                  !
    !phi_p =   phi_p(0,1)    phi_p(0,2)    phi_p(0,3)     phi_p(0,4)          phi_p(0,j) ----> pen�ltimo valor de phi                  !
    !             ...         ...            ...             ...                  ...                                                  !
    !         phi_p(p+1,1)  phi_p(p+1,2)   phi_p(p+1,3)   phi_p(p+1,4)        phi_p(p+1,i) ----> �ltimo valor de phi (phi_P)           !
    !----------------------------------------------------------------------------------------------------------------------------------!
   
    
     DO j = 0, SIZE(v)-1    
    
        DO l = p+1, 0, -1          
        
            IF (l .EQ. p+1) THEN 
                phi_p(p+1,j) = 2**(p+1)*ap(p)*v(j)
            ELSE 
                  phi_p(l,j) = 0.5*(ASIN( cp(l)/ap(l)*SIN(phi_p(l+1,j)) ) + phi_p(l+1,j))
            END IF
            
        
        END DO
     END DO
     !WRITE(*,*)'u(5) = ',u(5),'Phi_0 e Phi_1 s�o: ', phi_n(0,5),'e',phi_n(0,5),'n+1= ',n+1

    !************************************************************************************************************************************   
    ! ----------------------------------------------- 5) C�culo de sn, cn e dn -----------------------------------------------------------
    !************************************************************************************************************************************ 
    DO i = 0, SIZE(u)-1                   ! Loop que corre em todo vetor u

        IF (m == 0.) THEN                 ! Verifica��o de m == 0 ou 1
            
            snu(i) = SIN(u(i))
            cnu(i) = COS(u(i))
            dnu(i) = 1.
            
        ELSE IF (m == 1) THEN       
            
            snu(i) = TANH(u(i))
            cnu(i) = 1./COSH(u(i))
            dnu(i) = cnu(i)
            
        ELSE                             ! Verificou se que m n�o � 0 nem 1
     
            IF (u(i) == 0.) THEN         ! Verifica se u(i) == 0 ou u(i) ==  KK
            
                snu(i) = 0.
                cnu(i) = 1.
                dnu(i) = 1.
            
             ELSE IF (u(i) == KK) THEN
            
                snu(i) = 1.
                cnu(i) = 0.
                dnu(i) = kc
            
             ELSE

                snu(i) = SIN(phi_n(0,i))
                cnu(i) = COS(phi_n(0,i))
                dnu(i) = SQRT(1.-k**2*snu(i)**2)   !------> V�lido somente quando m e u s�o estritamente positivos
                !dnu(i) = COS(phi_n(0,i))/( COS(phi_n(1,i) - phi_n(0,i)) )
        
         END IF
        END IF
    END DO
 
    
    
    DO j = 0, SIZE(v)-1                 ! Loop que corre em todo vetor v

        IF (m == 0) THEN                ! Verifica��o de m == 0 ou 1
            
            snv(j) = SIN(v(j))
            cnv(j) = COS(v(j))
            dnv(j) = 1.
            
        ELSE IF (m == 1) THEN       
            
            snv(j) = TANH(v(j))
            cnv(j) = 1./COSH(v(j))
            dnv(j) = cnv(j)
            
        ELSE                            ! Verificou que m n�o � 0 nem 1
     
            IF (v(j) == 0) THEN         ! Verifica se v(j) == 0 ou v(j) ==  KKp
            
                snv(j) = 0.
                cnv(j) = 1.
                dnv(j) = 1.
            
             ELSE IF (v(j) == KKp) THEN
            
                snv(j) = 1            
                cnv(j) = 0
                dnv(j) = kc
            
             ELSE


            snv(j) = SIN(phi_p(0,j))
            cnv(j) = COS(phi_p(0,j))
            dnv(j) = SQRT(1.-kc**2*snv(j)**2)  !-------> V�lido somente quando m e v s�o estritamente positivos
            !dnv(j) = COS(phi_p(0,j))/( COS(phi_p(1,j) - phi_p(0,j)) )
             
             END IF
        END IF
    END DO
    
    
    !************************************************************************************************************************************   
    ! ------------------------------------------------ 6) C�culo de x e y ---------------------------------------------------------------
    !************************************************************************************************************************************ 
    
    ! j -->> ser� o n�mero de colunas, corresponde aos n�meros de v
    ! i -->> ser� o n�mero de linhas, corresponde aos n�meros de u
    
    DO i = 0, SIZE(u)-1     !linhas
        DO j = 0, SIZE(v)-1 !colunas
        
        x_contorno(j,i) = a_proporcionalidade*(cnu(i)*cnv(j) ) / ( cnv(j)**2 + m*snu(i)**2*snv(j)**2 )
        y_contorno(j,i) = a_proporcionalidade*(snu(i)*dnu(i)*snv(j)*dnv(j) )  / ( cnv(j)**2 + m*snu(i)**2*snv(j)**2 )
        
        IF ( IsNaN(x_contorno(j,i)) ) THEN  ! O algoritmo fonece um 0/0 quando u = KK e v = KKp, substitui-se por 0 para resolver esse problema
        x_contorno(j,i) = 0.
        ELSE
       END IF
        
        IF ( IsNaN(y_contorno(j,i)) ) THEN ! O algoritmo fonece um 0/0 quando u = KK e v = KKp, substitui-se por 0 para resolver esse problema
        y_contorno(j,i) = 0.
         ELSE
        END IF
        
        END DO
    END DO
    
 ! Salvando todos os contornos gerados na UNIT=100   

100 FORMAT(32(f14.8),/)
    WRITE(100,*) 'm = ',m,'# colunas: ',j, '# colunas: ',i, 'K = ', KK, 'KKp = ', KKp
    WRITE(100,*) 'Matriz x, linhas--> u, colunas--> v'
    WRITE(100,200) x_contorno
    WRITE(100,*) 'Matriz y, linhas--> u, colunas--> v' 
    WRITE(100,200) y_contorno


!************************************************************************************************************************************   
! ---------------------------------- 7) Escolhe um contorno �nico (fixa um v), com o m dado -----------------------------------------
!************************************************************************************************************************************     
 v=u ! J� salvei os contornos de interesse e posso fazer v = u. Por qu�? As matrizes de malha seguintes s�o compartilhadas ...
     ! ...pelo m�dulo dadosGerais entre Ellipk_geo e Cassini_geo. Em cassini trabalho com 'v', ent�o decidi trabalhar 
     ! ...com essa mesma letra aqui tamb�m para n�o dar confus�o no m�dulo, que est� todo em fun��o de v.
    coluna_escolhida = 16    ! Fixa o v, ou seja, o contorno com o m fixo para ser o molde da geometria
    x_curva_interna=x_contorno(coluna_escolhida,:)
    y_curva_interna=y_contorno(coluna_escolhida,:)
    
    !**ATEN��O: ESCREVENDO A CURVA INTERNA AGORA EM SENTIDO HOR�RIO PARA N�O TER DISCONTINUIDADE QUANDO UTILIZAR A SIMETRIA DE 1/2 DO LADO DIREITO
    x_curva_interna_auxiliar = x_curva_interna
    y_curva_interna_auxiliar = y_curva_interna
    
    k=SIZE(v)
        DO j = 0, SIZE(v) - 1
    k=k-1
    x_curva_interna(j) =  x_curva_interna_auxiliar(k)
    y_curva_interna(j) =  y_curva_interna_auxiliar(k)
    ENDDO

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
    ALLOCATE(X_Malha_FUEL(0:num_bic-1, 0:2*(SIZE(v) -1)))
    ALLOCATE(Y_Malha_FUEL(0:num_bic-1, 0:2*(SIZE(v) -1)))
    
    ! Coordenadas iniciais de X_MALHA_FUEL e Y_MALHA_FUEL    
    X_Malha_FUEL(0,:) = x_curva_interna(:) ! Toda a coluna 0 � prenchida pela bic�ncava definida pela curva interna (primeira bic�ncava conc�ntrica)
    Y_Malha_FUEL(0,:) = y_curva_interna(:) ! Toda a coluna 0 � prenchida pela bic�ncava definida pela curva interna (primeira bic�ncava conc�ntrica)
    
    ! N�o tenho a curva externa ainda, mas tenho as duas "quinas" externas, preciso deles porque uso o ponto da frente e o de tr�s das bic�ncavas interna e externa para o c�culo do coeficiente angular da reta secante.
    X_Malha_FUEL(num_bic-1, 0) = x_curva_interna(0) + df                    ! O primeiro x da bic. externa � apenas a soma do primeiro pt. da bic. interna + a espessura.
    Y_Malha_FUEL(num_bic-1, 0) = 0                                          ! O primeiro y da bic. externa � zero.
            
    X_Malha_FUEL(num_bic-1, SIZE(v) - 1) = 0                                ! O �ltimo x da bic. externa � zero.
    Y_Malha_FUEL(num_bic-1, SIZE(v) - 1) = y_curva_interna(SIZE(v) - 1) + df ! O �ltimo y da bic. externa � apenas a soma do �ltimo y da bic. interna + a espessura.
   
    intervalos = num_bic - 1                ! N�mero total de intervalos entre as bic�ncavas conc�ntricas.
    step    = df/(intervalos)               ! Passo entre os pontos das retas radiais (i. e. bic�ncavas), onde d � a dist�ncia entre a primeira e a �ltima bic�ncava (curva interna e externa).
    ! Loop para encontrar os pontos internos de x_ort e, ap�s isso, y_ort
    DO j = 1, SIZE(v) - 2
        DO g = 1, num_bic - 1                                        ! O 'g' � o �ndice da bic�ncava. A primeira � a curva_interna. A �ltima ser� calculada quando g = intervalos = num_bic-1 = �ltima bic�ncava.
  
    di = step*g                                                      ! Cada ponto da reta radial pertence a uma bic�ncava conc. Para mudar de bic. basta multiplicar g pelo passo 'step'.
              
    X_Malha_FUEL(g,0) = X_Malha_FUEL(0,0)+di
    Y_Malha_FUEL(g,SIZE(v)-1)=Y_Malha_FUEL(0,SIZE(v)-1)+di
              
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
    !2) Para continuar escrevendo a curva no sentido anti-hor�rio, vou espelhar tanto x quanto y;
    !3) O 1o QUADRANTE tem um pt a mais que � o x=0, ou seja, se o primeiro quadrante � SIZE(v) = 51, o total de pontos � 101 e o segundo tem 50 pts.        
   
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
    ALLOCATE(X_Malha_CLAD1(0:C1-1, 0:2*(SIZE(v) -1)))
    ALLOCATE(Y_Malha_CLAD1(0:C1-1, 0:2*(SIZE(v) -1)))

    ! N�o tenho a curva interna do CLAD ainda, mas tenho as duas "quinas" internas, preciso delas porque uso o ponto da frente e o de tr�s das bic�ncavas interna e externa para o c�culo do coeficiente angular da reta secante.
    X_Malha_CLAD1(0, 0) = X_malha_FUEL(0,0) - dc1                ! O primeiro x da bic. interna do CLAD � apenas a subtra��o do primeiro pt. da bic. interna do FUEL - a espessura do CLAD.
    Y_Malha_CLAD1(0, 0) = 0                                          ! O primeiro y da bic. interna do CLAD � zero.
            
    X_Malha_CLAD1(0, SIZE(v)-1) = 0                                  ! O �ltimo x da bic. interna do CLAD � zero.
    Y_Malha_CLAD1(0, SIZE(v)-1) = Y_malha_FUEL(0,SIZE(v) - 1) - dc1 ! O �ltimo y da bic. interna do CLAD � apenas a subtra��o do �ltimo y da bic. interna do FUEL - a espessura do CLAD.
    
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
    
    X_Malha_CLAD1(g,0)                = X_Malha_FUEL(0,0) - di       !calculando a reta radial que fica sobreposta ao eixo x
    Y_Malha_CLAD1(g,0)                 = 0                           !calculando a reta radial que fica sobreposta ao eixo x

    X_Malha_CLAD1(g,SIZE(v)-1)         = 0                           !calculando a reta radial que fica sobreposta ao eixo y
    Y_Malha_CLAD1(g,SIZE(v)-1) = Y_Malha_FUEL(0,SIZE(v)-1)-di        !calculando a reta radial que fica sobreposta ao eixo y
                            
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
        linha_radial = SIZE(v)-1
        
        DO j = SIZE(v), 2*(SIZE(v)-1)
            
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
    ALLOCATE(X_Malha_CLAD2(0:C2-1, 0:2*(SIZE(v) -1)))
    ALLOCATE(Y_Malha_CLAD2(0:C2-1, 0:2*(SIZE(v) -1)))
    
   ! Preenchendo a primeira bic�ncava do CLAD2 com o �ltima curva do FUEL
   !! X_Malha_CLAD2(0,:) = X_Malha_FUEL(num_bic-1,:) ! Toda a coluna 0 � prenchida pela bic�ncava definida pela curva externa do FUEL (primeira bic�ncava conc�ntrica do CLAD2)
    !!Y_Malha_CLAD2(0,:) = Y_Malha_FUEL(num_bic-1,:) ! Toda a coluna 0 � prenchida pela bic�ncava definida pela curva externa do FUEL (primeira bic�ncava conc�ntrica do CLAD2)

    ! N�o tenho a curva externa ainda, mas tenho as duas "quinas" externas, preciso deles porque uso o ponto da frente e o de tr�s das bic�ncavas interna e externa para o c�culo do coeficiente angular da reta secante.
    X_Malha_CLAD2(C2-1, 0) = X_Malha_FUEL(num_bic-1,0) + dc2               ! O primeiro x da bic. externa � apenas a soma do primeiro pt. da bic. externa do FUEL + a espessura do CLAD.
    Y_Malha_CLAD2(C2-1, 0) = 0                                      ! O primeiro y da bic. externa � zero.
            
    X_Malha_CLAD2(C2-1, SIZE(v)-1) = 0                                ! O �ltimo x da bic. externa � zero.
    Y_Malha_CLAD2(C2-1, SIZE(v)-1) = Y_Malha_FUEL(num_bic-1,SIZE(v) - 1) + dc2 ! O �ltimo y da bic. externa � apenas a soma do �ltimo y da bic. externa do FUEL + a espessura.
   
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
       
    X_Malha_CLAD2(g,0)         = X_Malha_FUEL(num_bic-1,0) + di          !calculando a reta radial que fica sobreposta ao eixo x
    Y_Malha_CLAD2(g,0)         = 0                               !calculando a reta radial que fica sobreposta ao eixo x

    X_Malha_CLAD2(g,SIZE(v)-1) = 0                               !calculando a reta radial que fica sobreposta ao eixo y
    Y_Malha_CLAD2(g,SIZE(v)-1) = Y_Malha_FUEL(num_bic-1,SIZE(v)-1)+di           !calculando a reta radial que fica sobreposta ao eixo y
                   
    m_incli  = (Y_Malha_FUEL(num_bic-1,j+1) - Y_Malha_FUEL(num_bic-1,j-1)) / ( X_Malha_FUEL(num_bic-1,j+1) - X_Malha_FUEL(num_bic-1,j-1))   ! Coeficiente angular da reta secante (expans�o em S�rie de Taylor).
    
    Bpp  = Y_Malha_FUEL(num_bic-1,j) + X_Malha_FUEL(num_bic-1,j)/m_incli                                                     ! C�culo do coeficiente linear da equa��o da curva y_ortogonal(j) = (-1/m)*x_ortogonal(j) + Bpp.
    
    A_bask = 1.0 + 1.0/m_incli**2                                                             ! � o termo "a" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    B_bask = -2.0*X_Malha_FUEL(num_bic-1,j) + 2*Y_Malha_FUEL(num_bic-1,j)/m_incli - 2*Bpp/m_incli                            ! � o termo "b" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
    C_bask = X_Malha_FUEL(num_bic-1,j)**2 + Y_Malha_FUEL(num_bic-1,j)**2 - 2*Y_Malha_FUEL(num_bic-1,j)*Bpp + Bpp**2 - di**2    ! � o termo "c" da equa��o de Bh�skara para achar o x ortogonal (na mesma reta radial).
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
    
    ! Definindo o n�mero de bic�ncavas e renomeando para melhor compreens�o nas subrotinas seguintes
    biconcavas = C1 + num_bic + C2 ! Total de bic�ncavas (em n�mero de pontos da parte interna para a externa)
    radiais    = 2*SIZE(v)-1       ! Total de retas radiais (em n�mero pontos longitudinais)
    
    F1      = C1+1                 ! Se C1=15, CLAD1 vai de 0:14 (=0:C1-1) e F1 = C1+1 pts. F1 � LIMITE, C1 � N�MERO DE PONTOS!
    F2      = F1 + num_bic-1       ! F2 � o LIMITE EXTERNO do FUEL, logo FUEL come�a em F1 (=C1+1=ao ponto 16) e anda num_bic pts (F1:num_bic-1), num_bic-1 porque o F1 EST� INCLUSO!
    Ifinal  = F2+1 + C2-1          ! J� o CLAD2 come�a no ponto a frente do final do FUEL e termina C2-1 pts depois, J� QUE O F2+1 EST� INLCUSO tenho que tirar um ponto (como no caso do ZERO...)
 
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

    !Escreve a geometria malhada em geometria.txt, j� formatado para plot em Scilab
200 FORMAT(32(f14.8),/)
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
    CLOSE(UNIT=100)
    
        END SUBROUTINE
