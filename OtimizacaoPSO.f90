SUBROUTINE PSO
    
!####################################################################################!
! Arquivo criado em:  19/11/2018                                                     !
! �ltima atualiza��o: 30/11/2018                                                     !
! Programador: Larissa                                                               !
! Descri��o:                                                                         !
! Subrotina que implementa o PSO cl�ssico para miniza��o da fun��o esfera            !
!####################################################################################!
    
IMPLICIT NONE
INTEGER, PARAMETER                :: DP = SELECTED_REAL_KIND(14)    ! Par�metro selecionado pra dupla precis�o do compilador
INTEGER                           :: sizeSeed                       ! Tamanho do vetor que armazena a seed do randomizador
INTEGER, ALLOCATABLE              :: new(:)                         ! Vetor que ir� armazenar a nova seed

INTEGER                           :: N_populacao, dimensoes         ! N�mero de indiv�duos na popula��o e dimens�es do problema (entrada pelo usu�rio)
INTEGER                           :: geracoes                       ! N�mero de gera��es (�pocas) do loop externo
INTEGER                           :: j ,i, k                        ! Contadores: j -> colunas (dimens�es), i -> linhas (indiv�duo), k -> gera��es
REAL(KIND=DP)                     :: r1, r2                         ! Par�metros rand�micos de 0 a 1 do c�lculo da velocidade do PSO
REAL(KIND=DP)                     :: max_velocidade, min_velocidade ! Limita��o da velocidade
REAL(KIND=DP)                     :: max_posicao, min_posicao       ! Limita��o do espa�o de busca
REAL(KIND=DP), ALLOCATABLE        :: fitness(:)                     ! Vetor que guarda a avalia��o da fitness do indiv�duo i, tamanho: N_populacao
REAL(KIND=DP), ALLOCATABLE        :: fitness_pbest(:)               ! Vetor que guarda a avalia��o da fitness da melhor posi��o do indiv�duo i, tamanho: N_populacao
REAL(KIND=DP)                     :: fitness_gbest          ! Escalares que guardam a avalia��o da fitness da melhor posi��o do enxame, sendo gbest a melhor posi��o
REAL(KIND=DP), ALLOCATABLE        :: x(:,:), v(:,:)                 ! Matrizes que guardam as posi��es e velocidades dos indiv�duos, tamanho: dimensoes x N_Popula��o 
REAL(KIND=DP), ALLOCATABLE        :: pbest(:,:), gbest(:)                     ! pbest: matriz que guarda o melhor individual em cada dimens�o, vetor que guarda os melhores globais 
REAL(KIND=DP), PARAMETER          :: c1 = 2., c2 = 2., w0 =0.9      ! Par�metros da velocidade, c1: peso do termo cognitivo, c2: peso do aprendizado social, w: in�rcia
REAL(KIND=DP), ALLOCATABLE        :: w(:)                           ! In�rcia

REAL, DIMENSION(10,8)             :: F
INTEGER                          :: L, M, CONTADOR
WRITE(*,*)'!#############################################################################!'
WRITE(*,*)'!Arquivo criado em:  19/11/2018                                               !'
WRITE(*,*)'!Ultima atualizacao: 30/11/2018                                               !'
WRITE(*,*)'!Programador: Larissa                                                         !'
WRITE(*,*)'!Descricao:                                                                   !'
WRITE(*,*)'!Subrotina que implementa o PSO classico para minizacao da funcao esfera      !'
WRITE(*,*)'!#############################################################################!'


!***************************************************** RANDOMIZADOR *************************************************************
CALL RANDOM_SEED ( )                           ! O processador inicializa a subrotina RANDOM_SEED de acordo com a hora e data

CALL RANDOM_SEED (SIZE = sizeSeed)             ! Recupera-se o tamanho da seed, que � ent�o armazenado em I

ALLOCATE (new(sizeSeed))                       ! Aloca��o do vetor 'new' que definir� a nova seed de acordo com o usu�rio

new = 7850                                     ! Definindo um valor novo de seed, para sempre partirmos do mesmo valor
CALL RANDOM_SEED (PUT=new(1:sizeSeed))         ! Calculando a nova seed de acordo com o vetor 'new', cujo valor � passado acima
!********************************************************************************************************************************


WRITE(*,*) 'Entre com o numero de individuos'
READ(*,*) N_populacao
WRITE(*,*) 'Entre com o numero de dimensoes'
READ(*,*) dimensoes

!********************************************* INICIALIZA��O DOS VETORES *******************************************************
ALLOCATE(x(dimensoes,N_populacao), v(dimensoes,N_populacao), pbest(dimensoes,N_populacao))
ALLOCATE(fitness(N_populacao), fitness_pbest(N_populacao), gbest(dimensoes))

CALL RANDOM_NUMBER(x)
CALL RANDOM_NUMBER(v)
CALL RANDOM_NUMBER(r1)
CALL RANDOM_NUMBER(r2)
x= x*100
v = 0.
!********************************************** DEFINI��O DOS PAR�METROS ******************************************************
geracoes      = 5
ALLOCATE(w(geracoes))

pbest         = 30.**50
gbest         = 1000.
fitness_gbest = 30.**50
fitness_pbest = 30.**50

max_posicao    = 10.
min_posicao    = -10.
max_velocidade = 10.
min_velocidade = -10.
k              = 0
!****************************************** C�LCULO DAS VELOCIDADES E POSI��ES ************************************************

DO k = 1, geracoes                  ! Comente se for usar o loop DO WHILE

!DO WHILE(abs(gbest) >  10.**-10)     ! Comente se for usar o par�metro: geracoes
!k = k+1                              ! Comente se for usar o par�metros: geracoes

w(k) = ((w0-0.4)*(geracoes - k))/(geracoes+0.4)   ! Recalculando a in�rcia ap�s um itera��o completa

    DO i = 1, N_populacao 
    
            fitness(i) = 0. 
            fitness(i) = sum(x(:,i)**2)   ! Estou achando a fitness multiplicando toda as colunas (dimens�es) da linha i (indiv�duo)

                ! Encontrando e armazenando a melhor fitness do indiv�duo i        
                IF (fitness(i) < fitness_pbest(i)) THEN
                    fitness_pbest(i) = fitness(i)
                     pbest(:,i) = x(:,i)
                ELSE
                END IF
                
                 
 
                ! Encontrando e armazenando a melhor fitness do grupo
                IF (fitness(i) < fitness_gbest) THEN
                    fitness_gbest  = fitness(i)
                    gbest(:) = x(:,i)

                ELSE 
                END IF
                
       
            DO j = 1, dimensoes
               
                ! Atualizando a velocidade
                v(j,i) = w(k)*v(j,i) + r1*c1*(pbest(j,i) - x(j,i)) + r2*c2*(gbest(j) - x(j,i))
                
                ! Limitando a velocidade
                IF (v(j,i) > max_velocidade) THEN 
                    v(j,i) = v(j,i) - max_velocidade 
                ELSEIF (v(j,i) < min_velocidade) THEN
                    v(j,i) = v(j,i) + min_velocidade  
                END IF
                
                x(j,i) = v(j,i) + x(j,i)
                
                ! Limitando o espa�o de busca
                IF (x(j,i) > max_posicao) THEN 
                    x(j,i) = x(j,i) - max_posicao
                ELSEIF (x(j,i) < min_posicao) THEN
                    x(j,i) = x(j,i) + min_posicao
                END IF
                
             ! Imprime na tela os resultados atuais   
            !WRITE(*,*) 'O gbest e: ', gbest 
            ! WRITE(*,*) 'Geracao: ', k
             
             
        
         END DO !(j)
         
    END DO !(i)

    
END DO !(k)

!**************************************************************************************************************************

   CONTADOR = 0
             DO L = 1, 10 ! COLUNA
                 DO M = 1, 8 ! LINHA
                     CONTADOR = CONTADOR + 1
                     F(L,M)  = CONTADOR 
                     
                 END DO
             END DO
                  WRITE(*,*) F(:,1) 


!WRITE(*,*) 'Convergiu em:', k, 'iteracoes'
!WRITE(*,*) 'A melhor solucao e:', fitness_gbest

     END SUBROUTINE