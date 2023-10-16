!"include 'lapack.f90' !Descomentar e adicionar a biblioteca MKL para testar o Lapack
    PROGRAM Principal
    USE dadosGerais

	IMPLICIT NONE

	INTEGER::opcao, geometria, i, j
 

	!Programa que calcula a malha e as distribuicoes de Temperaturas - Alvim - 08_08_206

	WRITE(*,*) 'Entre com a opcao'

	WRITE(*,*) ' 0 - Unidimensional, 1 - Bidimensional, 2 - PSO'

	READ(*,*) opcao

	SELECT CASE (opcao)

	CASE(0)

	CALL calculaTemperaturaRadial

    CASE(1)
    WRITE(*,*) 'Escolha a geometria:'
	WRITE(*,*) ' 0 - Eliptico, 1 - Cassini, 2- Nao gerar a malha (usa o txt ja salvo), 3 - Teste MKL, 4 - Le a malha feita em scilab como chute inicial'
	READ(*,*) geometria
    SELECT CASE (geometria)

    CASE(0)
        CALL Ellipk_geo
    
    CASE(1)
        CALL Cassini_geo
    
    CASE(2)
        
        OPEN(UNIT=777, FILE='ResultadoMalhaMemoria.txt', STATUS='OLD')
        READ(777,*) biconcavas, radiais, F1, F2
        ALLOCATE(X            (1:biconcavas, 1:radiais))
        ALLOCATE(Y            (1:biconcavas, 1:radiais))
        READ(777,*) ((X(i,j),i=1,biconcavas),j=1,radiais)
        READ(777,*) ((Y(i,j),i=1,biconcavas),j=1,radiais)
        CLOSE(UNIT=777)
        
        OPEN(UNIT=888, FILE = 'ResultadoCoeficientesMemoria.txt', status = 'unknown') ! Essa unidade salva o x e o y em um txt para não precisar gerar a malha em todas as rodadas
        ALLOCATE(alfa         (1:biconcavas, 1:radiais))
        ALLOCATE(beta         (1:biconcavas, 1:radiais))
        ALLOCATE(gama         (1:biconcavas, 1:radiais))
        ALLOCATE(Jacobiano    (1:biconcavas, 1:radiais))
        READ(888,*) ((alfa(i,j),i=1,biconcavas),j=1,radiais)
        READ(888,*) ((beta(i,j),i=1,biconcavas),j=1,radiais)
        READ(888,*) ((gama(i,j),i=1,biconcavas),j=1,radiais)
        READ(888,*) ((Jacobiano(i,j),i=1,biconcavas),j=1,radiais)
        CLOSE(UNIT=888)

        GO TO 1
    
    CASE(3)
        !Call MKL_test
        
    CASE(4)
        CALL LeDadosProblema
    END SELECT
     
    ! Chama a subrotina que gera a malha
	CALL Geramalha

	!Chama a subrotina que calcula as temperaturas nessa malha

1    CALL CalculaDistribuicaoTemperaturas

    !CASE(2)
        
        !CALL PSO
 

END SELECT
	END PROGRAM Principal