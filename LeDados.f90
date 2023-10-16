    Subroutine Ledados
    !L� os dados do problema: condi��es de contorno e pontos estimados de malha
    
	Use DadosGerais

	implicit none

	integer::i !n�mero de colunas, c�rculos conc�ntricos no dom�nio f�sico
    integer::j !n�mero de linhas, retas radiais que cortam os c�rculos no dom�nio f�sico

	open(unit=1,file='geometria.txt',status='old')
    
    !ALTERADO: Lembrar que no geometria.txt declaro primeiro o biconcavas (bic�ncavas conc�ntricas) e depois radiais (retas radiais) 
    read(1,*)  biconcavas, radiais
    WRITE(*,*) biconcavas, radiais

    read(1,*)   ((x(i,j),i=1,biconcavas),j=1,radiais)


	  read(1,*) ((y(i,j),i=1,biconcavas),j=1,radiais)


    close(unit=1)

	
	End subroutine Ledados

















