    Subroutine Ledados
    !Lê os dados do problema: condições de contorno e pontos estimados de malha
    
	Use DadosGerais

	implicit none

	integer::i !número de colunas, círculos concêntricos no domínio físico
    integer::j !número de linhas, retas radiais que cortam os círculos no domínio físico

	open(unit=1,file='geometria.txt',status='old')
    
    !ALTERADO: Lembrar que no geometria.txt declaro primeiro o biconcavas (bicôncavas concêntricas) e depois radiais (retas radiais) 
    read(1,*)  biconcavas, radiais
    WRITE(*,*) biconcavas, radiais

    read(1,*)   ((x(i,j),i=1,biconcavas),j=1,radiais)


	  read(1,*) ((y(i,j),i=1,biconcavas),j=1,radiais)


    close(unit=1)

	
	End subroutine Ledados

















