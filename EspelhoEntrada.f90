    subroutine EspelhoEntrada

	USE DadosGerais

	implicit none

	integer::i,j

	!imprime dados de entrada para simples conferência

	write(2,*) 'Pontos x(i,j), as colunas são o i-biconcavas concêntricas do domínio físico e as linhas são as j-retas radiais'

	write(2,1) ((x(i,j),i=1,biconcavas),j=1,radiais)

	write(2,*) 'Pontos y(i,j), as colunas são o i-biconcavas concêntricas do domínio físico e as linhas são as j-retas radiais '

	write(2,1) ((y(i,j),i=1,biconcavas),j=1,radiais)

1   format(32(E15.3),/)

	end subroutine EspelhoEntrada