    subroutine EscreveCoeficientes

	use DadosGerais

	implicit none 

	integer::i,j

	write(4,*) 'coeficientes alfa'

	write(4,1) ((i,j,alfa(i,j),i=1,biconcavas),j=1,radiais)

	write(4,*) 'coeficientes beta'

	write(4,1) ((i,j,beta(i,j),i=1,biconcavas),j=1,radiais)

	write(4,*) 'coeficientes gama'

	write(4,1) ((i,j,gama(i,j),i=1,biconcavas),j=1,radiais)

	write(4,*) 'Jacobianos'

	write(4,1) ((i,j,Jacobiano(i,j),i=1,biconcavas),j=1,radiais)
WRITE(*,*) 'esta calculando a rotina escreve coeficientes'

1	format(i10,i10,e19.8,1x,/)

	end subroutine EscreveCoeficientes