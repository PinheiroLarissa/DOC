    subroutine Escrevetempk(k)

	use dadosgerais

	use dadostermicos

	implicit none

	integer::k,i,j

	write(33,*) k


	write(33,200) (((i,j,temperatura(i,j)),i=1,radiais),j=1,biconcavas)



	200 Format(75(1x,i5,1x,i5,E15.8),/)

	end subroutine
