    subroutine EscreveTemperaturas

	use dadosGerais

	use DadosTermicos

	implicit none

	integer::i,j,k

	open(unit=13,file='Temperaturas.txt',status='unknown')

	write(13,*) 'temperaturas ='

	write(13,*) '      '

    write(13,200) (((temperatura(i,j)),i=1,biconcavas),j=1,radiais)
	!write(13,200) (((i,j,temperatura(i,j)),i=1,biconcavas),j=1,radiais)

    200 Format(32(F10.3),/)
	!200 Format(15(1x,i5,1x,i5,F10.3),/)

	close(unit=13)

	end subroutine EscreveTemperaturas

