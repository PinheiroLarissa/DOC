   Subroutine LedadosTermicos
    
    use DadosGerais
    
    use DadosTermicos
    
    implicit none
    
    integer::i,j
    
    open(unit=21,file='dadoscalor.txt',status='old')
    
    
    read(21,*)Kf, Kc, densidade_de_potencia, h_int, h_out, hg, T_int, T_out
    
    read(21,*)((Temperatura(i,j),i=1,biconcavas),j=1,radiais)
    
    close(unit=21)

	
    
  
    
    end subroutine LedadosTermicos