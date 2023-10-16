 subroutine EspelhoDadosTermicos
    
    use DadosGerais

    use DadosTermicos
    
    implicit none
    
    integer::i,j
    ALLOCATE(Temperatura   (1:biconcavas, 1:radiais), STAT = AllocateStatus)

    write(22,*) 'Kf,Kc, densidade_de_potencia,h_int,h_out,T_int,T_out'
    
    write(22,*) '                                      '
    
    write(22,100)  Kf, Kc, densidade_de_potencia, h_int, h_out, hg, T_int, T_out
    
    write(22,*) 'Temperaturas lidas'
    
    write(22,*) '                                      '

    write(22,1) ((temperatura(i,j),i=1,biconcavas),j=1,radiais)  
    
1   format(32(F13.4),1x,/)    
100 format(2F10.2,E10.2,5F10.3)       
    end subroutine EspelhoDadosTermicos 
    