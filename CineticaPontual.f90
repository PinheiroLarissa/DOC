SUBROUTINE Cinetica_Pontual
    USE dadosGerais
    USE DadosTermicos
    USE double_Precision
    
   
    
    
    DO n = 0, tfinal, dt
        
        dt = tfinal
        DO j=1,6
        sum_C(n) = frac(j) * C(j,n)
        END DO
        P(n+1)= ((react(n) - beta_delayed)/Lambda)*P(n) + sum_C))*dt + P(n)
        
        C(j,n+1) = 
        
    
    
    
    
    
    
    
    
    
    
    END SUBROUTINE Cinetica_Pontual