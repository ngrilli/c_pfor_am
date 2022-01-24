c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c
c Aug. 12th, 2021 - 1st working version
c
	module sliphardlaws
      implicit none
      contains
c      
c      
c      
      subroutine sliphard(sstate,gsum,gint,gdot,temp,sstate0,sXdist,
     &sstatedot)
      use globalvars, only: modelno,sliphard_param,Rgas,numstvar,
     &numslip, G
      implicit none
      real(8) gsum, gint, gdot, temp, sXdist
      real(8)	sstatedot(numstvar+1), sstate(numstvar), sstate0(numstvar)
c     generic variables
      real(8) tauc, x
c     variables for Voce hardening law      
      real(8) h0_1, ss_1, n_1
      real(8) hb_1
c     variables for isotropic and kinematic hardening with AR thermal softening 
      real(8) tauc0_2, h0_2, n_2, A_2, Q_2, d_2, h_2, hD_2
c     variables for isotropic and kinematic hardening (Code Aster model) 
      real(8) tauc0_3, b_3, Q_3, C_3, D_3    
c     variable used in GND model with slip gradients
      real(8) b_4, yc_4, K_4
      real(8) rhoSSD, rhotot      
c     variables for Voce hardening law with slip distance     
      real(8) h_5, hD_5, alpha_5, b_5, C_5, yc_5, K_5
      real(8) mKM_5
c      
c      
      sstatedot = 0.0d+0
c      
c      
c     Three parameter hardening law      
      if (modelno.eq.1) then
c          
c
c
c
c          
          tauc = sstate(1)
c          
          h0_1 = sliphard_param(2)
          ss_1 = sliphard_param(3)
          n_1 = sliphard_param(4)
c          
          hb_1 = h0_1*(dabs(1.0d+0-(tauc/ss_1)))**n_1
          
c         self-hardening          
		sstatedot(1) = hb_1*dabs(gdot)
c          
c          
c          write(6,*) 'inside sliphar'
c          write(6,*) 'selfhrd'
c         Other slip hardening laws come here!!!
c
c          write(6,*) 'sstatedot(1)', sstatedot(1)
c
c           write(6,*) 'sstate(1)', sstate(1)
c
      elseif (modelno.eq.2) then
c          
c    
c          
c          
          tauc = sstate(1)
          x = sstate(2)
c          
          tauc0_2 = sstate0(1)
          h0_2 = sliphard_param(2)
          n_2 = sliphard_param(3)
          Q_2 = sliphard_param(4)
          d_2 = sliphard_param(5)
          A_2 = sliphard_param(6)
          h_2 = sliphard_param(7)
          hD_2 = sliphard_param(8)
c          
c          
c         self-hardening          
		sstatedot(1) = h0_2*dabs(gdot)*(1.0d+0 + 
     &(h0_2*gsum/tauc0_2/n_2))**(n_2-1.0d+0)
c          
c         Thermal softening
          sstatedot(3) = -A_2 * (tauc**d_2) * dexp(-Q_2/Rgas/temp) 
c 
c         kinematic hardening
          sstatedot(2) = h_2*gdot - hD_2*x*dabs(gdot)
c
c
      elseif (modelno.eq.3) then
c          
c    
c          
c          
          tauc = sstate(1)
          x = sstate(2)
c          
          b_3 = sliphard_param(2)
c          Q_3 = sliphard_param(3)
          C_3 = sliphard_param(4)
          D_3 = sliphard_param(5)
c          
c          
c         self-hardening          
c		sstatedot(1) = Q_3 * (1. - dexp(-b_3*gint))
          sstatedot(1) = (1.0d+0 - dexp(-b_3*gint))
          
 
c         kinematic hardening
          sstatedot(2) = C_3*gdot - D_3*x*dabs(gdot)
          
          
          
          
          
      elseif (modelno.eq.4) then
          
    
          
          
          rhoSSD = sstate(1)
c          write(6,*) 'rhoSSD'
c          write(6,*) rhoSSD
          
          
          
c         No rates are defined for GND and backstress state variables
          
          rhotot = sstate(1) + dabs(sstate(2)) + dabs(sstate(3))
          
          b_4 = sliphard_param(3)
          K_4 = sliphard_param(5)
          yc_4 = sliphard_param(6)



          
c        SSD evolution          
         sstatedot(1)=dabs(gdot)/b_4*(K_4*dsqrt(rhotot) - 
     & 2.0d+0*yc_4*rhoSSD)
 
c          write(6,*) 'sstatedot(1)', sstatedot(1)
          
          
c         GND evolution is calculated at the end of the former time increment  

c         GND rate is set to zero for implicit integration (states = 2 & 3)
c         Backstress rate is set to zero for implicit integration (state=4)

          
                    
          
c     Three parameter hardening law with slip distance effect (size)
      elseif (modelno.eq.5) then
c          
c          
          rhoSSD = sstate(1)
          x = sstate(2)
c          
c
c
c          
c          
          b_5 = sliphard_param(3)
          K_5 = sliphard_param(4)
          yc_5 = sliphard_param(5)
          alpha_5 = sliphard_param(6)
          h_5 = sliphard_param(7)
          hD_5 = sliphard_param(8)          
          C_5 = sliphard_param(9)
c            
c
c          
c         Slip distance will be sent as an input to this subroutine
c         Slip is not actually the statedot, statedot is used just to avoid new
c         definition of input variables
c         Xdist is nonlocal slip distance, that has a different value for each IP
c
c          
c         The RHS term represents the Taylor hardening term: alpha*G*b/L scaled with a constant C          
          mKM_5 = K_5 * dsqrt(rhoSSD) - 2.0d+0 * yc_5 * rhoSSD
     & + C_5 * alpha_5 * G * b_5 / sXdist
c          
c         self-hardening          
		sstatedot(1) = mKM_5 * dabs(gdot) / b_5
c          
c          
c         kinematic hardening
          sstatedot(2) = h_5*gdot - hD_5*x*dabs(gdot)          
c          
c          
c          write(6,*) 'inside sliphar'
c          write(6,*) 'selfhrd'
c         Other slip hardening laws come here!!!
c
c          write(6,*) 'sstatedot(1)', sstatedot(1)
c
c           write(6,*) 'sstate(1)', sstate(1)          
c          
c          
      endif
c      
c      
c      
c      
c      
c      
c      
c      
c      
      return
      end subroutine sliphard
c      
c      
c      
c      
c      
c      
c      
c      
      end module sliphardlaws