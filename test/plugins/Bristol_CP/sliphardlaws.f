c Chris Allen
c Edward Horton
c Eralp Demir
c Aug. 12th, 2021 - 1st working version
c
	module sliphardlaws
      implicit none
      contains
      
      
      
      subroutine sliphard(sstate, gsum, gdot, temp, sstate0, sstatedot)
      use globalvars, only: modelno,sliphard_param,Rgas,numstvar,
     &numslip
      implicit none
      real(8) gsum, gdot, temp
      real(8)	sstatedot(numstvar+1), sstate(numstvar), sstate0(numstvar)
c     variables for Voce hardening law      
      real(8) h0_1, ss_1, n_1
      real(8) hb_1
c     variables for is0tropic and kinematic hardening        
      real(8) tauc0_2, h0_2, n_2, A_2, Q_2, d_2, h_2, hD_2
      real(8) tauc, x
      
      
      sstatedot = 0.0d+0
      
      
c     Three parameter hardening law      
      if (modelno.eq.1) then
          

          

          
          tauc = sstate(1)
          
          h0_1 = sliphard_param(2)
          ss_1 = sliphard_param(3)
          n_1 = sliphard_param(4)
          
          hb_1 = h0_1*((1.-(tauc/ss_1))**n_1)
          
c         self-hardening          
		sstatedot(1) = hb_1*dabs(gdot)
          
          
c          write(6,*) 'inside sliphar'
c          write(6,*) 'selfhrd'
c     Other slip hardening laws come here!!!

c          write(6,*) 'sstatedot(1)', sstatedot(1)

c           write(6,*) 'sstate(1)', sstate(1)

      elseif (modelno.eq.2) then
          
    
          
          
          tauc = sstate(1)
          x = sstate(2)
          
          tauc0_2 = sstate0(1)
          h0_2 = sliphard_param(2)
          n_2 = sliphard_param(3)
          Q_2 = sliphard_param(4)
          d_2 = sliphard_param(5)
          A_2 = sliphard_param(6)
          h_2 = sliphard_param(7)
          hD_2 = sliphard_param(8)
          
          
c         self-hardening          
      sstatedot(1) = h0_2*dabs(gdot)*(1.+
     + (h0_2*gsum/tauc0_2/n_2))**(n_2-1.)
          
c         Thermal softening
          sstatedot(3) = -A_2 * (tauc**d_2) * dexp(-Q_2/Rgas/temp) 
 
c         kinematic hardening
          sstatedot(2) = h_2*gdot - hd_2*x*dabs(gdot)
          


      endif
      
      
      
      
      
      
      
      
      
      return
      end subroutine sliphard
      
      
      
      
      
      
      
      
      end module sliphardlaws