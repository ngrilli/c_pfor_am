c Chris Allen
c Edward Horton
c Eralp Demir
c Aug. 12th, 2021 - 1st working version
c
	module slipratelaws
      implicit none
      contains
      
      
      
      subroutine sliprate(tau, sstate, temp, gdot, dgdot_dtau)
      use globalvars, only: modelno, sliprate_param, numstvar
      implicit none
      real(8) tau, tauc, x, temp, gdot, dgdot_dtau
c     model constants for slip/creep     
      real(8) gdot0_1, m_1, crit_1
c     model constants for slip&creep with backstress
      real(8) gdot0s_2, gdot0c_2, ms_2, mc_2, crit_2
      real(8)	sstate(numstvar)
      
c     Slip/Creep  with no kinematic hardening (x = 0), s it has no effect
      if (modelno.eq.1) then
          
c          write(6,*) 'sliprate'
c          write(6,*) sstate
          
          tauc = sstate(1)          
          
          
          gdot0_1 = sliprate_param(1)
          
          m_1 = sliprate_param(2)
          
   
      gdot=gdot0_1*((dabs(tau)/tauc)**(1.0d0/m_1))*dsign(1.0d0,tau)

      dgdot_dtau=gdot0_1/m_1*((dabs(tau)/tauc)**((1./m_1)-1.))/tauc
          



c     Slip and Creep with kinematic hardening
      elseif (modelno.eq.2) then
          
          

    
          
          tauc = sstate(1)
          x = sstate(2)
          
          
          gdot0s_2 = sliprate_param(1)
          
          ms_2 = sliprate_param(2)
          
          gdot0c_2 = sliprate_param(3)
          
          mc_2 = sliprate_param(4)          
          
   
      gdot=gdot0s_2*((dabs(tau-x)/tauc)**(1.0d0/ms_2))
     + *dsign(1.0d0,tau-x) +
     + gdot0c_2*((dabs(tau-x)/tauc)**(1.0d0/mc_2))*dsign(1.0d0,tau-x)

      dgdot_dtau=gdot0s_2/ms_2*((dabs(tau-x)/tauc)**(1./ms_2-1.))
     + /tauc + gdot0c_2/mc_2*((dabs(tau-x)/tauc)**(1./mc_2-1.))/tauc
              
              

c          write(6,*) 'gdot', gdot

c     Other slip rate laws come here!!!


      endif
          
          
      
      
      
      
      
      
      
      
      
      
      
      
      return
      end subroutine sliprate
      
      
      
      
      
      
      
      
      
      end module slipratelaws