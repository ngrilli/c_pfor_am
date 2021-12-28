c Chris Allen
c Edward Horton
c Eralp Demir
c Aug. 12th, 2021 - 1st working version
c
	module slipratelaws
      implicit none
      contains
      
      
      
      subroutine sliprate(tau, sstate, temp, gdot, dgdot_dtau)
      use globalvars, only: modelno, sliprate_param, sliphard_param, 
     &numstvar, G
      implicit none
      real(8) tau, tauc, x, temp, gdot, dgdot_dtau
      real(8)	sstate(numstvar)
c     model constants for slip/creep     
      real(8) gdot0_1, m_1, crit_1
c     model constants for slip&creep with backstress
      real(8) gdot0s_2, gdot0c_2, ms_2, mc_2, crit_2
c     model constants for slip with backstress (Code Aster model)
      real(8) K_3, n_3, pl_3
c     model constatns for GND using slip gradients
      real(8) gdot0_4, m_4, alpha_4, b_4, tauc0_4
      
      
c     Slip/Creep  with no kinematic hardening (x = 0), s it has no effect
      if (modelno.eq.1) then
          
c          write(6,*) 'sliprate'
c          write(6,*) sstate
          
          tauc = sstate(1)          
          
          
          gdot0_1 = sliprate_param(1)
          
          m_1 = sliprate_param(2)
          
   
      gdot=gdot0_1*((dabs(tau)/tauc)**(1./m_1))*dsign(1.0d+0,tau)

		    dgdot_dtau=gdot0_1/m_1*((dabs(tau)/tauc)**((1./m_1)-1.))/tauc
          



c     Slip and Creep with kinematic hardening
      elseif (modelno.eq.2) then
          
          

    
          
          tauc = sstate(1)
          x = sstate(2)
          
          
          gdot0s_2 = sliprate_param(1)
          
          ms_2 = sliprate_param(2)
          
          gdot0c_2 = sliprate_param(3)
          
          mc_2 = sliprate_param(4)          
          
   
      gdot=gdot0s_2*((dabs(tau-x)/tauc)**(1./ms_2))*dsign(1.0d+0,tau-x)+
     + gdot0c_2*((dabs(tau-x)/tauc)**(1./mc_2))*dsign(1.0d+0,tau-x)

          dgdot_dtau=gdot0s_2/ms_2*((dabs(tau-x)/tauc)**(1./ms_2-1.))
     &/tauc + gdot0c_2/mc_2*((dabs(tau-x)/tauc)**(1./mc_2-1.))/tauc
              
              

c     Slip with  kinematic hardening (Code Aster)
      elseif (modelno.eq.3) then
          
c          write(6,*) 'sliprate'
c          write(6,*) sstate
          
          tauc = sstate(1)          
          x = sstate(2)
          
          K_3 = sliprate_param(1)
          
          n_3 = sliprate_param(2)
          
          pl_3 = (dabs(tau-x)-tauc)/K_3
          
          
          if (pl_3.le.0.0) then
          
              gdot = 0.0
              
              dgdot_dtau = 0.0
              
          else
          
      gdot=(pl_3**n_3)*dsign(1.0d+0,tau-x)

		    dgdot_dtau = pl_3**(n_3-1.)/K_3/n_3
                  
          
          endif
          
   

c     Slip/Creep considering GNDs with slip gradients
      elseif (modelno.eq.4) then
          
c          write(6,*) 'sliprate'
c          write(6,*) sstate
          


          
          gdot0_4 = sliprate_param(1)
          
          m_4 = sliprate_param(2)
          
          
c         Geometric factor for slip resistance calculation
          alpha_4 = sliphard_param(4)
          

          
c         Burgers vector for slip resistance calculation
          b_4 = sliphard_param(3)
          
c         Initial slip resistance
          tauc0_4 = sliphard_param(1)
          
          
c         Calculate slip resistance          
          tauc = tauc0_4 + alpha_4 * G * b_4 * 
     & dsqrt( sstate(1) + dabs(sstate(2)) + dabs(sstate(3)) )
          
c          write(6,*) 'tauc', tauc
          
c         x is not a state but it will be derived from states
c         calculate backstress here!          
          x = sstate(4)
          
c          write(6,*) 'x', x
          
          
c          if (dabs(tau).gt.dabs(x)) then
      gdot=gdot0_4*((dabs(tau-x)/tauc)**(1./m_4))*dsign(1.0d+0,tau-x)

             dgdot_dtau=gdot0_4/m_4*((dabs(tau-x)/tauc)**((1./m_4)-1.))
     &/tauc 
c          else
c              gdot= 0.0
c              dgdot_dtau= 0.0
c          endif              
              

c          write(6,*) 'gdot', gdot

c     Other slip rate laws come here!!!


      endif
          
          
      
      
      
      
      
      
      
      
      
      
      
      
      return
      end subroutine sliprate
      
      
      
      
      
      
      
      
      
      end module slipratelaws