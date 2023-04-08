
c Edward Horton

c 
c Aug. 12th, 2021 - 1st working version
c
      module creepphasefielddamage
      
      implicit none
      contains
      
c     R5 Model:     https://doi.org/10.1179/mht.2004.007  
      subroutine R5(f_ep_c,sigma,Strain_vec,Old_Strain_vec,dt)
      
      use globalvars, only: temp0
      
      real(8) sigma(6),dt,Strain_vec(6),Old_Strain_vec(6)
      real(8) p_R5,q_R5, A,P,n,m
      real(8) sigmaMaxP,sigmaEq,sigmaH
      real(8) f_ep_c0,f_ep_c
      real(8) PI
      real(8) epsilon_c, epsilon_f, epsilonMaxP, epsilon_fu
      real(8) Old_epsilonMaxP,epsilon_c_dot
      real(8) dummy1,dummy2
      integer i
       
      
      
      p_R5 = 0.15d+0
      q_R5 = 1.25d+0
      
      !From: doi:10.1016/j.ijpvp.2004.09.003 - For 347 Not 316 so need changing...
      A = 0.2723d+0
      P = 16268d+0
      n=0.4896d+0
      m = 2.873d+0
      !Test Values for checking Invariant Stresses Subroutine   
      !sigma(1) = -19
      !sigma(2) = 4.6
      !sigma(3) = -8.3
      !sigma(4) = -4.7
      !sigma(5) = 6.45
      !sigma(6) = 11.8
      !roots should be: 11.6178, -25.3163, -9.0015 and sigma_vm should be: 32.058
      ! keep initial damage variable.
      f_ep_c0 = f_ep_c
      
      PI = 4.D0*DATAN(1.D0)
      
      !These need to be added as inputs... maybe
 
      
      ! Find stresses.
       call InvariantStresses(sigma,sigmaMaxP,sigmaEq,sigmaH)
        
      ! creep strain
      call InvariantStresses(Strain_vec,epsilonMaxP,dummy1,dummy2)
      call InvariantStresses(Old_Strain_vec,Old_epsilonMaxP,
     1dummy1,dummy2)
      
      
      epsilon_c = abs(epsilonMaxP)
         if (isnan(epsilon_c)) then
          epsilon_c = 0d+0
         endif
         
         
         ! creep rate:
         epsilon_c_dot = epsilonMaxP-Old_epsilonMaxP
         !If negative then no damage.
         if  (epsilon_c_dot .le. 0d+0) then
              epsilon_c_dot = 0d+0
         endif
      ! Uniaxial Failure Strain
          epsilon_fu = A*exp(P/temp0)*(epsilon_c_dot**n)
     1*(sigmaMaxP**(-1*m))
          !write(6,*) ' epsilon_fu',epsilon_fu
          !Inwrite(6,*) ' epsilon_c_dot',epsilon_c_dot
          if (epsilon_fu .lt. 1e-6) then
              epsilon_fu = 0.5
          endif
          
      ! Multiaxial Failure Strain
         epsilon_f =  epsilon_fu * 
     1exp(p_R5*(1d+0 - (sigmaMaxP/sigmaEq))
     2+q_R5*(0.5d+0 - 1.5d+0*(sigmaH/sigmaEq)))
         
              ! write(6,*) 'epsilon_f', epsilon_f

          f_ep_c = epsilon_f/(epsilon_f + epsilon_c)
          
          
c     Damage cannot go backwards so: 
      if (f_ep_c .gt. f_ep_c0) then
          f_ep_c = f_ep_c0
      end if
      ! write(6,*) 'f_ep_c0', f_ep_c0
      ! write(6,*) 'f_ep_c',f_ep_c
      return
      end subroutine R5
      
      
      
      
!     Calculating max principal stresses and von mises stresses.
      subroutine InvariantStresses(sigma,sigmaMaxP,sigmaEq,sigmaH)
      
      real(8) sigma(6)
      real(8) I1,I2,I3,J2
      real(8) sigmaMaxP,sigmaEq,sigmaH
      real(8) Rc,Qc,PI,roots(3), theta
      integer i
      
      PI = 4.D0*DATAN(1.D0)
          !Invariants:
              I1 = sigma(1) + sigma(2) + sigma(3)
          
              I2 = (sigma(1)*sigma(2)) + (sigma(2)*sigma(3))
     &+(sigma(3)*sigma(1)) - sigma(4)**2 - sigma(5)**2 - sigma(6)**2
          
              I3 = (sigma(1)*sigma(2)*sigma(3)) +(2*sigma(4)
     &*sigma(5)*sigma(6)) - (sigma(1)*sigma(6)**2) - 
     &(sigma(2)*sigma(5)**2) - (sigma(3)*sigma(4)**2)
              !write(6,*) 'I1:',I1
              !write(6,*) 'I2:',I2
              !write(6,*) 'I3:',I3
    !     Von Mises Stress:
              J2 = I2 - (I1**2d+0)/3d+0
              sigmaEq =(-3d+0*J2)**0.5
              
    !          write(6,*) 'Von Mises:',sigmaEq
    !     Max Principal Stress:
    !     Finding roots of cube: https://www.continuummechanics.org/principalstress.html

              Qc = (3d+0*I2 - I1**2)/9d+0 
              Rc = ((2d+0*I1**3d+0)-(9d+0*I1*I2)+(27d+0*I3))/54d+0
          
              theta = acos(Rc/(-Qc**3d+0)**0.5) 
!        write(6,*)'theta:',theta
          roots(1) = 2*(-Qc)**0.5*cos(theta/3d+0) + (1d+0/3d+0 * I1)
          roots(2) = 2*(-Qc)**0.5*cos((theta+2*PI)/3d+0) 
     &+ (1d+0/3d+0 * I1)
          roots(3) = 2*(-Qc)**0.5*cos((theta+4*PI)/3d+0) 
     &+ (1d+0/3d+0 * I1)
          
          sigmaMaxP = -1.0d+57
          do i=1,3
              if (roots(i) > sigmaMaxP) then
                  sigmaMaxP = roots(i)
              endif
          enddo
          ! ============= Hydrostatic stress:
                sigmaH = I1/3
      end subroutine InvariantStresses
      
      
      end module creepphasefielddamage
