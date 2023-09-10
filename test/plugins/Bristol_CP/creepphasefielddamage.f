
c Edward Horton

c 
c Aug. 12th, 2021 - 1st working version
c
      module creepphasefielddamage
        implicit none
        contains
    
    
c     Update Creep Damage
    
        subroutine update_creep_damage(F_t,F,S,f_ep_c,dt,ph_no)

        use globalsubs, only: convert3x3to6
        use globalvars, only: I3
        
        !Deformation Gradients
        real(8), intent(in) :: F_t(3,3), F(3,3)

        ! Cauchy stress in vector form
        real(8), intent(in) :: S(6)

        !Cauchy Tensors 
        real(8) ::A(3,3),A_t(3,3)
        
        !Green-Lagrange strain tensors from beginning and end of increment.
        real(8) :: E(3,3), E_t(3,3)
        real(8) :: E_vec(6), E_t_vec(6)
        !2nd Piola Kirchoff Stress
        real(8) :: dt,PK2(3,3),PK2_vec(6)
        !Creep degradation function
        real(8) :: f_ep_c

        integer :: ph_no


    ! Calculate strain from deformation gradient... Currently using plastic. use total?
        A = matmul(transpose(F),F)
        A_t = matmul(transpose(F_t),F_t)
        E = 0.5d+0*(A-I3)
        E_t = 0.5d+0*(A_t-I3)
        ! Convert strains to vector
        call convert3x3to6(E,E_vec)
        call convert3x3to6(E_t,E_t_vec)

        call R5_SMDE(f_ep_c,S,E_vec,E_t_vec,dt,ph_no)

        return

        end subroutine update_creep_damage
        
            
c ==========    R5 - SMDE Model:     https://doi.org/10.1179/mht.2004.007  
        subroutine R5_SMDE(f_ep_c,sigma,E_vec,E_vec_t,dt,ph_no)
        
        use globalvars, only: temp0,Rgas
        !Cauchy Stress
        real(8) sigma(6)
        !Invariants Stresses
        real(8) sigmaMaxP,sigmaEq,sigmaH
        !Green-Lagrange Strain:
        real(8) E_vec(6),E_vec_t(6)
        real(8) p_R5,q_R5, A,Q,n,m
        
        ! Creep degradation function
        real(8) f_ep_c0,f_ep_c
        !various strains
        real(8) epsilon_c, epsilon_f,  epsilon_fu
        !Used for calculating Von Mises creep rate.
        real(8) E_MaxP_t,E_MaxP,epsilon_c_dot,E_Eq, E_Eq_t
        
        real(8) dummy1,dummy2,dt

        real(8) dmg0,dmg,dmg_inc
        
        !real(8) failure_strain(2),e_upper
        integer i,ph_no
        
        
        ! ==== CONSTANTS ===== for 316
        p_R5 = 0.15
        q_R5 = 1.25
        !p_R5 = c_damage_param(ph_no,1)
        !q_R5 = c_damage_param(ph_no,2)
        
        !From: doi:https://doi.org/10.1179/mht.2004.007 (For 304)
        A = 70d+0
        Q = 51123.0
        n= 0.2971
        m = 1.2695
        !From NIMS Data - 316H @ 650C
        !A = c_damage_param(ph_no,3)
        !Q = c_damage_param(ph_no,4)
        !n = c_damage_param(ph_no,5)
        !m = c_damage_param(ph_no,6)

        ! ==== CALCULATIONS =====
        f_ep_c0 = f_ep_c
        dmg0 = 1-f_ep_c
    
        
        
        !==== Find stresses.
        call InvariantStresses(sigma,sigmaMaxP,sigmaEq,sigmaH)
        
        
        !==== creep strain
        
        call InvariantStresses(E_vec,E_MaxP,E_Eq,dummy2)
        call InvariantStresses(E_vec_t,E_MaxP_t,E_Eq_t,dummy2)
        !Von Mises creep strain...
        epsilon_c = E_Eq
        ! If creep strain is negative then there will be no damage added.
        if (epsilon_c .lt. 0d+0) then
            epsilon_c = 0d+0
        endif
        
        
        
        !==== Von Mises creep strain rate: https://doi.org/10.1179/mht.2004.007
        epsilon_c_dot = (E_Eq-E_Eq_t)/dt
        !If negative then no damage.
        if  (epsilon_c_dot .le. 0d+0) then
                epsilon_c_dot = 0d+0
        endif
        
        
        !==== Uniaxial Failure Strain: https://doi.org/10.1179/mht.2004.007
            epsilon_fu = A*exp(Q/(Rgas*temp0))*(epsilon_c_dot**n)
     1*(sigmaMaxP**(-1*m))


            if (epsilon_fu .lt. 1e-6) then
                epsilon_fu = 0.9
            endif
            
            
        !==== Multiaxial Failure Strain
        epsilon_f=  epsilon_fu * 
     1exp(p_R5*(1d+0 - (sigmaMaxP/sigmaEq))
     2+q_R5*(0.5d+0 - 1.5d+0*(sigmaH/sigmaEq)))
        
        

        !==== Creep degradation function
            dmg_inc = (epsilon_c_dot*dt)/epsilon_f
            dmg = dmg0+dmg_inc
            f_ep_c = 1-dmg
            !f_ep_c = epsilon_f/(epsilon_f + epsilon_c)
            
            
c     Damage cannot recover so stop the degradation function increasing: 
        if (f_ep_c .gt. f_ep_c0) then
            f_ep_c = f_ep_c0
        end if


        return
        end subroutine R5_SMDE
        
        
        
        
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