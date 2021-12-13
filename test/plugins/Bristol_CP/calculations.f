c Chris Allen
c Edward Horton
c Eralp Demir
c Aug. 12th, 2021 - 1st working version
c
	module calculations
      implicit none
	contains
c
c	This subroutine calculates the two main variables: Stress and Consistent tangent	
      subroutine calcs(F_t,F,t,dt,temp,inc,el_no,ip_no,
     + sigma,jacob,pnewdt,coords)
      use globalvars, only: global_Fp,global_Fp_t,global_Fe,
     + global_Fe_t,
     + global_state,global_state_t,I3,inc_old,elas66,
     + global_ori,njaco,innoitmax,ounoitmax,
     + global_sigma,global_jacob_t,global_jacob,TF,phaseID,numslip,
     + global_gammadot,global_S,global_S_t,t_old,ratio_lb,ratio_ub,
     + dgamma_s,global_gamma_t,global_gamma,global_sigma_t,
     + global_gamma_sum_t,global_gamma_sum,numstvar,thermo,temp0,
     + GSeffect,grainID,grainsize_init,global_state0,tstep_forw,
     + tstep_back,numel,numip
      use globalsubs, only: convert3x3to6,invert3x3,skew,misorientation,
     + polar
      use initialization, only: initialize_grainsize
	implicit none
c	Inputs
      integer inc
      real(8) F_t(3,3),F(3,3),t,dt,temp
c	Outputs
      real(8) sigma(6),jacob(6,6),pnewdt,coords(3)
c	Variables used within this subroutine
      real(8) Fp_t(3,3),Fe_t(3,3),tauc_t(numslip)
      real(8) Fp(3,3),Fe(3,3),tauc(numslip),Cauchy(3,3)
      real(8) S(6),S_t(6),det,C(3,3,numslip),epsdot
      real(8) Lp(3,3),R(3,3),U(3,3),gammadot(numslip)
      real(8) dgammadot_dtau(numslip)
      integer is,el_no,ip_no,gr_no
      real(8) g1(3,3),ang,ax(3),dg(3,3)
      real(8) ratio,dgamma_max,gsum_t,gsum
      integer sconv, jconv, i, j
      real(8) state(numslip,numstvar),state_t(numslip,numstvar)
      real(8) state0(numslip,numstvar)
      real(8)	sstate(numstvar),sstate_t(numstvar),sstate0(numstvar)
      
      
      
c     - The step update is indicated with the change of time
c     - Change in INC does not work here!!! Because, during sub time stepping,
c       INC also increases while time remains the same
	if ((t.gt.t_old).or.(inc.ne.inc_old)) then
          write(6,*) 'Debug: Updated the state variables'
		  write(6,*) 'Time'
		  write(6,*) t
		  write(6,*) t_old
		  write(6,*) inc
		  write(6,*) inc_old
          t_old=t
		inc_old=inc
		global_Fp_t=global_Fp
		global_Fe_t=global_Fe
          global_S_t=global_S
          global_gamma_t=global_gamma
          global_gamma_sum_t=global_gamma_sum
		global_state_t=global_state
          global_jacob_t=global_jacob
          global_sigma_t=global_sigma
              
      endif
      
c     If length scale calculation is "ON"
      if (GSeffect.eq.1) then
c         If not initialized
          if (grainsize_init(el_no,ip_no).eq.0) then

              
c             flag for initialization
              grainsize_init(el_no,ip_no)=1

          
c              write(6,*) 'element no.: ', el_no
          
c              write(6,*) 'IP no.: ', ip_no
          
              gr_no = grainID(el_no)
          
c              write(6,*) 'grain no.: ', gr_no
           
          
c              write(6,*) 'coordinates: ', coords
          
              call initialize_grainsize(el_no,ip_no,gr_no,coords)
          
          endif
      endif
         
      

      
c      write(6,*) 'numel',numel
c      write(6,*) 'el_no',el_no
c      write(6,*) 'ip_no',ip_no
c      write(6,*) 'state_t(3)', global_state_t(el_no,ip_no,3,1)
      
      
      
      
      

      
c     In case of mechanical solver only, use initial temperatures      
      if (thermo.eq.0) then
          temp = temp0    
      endif
      
      
      
      
      
c              


      
      
c     Isotropic material calculations          
      if (phaseID(el_no).eq.0) then
          
          
          
          
c          write(6,*) 'J2 CALCULATIONS'
          
          
          Fp_t=global_Fp_t(el_no,ip_no,:,:)
          Fe_t=global_Fe_t(el_no,ip_no,:,:)
          gsum_t=global_gamma_sum_t(el_no,ip_no)
          
          
          sstate0=global_state0(el_no,ip_no,1,:)
          
          sstate_t=global_state_t(el_no,ip_no,1,:)
          


          
                    
          
c	    Calculate stress and flow stress using J2 plasticity
          call J2_main(dt,F,Fp_t,Fe_t,sstate_t,gsum_t,temp,sstate0,
     &	Fp,Fe,epsdot,sstate,gsum,sigma,sconv)
     


c         Calculate the time factor
c           pnewdt=1.25
          pnewdt=ratio_ub
     
     
     	    if (sconv.eq.0) then
		    write(6,*) 'Inner/Outer loop during stress calculation
     &	did not converge!'
		    write(6,*) 'el_no'
		    write(6,*) el_no
		    write(6,*) 'ip_no'
		    write(6,*) ip_no
		    write(6,*) 'F'
		    write(6,*) F
		    write(6,*) 'F_t'
		    write(6,*) F_t
		    write(6,*) 'Fp_t'
		    write(6,*) Fp_t
 		    write(6,*) 'sigma'
		    write(6,*) sigma          
		    write(6,*) 'state_t'
		    write(6,*) state_t
		    write(6,*) 'dt'
		    write(6,*) dt
cc		    QUIT ABAQUS
c		    call xit



c             Introduce cut-back              
c              pnewdt=0.5
              pnewdt=tstep_back
              write(6,*) 'pnewdt',pnewdt
              
c             Do not calculate the jacobian
              jacob=global_jacob_t(el_no,ip_no,:,:)

c             Assign the old value of stress   
              sigma = global_sigma_t(el_no,ip_no,:)

c




c         Regular time stepping
          else

c              jacob=global_jacob_t(el_no,ip_no,:,:)

     
              gammadot = epsdot*TF
c             Strain increment
              dgamma_max = maxval(dabs(gammadot))*dt
      
      
c             Time factor
              ratio = dgamma_max / dgamma_s
      
      
c             calculate the factor
              if (ratio.lt.ratio_lb) then
              
c                  pnewdt = 1.1
                  pnewdt = tstep_forw
              
              elseif (ratio.gt.ratio_ub) then
              
c                  pnewdt = 0.75
                  pnewdt = ratio_lb
              
              elseif ((ratio.ge.ratio_lb).or.(ratio.le.ratio_ub)) then
              
c                  pnewdt = 1.0/ratio
                  pnewdt = 1.0
              
              endif              


   



     
   
          
     
          
          

c             JACOBIAN CALCULATION FOR ISOTROPIC PHASE     
c	        For a number of increments do the jacobian calculation
c	        Note, jacobian is needed at the first calculation
	        if (modulo(inc,njaco).eq.0) then
c		    Calculate the material tangent (using perturbation)

    
      call J2_jacobian(dt,F_t,F,Fe_t,Fp_t,sstate_t,gsum_t,temp,
     + sstate0,sigma,jacob,jconv)
              
              
              
              

c                 jacob=zeta66

c		        When the jacobian calculation did not converge, assign the old jacobian!
		        if (jconv.eq.0) then
c			        Jacobian 			
                      write(6,*) 'Jacobian has not converged - J2'
                      write(6,*) 'el_no',el_no
                      write(6,*) 'ip_no',ip_no
			        jacob=global_jacob_t(el_no,ip_no,:,:)
			        
                      
c                      pnewdt=0.5
                      pnewdt=tstep_back
                      write(6,*) 'pnewdt',pnewdt
                      
		        endif
              else
                  
c                 Note this also works when inc=1 since it is elasticity matrix
		        jacob=global_jacob_t(el_no,ip_no,:,:)

              endif
     
     
     
c          write(6,*) 'jacob', jacob
     
     
          endif

          
                
c         Update the state variables
          global_Fp(el_no,ip_no,:,:) = Fp
          global_Fe(el_no,ip_no,:,:) = Fe
          global_S(el_no,ip_no,:) = 0.0d+0
c         Store the important variables
c         For J2 model assign the same values for all slip sytems    
          do i=1,numslip
              global_state(el_no,ip_no,i,1:numstvar)=sstate(1:numstvar)
              global_gammadot(el_no,ip_no,i) = epsdot*TF
              global_gamma(el_no,ip_no,i) = gsum
          enddo
          
          global_gamma_sum(el_no,ip_no) = gsum 

      
      
      
      
c     Single crystal calculations           
      else
c          
c          write(6,*) 'SINGLE CRYSTAL CALCULATIONS'      
c
c	    Assign the globally stored variables
	    Fp_t = global_Fp_t(el_no,ip_no,:,:)
	    S_t = global_S_t(el_no,ip_no,:)
	    state_t = global_state_t(el_no,ip_no,:,:)
          state0 = global_state0(el_no,ip_no,:,:)
          gsum_t = global_gamma_sum_t(el_no,ip_no)
          
          
c
c
c	    Calculate stress and shear resistance
c	    Note: el_no and ip_no are needed to get the values of Schmid vectors and
c	    elasticity tensor from the global variables
      call SC_main(dt,F,Fp_t,S_t,state_t,gsum_t,temp,state0,
     + C,S,Lp,Fp,Fe,sigma,gammadot,dgammadot_dtau,state,
     + gsum,sconv)
c

      
c         Calculate the time factor
c          pnewdt=1.25
          pnewdt=ratio_ub
       


	    if (sconv.eq.0) then
		    write(6,*) 'Inner or outer loop during stress calculation
     &	    did not converge!'
		    write(6,*) 'inc'
		    write(6,*) inc          
		    write(6,*) 'el_no'
		    write(6,*) el_no
		    write(6,*) 'ip_no'
		    write(6,*) ip_no
		    write(6,*) 'F'
		    write(6,*) F
		    write(6,*) 'F_t'
		    write(6,*) F_t
		    write(6,*) 'Fp_t'
		    write(6,*) Fp_t
		    write(6,*) 'Fp'
		    write(6,*) Fp
		    write(6,*) 'S_t'
		    write(6,*) S_t              
		    write(6,*) 'sigma'
		    write(6,*) sigma         
		    write(6,*) 'state_t'
		    write(6,*) state_t
		    write(6,*) 'gammadot'
		    write(6,*) gammadot
		    write(6,*) 'dgammadot_dtau'
		    write(6,*) dgammadot_dtau              
		    write(6,*) 'dt'
		    write(6,*) dt
              
c             Introduce cut-back              
c              pnewdt=0.5
              pnewdt=tstep_back
              write(6,*) 'pnewdt',pnewdt
              
c             Do not calculate the jacobian
              jacob=global_jacob_t(el_no,ip_no,:,:)
              
c             Assign the old value of stress              
              sigma = global_sigma_t(el_no,ip_no,:)
              
              
c		    QUIT ABAQUS
c		    call xit


c         Regular time stepping
          else

c              jacob=global_jacob_t(el_no,ip_no,:,:)

     
c             Strain increment
              dgamma_max = maxval(dabs(gammadot))*dt
      
      
c             Time factor
              ratio = dgamma_max / dgamma_s
      
      
c             calculate the factor
              if (ratio.lt.ratio_lb) then
              
c                  pnewdt = 1.1
                  pnewdt = tstep_forw
              
              elseif (ratio.gt.ratio_ub) then
              
c                  pnewdt = 0.75
                  pnewdt = ratio_lb
              
              elseif ((ratio.ge.ratio_lb).or.(ratio.le.ratio_ub)) then
              
c                  pnewdt = 1.0/ratio
                  pnewdt = 1.0
              
              endif              
              
c              write(6,*) 'pnewdt',pnewdt
         
c	        For a number of increments do the jacobian calculation
c	        Note, jacobian is needed at the first calculation
      if (modulo(inc,njaco).eq.0) then
c		        Calculate the material tangent (using perturbation)
      call SC_jacobian(dt,F_t,F,S_t,Fp_t,state_t,
     + gsum_t,temp,state0,
     + sigma,jacob,jconv)
      
c		        When the jacobian calculation did not converge, assign the old jacobian!
		        if (jconv.eq.0) then
                      write(6,*) 'Jacobian has not converged - CP!'
		            write(6,*) 'inc'
		            write(6,*) inc                       
                      write(6,*) 'el_no',el_no
                      write(6,*) 'ip_no',ip_no
                  
c                      pnewdt=0.5
                      pnewdt=tstep_back
                      write(6,*) 'pnewdt',pnewdt
                  
c			        Jacobian 			
			        jacob=global_jacob_t(el_no,ip_no,:,:)
                  
c                     do i=1,6
c                         write(6,*) (jacob(i,j),j=1,6)
c                     enddo
                   

		        endif
              else
c		        Note this also works when inc=1 since it is elasticity matrix
		        jacob=global_jacob_t(el_no,ip_no,:,:)

              endif


	    endif
c
!	    if (ounoit.ge.ounoitmax) then
!		    write(6,*) 'Outer loop during stress calculation
!     &	    did not converge!'
!		    write(6,*) 'inc'
!		    write(6,*) inc                
!		    write(6,*) 'el_no'
!		    write(6,*) el_no
!		    write(6,*) 'ip_no'
!		    write(6,*) ip_no
!		    write(6,*) 'inner no it'
!		    write(6,*) innoit
!		    write(6,*) 'outer no it'
!		    write(6,*) ounoit
!		    write(6,*) 'F'
!		    write(6,*) F
!		    write(6,*) 'F_t'
!		    write(6,*) F_t
!		    write(6,*) 'Fp_t'
!		    write(6,*) Fp_t
!		    write(6,*) 'Fp'
!		    write(6,*) Fp              
!		    write(6,*) 'sigma'
!		    write(6,*) sigma              
!		    write(6,*) 'tauc_t'
!		    write(6,*) tauc_t
!		    write(6,*) 'gammadot'
!		    write(6,*) gammadot
!		    write(6,*) 'dgammadot_dtau'
!		    write(6,*) dgammadot_dtau    
!		    write(6,*) 'dt'
!		    write(6,*) dt
!              pnewdt=0.5
!c		    QUIT ABAQUS
!c		    call xit
!          endif
      
          
          
          
c	    Store the important variables
	    global_Fp(el_no,ip_no,:,:) = Fp
	    global_Fe(el_no,ip_no,:,:) = Fe
          global_S(el_no,ip_no,:) = S
	    global_state(el_no,ip_no,:,:) = state
          global_gammadot(el_no,ip_no,:) = gammadot 
          global_gamma(el_no,ip_no,:) = global_gamma_t(el_no,ip_no,:) +
     &    gammadot*dt
          global_gamma_sum(el_no,ip_no) = gsum
          
          
c

          

          
   
          
          
          
          
      endif
      
      
 
      
      
          
c     Store the important results    
      global_sigma(el_no,ip_no,:)=sigma
      
c     Assign the value of jacobian even if there is no convergence               
      global_jacob(el_no,ip_no,:,:)=jacob      
           
          
      !endif
          
      

      
      
      
      
      

c      if (inc.eq.1) then
c          
c          write(6,*) 'jacob', jacob
c          write(6,*) 'elas66', elas66
c          
c      endif
      !
      !
      !if (inc.eq.10) then
      !    
      !    write(6,*) 'exit'
      !    call xit
      !    
      !endif
      
      
      return
      end subroutine calcs
c
c
c

c      
c      
c
      
c      
c      
c     Main routine for jacobian calculation for Martensite 
c	This subroutine calculates consistent tangent
      subroutine J2_jacobian(dt,F_t,F,Fe_t,Fp_t,sstate_t,
     + gsum_t,temp,sstate0,
     + Cauchy_vec,jacob,jconv)
	 
      use globalvars, only : deps,innoitmax,ounoitmax,numslip,numstvar
      use globalsubs, only : convert3x3to6,convert6to3x3
      implicit none
c	Inputs
      real(8) F_t(3,3),F(3,3),dt,gsum_t,temp
      real(8) Fp_t(3,3),Fe_t(3,3),Cauchy_vec(6)
c	Outputs
	real(8) jacob(6,6)
	integer jconv
c	Variables used within this subroutine
      real(8) Fp(3,3),Fe(3,3),epsdot,sigc
      real(8) Fper(3,3),dFrel_vec(6),dFrel(3,3),Cauchy_per_vec(6),gsum
      integer i,iloop,oloop,sconv
      real(8) sstate_t(numstvar),sstate(numstvar),sstate0(numstvar)
c
c	Assign the convergent behavior
	jconv=1
      
cc	Cauchy stress
c	call convert6to3x3(Cauchy_t_vec,Cauchy_t)      
c
c	Increment 6 components of relative deformation gradient
	do i=1,6
c		Component-wise pertubation
		dFrel_vec=0.0
          
c          This doubles the shear terms in the jacobian XXXX which is WRONG!!!          
!          dFrel_vec(i)=deps
          

c         This is TRUE as in Kalididi's study - gives "G" as the shear term
		if (i.le.3) then
              dFrel_vec(i)=deps
          else   
c		Note it is not deps/2 since during conversion only one component is considered
		    dFrel_vec(i)=deps/2.0
          endif
          
c		Convert the vector to a matrix
		call convert6to3x3(dFrel_vec,dFrel)
		Fper=F+matmul(dFrel,F_t)
c		Call the calculation procedure  
           call J2_main(dt,Fper,Fp_t,Fe_t,sstate_t,gsum_t,temp,sstate0,
     &	Fp,Fe,epsdot,sstate,gsum,Cauchy_per_vec,sconv)
           
           if (sconv.eq.0) jconv=0

		jacob(1:6,i)=(Cauchy_per_vec-Cauchy_vec)/deps
      enddo


c     Make it symmetric      
      jacob=(transpose(jacob)+jacob)/2.0

	return
	end subroutine J2_jacobian
c
c
      
      
      
      
      
      
      
c      
c
c     Main routine for J2-plasticity
      subroutine J2_main(dt,F,Fp_t,Fe_t,sstate_t,gsum_t,temp,sstate0,
     &Fp,Fe,epsdot,sstate,gsum,Cauchy_vec,sconv)
      use globalsubs, only: invert3x3, convert3x3to6, trace, normmat, 
     &convert6to3x3, determinant
      use globalvars, only: elas66_iso, G, E, nu, TF, modelno,
     &innertol, outertol, innoitmax, ounoitmax, I3, numstvar, dS_cr
      use slipratelaws, only: sliprate
      use sliphardlaws, only: sliphard
	implicit none
c	Inputs
      real(8) F(3,3),dt, temp
	real(8) Fp_t(3,3),Fe_t(3,3),gsum_t
c	Outputs
	real(8) Fp(3,3),Fe(3,3),epsdot,sigc,Cauchy_vec(6),gsum
      integer sconv
c	Vairables used    

      real(8) invFp_t(3,3),GLS(3,3),GLSvec(6),det
	real(8) GLStr(3,3),GLStrvec(6)
	real(8) Cauchytrvec(6),Cauchytr(3,3),Cauchytrdev(3,3)
      real(8) hyd,sigmatr,Cauchy(3,3)
	real(8) Nptr(3,3),Cauchydev(3,3),sigma,dR1
	real(8) gdot,dgdot_dtau,depsdot_dsigma
	real(8) R1,xx1,xx2,Dp(3,3),invFp(3,3),detFp
      real(8) tau,tautr,dsigma
      real(8)	sstate_t(numstvar),sstate(numstvar),sstate0(numstvar)
      real(8) R2(numstvar),sstatedot(numstvar+1)
      integer iloop,oloop,i
      logical notnum
c	real(8) xx1_old,sigma_old
 
c      
	call invert3x3(Fp_t,invFp_t,det)
c	Trial elastic deformation	
	Fe=matmul(F,invFp_t)     


      
      
      
c	Trial elastic stretch
	GLStr=(matmul(transpose(Fe),Fe)-I3)/2.0
      
      
      
     
c	Vectorize strain
	call convert3x3to6(GLStr,GLStrvec)     
     

c	Vectorized trial stress
	Cauchytrvec=matmul(elas66_iso,GLStrvec)      
      
c	Trial stress
	call convert6to3x3(Cauchytrvec,Cauchytr)
c	Hydrostatic component of the trial stress
	call trace(Cauchytr,3,hyd)
	hyd=hyd/3.0
c	Deviatoric component of the trial stress
	Cauchytrdev=Cauchytr-(I3*hyd)
c	Equivalent trial stress
	call normmat(Cauchytrdev,3,sigmatr)
	sigmatr=dsqrt(3.0d+0/2.0d+0)*sigmatr
      
      
cc     Convert to slip system values
c      tautr = sigmatr/TF
      
      
c	Plastic flow direction (if statement is to correct sigmabar=0 situation which give rise to 1/0 problem)
	if (sigmatr.ne.0.) then
		Nptr=3.0/2.0*Cauchytrdev/sigmatr
	else
		Nptr=0.0
	endif

      
      
!c     Convert to slip system values      
!      tauc=sigc/TF
!      tauc_t=sigc_t/TF
!      
      
      
      
      
c     COMPUTE FORMER STRESS      
c	Former elastic stretch
	GLS=(matmul(transpose(Fe_t),Fe_t)-I3)/2.0     
c	Vectorize strain
	call convert3x3to6(GLS,GLSvec)    
c	Vectorized former stress
	Cauchy_vec=matmul(elas66_iso,GLSvec)     
c     3x3 former stress tensor
c	Initial guess for the stress      
	call convert6to3x3(Cauchy_vec,Cauchy)      
      


c	Deviatoric stress
	call trace(Cauchy,3,hyd)
	hyd=hyd/3.0
	Cauchydev=Cauchy-(I3*hyd)
c	Equivalent stress
	call normmat(Cauchydev,3,sigma)
	sigma=dsqrt(3.0d+0/2.0d+0)*sigma


      
cc     Convert to slip system values      
c      tau = sigma/TF
      
c	Assign the initial value for the parameters before starting the loop
c	Initial guess for the shear resistance
	sstate=sstate_t      
      
c	OUTER loop starts here
	do oloop=1,ounoitmax
c		INNER loop starts here
cc		Assign a very large value
c		xx1_old=1.0d+10
		do iloop=1,innoitmax
              
              
              tau = sigma / TF
              
c			Plastic strain rate
              call sliprate(tau, sstate, temp, gdot, dgdot_dtau)
              
              depsdot_dsigma = dgdot_dtau / TF**2.0
c              write(6,*) 'sstate',sstate


c              write(6,*) 'iloop',iloop
c              write(6,*) 'gdot',gdot
c              write(6,*) 'tau',tau
              
			epsdot = gdot / TF
c			
c			if (gammadot.gt.1.0) then
c				write(6,*) '*******************'
c				write(6,*) 'Shear rates flying!'
c				write(6,*) '*******************'
c				call quit(9100)
c			endif
c			Derivative of shear rate with respect to the stress
c			depsdot_dsigma=dgdot_dtau/(TF**2.)

c			Residual
			R1=sigma-sigmatr+(dt*3.0*G*epsdot)
c              write(6,*) 'R1',R1
              xx1=dabs(R1)
cc			Relative norm of the residual
c			if (taubar.eq.0) then 
c				xx1=0.0
c			else
c				
c			endif
c			Check the residual
c			if residual is smaller than the tolerance EXIT
			if (xx1.lt.innertol) then
				exit
c			if the residual have a converging behavior
              endif
c             Calculate tangent for N-R iteration
              dR1=1.0+(dt*3.0*G*depsdot_dsigma)
              
c             Stress increment
              dsigma = -R1/dR1
              
c             If the stress increment is larger than the critical value              
              if (dabs(dsigma).gt.dS_cr) then
                  dsigma = dsign(1.0d0,dsigma)*dS_cr
              endif
              
cc            Store the old value of stress to update the guess
c             sigma_old=sigma
c             Stress after the iteration
              sigma=sigma + dsigma
cc            Assign the old value of the norm for the next iteration
c             xx1_old=xx1
              
		enddo
c
c		End of INNER loop
          gsum = gsum_t + gdot * dt

          
c         Strain hardening
          call sliphard(sstate, gsum, gdot, temp, sstate0, sstatedot)
              
        
          
          
c         Slip hardening law of Dylan 
          if (modelno.eq.2) then

                            
              sstatedot(1) = sstatedot(1) - sstatedot(3)
              
          
          endif          
          
          
          
          
c          write(6,*) 'oloop',oloop
c          write(6,*) 'gdot',gdot
c          write(6,*) 'sstate',sstate
c          write(6,*) 'sstatedot',sstatedot

c		Increment in shear resistance
c		Residual in the shear resistance
      R2=sstate(1:numstvar)-sstate_t(1:numstvar)-
     + sstatedot(1:numstvar)*dt

c		Absolute tolerance
		xx2=maxval(dabs(R2))
          
c		Check the tolerance
		if (xx2.lt.outertol) exit
          
c		Update the shear resistance
		sstate=sstate_t(1:numstvar)+sstatedot(1:numstvar)*dt
          
	enddo
c
c

cc     Convert slip system values to aggregate
c      epsdot = gdot/TF
c      sigma = tau*TF
c      sigc = tauc*TF




c	Plastic stretch tensor
	Dp=epsdot*Nptr
      
c	Plastic part of the deformation gradient
	Fp=matmul(((Dp*dt)+I3),Fp_t)
      
c     Find the determinant and scale it
      call determinant(Fp,detFp)
      
c     Plastic part of the deformation gradient      
      Fp = Fp / detFp**(1./3.)
      
c     Invert Fp
	call invert3x3(Fp,invFp,detFp)      
      
c	Elastic part of the deformation gradient
	Fe=matmul(F,invFp)
      
c	Elastic stretch
	GLS=(matmul(transpose(Fe),Fe)-I3)/2.0d+0


c	Vectorize strain
	call convert3x3to6(GLS,GLSvec)

      
c	Vectorized stress
	Cauchy_vec=matmul(elas66_iso,GLSvec)

      

      sconv=1d+0     
c     Check if the stress value converged         
      do i=1,6
          notnum = isnan(Cauchy_vec(i))
          if (notnum) sconv=0d+0
      enddo
      
c     Check if the stress value converged       
      do i=1,6
          
          if (abs(Cauchy_vec(i)).gt.1.0d+50) sconv=0d+0
          
      enddo
      
      
c     Check if the slip rates are infinite
      if (dabs(epsdot).gt.1.0d+50) sconv=0d+0
      
      
      
c     Check for the number of iterations      
      if (sconv.eq.1) then
          if (iloop.ge.innoitmax) sconv=0d+0
     
          
          if (oloop.ge.ounoitmax) sconv=0d+0
          
      endif
      
            
      
      
      return
      end subroutine J2_main
c     
c     
c      
c      
c            
      
      
      

c	This subroutine calculates consistent tangent
      subroutine SC_jacobian(dt,F_t,F,S_vec_t,Fp_t,state_t,
     + gsum_t,temp,state0,
     + Cauchy_vec,jacob,jconv)
	 
      use globalvars, only : deps,innoitmax,ounoitmax,numslip,numstvar
      use globalsubs, only : convert6to3x3, determinant
      implicit none
c	Inputs
      real(8) F_t(3,3),F(3,3),S_vec_t(6),Cauchy_vec(6),dt
	real(8) Fp_t(3,3),Cauchy(3,3),gsum_t,temp
c	Outputs
	real(8) jacob(6,6)
	integer jconv
c	Variables used within this subroutine
	real(8) detF,Cauchy_per(3,3)
	real(8) F_per(3,3),dFrel_vec(6),dFrel(3,3),Cauchy_per_vec(6)
	real(8) S_per_vec(6),Fe_per(3,3)
	real(8) invFp_per(3,3),dummy1(3,3,numslip),dummy2(6),dummy3(3,3)
      real(8) dummy4(3,3),dummy5(3,3)
	real(8) dummy6(numslip),dummy7(numslip),dummy8,dummy9,dummy10
	integer i,j,sconv
      real(8)	state(numslip,numstvar), state_t(numslip,numstvar)
      real(8) state0(numslip,numstvar)
      
cc     Determinant of DG
c      call determinant(F,detF)
            
      
c	Assign the convergent behavior
	jconv=1
c
c	Increment 6 components of relative deformation gradient
	do i=1,6
c		Component-wise pertubation
		dFrel_vec=0.0
		!dFrel_vec(i)=deps
		if (i.le.3) then
              dFrel_vec(i)=deps
          else   
c		Note it is not deps/2 since during conversion only one component is considered
		    dFrel_vec(i)=deps/2.0
          endif

c		Convert the vector to a matrix
		call convert6to3x3(dFrel_vec,dFrel)
		F_per=F+matmul(dFrel,F_t)
c		Call the calculation procedure
		call SC_main(dt,F_per,Fp_t,S_vec_t,state_t,gsum_t,temp,state0,
     &dummy1,dummy2,dummy3,dummy4,dummy5,Cauchy_per_vec,dummy6,
     &dummy7,state,dummy8,sconv)
c
          if (sconv.eq.0) jconv=0

c
c		Assignment of jacobian components
		jacob(1:6,i)=(Cauchy_per_vec-Cauchy_vec)/deps
      enddo
      
      
       
c     Make it symmetric      
      jacob=(transpose(jacob)+jacob)/2.0
      
      
c      do i=1,6
c      write(6,*) 'jacob', (jacob(i,j), j=1,6)
c      enddo

      
      
	return
      end subroutine SC_jacobian
c            
      
      
   
      
c
c	main subroutine that calculates the stress and correct Fp through
c	Dawson-Maniatty intergration method
c	This contains all the sub-routines called during the calculation
c	INPUTS: F_T(3,3), Fe_t0(3,3), Fp_t0(3,3), tauc_t0(12), dt
c	OUPUTS: invFp_T(3,3), T_T_vec(6), gammadot(12), dgammadot_dtau(12),
c			 Lp(3,3), tauc(12), initno, ouitno
c	USES:	scale, innertol, outertol, innoitmax, ounoitmax
	subroutine SC_main(dt,F,Fp_t,S_vec_t,state_t,gsum_t,temp,state0,
     &C,S_vec,Lp,Fp,Fe,Cauchy_vec,gammadot,dgammadot_dtau,state,gsum,
     &sconv)
c      
c      
	use globalvars, only : ounoitmax, innoitmax, innertol, outertol,
     &numslip, Schmid, SchmidT, I3, elas66, I6, dS_cr, numstvar
	use globalsubs, only : invert3x3, determinant,convert3x3to6,
     &convert6to3x3, inverse, invertnxn
	implicit none
c	Input variable declarations
      real(8) dt,F(3,3),Fp_t(3,3),S_vec_t(6),gsum_t,temp
c	Output variable declarations
      real(8) C(3,3,numslip),S_vec(6),Lp(3,3)
      real(8) Fp(3,3),Fe(3,3),Cauchy(3,3)
      real(8) Cauchy_vec(6), gammadot(numslip),dgammadot_dtau(numslip)
      real(8) gsum
      integer ounoit,innoit,i,j,is,sconv,count
c	Variables used within this subroutine
      real(8) detFp_t,invFp_t(3,3),A(3,3)
	  real(8) B(3,3,numslip),B_vec(6,numslip)
      real(8) E_tr(3,3),E_vec_tr(6),S_vec_tr(6),C_vec(6,numslip)
      real(8) sumG(6), tau(numslip),inres
      real(8) dG(6,6),dS_vec(6),invdG(6,6),detdG,G_vec(6)
	  real(8) outres,invFp(3,3),detFp, state0(numslip,numstvar)
      real(8) Ee(3,3), Ee_vec(6), Rstate_vec(numslip*numstvar)
      real(8)	state(numslip,numstvar), state_t(numslip,numstvar)
      real(8) statedot(numslip,numstvar), Rstate(numslip,numstvar)
      logical notnum
      
      
c	Calculation of known quantities
      call invert3x3(Fp_t,invFp_t,detFp_t)

      A=matmul(transpose(invFp_t),
     + matmul(matmul(transpose(F),F),invFp_t))
      B=0.
      B_vec=0.
      do is=1,numslip
          B(:,:,is)=matmul(A,Schmid(is,:,:)) + matmul(SchmidT(is,:,:),A)
          call convert3x3to6(B(:,:,is),B_vec(:,is))
      enddo
      
c     Trial strain
      E_tr=0.5*(A-I3)
          
      
c     Vectorization of trial strain
      call convert3x3to6(E_tr,E_vec_tr)
      
      
c     Trial stress
      S_vec_tr = matmul(elas66,E_vec_tr)
      
      
c     Calculation of constant C
	C=0.
      C_vec=0.
      
      do is=1,numslip
          
          C_vec(:,is) = 0.5*matmul(elas66,B_vec(:,is))
          call convert6to3x3(C_vec(:,is),C(:,:,is))
          
      enddo
      
      
      
c	Assign initial variables	
	state=state_t
      S_vec=S_vec_t
      

      
c	Outer loop starts HERE!
	do ounoit=1,ounoitmax
cc		Assign a high value for the residual (Necessary for scale-back algorithm)
c		inres=1.0d+10
c		Inner loop starts HERE!
		do innoit=1,innoitmax


c			Constitutive calculations
			call constitutive(S_vec,C_vec,state,temp,
     &                tau,gammadot,dgammadot_dtau,Lp,sumG)
              
c              
c	        Residual
              G_vec = S_vec - S_vec_tr + sumG*dt
      

c              write(6,*) 'state',state 
               
              
              
c              write(6,*) 'tau', tau
c			write(3,*) 'R_vec'
c			write(3,*) R_vec	
c			Assign old value of residual 
c			inres0=inres
			inres=maxval(dabs(G_vec))

c              write(6,*) 'innoit',innoit
c              write(6,*) 'gammadot',gammadot
c              write(6,*) 'dgammadot_dtau',dgammadot_dtau
c              write(6,*) 'inres',inres
              
c			If the residual is smaller than the tolerance exit the loop
			if (inres.lt.innertol) then
c				write(3,*) innoit
c				write(3,*) 'converged'
				exit
c			
			else
c				write(3,*) innoit
				call tangent(C_vec,dgammadot_dtau,dt,dG)
c

                  
                  
                  
				
c                 Inverse of the tangent for NR iteration
c                  call  invertnxn(dG,6,invdG,detdG)
                  call inverse(dG,invdG,6)
c                  write(6,*) 'detdG',detdG
 
                  
                  
c				Stress increment
                  dS_vec = matmul(-invdG,G_vec)
                  
                  
                  
c			write(3,*) 'dS_vec'
c			write(3,*) dS_vec                  
                  
                  
c                 Stress component checks
                  do i=1,6
                      
                      if (dabs(dS_vec(i)).gt.dS_cr) then
                          
                          dS_vec(i) = dsign(1.0d0,dS_vec(i))*dS_cr
                          
                      endif
                  
                      
                  enddo
                  
c			write(3,*) 'dS_vec'
c			write(3,*) dS_vec                


                  
                  
                  S_vec = S_vec + dS_vec
                  
c
c				write(3,*) 'dFp_T'
c				write(3,*) dFp_T
			endif
          enddo
          
          
          gsum = gsum_t + sum(dabs(gammadot))*dt
                   

  
c		write(3,*) 'innoit:  ',innoit
c		END of INNER iteration loop	
c		Calculate amount of slip system hardening	
		call hardening(gammadot, state, gsum, temp, state0, statedot)
c          write(6,*) 'statedot',statedot
c		Residual increments of slip resistance
          count=0
          Rstate_vec = 0.0d+0
          do j=1,numstvar
              do i=1,numslip
                  Rstate(i,j)=state(i,j)-state_t(i,j)-statedot(i,j)*dt
c                 Vectorize state variables
                  count=count+1
                  Rstate_vec(count) = Rstate(i,j)
              enddo
          enddo
          
		
c		Relative maximum change in slip system resistivity
          outres=maxval(dabs(Rstate_vec))
c          write(6,*) 'Rstate_vec'
c          write(6,*) Rstate_vec
c          write(6,*) 'ounoit',ounoit
c          write(6,*) 'outres',outres
c          write(6,*) 'gammadot',gammadot
c          write(6,*) 'tauc',tauc
c          write(6,*) 'dtauc',tauc
          
c		Check if the change of hardening is within the tolerances
		if (outres.lt.outertol) exit
c		If not accept the result and continue the implicit iteration
          do j=1,numstvar
              do i=1,numslip
                  state(i,j)=state_t(i,j)+statedot(i,j)*dt
              enddo
          enddo
c		state=state_t+statedot*dt
          

          !write(6,*) 'state_t'
          !write(6,*) state_t
          !
          !write(6,*) 'state'
          !write(6,*) state
          
      enddo
c	END of OUTER iteration loop	
c	write(3,*) 'ounoit:  ',ounoit


      
      
cc     Constitutive calculations
c	call constitutive(C,invFp,tauc,el_no,ip_no,
c     &				Ce,PK2_vec,tau,gammadot,dgammadot_dtau,Lp)

     
c     Calculate plastic deformation gradient
      Fp = matmul((I3 + Lp*dt),Fp_t)
      
c     Find the determinant and scale it
      call determinant(Fp,detFp)
      
c     Plastic part of the deformation gradient      
      Fp = Fp / detFp**(1./3.)
      
c     Invert Fp
	call invert3x3(Fp,invFp,detFp)
      
c	Calculate the elastic part of the deformation
	Fe=matmul(F,invFp)
     
      
      
c     Calculate elastic strains
      Ee = (matmul(transpose(Fe),Fe)-I3)/2.0
      
c     Vectorize strains
      call convert3x3to6(Ee, Ee_vec)
      
c     Calculate the stresses
      S_vec = matmul(elas66,Ee_vec)
      
      
     
c	Calculate Cauchy stress
	call cauchystress(S_vec,Fe,Cauchy,Cauchy_vec)
      
      
c     Set flag for convergence
      sconv=1d+0
      
c     Check if the stress value converged         
      do i=1,6
          notnum = isnan(Cauchy_vec(i))
          if (notnum) sconv=0d+0
      enddo
      
c     Check if the stress value converged       
      do i=1,6
          
          if (abs(Cauchy_vec(i)).gt.1.0d+50) sconv=0d+0
          
      enddo
      
      
c     Check if the slip rates are infinite
      do i=1,numslip
          
          if (dabs(gammadot(is)).gt.1.0d+50) sconv=0d+0
          
      enddo
      
      
c     Check for the number of iterations      
      if (sconv.eq.1) then
          
          if (innoit.ge.innoitmax) sconv=0d+0

          
          if (ounoit.ge.ounoitmax) sconv=0d+0

          
      endif
      
      


	return
	end subroutine SC_main





c	This subroutine includes the constitutive calculations
c	INPUTS: C(3,3), invFp_T(3,3), tauc(12)
c	OUTPUTS:Ce(3,3), T_T_vec(6), tau(12), gammadot(12), dgammadot_dtau(12), Lp(3,3)
c	USES: Schmid(12,3,3), Schmid_vec(12,6), zeta66(6,6), gammadot0, mm, I3(3,3)
      subroutine constitutive(S_vec,C_vec,state,temp,
     + tau,gammadot,dgammadot_dtau,Lp,sumG)
      
	  use globalvars, only : I3,elas66,Schmid
      use globalvars, only : Schmid_vec,numslip,numstvar
      use globalsubs, only : convert3x3to6
      use slipratelaws, only: sliprate
      implicit none
c	Input variable declarations
      real(8) S_vec(6),C_vec(6,numslip),temp
c	Output variable declarations
      real(8) tau(numslip),gammadot(numslip)
      real(8) dgammadot_dtau(numslip),Lp(3,3)
      real(8) sumG(6)
c	Variables used within the code
      real(8) E(3,3),E_vec(6)
      integer xx,i
      real(8)	state(numslip,numstvar)
c	ASSIGNMENT OF GLOBAL VARIABLES


	sumG =0.0
      


c	Calculation of plastic part of the velocity gradient
	Lp=0.0
      
	do xx=1,numslip
c		Calculate resolved shear stress	
		tau(xx)=0.0
		do i=1,6
			tau(xx)=tau(xx)+(Schmid_vec(xx,i)*S_vec(i))
          enddo
		

          
c         Calculate slip rates
          call sliprate(tau(xx),state(xx,1:numstvar),temp,
     &            gammadot(xx),dgammadot_dtau(xx))
          

c         Calculate sumG
          sumG = sumG + gammadot(xx)*C_vec(:,xx)
          
c         Plastic part of the velocity gradient          
		Lp=Lp+(gammadot(xx)*Schmid(xx,1:3,1:3))
	enddo


	return
	end subroutine constitutive



	






c	This subroutine calculates the tangent and the increment in plastic part of the deformation gradient
c	for the Newton-Raphson scheme
c	INPUTS:		R_vec(9), Ce(3,3), invFp_T(3,3), detFp_T, Fp_t0(3,3), dgammadot_dtau(12), dt
c	OUTPUTS:	dFp_T(3,3)
c	USES:		Schmid(12,3,3), zeta3333(3,3,3,3), lambda_p, I3(3,3)
	subroutine tangent(C_vec,dgammadot_dtau,dt,dG)
	use globalvars, only: Schmid_vec,numslip,I6
	implicit none
c	Input variable declarations
	real(8) C_vec(6,numslip), dgammadot_dtau(numslip),dt
c	Output variable declaration
	real(8) sumdG(6,6), dG(6,6)
c	Variables used within this subroutine
	integer i,j,xx


      

c	Calculation tangent
	sumdG=0.0

	do i=1,6
		do j=1,6
              do xx=1,numslip 
                sumdG(i,j) = sumdG(i,j)+C_vec(i,xx)*Schmid_vec(xx,j)*
     &                            dgammadot_dtau(xx)
              enddo
		enddo
	enddo
	
	

c     Tangent of the residual
      dG = sumdG*dt + I6



	return
	end subroutine tangent





c	This subroutine calculates the amount of hardening for a given shear rate
c	INPUTS:		gammadot(12), tauc(12)
c	OUTPUTS:	dtauc(12)
c	USES:		h0, ss, a, hardmat
      subroutine hardening(gammadot, state, gsum, temp, 
     + state0, statedot)
	
      use globalvars, only: intmat, numslip, modelno, numstvar
      use sliphardlaws, only: sliphard
      implicit none
c	Input variable declarations
      real(8) gammadot(numslip), gsum, temp
c	Input variable declarations
c	Variables used within this subroutine
      real(8) taucdot(numslip), tothard
      integer xx, i, j
      real(8)	state(numslip,numstvar),state0(numslip,numstvar)
      real(8)	statedot(numslip,numstvar)
      real(8) statedot_(numslip,numstvar+1)

      
      statedot = 0.

	do xx=1,numslip
		
          
          
      

          call sliphard(state(xx,1:numstvar), gsum, gammadot(xx), temp,
     &            state0(xx,1:numstvar),statedot_(xx,1:numstvar+1))
          
           
          
           
      enddo

      
      statedot = statedot_(1:numslip,1:numstvar)
      

      
      tothard = 0.
c     Slip hardening law of Dylan - relies on the cumulative sum of hardening rates
      if (modelno.eq.2) then

          do xx=1,numslip
              tothard = tothard + statedot_(xx,1)
          enddo
          
          
          do xx=1,numslip
              
              statedot(xx,1) = tothard - statedot_(xx,3)
              
          enddo

         
          
      endif
      

    
      
c     State variable depends on the model so, an if statement is placed!     
c     Slip interaction matrix / latent hardening effects
c     Applies to the 1st state variable only! (i.e. tauc, rho, etc.)
c     An if statement is placed specific to the model due to the state variables
      if (modelno.eq.1) then
      

          taucdot =  statedot(:,1)
      
          taucdot = matmul(intmat,taucdot)
      
          statedot(:,1) = taucdot

           
          
      elseif (modelno.eq.2) then
          
          taucdot =  statedot(:,1)
      
          taucdot = matmul(intmat,taucdot)
      
          statedot(:,1) = taucdot 
          
          
      endif
      
      
      
      
	return
	end subroutine hardening


c	This subroutine calculates the Cauchy stress 
c	INPUTS:	T_T_vec(6), invFp_T(3,3)
c	OUTPUTS: Cauchy(3,3), Cauchy_ve(6)
	subroutine cauchystress(PK2_vec,Fe,Cauchy,Cauchy_vec)
	use globalsubs, only: convert6to3x3,determinant,convert3x3to6
	implicit none
c	Input variable declarations
	real(8) PK2_vec(6), Fe(3,3)
c	Input variable declarations
	real(8) Cauchy(3,3),Cauchy_vec(6)
c	Variables used within this subroutine
	real(8) detFe, PK2(3,3)
c	 2nd PK stress
	call convert6to3x3(PK2_vec,PK2)
c	Determinant of the elastic part of the overall deformation
	call determinant(Fe,detFe)
c	Cauchy stress
	Cauchy=matmul(matmul(Fe,PK2),transpose(Fe))/detFe
c	Vectorize Cauchy stress
	call convert3x3to6(Cauchy,Cauchy_vec)
      
	return
	end subroutine cauchystress






	end module calculations