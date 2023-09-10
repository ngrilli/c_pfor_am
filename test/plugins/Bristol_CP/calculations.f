c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c Nicolò Grilli
c
c Aug. 12th, 2021 - 1st working version
c
	module calculations
      implicit none
	contains
c
c	This subroutine calculates the two main variables: Stress and Consistent tangent
      subroutine calcs(F_t,F,t,dt,temp,inc,el_no,ip_no,
     &sigma,jacob,pnewdt,coords) ! sigma_damaged must be used for sigma
	use globalvars, only: global_Fp,global_Fp_t,global_Fe,global_Fe_t,
     &global_state,global_state_t,I3,inc_old,elas66,global_coords,
     &global_ori,njaco,innoitmax,ounoitmax,global_gammadot_t,smallnum,
     &global_sigma,global_jacob_t,global_jacob,TF,phaseID,numslip,
     &global_gammadot,global_S,global_S_t,t_old,ratio_lb,ratio_ub,
     &dgamma_s,global_gamma_t,global_gamma,global_sigma_t,GNDeffect,
     &global_gamma_sum_t,global_gamma_sum,numstvar,thermo,temp0,
     &GSeffect,grainID,grainsize_init,global_state0,tstep_forw,GND_init,
     &tstep_back,numel,numip,global_Fr0,tres,resdef,mtdjaco,coords_init,
     &grainmorph,global_damage,global_F_pos,global_F_neg,global_pk2_pos,
     &phasefielddamage,phaseind,maxnumslip,global_Wp,global_Wp_t,
     &global_f_ep_c,global_f_ep_c_t,creepphasefieldflag,
     &global_sigma_damaged,global_sigma_damaged_t,
     &global_S_damaged,global_S_damaged_t
     
      use initialization, only: initialize_grainsize,
     &initialize_gndslipgradel
      use phasefieldfracture, only: update_plastic_work
      use creepphasefielddamage, only: update_creep_damage
      
      implicit none
c	Inputs
      integer inc
      real(8) F_t(3,3),F(3,3),t,dt,temp
c	Outputs
	real(8) sigma(6),jacob(6,6),pnewdt,coords(3)
      real(8) sigma_damaged(6)
c	Variables used within this subroutine
	real(8) Fp_t(3,3),Fe_t(3,3),tauc_t(maxnumslip), Fr(3,3),Fr_t(3,3)
	real(8) Fp(3,3),Fe(3,3),tauc(maxnumslip),Cauchy(3,3)
	real(8) S(6),S_t(6),det,C(3,3,maxnumslip)
      real(8) S_damaged(6),S_damaged_t(6)
	real(8) Lp(3,3),R(3,3),U(3,3),gammadot(maxnumslip)
	real(8) dgammadot_dtau(maxnumslip)
	integer is,el_no,ip_no,gr_no,ph_no
	real(8) g1(3,3),ang,ax(3),dg(3,3),gint(maxnumslip)
      real(8) ratio,dgamma_max,gsum_t,gsum,gint_t(maxnumslip)
      integer sconv, jconv, i, j, flag, nss
c     SC-model inputs      
      real(8)	state(maxnumslip,numstvar),state_t(maxnumslip,numstvar),
     &state0(maxnumslip,numstvar),Xdist(maxnumslip)
c     J2-model inputs
      real(8)	sstate(numstvar),sstate_t(numstvar),sstate0(numstvar),
     &sXdist,sgint_t,sgint,sgsum,sgsum_t,Np(3,3),eta,depsdot_dsigma,
     &epsdot 
c     Variables related with phase field damage
      real(8) F_pos, F_neg, dam, pk2_pos_mat(3,3)
      integer damflag
      real(8) Wp, Wp_t
c     Variables related to Ed's Creep Damage Models      
      real(8) f_ep_c
      
	  ! hard coded flag to activate debug messages when
	  ! the NR loop does not converge
      integer nonconv_debug_messages
	  nonconv_debug_messages = 0

c     - The step update is indicated with the change of time
c     - Change in INC does not work here!!! Because, during sub time stepping,
c       INC also increases while time remains the same
c	if ((t.gt.t_old).or.(inc.ne.inc_old)) then
      if ((t.gt.t_old(el_no,ip_no)).and.(inc.ne.inc_old(el_no,ip_no)))
     &then
c         write(6,*) 'Updated the state variables'
         t_old(el_no,ip_no)=t
         inc_old(el_no,ip_no)=inc
         global_Fp_t(el_no,ip_no,:,:)=global_Fp(el_no,ip_no,:,:)
         global_Fe_t(el_no,ip_no,:,:)=global_Fe(el_no,ip_no,:,:)
         global_S_t(el_no,ip_no,:)=global_S(el_no,ip_no,:)
      global_S_damaged_t(el_no,ip_no,:)=global_S_damaged(el_no,ip_no,:) 
         global_gammadot_t(el_no,ip_no,:)=global_gammadot(el_no,ip_no,:)
         global_gamma_t(el_no,ip_no,:)=global_gamma(el_no,ip_no,:)
         global_gamma_sum_t(el_no,ip_no)=global_gamma_sum(el_no,ip_no)
         global_state_t(el_no,ip_no,:,:)=global_state(el_no,ip_no,:,:)
         global_jacob_t(el_no,ip_no,:,:)=global_jacob(el_no,ip_no,:,:)
         global_sigma_t(el_no,ip_no,:)=global_sigma(el_no,ip_no,:)
      global_sigma_damaged_t(el_no,ip_no,:) =
     & global_sigma_damaged(el_no,ip_no,:)
         global_Wp_t(el_no,ip_no) = global_Wp(el_no,ip_no)
         global_f_ep_c_t(el_no,ip_no) = global_f_ep_c(el_no,ip_no)
      endif
	  

c     Grain ID     
      gr_no = grainID(el_no)
      
c     Phase index (not phase ID)
      ph_no = phaseind(el_no)         

      
c     Number of slip systems of the current phase
      nss = numslip(ph_no)
	 
      
         
      
c     If length scale calculation is "ON" - Dylan's version
      if ((GSeffect.eq.1d+0).or.(GSeffect.eq.2d+0)) then
c         If not initialized
          if (grainsize_init(el_no,ip_no).eq.0d+0) then

              
c             flag for initialization
              grainsize_init(el_no,ip_no)=1d+0

          
c              write(6,*) 'element no.: ', el_no
          
c              write(6,*) 'IP no.: ', ip_no
          

          
c              write(6,*) 'grain no.: ', gr_no
           
          
c              write(6,*) 'coordinates: ', coords
          
              call initialize_grainsize(el_no,ip_no,gr_no,ph_no,coords)
          
          endif
      endif
         
       

c     If not initialized
      if (coords_init(el_no,ip_no).eq.0d+0) then

              
c         flag for initialization
          coords_init(el_no,ip_no) = 1d+0

          
c          write(6,*) 'element no.: ', el_no
          
c          write(6,*) 'IP no.: ', ip_no
          
c          write(6,*) 'coordinates: ', coords
          
          global_coords(el_no,ip_no,1:3) = coords
          
      endif      
      


c     Not yet initialized (done only once)
      if (GND_init.eq.0d+0) then
          
          
c         Once the calculations are complete
          flag=0d+0
          do i=1, numel
              do j=1,numip
                  flag = flag + coords_init(i,j)
              enddo
          enddo
          
c         Once the coordinates are computed          
          if (flag.eq.numel*numip) then
              
              GND_init = 1d+0
                    
c             Initialize GND calculation from slip gradients
              if ((GNDeffect.eq.1d+0).or.(GNDeffect.eq.2d+0)) then

c                 initialize calculations for GND mapping after all the element information is complete!
c                 This is done ONCE!
                  call initialize_gndslipgradel
                  
            
              endif
              
          endif
          

          
          
          
      endif      
      

      


      
c     In case of mechanical solver only, use initial temperatures      
      if (thermo.eq.0d+0) then
          temp = temp0    
      endif
      
      
      

      
c              


      
c     J2-model      
c     Isotropic material calculations          
      if (phaseID(el_no).eq.0d+0) then
          
          
          
          
c          write(6,*) 'J2 CALCULATIONS'
          
          
          Fp_t=global_Fp_t(el_no,ip_no,:,:)
          Fe_t=global_Fe_t(el_no,ip_no,:,:)
          S_t=global_S_t(el_no,ip_no,:)
          sgsum_t=global_gamma_sum_t(el_no,ip_no)
          sgint_t=global_gamma_t(el_no,ip_no,1)
          sstate0=global_state0(el_no,ip_no,1,:)       
          sstate_t=global_state_t(el_no,ip_no,1,:)
c         Average slip distance          
          sXdist=sum(grainmorph(el_no,ip_no,:))/nss
        
                    
          
c	    Calculate stress and flow stress using J2 plasticity
          call J2_main(dt,F,Fp_t,Fe_t,S_t,sstate_t,sgsum_t,sgint_t,
     &    temp,sstate0,sXdist,ph_no,Fp,Fe,S,Np,eta,epsdot,
     &    depsdot_dsigma,sstate,sgsum,sgint,sigma,sconv)
 
              
c         Calculate the time factor
c           pnewdt=1.25
          pnewdt=ratio_ub
     
c         If not converged     
     	    if (sconv.eq.0d+0) then
			if (nonconv_debug_messages.eq.1) then
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
              write(6,*) 'Fp'
		    write(6,*) Fp
 		    write(6,*) 'sigma'
		    write(6,*) sigma
 		    write(6,*) 'epsdot'
		    write(6,*) epsdot          
		    write(6,*) 'sstate_t'
		    write(6,*) sstate_t
		    write(6,*) 'sstate'
		    write(6,*) sstate              
		    write(6,*) 'dt'
		    write(6,*) dt
			end if
cc		    QUIT ABAQUS
c		    call xit



c             Introduce cut-back              
c              pnewdt=0.5
              pnewdt=tstep_back
			  if (nonconv_debug_messages.eq.1) then
              write(6,*) 'pnewdt',pnewdt
			  end if
              
c             Do not calculate the jacobian - assign former value
              jacob=global_jacob_t(el_no,ip_no,:,:)

c             Assign the old value of stress   
              sigma = global_sigma_t(el_no,ip_no,:)

c
c             Assign former values if not converge
              Fp = Fp_t
              Fe = Fe_t          
              S = S_t
              epsdot = global_gammadot_t(el_no,ip_no,1)/TF
              sstate = sstate_t             
              sgsum = sgsum_t
              sgint = sgint_t
              




c         If converged - regular time stepping
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
                  pnewdt = 1.0d+0
              
              endif              


   


          

c             JACOBIAN CALCULATION FOR ISOTROPIC PHASE     
c	        For a number of increments do the jacobian calculation
c	        Note, jacobian is needed at the first calculation
	        if (modulo(inc,njaco).eq.0) then



c                 Calculate the material tangent (using perturbation)
                  if (mtdjaco.eq.1d+0) then
				  
                      call J2_jacobian_per(dt,F,F_t,Fe_t,Fp_t,S_t,
     &                sstate_t,sgsum_t,sgint_t,temp,sstate0,sXdist,
     &                ph_no,sigma,jacob,jconv)


                      
c                 Calculate the material tangent (using analytical tangent)                      
                  elseif (mtdjaco.eq.2d+0) then
                      

                      
                      call J2_jacobian_ana(dt,eta,Np,depsdot_dsigma,
     &                ph_no,jacob,jconv)
                      

                      
c                 Calculate the material tangent (using elasticity)                      
                  elseif (mtdjaco.eq.3d+0) then
                      
                      jconv = 1d+0
                      jacob = global_jacob_t(el_no,ip_no,:,:)   
                      
                      
                  endif


c		        When the jacobian calculation did not converge, assign the old jacobian!
		        if (jconv.eq.0d+0) then
c			        Jacobian 	
                if (nonconv_debug_messages.eq.1) then		
                      write(6,*) 'Jacobian has not converged - J2'
                      write(6,*) 'el_no',el_no
                      write(6,*) 'ip_no',ip_no
                end if
			        jacob=global_jacob_t(el_no,ip_no,:,:)
			        
                      
c                      pnewdt=0.5
                      pnewdt=tstep_back
					  if (nonconv_debug_messages.eq.1) then
                      write(6,*) 'pnewdt',pnewdt
					  end if
                      
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
          global_S(el_no,ip_no,:) = S
c         Store the important variables
c         For J2 model uses the 1st slip sytem only
          do is=1,nss
              global_state(el_no,ip_no,is,1:numstvar)=sstate(1:numstvar)
              global_gammadot(el_no,ip_no,is) = epsdot*TF
              global_gamma(el_no,ip_no,is) = sgint
          enddo
          
          global_gamma_sum(el_no,ip_no) = sgsum 

      
      
      
      
c     Single crystal calculations
      else
	  

c	    Assign the globally stored variables
        Fe_t = global_Fe_t(el_no,ip_no,:,:)
	    Fp_t = global_Fp_t(el_no,ip_no,:,:)
	    S_t = global_S_t(el_no,ip_no,:)
        S_damaged_t = global_S_damaged_t(el_no,ip_no,:)
	    state_t = global_state_t(el_no,ip_no,:,:)
        state0 = global_state0(el_no,ip_no,:,:)
        gint_t = global_gamma_t(el_no,ip_no,:)
        gsum_t = global_gamma_sum_t(el_no,ip_no)
        Xdist = grainmorph(el_no,ip_no,:)
        Wp_t = global_Wp_t(el_no,ip_no)
        f_ep_c = global_f_ep_c_t(el_no,ip_no)

          ! assign local damage phase field variable
          dam = global_damage(el_no,ip_no)
          damflag = phasefielddamage

c         Fr is scaled with time: i.e. tres = 1 seconds
          if (resdef.eq.1) then
              if (t.le.tres) then
                  Fr= (global_Fr0(el_no,ip_no,:,:)-I3)*(t+dt)/tres + I3
                  Fr_t = (global_Fr0(el_no,ip_no,:,:)-I3)*t/tres + I3
              else
                  Fr = global_Fr0(el_no,ip_no,:,:)
                  Fr_t = global_Fr0(el_no,ip_no,:,:)
              endif
          else
                  Fr = I3
                  Fr_t = I3
          endif
		  

c
c	    Calculate stress and shear resistance
c	    Note: el_no and ip_no are needed to get the values of Schmid vectors and
c	    elasticity tensor from the global variables

          call SC_main(dt,F,Fp_t,Fr,S_t,S_damaged_t,state_t,gsum_t,
     & gint_t,temp,
     & state0,Xdist,dam,damflag,ph_no,C,S,Lp,Fp,Fe,sigma,gammadot,
     & dgammadot_dtau,state,gsum,gint,F_pos,F_neg,pk2_pos_mat,sconv,
     & Wp_t)



c         Calculate the time factor
c          pnewdt=1.25
          pnewdt=ratio_ub






	    if (sconv.eq.0d+0) then
		
          if (nonconv_debug_messages .eq. 1) then
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
		    write(6,*) 'Fe_t'
		    write(6,*) Fe_t
		    write(6,*) 'Fe'
		    write(6,*) Fe
		    write(6,*) 'S_t'
		    write(6,*) S_t
		    write(6,*) 'sigma'
		    write(6,*) sigma
		    write(6,*) 'state_t'
		    write(6,*) state_t
		    write(6,*) 'state'
		    write(6,*) state
		    write(6,*) 'gammadot'
		    write(6,*) gammadot
		    write(6,*) 'dgammadot_dtau'
		    write(6,*) dgammadot_dtau
		    write(6,*) 'dt'
		    write(6,*) dt
        endif

c             Introduce cut-back
c              pnewdt=0.5
              pnewdt=tstep_back
			  if (nonconv_debug_messages.eq.1) then
              write(6,*) 'pnewdt',pnewdt
			  end if

c             Do not calculate the jacobian
              jacob=global_jacob_t(el_no,ip_no,:,:)

c             Assign the old value of stress
              sigma = global_sigma_t(el_no,ip_no,:)



c             Assign the former values
              Fe = Fe_t
              Fp = Fp_t
              S = S_t
              gammadot = global_gammadot_t(el_no,ip_no,:)
              state=state_t
              gsum=gsum_t
              gint=gint_t
							Wp = Wp_t




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
                  pnewdt = 1.0d+0

              endif

c              write(6,*) 'pnewdt',pnewdt




c             For a number of increments do the jacobian calculation
c             Note, jacobian is needed at the first calculation
              if (modulo(inc,njaco).eq.0d+0) then
c                 Calculate the material tangent (using perturbation)
                  if (mtdjaco.eq.1d+0) then

      call SC_jacobian_per(dt,F_t,F,S_t,Fp_t,Fr,state_t,
     & gsum_t,gint_t,temp,state0,Xdist,dam,damflag,ph_no,
     & sigma,jacob,jconv)


c                 Calculate the material tangent (using analytical tangent)
                  elseif (mtdjaco.eq.2d+0) then



                      call SC_jacobian_ana(dt,F,Fe,Fr,S,F_t,Fe_t,
     & Fr_t,gammadot,dgammadot_dtau,C,ph_no,jacob,jconv,
     & dam,damflag)



c                 Calculate the material tangent (using elasticity)
                  elseif (mtdjaco.eq.3d+0) then

                      jconv = 1d+0
                      jacob = global_jacob_t(el_no,ip_no,:,:)


                  endif


c                 When the jacobian calculation did not converge, assign the old jacobian!
                  if (jconv.eq.0) then
				    if (nonconv_debug_messages.eq.1) then
                      write(6,*) 'Jacobian has not converged - CP!'
                      write(6,*) 'inc'
                      write(6,*) inc
                      write(6,*) 'el_no',el_no
                      write(6,*) 'ip_no',ip_no
					end if

c                     pnewdt=0.5
                      pnewdt=tstep_back
					  if (nonconv_debug_messages.eq.1) then
                      write(6,*) 'pnewdt',pnewdt
					  end if

c                     Jacobian
                      jacob=global_jacob_t(el_no,ip_no,:,:)

c                     do i=1,6
c                         write(6,*) (jacob(i,j),j=1,6)
c                     enddo



                  endif

              else

c                 Note this also works when inc=1 since it is elasticity matrix
                  jacob=global_jacob_t(el_no,ip_no,:,:)

              endif


c
            if (phasefielddamage.eq.1d+0) then

              call update_plastic_work(Fp,Fp_t,Fe,S_damaged_t,Wp_t,Wp)

            endif

            if (creepphasefieldflag .eq. 1d+0) then
              !sigma is the cauchy stress.
              call update_creep_damage(Fp_t,Fp,sigma,f_ep_c,dt,ph_no)
            endif

          endif ! end of Regular time stepping


c	    Store the important variables
          global_Fp(el_no,ip_no,:,:) = Fp
          global_Fe(el_no,ip_no,:,:) = Fe
          global_S(el_no,ip_no,:) = S
          global_S_damaged(el_no,ip_no,:) = S_damaged
          global_state(el_no,ip_no,1:nss,1:numstvar) = state
          global_gammadot(el_no,ip_no,1:nss) = gammadot
          global_gamma(el_no,ip_no,1:nss) = gint
          global_gamma_sum(el_no,ip_no) = gsum

	  ! assign to global variables
          global_F_pos(el_no,ip_no) = F_pos
          global_F_neg(el_no,ip_no) = F_neg
	    global_pk2_pos(el_no,ip_no,1:3,1:3) =
     & pk2_pos_mat(1:3,1:3)

        ! Store updated plastic work to global variables
        global_Wp(el_no,ip_no) = Wp
        
        ! Store creep damage to global variables
        global_f_ep_c(el_no,ip_no) = f_ep_c

      endif





c     Store the important results
      global_sigma(el_no,ip_no,:)=sigma
      global_sigma_damaged(el_no,ip_no,:)=sigma_damaged

c     Assign the value of jacobian even if there is no convergence
      global_jacob(el_no,ip_no,:,:)=jacob

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
      subroutine J2_jacobian_per(dt,F,F_t,Fe_t,Fp_t,S_vec_t,sstate_t,
     & sgsum_t,sgint_t,temp,sstate0,sXdist,ph_no,Cauchy_vec,
     & jacob,jconv)

	use globalvars, only : deps,numstvar
	use globalsubs, only : convert6to3x3
	implicit none
c	Inputs
      real(8) F_t(3,3),F(3,3),dt,sgsum_t,sgint_t,temp,sXdist
	real(8) Fp_t(3,3),Fe_t(3,3),Cauchy_vec(6),sstate0(numstvar)
      real(8) sstate_t(numstvar), S_vec_t(6)
      integer ph_no
c	Outputs
	real(8) jacob(6,6)
	integer jconv
c	Variables used within this subroutine
	real(8) Fp(3,3),Fe(3,3),S(6),Np(3,3),epsdot,depsdot_dsigma
	real(8) Fper(3,3),dFrel_vec(6),dFrel(3,3),Cauchy_per_vec(6),sgsum
      real(8) sstate(numstvar),sgint,eta
	integer i,sconv
c	Assign the convergent behavior
	jconv=1d+0
      
cc	Cauchy stress
c	call convert6to3x3(Cauchy_t_vec,Cauchy_t)      
c
c	Increment 6 components of relative deformation gradient
	do i=1,6
c		Component-wise pertubation
		dFrel_vec=0.0d+0
          
         
          dFrel_vec(i)=deps
cc         This is TRUE as in Kalididi's study - gives "G" as the shear term
c		if (i.le.3) then
c              dFrel_vec(i)=deps
c          else   
cc		Note it is not deps/2 since during conversion only one component is considered
c		    dFrel_vec(i)=deps/2.0d+0
c          endif
          
c		Convert the vector to a matrix
		call convert6to3x3(dFrel_vec,dFrel)
		Fper=F+matmul(dFrel,F_t)
c		Call the calculation procedure  
          call J2_main(dt,Fper,Fp_t,Fe_t,S_vec_t,sstate_t,sgsum_t,
     &sgint_t,temp,sstate0,sXdist,ph_no,Fp,Fe,S,Np,eta,epsdot,
     &depsdot_dsigma,sstate,sgsum,sgint,Cauchy_per_vec,sconv)

          
           if (sconv.eq.0) jconv=0d+0

		jacob(1:6,i)=(Cauchy_per_vec-Cauchy_vec)/deps
      enddo


c     Make it symmetric      
      jacob=(transpose(jacob)+jacob)/2.0d+0

	return
	end subroutine J2_jacobian_per
c
c
      
      
      
c      
c      
c     Main routine for jacobian calculation for Martensite 
c	This subroutine calculates consistent tangent
      subroutine J2_jacobian_ana(dt,eta,Np,depsdot_dsigma,
     & ph_no,jacob,jconv)
c
	use globalvars, only : I3, G, kappa
	use globalsubs, only : convert3x3x3x3to6x6
	implicit none
c	Inputs
      real(8) dt, eta, Np(3,3),depsdot_dsigma
      integer ph_no
c	Outputs
	real(8) jacob(6,6)
	integer jconv
c	Variables used within this subroutine
	real(8) jacob4(3,3,3,3),dsigma_dsigmatr
	integer i,j,k,l

c
c	Assign the convergent behavior
	jconv=1d+0
      
      
cc     Dai thesis      
cc     dsigma_dsigmatr
c      dsigma_dsigmatr = 1.0d+0/(1.0d+0 + 3.0d+0*G*dt*depsdot_dsigma)
      
      
      jacob4=0.0d+0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      
                 
c     F. Dunne book
                      jacob4(i,j,k,l)=2.0*eta*G(ph_no)*I3(i,k)*I3(j,l)+
     &3.0*G(ph_no)*(1.0/(1.0 + 3.0*G(ph_no)*dt*depsdot_dsigma)-eta)*
     &Np(i,j)*Np(k,l) + (kappa(ph_no)-2.0*G(ph_no)/3.0)*I3(i,j)*I3(k,l)
      
                      
cc     Dai thesis                               
c                      jacob4(i,j,k,l)= 2.0d+0*eta*G(ph_no)*I3(i,k)*I3(j,l) + 
c     &                (kappa-2.0d+0/3.0d+0*eta*G)*I3(i,j)*I3(k,l) -
c     &                3.0d+0*G*eta*Np(i,j)*Np(k,l) +
c     &                3.0d+0*dsigma_dsigmatr*Np(i,j)*Np(k,l)
          
                      
                  enddo
              enddo
          enddo
      enddo
      
      
      
      call convert3x3x3x3to6x6(jacob4,jacob)
      
      

c     Make it symmetric      
      jacob=(transpose(jacob)+jacob)/2.0d+0

	return
	end subroutine J2_jacobian_ana
c
c      
      
      
      
c      
c
c     Main routine for J2-plasticity
      subroutine J2_main(dt,F,Fp_t,Fe_t,S_t,sstate_t,sgsum_t,sgint_t,
     &temp,sstate0,sXdist,ph_no,Fp,Fe,S_vec,Np,eta,epsdot,
     &depsdot_dsigma,sstate,sgsum,sgint,Cauchy_vec,sconv)
      use globalsubs, only: invert3x3, convert3x3to6, trace, normmat, 
     &convert6to3x3, determinant, polar
      use globalvars, only: elas66_iso, G, E, nu, TF, modelno, largenum,
     &innertol, outertol, innoitmax, ounoitmax, I3, numstvar, dS_cr,
     &outerabstol 
      use slipratelaws, only: sliprate
      use sliphardlaws, only: sliphard
	implicit none
c	Inputs
      real(8) F(3,3),dt, temp, S_t(6)
	real(8) Fp_t(3,3),Fe_t(3,3),sgsum_t,sgint_t,sXdist
      real(8) sstate_t(numstvar),sstate0(numstvar)
      integer ph_no
c	Outputs
	real(8) Fp(3,3),Fe(3,3),epsdot,sigc,Cauchy_vec(6),sgsum,sgint
      real(8) Np(3,3), eta, sstate(numstvar),depsdot_dsigma,S_vec(6)
      integer sconv
c	Vairables used    

      real(8) invFp_t(3,3),detFp_t,Fetr(3,3)
	real(8) Eetr(3,3),Eetr_vec(6),Ee(3,3),Ee_vec(6)
	real(8) Str_vec(6),Str(3,3),Strdev(3,3)
      real(8) hyd,sigmatr,Cauchy(3,3)
	real(8) sigma,dRi,S(3,3),Sdev(3,3)
	real(8) gdot,dgdot_dtau
	real(8) Ri,nRi,Dp(3,3),invFp(3,3),detFp
      real(8) tau,dsigma,Re(3,3),Ue(3,3)
      real(8) Ro(numstvar),sstatedot(numstvar+1)
      integer iloop,oloop,i
      logical notnum, conv
c	real(8) xx1_old,sigma_old
 
c      
	call invert3x3(Fp_t,invFp_t,detFp_t)
cc	Trial elastic deformation	
	Fetr=matmul(F,invFp_t)     


      
      
c	Trial elastic stretch
	Eetr=(matmul(transpose(Fetr),Fetr)-I3)/2.0d+0
      
      
      
     
c	Vectorize strain
	call convert3x3to6(Eetr,Eetr_vec)     
     
c	Vectorized trial stress
	Str_vec=matmul(elas66_iso(ph_no,:,:),Eetr_vec)      
      
c	Trial stress
	call convert6to3x3(Str_vec,Str)
c	Hydrostatic component of the trial stress
	call trace(Str,3,hyd)
	hyd=hyd/3.0d+0
c	Deviatoric component of the trial stress
	Strdev=Str-(I3*hyd)
c	Equivalent trial stress
	call normmat(Strdev,3,sigmatr)
	sigmatr=dsqrt(3.0d+0/2.0d+0)*sigmatr
      
      
      
      

      
      
c	Plastic flow direction (if statement is to correct sigmabar=0 situation which give rise to 1/0 problem)
	if (sigmatr.ne.0.0d+0) then
		Np=dsqrt(3.0d+0/2.0d+0)*Strdev/sigmatr
	else
		Np=0.0d+0
	endif

       
      
   
c     Use the stress at the former time step as the initial guess
      S_vec=S_t
c	Trial stress
	call convert6to3x3(S_vec,S)
c	Hydrostatic component of the trial stress
	call trace(S,3,hyd)
	hyd=hyd/3.0d+0
c	Deviatoric component of the trial stress
	Sdev=S-(I3*hyd)
c	Equivalent trial stress
	call normmat(Sdev,3,sigma)
	sigma=dsqrt(3.0d+0/2.0d+0)*sigma
      

      

      
c	Assign the initial value for the parameters before starting the loop
c	Initial guess for the shear resistance
	sstate=sstate_t      
      
c	OUTER loop starts here
	do oloop=1,ounoitmax
c		INNER loop starts here
cc		Assign a very large value
c		xx1_old=1.0d+10
		do iloop=1,innoitmax

c              write(6,*) 'iloop',iloop
c              write(6,*) 'sigma',sigma
              
              tau = sigma / TF
              
c			Plastic strain rate
              call sliprate(tau,sstate,temp,sXdist,ph_no,
     &        gdot,dgdot_dtau)
              
c             Convert slip rates to strain rates uisng Taylor Factor              
              epsdot = gdot / TF
              depsdot_dsigma = dgdot_dtau / TF**2.0d+0
c              write(6,*) 'sstate',sstate



c              write(6,*) 'gdot',gdot
c              write(6,*) 'tau',tau

           
      
c      write(6,*) 'sigma_before', sigma              

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
			Ri=sigma-sigmatr+(dt*3.0d+0*G(ph_no)*epsdot)
c              write(6,*) 'Ri',Ri
              nRi=dabs(Ri)
cc			Relative norm of the residual
c			if (taubar.eq.0) then 
c				xx1=0.0
c			else
c				
c			endif
c			Check the residual
c			if residual is smaller than the tolerance EXIT
			if (nRi.lt.innertol) then
				exit
c			if the residual have a converging behavior
              endif
c             Calculate tangent for N-R iteration
              dRi = 1.0d+0 + (dt * 3.0d+0 * G(ph_no) * depsdot_dsigma)
              
c             Stress increment
              dsigma = -Ri/dRi
              
c             If the stress increment is larger than the critical value              
              if (dabs(dsigma).gt.dS_cr(ph_no)) then
                  dsigma = dsign(1.0d+0,dsigma)*dS_cr(ph_no)
              endif
              
cc            Store the old value of stress to update the guess
c             sigma_old=sigma
c             Stress after the iteration
              sigma=sigma + dsigma
cc            Assign the old value of the norm for the next iteration
c             xx1_old=xx1
c              write(6,*) 'sigma_after', sigma




              
		enddo
c
c		End of INNER loop
          sgsum = sgsum_t + gdot * dt
          
          sgint = sgint_t + gdot * dt
          
          
c         Strain hardening
          call sliphard(sstate,sgsum,sgint,gdot,temp,
     &sstate0,sXdist,ph_no,sstatedot)
          
c         Slip hardening law of Dylan 
          if (modelno.eq.2d+0) then
              sstatedot(1) = sstatedot(1) - sstatedot(3)
          endif          
          
          
c          write(6,*) 'statedot',statedot
c         Assign convergence flag
          conv = .true.
c		Residual increments of slip resistance
          do i=1,numstvar
            
              Ro(i)=sstate(i)-sstate_t(i)-sstatedot(i)*dt
              
c             Check if it is within the tolerance
              if (dabs(Ro(i)).gt.outerabstol(ph_no,i)) then
                  conv = .false.
              endif
              
          enddo          
          
c          write(6,*) 'oloop',oloop
c          write(6,*) 'gdot',gdot
c          write(6,*) 'sstate',sstate
c          write(6,*) 'sstatedot',sstatedot


          
c		Check the tolerance
		if (conv) exit
          
c		Update the shear resistance
          do i=1,numstvar
		    sstate(i)=sstate_t(i)+sstatedot(i)*dt
          enddo
          
          
	enddo
c
c

cc     Convert slip system values to aggregate
c      epsdot = gdot/TF
c      sigma = tau*TF
c      sigc = tauc*TF




c	Plastic stretch tensor
	Dp=dsqrt(3.0d+0/2.0d+0)*epsdot*Np
      
c	Plastic part of the deformation gradient
	Fp=matmul(Dp*dt + I3,Fp_t)
      
c     Find the determinant and scale it
      call determinant(Fp,detFp)
      
c     Plastic part of the deformation gradient      
      Fp = Fp / detFp**(1.0d+0/3.0d+0)
      
c     Invert Fp
	call invert3x3(Fp,invFp,detFp)
      
c	Elastic part of the deformation gradient
	Fe=matmul(F,invFp)
      
c	Elastic stretch
	Ee=(matmul(transpose(Fe),Fe)-I3)/2.0d+0


c	Vectorize strain
	call convert3x3to6(Ee,Ee_vec)

      
c	Vectorized stress - 2PK stress
	S_vec=matmul(elas66_iso(ph_no,:,:),Ee_vec)


      
      
      eta = sigma/sigmatr      
c     This can happen initially      
      if (isnan(eta)) then
          eta=1.0d+0
      endif      
      
c      S = eta*Strdev + I3*hyd
c      
cc	Vectorize 2nd PK stress
c	call convert3x3to6(S,S_vec)



cc     Polar decomposition
c      call polar(Fe,Re,Ue)


c      
c     Calculate Cauchy stress for large deformations
      call cauchystress(S_vec,Fe,Cauchy,Cauchy_vec)
      
      

      

      
      

      
c     Flag for convergence      
      sconv=1d+0     
c     Check if the stress value converged         
      do i=1,6
          notnum = isnan(Cauchy_vec(i))
          if (notnum) sconv=0d+0
      enddo
      
c     Check if the stress value converged       
      do i=1,6
          
          if (abs(Cauchy_vec(i)).gt.largenum) sconv=0d+0
          
      enddo
      
      
c     Check if the slip rates are infinite
      if (dabs(epsdot).gt.largenum) sconv=0d+0
      
      
      
c     Check for the number of iterations      
      if (sconv.eq.1d+0) then
          if (iloop.eq.innoitmax) sconv=0d+0
     
          
          if (oloop.eq.ounoitmax) sconv=0d+0
          
      endif
      
            
      
      
      return
      end subroutine J2_main

c
c

c
c
c
c
c
c
c
c
c
c	This subroutine calculates consistent tangent
      subroutine SC_jacobian_per(dt,F_t,F,S_vec_t,
     & Fp_t,Fr,state_t,gsum_t,gint_t,temp,state0,
     & Xdist,dam,damflag,ph_no,Cauchy_vec,jacob,jconv)

	use globalvars, only : deps,innoitmax,ounoitmax,numslip,numstvar,
     &maxnumslip 
	use globalsubs, only : convert6to3x3, determinant
	implicit none
c	Inputs
      real(8) F_t(3,3),F(3,3),S_vec_t(6),Cauchy_vec(6),dt,Fr(3,3)
	real(8) Fp_t(3,3),Cauchy(3,3),gsum_t,gint_t(maxnumslip),temp
      real(8) S_damaged_vec_t(6)
      real(8) Xdist(maxnumslip), dam
      integer damflag, ph_no
c	Outputs
	real(8) jacob(6,6)
	integer jconv
c	Variables used within this subroutine
      real(8) detF,Cauchy_per(3,3)
      real(8) F_per(3,3),dFrel_vec(6),dFrel(3,3),Cauchy_per_vec(6)
      real(8) S_per_vec(6),Fe_per(3,3)
      real(8) invFp_per(3,3),dummy1(3,3,maxnumslip),dummy2(6)
	  real(8) dummy3(3,3)
      real(8) dummy4(3,3),dummy5(3,3),dummy9(maxnumslip)
      real(8) dummy6(maxnumslip),dummy7(maxnumslip),dummy8
      real(8) dummy10, dummy11, dummy12(3,3), dummy13
      integer i,j,sconv,nss
      real(8)	state(maxnumslip,numstvar), state_t(maxnumslip,numstvar)
      real(8) state0(maxnumslip,numstvar)

cc     Determinant of DG
c      call determinant(F,detF)
      
c     Number of slip systems of the current phase
      nss = numslip(ph_no)

      
c	Assign the convergent behavior
	jconv=1d+0
c
c	Increment 6 components of relative deformation gradient
	do i=1,6
c		Component-wise pertubation
		dFrel_vec=0.0d+0

		if (i.le.3) then
              dFrel_vec(i)=deps
          else   
c		Note it is not deps/2 since during conversion only one component is considered
		    dFrel_vec(i)=deps/2.0d+0
          endif

c		Convert the vector to a matrix
		call convert6to3x3(dFrel_vec,dFrel)
		F_per=F+matmul(dFrel,F_t)
c		Call the calculation procedure

		call SC_main(dt,F_per,Fp_t,Fr,S_vec_t,S_damaged_vec_t,
     &    state_t,gsum_t,gint_t,
     &    temp,state0,Xdist,dam,damflag,ph_no,dummy1,dummy2,dummy3,
     &    dummy4,dummy5,Cauchy_per_vec,dummy6,dummy7,state,dummy8,
     &    dummy9,dummy10,dummy11,dummy12,sconv,dummy13)

c
          if (sconv.eq.0d+0) jconv=0d+0

c
c		Assignment of jacobian components
		jacob(1:6,i)=(Cauchy_per_vec-Cauchy_vec)/deps
      enddo



c     Make it symmetric
      jacob=(transpose(jacob)+jacob)/2.0d+0


c      do i=1,6
c      write(6,*) 'jacob', (jacob(i,j), j=1,6)
c      enddo



	return
      end subroutine SC_jacobian_per
c
c
c

c	This subroutine calculates consistent tangent
      subroutine SC_jacobian_ana(dt,F,Fe,Fr,T_vec,F_t,Fe_t,Fr_t,
     & gammadot,dgammadot_dtau,C,ph_no,jacob,jconv,dam,damflag)
      use globalvars, only: numslip, Schmid,
     & elas3333, I3333, I6, largenum, maxnumslip
      use globalsubs, only: invert3x3, convert6to3x3, polar,
     & convert3x3x3x3to9x9, invertnxn, convert9x9to3x3x3x3,
     & convert3x3x3x3to6x6
      use globalsubs, only: determinant
	  
      implicit none
c	Inputs
      real(8) dt,F(3,3),Fe(3,3),Fr(3,3),T_vec(6),F_t(3,3),Fe_t(3,3)
      real(8) Fr_t(3,3),gammadot(maxnumslip)
	  real(8) dgammadot_dtau(maxnumslip)
      real(8) C(3,3,maxnumslip)
      integer ph_no
c	Outputs
      real(8) jacob(6,6)
      integer jconv
c	Variables used within this subroutine
      real(8) detF_t,invF_t(3,3),dF(3,3),T(3,3),dR(3,3),dU(3,3)
      real(8) Fer(3,3),Fer_t(3,3), detFer, invFer(3,3), L4(3,3,3,3)
      real(8) sum, D4(3,3,3,3), G4(maxnumslip,3,3,3,3)
      real(8) B2(maxnumslip,3,3), K4(3,3,3,3), K99(9,9), invK99(9,9)
      real(8) invK4(3,3,3,3), Q4(3,3,3,3), sum_, R2(maxnumslip,3,3)
      real(8) S4(3,3,3,3), W4(3,3,3,3), J4(maxnumslip,3,3,3,3)
      integer i,j,k,l,m,n,q,p,is,nss
      logical notnum
      integer damflag
      real(8) dam
      real(8) Je
	  
c     Number of slip systems of the current phase
      nss = numslip(ph_no)


c	Assign the convergent behavior
	jconv=1d+0



c     Inverse of deformation gradient at former time step
      call invert3x3(F_t,invF_t,detF_t)



c     Relative deformation gradient
      dF=matmul(F,invF_t)

c     2nd PK stress tensor
      call convert6to3x3(T_vec,T)



c     Polar decomposition of relative deformation
      call polar(dF,dR,dU)

      !if (sing.eq.0) then




c     Fe with Fr (considering residual deformations)
      Fer = matmul(Fe,Fr)

c     Fe with Fr (considering residual deformations)
      Fer_t = matmul(Fe_t,Fr_t)



      call invert3x3(Fer,invFer,detFer)



c     Step-1. Calculation of L
      L4 = 0.0d+0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      sum=0.0d+0
					  do m=1,3
						sum = sum + Fer_t(k,i) * dU(l,m) * Fer_t(m,j)
     & + Fer_t(m,i) * dU(m,k) * Fer_t(l,j)

                      enddo

                      L4(i,j,k,l) = sum
                  enddo
              enddo
          enddo
      enddo

      call determinant(Fe, Je)

c     Step-3. Calculation of D
      D4 = 0.0d+0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      sum=0.0d+0
                      do m=1,3
                          do n=1,3

      if (damflag.eq.1d+0 .and. Je>=1d+0) then ! remove this damage part

      sum = sum + 0.5d+0 * elas3333(ph_no,i,j,m,n) 
     & * L4(m,n,k,l) * (1.0 - dam) * (1.0 - dam)

      else

      sum = sum + 0.5d+0*elas3333(ph_no,i,j,m,n)
     & * L4(m,n,k,l)

      endif

                          enddo
                      enddo

                      D4(i,j,k,l) = sum

                  enddo
              enddo
          enddo
      enddo




c     Step-4.
c     (a) Calculation of G
      G4 = 0.0d+0
      do is=1,nss
          do m=1,3
              do n=1,3
                  do k=1,3
                      do l=1,3

                          sum=0.0d+0
                          do p=1,3

                              sum = sum + L4(m,p,k,l)*
     &Schmid(ph_no,is,p,n) + Schmid(ph_no,is,p,m)*L4(p,n,k,l)

                          enddo

                          G4(is,m,n,k,l) = sum

                      enddo
                  enddo
              enddo
          enddo
      enddo


c     (b) Calculation of J
      J4 = 0.0d+0
      do is=1,nss
          do i=1,3
              do j=1,3
                  do k=1,3
                      do l=1,3


                          sum=0.0d+0
                          do m=1,3
                              do n=1,3

      if (damflag.eq.1d+0 .and. Je>=1d+0) then

      sum = sum + 0.5d+0
     & *elas3333(ph_no,i,j,m,n)*G4(is,m,n,k,l)*(1.0-dam)*(1.0-dam)

      else

      sum = sum + 0.5d+0
     & * elas3333(ph_no,i,j,m,n) * G4(is,m,n,k,l)

      endif

                              enddo
                          enddo

                          J4(is,i,j,k,l) = sum

                      enddo
                  enddo
              enddo
          enddo
      enddo



c     Step-5. Calculation of B
      B2 = 0.0d+0
      do is=1,nss
          do i=1,3
              do j=1,3

                  B2(is,i,j)=dgammadot_dtau(is)*dt*Schmid(ph_no,is,i,j)

              enddo
          enddo
      enddo


c     Step-6. Calculation of Q
c     (a) Calculate K
      K4=0.0d+0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3


                      sum = 0.0d+0
                      do is =1,nss

                          sum = sum + C(i,j,is) * B2(is,k,l)

                      enddo

                      K4(i,j,k,l) = I3333(i,j,k,l) + sum


                  enddo
              enddo
          enddo
      enddo


c     Convert 4th rank tensor "K" to 9x9 matrix
      call convert3x3x3x3to9x9(K4,K99)



c     Invert the matrix
      call invertnxn(K99,invK99,9)



c     Convert to a 4th rank tensor
      call convert9x9to3x3x3x3(invK99,invK4)


c     (b) Calculate Q
      Q4=0.0d+0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3


                      sum = 0.0d+0

                      do m=1,3
                          do n=1,3

                              sum_ = 0.0d+0

                              do is=1,nss

                                  sum_ = sum_ + invK4(i,j,m,n)
     & * J4(is,m,n,k,l) * gammadot(is) * dt

                              enddo

                              sum = sum + invK4(i,j,m,n)
     & * D4(m,n,k,l) - sum_


                          enddo
                      enddo

                      Q4(i,j,k,l) = sum

                  enddo
              enddo
          enddo
      enddo


c     Step-7. Calculation of R and S
c     (a) Calcualtion of R

      R2 = 0.0d+0
      do is=1,nss
          do i=1,3
              do j=1,3
                  sum=0.0d+0

                  do k=1,3
                      do l=1,3

                          sum = sum + B2(is,k,l) * Q4(k,l,i,j)

                      enddo
                  enddo

                  R2(is,i,j) = sum

              enddo
          enddo
      enddo

c     (b) Calculation of S
      S4 = 0.0d+0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3

                      sum = 0.0d+0

                      do p=1,3

                          do is=1,nss

                              sum = sum + dR(i,k) * Fer_t(l,p) 
     &*gammadot(is)*dt*Schmid(ph_no,is,p,j)



                          enddo
                      enddo


                      sum_=0.0d+0
                      do m=1,3
                          do n=1,3
                              do p=1,3
                                  do is=1,nss

                                      sum_= sum_ + dR(i,m) * dU(m,n)
     &*Fer_t(n,p)*R2(is,k,l)*Schmid(ph_no,is,p,j)

                                  enddo
                              enddo
                          enddo
                      enddo


                      S4(i,j,k,l)= dR(i,k) * Fer_t(l,j) - sum - sum_

                  enddo
              enddo
          enddo
      enddo






      W4 = 0.0d+0
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3


                      sum_=0.0d+0
                      do q=1,3
                          do p=1,3

                              sum_= sum_ + S4(p,q,k,l) * invFer(q,p)

                          enddo
                      enddo



                      sum=0.0d+0
                      do m=1,3
                          do n=1,3

                              sum = sum +
     & S4(i,m,k,l) * T(m,n) * Fer(j,n) +
     & Fer(i,m) * Q4(m,n,k,l) * Fer(j,n) +
     & Fer(i,m) * T(m,n) * S4(j,n,k,l) -
     & Fer(i,m) * T(m,n) * Fer(j,n) * sum_

                          enddo
                      enddo


                      W4(i,j,k,l) = sum / detFer



                  enddo
              enddo
          enddo
      enddo


      call convert3x3x3x3to6x6(W4,jacob)



c     Make it symmetric
      jacob=(transpose(jacob)+jacob)/2.0d+0



c     Check if the jacob value converged or not
      do i=1,6
          do j=1,6
              notnum = isnan(jacob(i,j))
              if (notnum) then
                  jconv=0d+0
                  jacob = I6
              endif
          enddo
      enddo

c     Check if the jacob value converged or not
      do i=1,6
          do j=1,6
              if (abs(jacob(i,j)).gt.largenum) then
                  jconv=0d+0
                  jacob=I6
              endif
          enddo
      enddo



c




c      write(6,*) 'jacob'
c      do i=1,6
c      write(6,*) (jacob(i,j), j=1,6)
c      enddo



	return
      end subroutine SC_jacobian_ana
c
c
c
c
c
c	main subroutine that calculates the stress and correct Fp through
c	Kalidindi-Anand intergration method
c	This contains all the sub-routines called during the calculation
c	INPUTS: F_T(3,3), Fe_t0(3,3), Fp_t0(3,3), tauc_t0(12), dt
c	OUPUTS: invFp_T(3,3), T_T_vec(6), gammadot(12), dgammadot_dtau(12),
c			 Lp(3,3), tauc(12), initno, ouitno
c	USES:	scale, innertol, outertol, innoitmax, ounoitmax

      subroutine SC_main(dt,F,Fp_t,Fr,S_vec_t,S_damaged_vec_t,state_t,
     & gsum_t,gint_t,temp,state0,Xdist,dam,damflag,ph_no,C,
     & S_vec,Lp,Fp,Fe,Cauchy_vec,gammadot,dgammadot_dtau,
     & state,gsum,gint,F_pos,F_neg,pk2_pos_mat,sconv,Wp)
c
c
c
      use globalvars, only : ounoitmax, innoitmax, 
     &innertol, outerabstol, I6, dS_cr, numstvar, largenum,
     &Schmid, SchmidT, I3, elas66, numslip, maxnumslip
c
c
c
      use phasefieldfracture, only : computeStrainVolumetric
c
      use globalsubs, only : invert3x3, determinant,convert3x3to6,
     &convert6to3x3, invertnxn
      
      implicit none
c	Input variable declarations
      real(8) dt,F(3,3),Fp_t(3,3),S_vec_t(6),gsum_t,temp
      real(8) S_damaged_vec_t(6)
      real(8) Fr(3,3),Xdist(maxnumslip),dam,gint_t(maxnumslip)
      integer damflag, ph_no
c	Output variable declarations
      real(8) C(3,3,maxnumslip),S_vec(6),Lp(3,3)
      real(8) Fp(3,3),Fe(3,3),Cauchy(3,3),gsum
      real(8) Cauchy_vec(6), gammadot(maxnumslip)
      real(8) gint(maxnumslip),dgammadot_dtau(maxnumslip)
c     Output variables related with phase field damage
      real(8) F_pos, F_neg, pk2_pos_mat(3,3)

c	Variables used within this subroutine
      real(8) detFp_t,invFp_t(3,3),A(3,3)
      real(8) B(3,3,maxnumslip),B_vec(6,maxnumslip)
      real(8) E_tr(3,3),E_vec_tr(6),S_vec_tr(6),C_vec(6,maxnumslip)
      real(8) sumG(6), tau(maxnumslip),inres,invFr(3,3),detFr,Fsum(3,3)
      real(8) dG(6,6),dS_vec(6),invdG(6,6),detdG,G_vec(6)
      real(8) outres,invFp(3,3),detFp,state0(maxnumslip,numstvar)
      real(8) Ee(3,3),Ee_vec(6)
      real(8)	state(maxnumslip,numstvar),state_t(maxnumslip,numstvar)
      real(8) statedot(maxnumslip,numstvar),Rstate(maxnumslip,numstvar)
      logical notnum, conv
      integer ounoit,innoit,i,j,is,sconv,nss
      real(8) Je
			real(8) Wp
      integer nonconv_debug_messages
	  
	  ! hard coded flag to activate debug messages when
	  ! the NR loop does not converge
      nonconv_debug_messages = 0

c     Number of slip systems      
      nss = numslip(ph_no)   

c	Calculation of known quantities
	call invert3x3(Fp_t,invFp_t,detFp_t)

c	Calculation of known quantities
	call invert3x3(Fr,invFr,detFr)

c     First calculate the original one
      A = matmul(transpose(invFp_t),
     & matmul(matmul(transpose(F),F),invFp_t))

c     Modified for residual deformation
      A = matmul(transpose(invFr),matmul(A,invFr))


      B=0.0d+0
      B_vec=0.0d+0
      do is=1,nss
          B(:,:,is)=matmul(A,Schmid(ph_no,is,:,:)) +
     &    matmul(SchmidT(ph_no,is,:,:),A)
          call convert3x3to6(B(:,:,is),B_vec(:,is))
      enddo

c     Trial strain
      E_tr=0.5d+0*(A-I3)


c     Vectorization of trial strain
      call convert3x3to6(E_tr,E_vec_tr)

c      write(6,*) 'F'
c      write(6,*) F


c      write(6,*) 'E_vec_tr'
c      write(6,*) E_vec_tr

c     Trial stress

      !if (damflag.eq.0d+0) then

      S_vec_tr = matmul(elas66(ph_no,:,:),E_vec_tr)
      F_neg=0.0d+0
      F_pos=0.0d+0
      pk2_pos_mat=0.0d+0

      !else ! phase field damage model

	  ! calculate original elastic deformation gradient
	  ! and modify it for residual deformation
      !Fe = matmul(F,invFp_t)
      !Fe = matmul(Fe,invFr)
      !call determinant(Fe, Je)

      !call computeStrainVolumetric(ph_no,E_tr,A,Fe,dam,
      !&S_vec_tr,F_pos,F_neg,pk2_pos_mat,Wp)

      !end if

c     Calculation of constant C
	C=0.0d+0
      C_vec=0.0d+0

      do is=1,nss

	    ! include damage in the NR Jacobian
        if (damflag.eq.1d+0 .and. Je>=1d+0) then

      C_vec(:,is) = 0.5d+0*(1.0-dam)*(1.0-dam)
     & *matmul(elas66(ph_no,:,:),B_vec(:,is))

        else

          C_vec(:,is) = 0.5d+0*matmul(elas66(ph_no,:,:),B_vec(:,is))

        endif

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
			call constitutive(S_vec,C_vec,state,temp,Xdist,
     &        ph_no,tau,gammadot,dgammadot_dtau,Lp,sumG)

c
c	        Residual
              G_vec = S_vec - S_vec_tr + sumG*dt


c              write(6,*) 'state',state
      do i=1,6
        if (isnan(S_vec(i))) then
		if (nonconv_debug_messages.eq.1) then
        write(6,*) 'state', state
        write(6,*) 'S_vec', S_vec
        write(6,*) 'gammadot', gammadot
        write(6,*) 'tau', tau
		end if
        endif
      enddo


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
				call tangent(C_vec,dgammadot_dtau,dt,ph_no,dG)
c





c                 Inverse of the tangent for NR iteration
c                  call  invertnxn(dG,6,invdG,detdG)
                  call invertnxn(dG,invdG,6)
c                  write(6,*) 'detdG',detdG



c				Stress increment
                  dS_vec = matmul(-invdG,G_vec)



c			write(3,*) 'dS_vec'
c			write(3,*) dS_vec


c                 Stress component checks
                  do i=1,6

                      if (dabs(dS_vec(i)).gt.dS_cr(ph_no)) then
                          
                      dS_vec(i) = dsign(1.0d+0,dS_vec(i))*dS_cr(ph_no)
                          
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

c         cumulative slip
          gsum = gsum_t + sum(dabs(gammadot))*dt

c         time integrated value of slip
          gint = gint_t + dabs(gammadot)*dt



c		write(3,*) 'innoit:  ',innoit
c		END of INNER iteration loop
c		Calculate amount of slip system hardening

          call hardening(gammadot,state0,state_t,
     &    state,gsum,gint,temp,Xdist,ph_no,dt,
     &    statedot)

c          write(6,*) 'statedot',statedot
c         Assign convergence flag
          conv = .true.
c		Residual increments of slip resistance
          do j=1,numstvar
              do i=1,nss
                  Rstate(i,j)=state(i,j)-state_t(i,j)-statedot(i,j)*dt
                  
c                 Check if it is within the tolerance
                  if (dabs(Rstate(i,j)).gt.outerabstol(ph_no,j)) then
                      conv = .false.
                  endif
                  
c                 Vectorize state variables
c                  count=count+1
c                  Rstate_vec(count) = Rstate(i,j)
              enddo
          enddo


c		Relative maximum change in slip system resistivity
c          outres=maxval(dabs(Rstate_vec))
c          write(6,*) 'Rstate_vec'
c          write(6,*) Rstate_vec
c          write(6,*) 'ounoit',ounoit
c          write(6,*) 'outres',outres
c          write(6,*) 'gammadot',gammadot
c          write(6,*) 'tauc',tauc
c          write(6,*) 'dtauc',tauc

c		Check if the change of hardening is within the tolerances
		if (conv) exit
c		If not accept the result and continue the implicit iteration
          do j=1,numstvar
              do i=1,nss
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
      Fp = Fp / detFp**(1.0d+0/3.0d+0)

c     Invert Fp
      call invert3x3(Fp,invFp,detFp)

c     Modifed for residual deformation
      Fe=matmul(F,invFp)

c     Calculate inverse of the residual deformation



c	Calculate the elastic part of the deformation
      Fe=matmul(Fe,invFr)

c     A is needed by the damage model

      A = matmul(transpose(Fe),Fe)

c     Calculate elastic strains
      Ee = (A-I3)/2.0d+0

c     Vectorize strains
      call convert3x3to6(Ee, Ee_vec)

c     Calculate the stresses

      if (damflag.eq.0d+0) then

        S_vec = matmul(elas66(ph_no,:,:),Ee_vec)

      else ! phase field damage model

          call computeStrainVolumetric(ph_no,Ee,A,Fe,dam,
     &S_vec,F_pos,F_neg,pk2_pos_mat,Wp)

      end if

c     Modifed for residual deformation
      Fsum = matmul(Fe,Fr)

c     Modifed for residual deformation
c	Calculate Cauchy stress
	call cauchystress(S_vec,Fsum,Cauchy,Cauchy_vec)


c     Set flag for convergence
      sconv=1d+0

c     Check if the stress value converged
      do i=1,6
          notnum = isnan(Cauchy_vec(i))
          if (notnum) then

              sconv=0d+0

              if (nonconv_debug_messages.eq.1) then
              write(6,*) 'Cauchy stress has NaN!'
			  end if

          endif

      enddo

c     Check if the stress value converged
      do i=1,6

          if (dabs(Cauchy_vec(i)).gt.largenum) then

              sconv=0d+0

              if (nonconv_debug_messages.eq.1) then
              write(6,*) 'Cauchy stress overshoots!'
			  end if

          endif


      enddo


c     Check if the slip rates are infinite
      do is=1,nss

          if (dabs(gammadot(is)).gt.largenum) then

              sconv=0d+0

              if (nonconv_debug_messages.eq.1) then
              write(6,*) 'slip rates overhoot!'
			  end if

          endif


      enddo


c     Check for the number of iterations
      if (sconv.eq.1) then

          if (innoit.eq.innoitmax) then
              sconv=0d+0
			  if (nonconv_debug_messages.eq.1) then
              write(6,*) 'inner loop diverges!'
			  end if
          endif


          if (ounoit.eq.ounoitmax) then
              sconv=0d+0
			  if (nonconv_debug_messages.eq.1) then
              write(6,*) 'outer loop diverges!'
			  end if
          endif


      endif




	return
	end subroutine SC_main





c	This subroutine includes the constitutive calculations
c	INPUTS: C(3,3), invFp_T(3,3), tauc(12)
c	OUTPUTS:Ce(3,3), T_T_vec(6), tau(12), gammadot(12), dgammadot_dtau(12), Lp(3,3)
c	USES: Schmid(12,3,3), Schmid_vec(12,6), zeta66(6,6), gammadot0, mm, I3(3,3)
	subroutine constitutive(S_vec,C_vec,state,temp,Xdist,
     &ph_no,tau,gammadot,dgammadot_dtau,Lp,sumG)
      use globalvars, only : I3,elas66,Schmid,
     &Schmid_vec,numslip,numstvar,maxnumslip,
     &modelno,intmat,intmat1,intmat2
	use globalsubs, only : convert3x3to6
      use slipratelaws, only: sliprate
	implicit none
c	Input variable declarations
	real(8) S_vec(6),C_vec(6,maxnumslip),temp,Xdist(maxnumslip)
      integer ph_no
c	Output variable declarations
	real(8) tau(maxnumslip),gammadot(maxnumslip)
	real(8) dgammadot_dtau(maxnumslip),Lp(3,3)
	real(8) sumG(6)
c	Variables used within the code
	real(8) E(3,3),E_vec(6)
      real(8)	state(maxnumslip,numstvar), state_(maxnumslip,numstvar)
      integer is,i,nss
c	ASSIGNMENT OF GLOBAL VARIABLES

c     Number of slip systems
      nss=numslip(ph_no)


	sumG = 0.0d+0
      
c     No changes for models 1-2-3      
      if (modelno.eq.1d+0) then
          
          state_=state
          
      elseif (modelno.eq.2d+0) then
          
          state_=state
          
      elseif (modelno.eq.3d+0) then
      
          state_=state
      
      
c     If model 4 or 5 apply interaction matrix
c     Dislocation density based model
      elseif (modelno.eq.4d+0) then
          
      
          state_(1:nss,1)=matmul(intmat(ph_no,:,:),
     &    state(1:nss,1))
          
          state_(1:nss,2)=matmul(intmat(ph_no,:,:),
     &    dabs(state(1:nss,2)))
          
          state_(1:nss,3)=matmul(intmat(ph_no,:,:),
     &    dabs(state(1:nss,3)))
          
          state_(:,4) = state(:,4)
          
          
          
          
          
      elseif (modelno.eq.5d+0) then
          
          state_(1:nss,1) = matmul(intmat(ph_no,:,:),
     &    state(1:nss,1))
          
          state_(:,2) = state(:,2)
          
        
c Irradiation model for tauc          	
      elseif (modelno.eq.6d+0) then	
          	
      	
          state_(1:nss,1) = matmul(intmat1(ph_no,:,:),
     &    state(1:nss,1))	
          	
          state_(1:nss,2) = matmul(intmat1(ph_no,:,:),
     &    dabs(state(1:nss,2)))	
          	
          state_(1:nss,3) = matmul(intmat1(ph_no,:,:),
     &    dabs(state(1:nss,3)))	
          	
          state_(:,4) = state(:,4)	
          	
          state_(1:nss,5) = matmul(intmat2(ph_no,:,:),
     &    dabs(state(1:nss,5)))
               
          
          
          
      endif
      
      
      
      
      
      

c	Calculation of plastic part of the velocity gradient
	Lp=0.0d+0
      gammadot=0.0d+0
      dgammadot_dtau=0.0d+0
      
	do is=1,nss
          
c		Calculate resolved shear stress	
		tau(is)=0.0d+0
		do i=1,6
			tau(is)=tau(is)+(Schmid_vec(ph_no,is,i)*S_vec(i))
          enddo
		

          
c         Calculate slip rates
          call sliprate(tau(is),state_(is,1:numstvar),temp,Xdist(is),
     &    ph_no,gammadot(is),dgammadot_dtau(is))
          
c Michael will check      	
c          do plane=1,3	
c              gdotloop=gdotloop+gammadot(plane)	
c          enddo
          
c         Calculate sumG
          sumG = sumG + gammadot(is)*C_vec(:,is)
          
c         Plastic part of the velocity gradient          
		Lp=Lp+(gammadot(is)*Schmid(ph_no,is,1:3,1:3))
      enddo

      


	return
	end subroutine constitutive



	






c	This subroutine calculates the tangent and the increment in plastic part of the deformation gradient
c	for the Newton-Raphson scheme
c	INPUTS:		R_vec(9), Ce(3,3), invFp_T(3,3), detFp_T, Fp_t0(3,3), dgammadot_dtau(12), dt
c	OUTPUTS:	dFp_T(3,3)
c	USES:		Schmid(12,3,3), zeta3333(3,3,3,3), lambda_p, I3(3,3)
	subroutine tangent(C_vec,dgammadot_dtau,dt,ph_no,dG)
	use globalvars, only: Schmid_vec,numslip,I6,maxnumslip
	implicit none
c	Input variable declarations
	real(8) C_vec(6,maxnumslip),dgammadot_dtau(maxnumslip),dt
c	Output variable declaration
	real(8) sumdG(6,6),dG(6,6)
c	Variables used within this subroutine
	integer i,j,is,ph_no,nss

c     Number of slip systems for this phasew
      nss = numslip(ph_no)      

c	Calculation tangent
	sumdG=0.0d+0

	do i=1,6
		do j=1,6
              do is=1,nss
                sumdG(i,j) = sumdG(i,j)+C_vec(i,is)*
     &Schmid_vec(ph_no,is,j)*dgammadot_dtau(is)
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
	subroutine hardening(gammadot,state0,state_t,state,gsum,gint,temp,
     &Xdist,ph_no,dt,statedot)
	use globalvars, only: intmat, numslip, modelno, numstvar,
     &sliphard_param,maxnumslip
      use sliphardlaws, only: sliphard
	implicit none
c	Input variable declarations
	real(8) gammadot(maxnumslip), gsum, gint(maxnumslip)
      real(8) temp, Xdist(maxnumslip), dt
      integer ph_no
c	Input variable declarations
c	Variables used within this subroutine
      real(8) taucdot(maxnumslip), tothard, tauc(maxnumslip)
	integer is, i, j, nss
      real(8)	state0(maxnumslip,numstvar), state_t(maxnumslip,numstvar)
      real(8)	state(maxnumslip,numstvar), statedot(maxnumslip,numstvar)
      real(8) statedot_(maxnumslip,numstvar+1)
      real(8) Q_3, tauc0_3, rhoSSDdot(maxnumslip)

      
      statedot = 0.0d+0

      nss = numslip(ph_no)
      
	do is=1,nss

          call sliphard(state(is,1:numstvar),gsum,gint(is),gammadot(is),
     &    temp,state0(is,1:numstvar),Xdist(is),ph_no,
     &    statedot_(is,1:numstvar+1))
          
      enddo

      
      statedot = statedot_(1:nss,1:numstvar)
      

      
      tothard = 0.0d+0
c     Slip hardening law of Dylan - relies on the cumulative sum of hardening rates
      if (modelno.eq.2d+0) then

          do is=1,nss
              
              tothard = tothard + statedot_(is,1)
              
          enddo
          
          
          do is=1,nss
              
              statedot(is,1) = tothard - statedot_(is,3)
              
          enddo

         
          
      endif
      

    
      
c     State variable depends on the model so, an if statement is placed!     
c     Slip interaction matrix / latent hardening effects
c     Applies to the 1st state variable only! (i.e. tauc, rho, etc.)
c     An if statement is placed specific to the model due to the state variables
      if (modelno.eq.1d+0) then
      

          taucdot =  statedot(:,1)
      
          taucdot(1:nss) = matmul(intmat(ph_no,:,:),
     &    taucdot(1:nss))
      
          statedot(:,1) = taucdot
      

           
          
      elseif (modelno.eq.2d+0) then
          
          taucdot =  statedot(:,1)
      
          taucdot(1:nss) = matmul(intmat(ph_no,:,:),
     &    taucdot(1:nss))
      
          statedot(:,1) = taucdot

          
          
      elseif (modelno.eq.3d+0) then
          
          taucdot =  statedot(:,1)
          
c         Hardening parameter of Code Aster - MFRONT          
          Q_3 = sliphard_param(ph_no,3)
          

          taucdot(1:nss) = Q_3*matmul(intmat(ph_no,:,:),
     &taucdot(1:nss))
          
          
c         Calculate the rate to fool the hardening integraion    
          statedot(:,1) = taucdot
          
 
c     Interaction matrices are applied at the slip rate laws
c     No changes are imposed for dislocation density based models; models 4-5
      elseif (modelno.eq.4d+0) then
          
c          rhoSSDdot = statedot(:,1)
                
c          rhoSSDdot = matmul(intmat,rhoSSDdot)
          
c          statedot(:,1) = rhoSSDdot
      
         
c          write(6,*) rhoSSDdot
          

      elseif (modelno.eq.5d+0) then
      

c          taucdot =  statedot(:,1)
      
c          taucdot = matmul(intmat,taucdot)
      
c          statedot(:,1) = taucdot
          
      elseif (modelno.eq.6d+0) then	
          	
c          rhoSSDdot = statedot(:,1)	
                	
c          rhoSSDdot = matmul(intmat,rhoSSDdot)	
          	
c          statedot(:,1) = rhoSSDdot	
      	
         	
c          write(6,*) rhoSSDdot          


          
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
c
c      
c      
c      
c     This subroutine calculates the misorientation with respect to the initial orientation
      subroutine calculate_misorientation(el_no,ip_no,mis)
      use globalvars, only: global_ori, global_Fe, phaseID
      use globalsubs, only: polar, misorientation
      implicit none
      integer el_no, ip_no, typ
      real(8) g1(3,3), g2(3,3), dg(3,3), mis, ax(3), U(3,3), Fe(3,3)
      
      

c     FCC - cubic symmetry operators              
      if (phaseID(el_no).eq.1d+0) then
          typ = 1d+0
c     BCC - cubic symmetry operators                     
      elseif (phaseID(el_no).eq.2d+0) then
          typ = 1d+0
c    BCT- body centered tetragonal symmetry operators
      elseif (phaseID(el_no).eq.3d+0) then
          typ = 2d+0      
c     HCP - hexagonal symmetry operators
      elseif (phaseID(el_no).eq.3d+0) then
          typ = 3d+0
c     Isotropic
      elseif (phaseID(el_no).eq.0d+0) then
          typ = 0d+0
      endif
c              
c       
c     If the phase is non-isotropic
      if (typ.gt.0d+0) then
c         Initial orientation
          g1=global_ori(el_no,ip_no,:,:)
c         Polar decomposition of the elastic part of the deformation gradient
          Fe = global_Fe(el_no,ip_no,:,:)
          call polar(Fe,g2,U)
c
c         Calcualte misorientation using the subroutine
          call misorientation(g1,g2,typ,ax,mis,dg)
          
          
      else  
          mis=0.0d+0
      endif
      
          
c      
      end subroutine calculate_misorientation 
c
c
c
	end module calculations
