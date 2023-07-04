c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c Nicolò Grilli
c
c
c Kaldindi-Anand semi-implicit two-level time integration scheme is implemented by Eralp Demir.
c
c
c
c
      include "globalvars.f"
      include "globalsubs.f"
      include "phasefieldfracture.f"
      include "creepphasefielddamage.f"
      include "lengthscale.f"
      include "gndslipgrad.f"
      include "initialization.f"
      include "slipratelaws.f"
      include "sliphardlaws.f"
      include "calculations.f"
c
c
c
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
c     Subroutine for initialization
      use initialization, only : initialize_all
      use globalvars, only: foldername, GND_init, GNDeffect
      use gndslipgrad, only: calculategnds, calculatebackstress1,
     &calculatebackstress2

      implicit none
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
      integer,                        intent(in) ::
     & LOP,
     & LRESTART,
     & KSTEP,
     & KINC
      real(8), dimension(1),          intent(in) ::
     & DTIME
      real(8), dimension(2),          intent(in) ::
     & TIME

      ! variables required to read the current working directory
      integer*4 :: status_getcwd
	  character(len=255) :: foldername255
c
c
c
c
c
c     At the start of the analysis (only ONCE!)
      if (LOP.eq.0 .or. LOP.eq.4) then
c
c
      status_getcwd = getcwd(foldername255)
      if ( status_getcwd .ne. 0 ) stop
     & 'Error: simulation files must be in the cwd'
	  foldername = trim(foldername255)

          write(6,*) 'initialization has started!'
          write(6,*) '********************************'
          call initialize_all(foldername)
          write(6,*) '********************************'
          write(6,*) 'initialization has ended!'
c
      endif
c
c
c
c
c     At the end of each increment
c     Update and calculate GNDs (nonlocal calculations)
c     This is done at the end of calculations.
c     GNDs that belong to the PREVOUS time step are used.
c     Initially GNDs are assumed to have "0" values.
      if (LOP.eq.2) then
c
c          write(6,*)'LOP=2'
cc         Reset the flag for cutback
c          GND_cutback=0d+0
c         GND calculations - if selected
          if (GND_init.eq.1d+0) then
c
              if (GNDeffect.eq.1d+0) then
c
c                 At the end of each increment
c                 Update and calculate GNDs (nonlocal calculations)
c                 This is done at the end of calculations.
c                 GNDs that belong to the PREVOUS time step are used.
c                 Initially GNDs are assumed to have "0" values.

c
c
c                 Calculate GNDs
                  call calculategnds(DTIME(1))
c
c                 Calculate Backstress from GND gradients
                  call calculatebackstress1
c
c
                  write(6,*) 'end of increment: ', KINC
                  write(6,*) 'GND & backstress calculations completed!'
c
c
c
c
c             GND calculations - if selected
              elseif (GNDeffect.eq.2d+0) then
c
c
c
c
c                 Calculate GNDs
                  call calculategnds(DTIME(1))
c
c                 Calculate Backstress from GND gradients
                  call calculatebackstress2
c
c
                  write(6,*) 'end of increment: ', KINC
                  write(6,*) 'GND & backstress calculations completed!'
c
c
c

              endif
c
          endif
c
      endif
c
c
c
      RETURN
      END
c
c
c
c
c
c
c
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c
      use calculations, only: calcs

      use globalvars, only : phasefielddamage, creep_param
      use globalvars, only : creepphasefieldflag

      use phasefieldfracture, only : moose_interface_input
      use phasefieldfracture, only : moose_interface_output

c
c
c
      implicit none
      integer,                        intent(in) ::
     & nDi,       !< Number of direct stress components at this point
     & nShr,      !< Number of engineering shear stress components
     & nTens,     !< Size of the stress / strain components (NDI + NSHR)
     & nStatV,    !< Number of solution-dependent state variables
     & nProps,    !< User-defined number of material constants
     & noEl,      !< element number
     & nPt,       !< integration point number
     & layer,     !< layer number (shell elements etc.)
     & kSpt,      !< section point within the current layer
     & kStep,     !< step number
     & kInc       !< increment number
      character(len=80),              intent(in) ::
     & cmname     !< uses-specified material name, left justified
      real(8),                        intent(in) ::
     & DTIME,
     & TEMP,
     & DTEMP,
     & CELENT
      real(8), dimension(1),          intent(in) ::
     & PREDEF,
     & DPRED
      real(8), dimension(2),          intent(in) ::
     & TIME       !< step time/total time at begin, of the current inc.
      real(8), dimension(3),          intent(in) ::
     & COORDS
      real(8), dimension(nTens),       intent(in) ::
     & STRAN,     !< total strains at beginning of the increment
     & DSTRAN     !< strain increments
      real(8), dimension(nProps),      intent(in) ::
     & PROPS
      real(8), dimension(3,3),         intent(in) ::
     & DROT,      !< rotation increment matrix
     & DFGRD0,    !< F at beginning of increment
     & DFGRD1     !< F at end of increment
      real(8),                         intent(inout) ::
     & PNEWDT,    !< ratio of suggested new time increment
     & SSE,       !< specific elastic strain engergy
     & SPD,       !< specific plastic dissipation
     & SCD,       !< specific creep dissipation
     & RPL,       !< volumetric heat generation
     & DRPLDT     !< varation of RPL with respect to the temperature
      real(8), dimension(nTens),       intent(inout) ::
     & STRESS     !< stress tensor at the beginning of the increment
      real(8), dimension(nStatV),      intent(inout) ::
     & STATEV     !< solution-dependent state variables
      real(8), dimension(nTens),       intent(out) ::
     & DDSDDT,
     & DRPLDE
      real(8), dimension(nTens,nTens), intent(out) ::
     & DDSDDE
c
      real(8) sigma(6), jacobi(6,6)
c
c     To avoid warning messages set the following to zero
      DDSDDT = 0.0d+0
      DRPLDE = 0.0d+0
c
c
c
c
c
c     Phase field damage model is added by Nicolò Grilli.
c     Get information for moose
c     phase field damage model
c     through state variables
c     9 state variables must be declared
c     for the phase field damage model
c
      if (phasefielddamage.eq.1d+0 .or. 
     & creepphasefieldflag.eq.1d+0) then
c
          call moose_interface_input(NOEL,NPT,STATEV,NSTATV)
c
      end if
c
c
c
c     Perform crystal plasticity/j2 plasticity calculations
      call calcs(DFGRD0,DFGRD1,TIME(2),DTIME,TEMP,KINC,NOEL,NPT,
     &            sigma,jacobi,PNEWDT,COORDS)
c
c

c     3D
      if (NTENS.eq.6) then

          STRESS = sigma
          DDSDDE = jacobi

c     Plane Strain
      elseif (NTENS.eq.4) then

          STRESS=0.0
          STRESS(1:4) = sigma(1:4)
          DDSDDE=0.0
          DDSDDE(1,1) = jacobi(1,1)
          DDSDDE(1,2) = jacobi(1,2)
          DDSDDE(2,1) = jacobi(2,1)
          DDSDDE(2,2) = jacobi(2,2)
          DDSDDE(3,1) = jacobi(3,1)
          DDSDDE(3,2) = jacobi(3,2)
          DDSDDE(4,4) = jacobi(4,4)

c     Plane Stress
      elseif (NTENS.eq.3) then

          STRESS=0.0
          STRESS(1) = sigma(1)
          STRESS(2) = sigma(2)
          STRESS(3) = sigma(4)
          DDSDDE=0.0
          DDSDDE(1,1) = jacobi(1,1)
          DDSDDE(1,2) = jacobi(1,2)
          DDSDDE(2,1) = jacobi(2,1)
          DDSDDE(2,2) = jacobi(2,2)
          DDSDDE(3,3) = jacobi(4,4)


      endif


c
c
c     Phase field damage model is added by Nicolò Grilli.
c     Send information for moose
c     Phase field damage model
c
      if (phasefielddamage.eq.1d+0 .or.
     & creepphasefieldflag.eq.1d+0) then

          call moose_interface_output(NOEL,NPT,STATEV,NSTATV)

      end if
c
c
c
c
c      write(6,*) 'KINC',KINC
c      write(6,*) 'NOEL',NOEL
c      write(6,*) 'NPT',NPT
c      write(6,*) 'KSTEP',KSTEP
c      write(6,*) 'DFGRD0',DFGRD0
c      write(6,*) 'DFGRD1',DFGRD1
c      write(6,*) 'STRESS',STRESS
c      write(6,*) 'DDSDDE',DDSDDE
c
c
c
c
      RETURN
      END
c
c
c
c
c
c
c
c
c
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
      use globalvars, only: global_gammadot,
     & global_state, numslip, numstvar, output_vars, global_gamma_sum,
     & global_sigma, numel, numip, foldername, grainmorph, maxnumslip,
     & phaseind, global_f_ep_c
      use globalsubs, only: convert6to3x3, vonmises_stress
      use calculations, only: calculate_misorientation
      implicit none
c
c
      integer,                        intent(in) ::
     & NUVARM,
     & NOEL,
     & NPT,
     & LAYER,
     & KSPT,
     & KSTEP,
     & KINC,
     & NDI,
     & NSHR
      real(8), dimension(3,3),        intent(in) ::
     & DIRECT,
     & T
      real(8), dimension(2),          intent(in) ::
     & TIME
      real(8), dimension(1),          intent(in) ::
     & DTIME
      character(len=80),              intent(in) ::
     & CMNAME,
     & ORNAME
      real(8), dimension(3),          intent(in) ::
     & COORD
      real(8), dimension(*),          intent(in) ::
     & JMAC,
     & JMATYP,
     & MATLAYO,
     & LACCFLA
      real(8), dimension(NUVARM),     intent(out) ::
     & UVAR

c
c
c
c
      integer i, j, k, varno, noph, nss
      real(8) mis
c      real(8) sigma_av(7), sigma33(3,3), evm
c
c
c     The dimensions of the variables FLGRAY, ARRAY and JARRAY
c     must be set equal to or greater than 15.
c
c
c
c
      varno=1d+0
c
c     Phase index of the element
      noph=phaseind(NOEL)
c
c     Number of slip systems
      nss=numslip(noph)
c
c     Misorientation with respect to the initial orientation
c     Misorientation is calculated in case it is requested since it is tedious!
      if (output_vars(1).eq.1d+0) then
c
c
      call calculate_misorientation(NOEL,NPT,mis)
c
c      UVAR(varno) = mis
      UVAR(varno) = global_f_ep_c(NOEL,NPT)
c
      varno=varno + 1d+0
c
c
      endif
c
c
c     Cumlative slip
      if (output_vars(2) .eq. 1d+0) then
        UVAR(varno) = global_gamma_sum(NOEL,NPT)
        varno=varno + 1d+0
      endif
c
c
c
c
c
c     Average state variables
c     Note that the number of outpus could be 1 or 2 depending
c     on the number of state variables
      if (output_vars(3).eq.1d+0) then
c
c         Loop through the number of state variables
      do i=1,numstvar
c
c             Calculate the average of the state variable
        UVAR(varno)=sum(dabs(global_state(NOEL,NPT,1:nss,i)))/nss
        varno=varno + 1d+0
c
      enddo
c
c
      endif
c
c
c     Slip rates
      if (output_vars(4).eq.1d+0) then
        do i=1,maxnumslip
          UVAR(varno) = global_gammadot(NOEL,NPT,i)
          varno=varno + 1d+0
        enddo
      endif
c
c
c     State variables per slip system
c     Note that the number of outpus could be 1 or 2 depending
c     on the number of state variables
      if (output_vars(5).eq.1d+0) then
        do i=1,numstvar
          do j=1,maxnumslip
            UVAR(varno)= dabs(global_state(NOEL,NPT,j,i))
            varno=varno + 1d+0
          enddo
        enddo
      endif
c
c
c
      RETURN
      END
	  
