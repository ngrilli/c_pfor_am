c Chris Allen
c Edward Horton
c Eralp Demir
c Aug. 12th, 2021 - 1st working version
c
      
      include "globalvars.f"
      include "globalsubs.f"
      include "lengthscale.f"
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
      use globalvars, only: foldername
c
      implicit none
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
c      
c
c
c      
c      
      foldername= "/home/nicolo/projects/c_pfor_am/test/tests/umat"
c      
c
c      
c
c
c     At the start of the analysis (only ONCE!)
      if (LOP.eq.0) then
          write(6,*) 'initialization has started!'
           write(6,*) '********************************'
          call initialize_all(foldername)
          write(6,*) '********************************'
          write(6,*) 'initialization has ended!'
      endif
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
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c
      use calculations, only: calcs
c
c
c
      implicit none
      integer,                        intent(in) :: 
     & nDi,       !< Number of direct stress components at this point
     & nShr,      !< Number of engineering shear stress components at this point
     & nTens,     !< Size of the stress or strain component array (NDI + NSHR)
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
c     To avoid warning messages set the following to zero
      DDSDDT = 0.0d+0
      DRPLDE = 0.0d+0
c
c 
c      
c      
c      integer :: defaultNumThreadsInt                 !< default value set by Abaqus
c      include "omp_lib.h"
c      defaultNumThreadsInt = omp_get_num_threads()    ! remember number of threads set by Marc
c      call omp_set_num_threads(defaultNumThreadsInt)  ! set number of threads for parallel execution set by DAMASK_NUM_THREADS      
c      
c
c      write(6,*) 'KINC',KINC
c      write(6,*) 'NOEL',NOEL
c      write(6,*) 'NPT',NPT
c      write(6,*) 'KSTEP',KSTEP
c      write(6,*) 'DFGRD0',DFGRD0
c      write(6,*) 'DFGRD1',DFGRD1
c      write(6,*) 'STRESS',STRESS
c      
c
c
c     Perform all the calculations     
      call calcs(DFGRD0,DFGRD1,TIME(2),DTIME,TEMP,KINC,NOEL+1,NPT+1,
     &            STRESS,DDSDDE,PNEWDT,COORDS)

c      
c
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
      use globalvars, only: global_ori, global_Fe, global_gammadot,
     & global_state, numslip, numstvar, output_vars, global_gamma_sum,
     & mattyp, global_sigma, numel, numip, foldername
      use globalsubs, only: misorientation, polar, convert6to3x3,
     & vonmises_stress
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
      integer i, j, k, typ, varno
      real(8) U(3,3)
c      real(8) sigma_av(7), sigma33(3,3), evm
      real(8) g1(3,3), g2(3,3), dg(3,3), mis, ax(3)
c
c     The dimensions of the variables FLGRAY, ARRAY and JARRAY
c     must be set equal to or greater than 15.
c
c
c --- ED HORTON EDIT ---
c
      varno=1d+0
c
c     
c     Misorientation with respect to the initial orientation
c     Misorientation is calculated in case it is requested since it is tedious!
      if (output_vars(1) .eq. 1d+0) then
c
c          
c         Set misorientation to zero for isotropic material type
          if (mattyp.eq.0) then
c              
              mis = 0.0d+0
c              
c         Calculate misorientations if material is not isotropic              
          else
c              
c             FCC - cubic symmetry operators              
              if (mattyp.eq.1d+0) then
                  typ = 1d+0
c             BCC - cubic symmetry operators                     
              elseif (mattyp.eq.2d+0) then
                  typ = 1d+0
c             HCP - hexagonal symmetry operators
              elseif (mattyp.eq.3d+0) then
                  typ = 3d+0
              endif
c              
c              
c	        Initial orientation
	        g1=global_ori(NOEL,NPT,:,:)
c	        Polar decomposition of the elastic part of the deformation gradient
	        call polar(global_Fe(NOEL,NPT,:,:),g2,U)
!c	        If singularity does no exist calculate
!	        if (sing.eq.0d+0)	then
c		        Calcualte misorientation using the subroutine
		        call misorientation(g1,g2,typ,ax,mis,dg)
         !     else
         !         mis=0.0d+0
	        !endif                  
c                  
          endif
c          
c
          UVAR(varno) = mis
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
      if (output_vars(3) .eq. 1d+0) then
c          
c         Loop through the number of state variables
          do i=1,numstvar
c              
c             Calculate the average of the state variable    
              UVAR(varno) = sum(global_state(NOEL,NPT,:,i))/numslip
              varno=varno + 1d+0
c                  
          enddo
c          
c      
      endif
c
c
c     Slip rates
      if (output_vars(4) .eq. 1d+0) then
            do i=1,numslip
                  UVAR(varno) = global_gammadot(NOEL,NPT,i)
                  varno=varno + 1d+0
            enddo
      endif
c
c      
c     State variables per slip system
c     Note that the number of outpus could be 1 or 2 depending
c     on the number of state variables
      if (output_vars(5) .eq. 1d+0) then
            do i=1,numstvar
                  do j=1,numslip
                        UVAR(varno)= global_state(NOEL,NPT,j,i)
                        varno=varno + 1d+0
                  enddo
            enddo
      endif
c ---ED HORTON EDIT END ---
c          
c
c
!c     Add the average of the stresses over all of the elements and IPs
!c     Do this calculation only ONCE!
!      if (NOEL.eq.1) then
!          if (NPT.eq.1) then
!              sigma_av=0.0
!              do i=1,numel
!                  do j=1,numip
!c                      
!c                     Calculate the Von-Mises stress                      
!                      call convert6to3x3(global_sigma(i,j,:),sigma33)
!                      call vonmises_stress(sigma33,evm)
!                      sigma_av(7) = sigma_av(7) + evm
!c                      
!                      do k =1,6
!                          sigma_av(k)=sigma_av(k) + global_sigma(i,j,k)
!                      enddo
!                  enddo
!              enddo
!              sigma_av = sigma_av / numel / numip
!c
!c
!c      
!c
!c     Print the results to a separate .txt file
!900           format(I4, 1X, 8(F7.2, 1X), /)
!              open(1000, file = foldername //
!     &'/avCauchyStress.txt', action='readwrite', status='old')
!c      
!c             read dummy lines      
!              read(1000,*)
!              do i=1,KINC-1
!                  read(1000,*)
!              enddo
!c      
!c      
!c             Write the following in the same order as:
!c             increment number, time, average stress components (11-22-33-12-13-23), and average Von-Mises stress      
!              write(1000,900) KINC, TIME(2), sigma_av(1:7)
!c      
!              close(1000)
!c      
!          endif
!      endif
c
      RETURN
      END