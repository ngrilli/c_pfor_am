!     Chris Allen
!     Edward Horton
!     Eralp Demir
!     Aug. 12th, 2021 - 1st working version
!
!     This module contains all the global variables stored in the code
      module globalvars    
      implicit none

!	MESH CONSTANTS
!	There must be only SINGLE type of element in the mesh!
!	Number of integration points per element
!	Total number of elements in the mesh
	integer, public :: numel
!	Number of integration points per element
	integer, public :: numip
!     Number of grains in the mesh      
	integer, public :: numgrain

!     FOLDER PARAMETERS
      character(len=:), allocatable, public   :: foldername
!
!	IMPORTANT VARIABLES THAT WERE STORED AT THE TIME INCREMENT 't'
!	Plastic part of the deformation gradient at 't'
      real(8), allocatable, public ::  global_Fp(:,:,:,:)
      real(8), allocatable, public ::  global_Fp_t(:,:,:,:)
!	Elastic part of the deformation gradient at 't'
      real(8), allocatable, public ::  global_Fe(:,:,:,:)
      real(8), allocatable, public ::  global_Fe_t(:,:,:,:)
!	Slip resistance at 't'
      real(8), allocatable, public ::  global_state0(:,:,:,:)
      real(8), allocatable, public ::  global_state(:,:,:,:)
      real(8), allocatable, public ::  global_state_t(:,:,:,:)     
!     2nd Piola-Kirchoff stress at 't' (vectorized)
      real(8), allocatable, public ::  global_S(:,:,:)
      real(8), allocatable, public ::  global_S_t(:,:,:)      
!	Global orientation matrix
      real(8), allocatable, public ::  global_ori(:,:,:,:)
!	Global total slip per slip system (-) --- output
      real(8), allocatable, public ::	 global_gamma(:,:,:)
      real(8), allocatable, public ::	 global_gamma_t(:,:,:)
!	Global total overall slip (-) --- output
      real(8), allocatable, public ::	 global_gamma_sum(:,:)
      real(8), allocatable, public ::	 global_gamma_sum_t(:,:) 
!	Global total sliprate per slip system (-) --- output
      real(8), allocatable, public ::	 global_gammadot(:,:,:)            
!	Global old jacobian matrix
      real(8), allocatable, public ::	 global_jacob_t(:,:,:,:)      
!	Global current jacobian matrix
      real(8), allocatable, public ::	 global_jacob(:,:,:,:)         
!	Global current stress vector - Cauchy
      real(8), allocatable, public ::	 global_sigma(:,:,:)
      real(8), allocatable, public ::	 global_sigma_t(:,:,:)  
!	Constants used everywhere
!	Constant pi
      real(8), parameter, public :: pi=3.141592654d+0
!	Taylor factor for a isotropic polycrytal aggregte
      real(8), parameter, public :: TF = 3.1d+0
!   Universal gas contant [J/mol/K]
      real(8), parameter, public :: Rgas = 8.31d+0
!   Number of slip systems
      integer, public :: numslip
!	Schmid tensor
      real(8)	, allocatable, public :: Schmid(:,:,:)
!	Climb tensor
      real(8)	, allocatable, public :: Climb(:,:,:)      
!	Schmid tensor - transpose
      real(8)	, allocatable, public :: SchmidT(:,:,:)      
!	Vectorized Schmid tensor
      real(8)	, allocatable, public :: Schmid_vec(:,:)
!	Vectorized Climb tensor
      real(8)	, allocatable, public :: Climb_vec(:,:)      
!	Dyadic product of Schmid tensor with Schmid tensor
      real(8)	, allocatable, public :: SchxSch(:,:,:,:,:)
!	Normalized slip directions
      real(8)	, allocatable, public :: b_slip(:,:)
!	Normalized slip plane normals
      real(8)	, allocatable, public :: n_slip(:,:)
!	Vector product of Burgers vector with slip plane normals
      real(8)	, allocatable, public :: l_slip(:,:)
!
!	CRYSTAL ORIENTATION ANGLES
   
!	Euler angles (deg)
	real(8)	, allocatable, public :: Euler(:,:)
!      
!	CONSTITUTIVE PARAMETERS
!         
!      
!     Flag to indicate the existing material types
!     FCC/BCC/HCP (1/2/3)
      integer, public ::	mattyp
!      
!	Total number of state variables
	integer, public :: numstvar
!      
!	Thermally-coupled or mechanical problem
      integer, public :: thermo
!      
!	Temperature (K) - initial temperature
      real(8), public :: temp0
!
!	Temperature dependent properties
!     "0": no temperature dependence
!     "1": linear temperature depedence
      integer, public :: tempdep
!      
!      
!     Flag to indicate material type
      integer, allocatable, public ::  phaseID(:)
!      
!     Flag to indicate grain number
      integer, allocatable, public ::  grainID(:)      
!
!	Elasticity constants (MPa)
      real(8), allocatable, public ::  elas_param(:)
      
!
!	Global variables that store the single crystal elasticity matrices
!	Elasticity tensor at the reference frame 
	real(8), public :: elas3333(3,3,3,3)
      
!	Work-conjugate elasticity matrix to stress at the reference frame
	real(8), public :: elas66(6,6)
!      
!	Work-conjugate elasticity matrix for isotropic phase
	real(8), public :: elas66_iso(6,6)
!      
!	modelno: 1/2/3/4/5
	integer, public :: modelno
      
!	interno: 0/1/2
	integer, public :: interno
      
!	Strain hardening interaction: 0 (OFF) /1 (ON)
	integer, public :: GSeffect
      
!	Strain rate parameters (1:param)
	real(8), allocatable, public :: sliprate_param(:)
!      
!	Strain hardening parameters (1:param)
      real(8), allocatable, public :: sliphard_param(:)
!      
!	Strain hardening parameters (1:param)
      real(8), allocatable, public :: slipint_param(:)
!      
!	Hardening interaction matrix
	real(8), allocatable, public :: intmat(:,:)
      
!	Length scale parameters (1:2)
      real(8), public :: grainsize_param(2)
      integer, allocatable, public :: grainsize_init(:,:)
!     grainsize(2): "k", linear hardening coefficient
!     grainsize(3): "c", hardening exponent for grain boundary mismatch
!
!     Global variables related with length scale calculation
!      real(8), allocatable, public :: oriensh(:,:)
      real(8), allocatable, public :: grainori(:,:)
      real(8), allocatable, public :: nodex(:,:)
      real(8), allocatable, public :: nodey(:,:)
      real(8), allocatable, public :: nodez(:,:)
      real(8), allocatable, public :: boundgrain(:,:)
      real(8), allocatable, public :: elcent(:,:)
      integer, public :: nodeout


      

!     Equivalent isotropic Young's modulus
      real(8), public :: E
      
!     Equivalent isotropic Poisson's ratio
      real(8), public :: nu

!     Equivalent isotropic Shear modulus
      real(8), public :: G
!      

! --- ED HORTON EDIT
!	OUTPUT VARIABLE CHECKS
!     1: misorientation angle
!     2: cumulative slip
!     3: average of state variables over slip systems
!     4: slip rates per slip system
!     5: state variables per slip system
	integer, public :: output_vars(5)
! --- ED HORTON EDIT END          
      
      
!
!	NUMERICAL CONSTANTS

!	Inner loop exit tolerance (ABSOLUTE)
	real(8), public :: innertol
!	Outer loop exit tolerance (RELATIVE)
	real(8), public :: outertol
!	Maximum number of iterations allowed for inner loop
	integer, public :: innoitmax
!	Maximum number of iterations allowed for outer loop
	integer, public :: ounoitmax
!	Constant for the jacobian calculation
	real(8), public :: deps
!	Method of jacobian calculation (1: perturbation, 2: analytical)
	integer, public :: mtdjaco      
!	Number of increments skipped for the jacobian calculation
	integer, public :: njaco
      
!     Critical threshold value for stress-update algorithm
      real(8), public :: dS_cr
!      
!     Factor for critical threshold stress calculation
      real(8), public :: dSratio_cr

!     Specified amount of slip for time-stepping algorithm
      real(8), public :: dgamma_s
      
!     Upper bound for time-stepping
      real(8), public :: ratio_ub


!     Lower bound for time-stepping
      real(8), public :: ratio_lb      
      
      
!     Forward fraction of time
      real(8), public :: tstep_forw


!     Backward fraction of time
      real(8), public :: tstep_back      
      
      
!
!
!	GLOBAL FLAGS

!	Increment number that is stored
	integer, public :: inc_old
	data    inc_old     /0/
	real(8), public :: t_old
	data    t_old     /0.0d+0/
      
!
!
!
!     Global variables used within the global subroutines

!	Order of tensor components (not necessarily symmetric)
!	11-22-33-12-21-13-31-23-32
	integer, public :: order(9,2)
	data order(1,:)		/1, 1/
	data order(2,:)		/1, 2/
	data order(3,:)		/1, 3/
	data order(4,:)		/2, 1/
	data order(5,:)		/2, 2/
	data order(6,:)		/2, 3/
	data order(7,:)		/3, 1/
	data order(8,:)		/3, 2/
	data order(9,:)		/3, 3/
!
!
!
!	Identity tensors used
	real(8), public :: I3(3,3)
	real(8), public :: I6(6,6)
	real(8), public :: I9(9,9)
	real(8), public :: eijk(3,3,3)
!	
!
!	BH stress states
	real(8), public :: BHstress(28,6)
      data BHstress(1,:) /1,-1,0,0,0,0/
	data BHstress(2,:) /0,1,-1,0,0,0/
	data BHstress(3,:) /-1,0,1,0,0,0/
	data BHstress(4,:) /0,0,0,1,0,0/
	data BHstress(5,:) /0,0,0,0,1,0/
	data BHstress(6,:) /0,0,0,0,0,1/
	data BHstress(7,:) /0.5,-1,0.5,0,0.5,0/
	data BHstress(8,:) /0.5,-1,0.5,0,-0.5,0/
	data BHstress(9,:) /-1,0.5,0.5,0.5,0,0/
	data BHstress(10,:) /-1,0.5,0.5,-0.5,0,0/
	data BHstress(11,:) /0.5,0.5,-1,0,0,0.5/
	data BHstress(12,:) /0.5,0.5,-1,0,0,-0.5/
	data BHstress(13,:) /0.5,0,-0.5,0.5,0,0.5/
	data BHstress(14,:) /0.5,0,-0.5,-0.5,0,0.5/
	data BHstress(15,:) / 0.5,0,-0.5,0.5,0,-0.5/
	data BHstress(16,:) /0.5,0,-0.5,-0.5,0,-0.5/
	data BHstress(17,:) /0,-0.5,0.5,0,0.5,0.5/
	data BHstress(18,:) /0,-0.5,0.5,0,-0.5,0.5/
	data BHstress(19,:) /0,-0.5,0.5,0,0.5,-0.5/
	data BHstress(20,:) /0,-0.5,0.5,0,-0.5,-0.5/
	data BHstress(21,:) /-0.5,0.5,0,0.5,0.5,0/
	data BHstress(22,:) /-0.5,0.5,0,-0.5,0.5,0/
	data BHstress(23,:) /-0.5,0.5,0,0.5,-0.5,0/
	data BHstress(24,:) /-0.5,0.5,0,-0.5,-0.5,0/
	data BHstress(25,:) /0,0,0,0.5,0.5,-0.5/
	data BHstress(26,:) /0,0,0,0.5,-0.5,0.5/
	data BHstress(27,:) /0,0,0,-0.5,0.5,0.5/
	data BHstress(28,:) /0,0,0,0.5,0.5,0.5/
!
!
!	Symmetry operators for cubic crystals
!     These are necessary for misorientation calculation
!     4-3-2 symmetry
	real(8), public :: cubsym(24,3,3)
	data cubsym(1,1,:) /1, 0, 0/
	data cubsym(1,2,:) /0, 1, 0/
	data cubsym(1,3,:) /0, 0, 1/

	data cubsym(2,1,:) /0, 0, 1/
	data cubsym(2,2,:) /1, 0, 0/
	data cubsym(2,3,:) /0, 1, 0/

	data cubsym(3,1,:) /0, 1, 0/
	data cubsym(3,2,:) /0, 0, 1/
	data cubsym(3,3,:) /1, 0, 0/

	data cubsym(4,1,:) /0, -1, 0/
	data cubsym(4,2,:) /0, 0, 1/
	data cubsym(4,3,:) /-1, 0, 0/

	data cubsym(5,1,:) /0, -1, 0/
	data cubsym(5,2,:) /0, 0, -1/
	data cubsym(5,3,:) /1, 0, 0/

	data cubsym(6,1,:) /0, 1, 0/
	data cubsym(6,2,:) /0, 0, -1/
	data cubsym(6,3,:) /-1, 0, 0/

	data cubsym(7,1,:) /0, 0, -1/
	data cubsym(7,2,:) /1, 0, 0/
	data cubsym(7,3,:) /0, -1, 0/

	data cubsym(8,1,:) /0, 0, -1/
	data cubsym(8,2,:) /-1, 0, 0/
	data cubsym(8,3,:) /0, 1, 0/

	data cubsym(9,1,:) /0, 0, 1/
	data cubsym(9,2,:) /-1, 0, 0/
	data cubsym(9,3,:) /0, -1, 0/

	data cubsym(10,1,:) /-1, 0, 0/
	data cubsym(10,2,:) /0, 1, 0/
	data cubsym(10,3,:) /0, 0, -1/

	data cubsym(11,1,:) /-1, 0, 0/
	data cubsym(11,2,:) /0, -1, 0/
	data cubsym(11,3,:) /0, 0, 1/

	data cubsym(12,1,:) /1, 0, 0/
	data cubsym(12,2,:) /0, -1, 0/
	data cubsym(12,3,:) /0, 0, -1/

	data cubsym(13,1,:) /0, 0, -1/
	data cubsym(13,2,:) /0, -1, 0/
	data cubsym(13,3,:) /-1, 0, 0/

	data cubsym(14,1,:) /0, 0, 1/
	data cubsym(14,2,:) /0, -1, 0/
	data cubsym(14,3,:) /1, 0, 0/

	data cubsym(15,1,:) /0, 0, 1/
	data cubsym(15,2,:) /0, 1, 0/
	data cubsym(15,3,:) /-1, 0, 0/

	data cubsym(16,1,:) /0, 0, -1/
	data cubsym(16,2,:) /0, 1, 0/
	data cubsym(16,3,:) /1, 0, 0/

	data cubsym(17,1,:) /-1, 0, 0/
	data cubsym(17,2,:) /0, 0, -1/
	data cubsym(17,3,:) /0, -1, 0/

	data cubsym(18,1,:) /1, 0, 0/
	data cubsym(18,2,:) /0, 0, -1/
	data cubsym(18,3,:) /0, 1, 0/

	data cubsym(19,1,:) /1, 0, 0/
	data cubsym(19,2,:) /0, 0, 1/
	data cubsym(19,3,:) /0, -1, 0/

	data cubsym(20,1,:) /-1, 0, 0/
	data cubsym(20,2,:) /0, 0, 1/
	data cubsym(20,3,:) /0, 1, 0/

	data cubsym(21,1,:) /0, -1, 0/
	data cubsym(21,2,:) /-1, 0, 0/
	data cubsym(21,3,:) /0, 0, -1/

	data cubsym(22,1,:) /0, 1, 0/
	data cubsym(22,2,:) /1, 0, 0/
	data cubsym(22,3,:) /0, 0, -1/

	data cubsym(23,1,:) /0, 1, 0/
	data cubsym(23,2,:) /1, 0, 0/
	data cubsym(23,3,:) /0, 0, -1/

	data cubsym(24,1,:) /0, -1, 0/
	data cubsym(24,2,:) /1, 0, 0/
	data cubsym(24,3,:) /0, 0, 1/


!	Symmetry operators for tetragonal materials
!	4-2 symmetry
	real(8), public :: tetsym(8,3,3)
	data tetsym(1,1,:) /1, 0, 0/
	data tetsym(1,2,:) /0, 1, 0/
	data tetsym(1,3,:) /0, 0, 1/

	data tetsym(2,1,:) /-1, 0, 0/
	data tetsym(2,2,:) /0, 1, 0/
	data tetsym(2,3,:) /0, 0, -1/

	data tetsym(3,1,:) /1, 0, 0/
	data tetsym(3,2,:) /0, -1, 0/
	data tetsym(3,3,:) /0, 0, -1/

	data tetsym(4,1,:) /-1, 0, 0/
	data tetsym(4,2,:) /0, -1, 0/
	data tetsym(4,3,:) /0, 0, 1/

	data tetsym(5,1,:) /0, 1, 0/
	data tetsym(5,2,:) /-1, 0, 0/
	data tetsym(5,3,:) /0, 0, 1/

	data tetsym(6,1,:) /0, -1, 0/
	data tetsym(6,2,:) /1, 0, 0/
	data tetsym(6,3,:) /0, 0, 1/

	data tetsym(7,1,:) /0, 1, 0/
	data tetsym(7,2,:) /1, 0, 0/
	data tetsym(7,3,:) /0, 0, -1/

	data tetsym(8,1,:) /0, -1, 0/
	data tetsym(8,2,:) /-1, 0, 0/
	data tetsym(8,3,:) /0, 0, -1/

!	Symmetry operators for hexagonal materials
!	Constant for hexagonal materials
	real(8), parameter :: hex = 0.866025403
	real(8), public :: hcpsym(12,3,3)
!	6-2 symmetry
	data hcpsym(1,1,:) /1, 0, 0/
	data hcpsym(1,2,:) /0, 1, 0/
	data hcpsym(1,3,:) /0, 0, 1/
	data hcpsym(2,1,:) /-0.5, 0.866025403, 0/
	data hcpsym(2,2,:) /-0.866025403, -0.5 , 0/
	data hcpsym(2,3,:) /0, 0, 1/
	data hcpsym(3,1,:) /-0.5, -0.866025403, 0/
	data hcpsym(3,2,:) /0.866025403, -0.5 , 0/
	data hcpsym(3,3,:) /0, 0, 1/
	data hcpsym(4,1,:) /0.5, 0.866025403, 0/
	data hcpsym(4,2,:) /-0.866025403, 0.5 , 0/
	data hcpsym(4,3,:) /0, 0, 1/	
	data hcpsym(5,1,:) /-1, 0, 0/
	data hcpsym(5,2,:) /0, -1 , 0/
	data hcpsym(5,3,:) /0, 0, 1/
	data hcpsym(6,1,:) /0.5, -0.866025403, 0/
	data hcpsym(6,2,:) /0.866025403, 0.5 , 0/
	data hcpsym(6,3,:) /0, 0, 1/
	data hcpsym(7,1,:) /-0.5, -0.866025403, 0/
	data hcpsym(7,2,:) /-0.866025403, 0.5 , 0/
	data hcpsym(7,3,:) /0, 0, -1/
	data hcpsym(8,1,:) /1, 0, 0/
	data hcpsym(8,2,:) /0, -1 , 0/
	data hcpsym(8,3,:) /0, 0, -1/
	data hcpsym(9,1,:) /-0.5, 0.866025403, 0/
	data hcpsym(9,2,:) /0.866025403, 0.5 , 0/
	data hcpsym(9,3,:) /0, 0, -1/
	data hcpsym(10,1,:) /0.5, 0.866025403, 0/
	data hcpsym(10,2,:) /0.866025403, -0.5 , 0/
	data hcpsym(10,3,:) /0, 0, -1/
	data hcpsym(11,1,:) /-1, 0, 0/
	data hcpsym(11,2,:) /0, 1 , 0/
	data hcpsym(11,3,:) /0, 0, -1/
	data hcpsym(12,1,:) /0.5, -0.866025403, 0/
	data hcpsym(12,2,:) /-0.866025403, -0.5 , 0/
	data hcpsym(12,3,:) /0, 0, -1/


      end module globalvars
