c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c Nicolò Grilli
c
c
c
c     This module contains all the global variables stored in the code
      module globalvars
      implicit none

c
c     ______________________________________________________
C	MESH CONSTANTS
c	There must be only SINGLE type of element in the mesh!
c	Number of integration points per element
c	Total number of elements in the mesh
	integer, public :: numel
c	Number of integration points per element
      integer, public :: numip
c     Number of nodes per element
      integer, public :: nnpe
c     Number of grains in the mesh (computed from materials.dat file)
      integer, public :: numgrain
c     ______________________________________________________
c
c
c
c     ______________________________________________________
c     FOLDER PARAMETERS
      character(len=:), allocatable, public   :: foldername
c     ______________________________________________________
c
C	IMPORTANT VARIABLES THAT WERE STORED AT THE TIME INCREMENT 't'
c	Plastic part of the deformation gradient at 't'
	real(8), allocatable, public ::  global_Fp(:,:,:,:)
	real(8), allocatable, public ::  global_Fp_t(:,:,:,:)
c	Elastic part of the deformation gradient at 't'
	real(8), allocatable, public ::  global_Fe(:,:,:,:)
	real(8), allocatable, public ::  global_Fe_t(:,:,:,:)
c	Residual deformations (inverse) defined at the beginning of deformation
	real(8), allocatable, public ::  global_Fr0(:,:,:,:)
c	Slip resistance at 't'
      real(8), allocatable, public ::  global_state0(:,:,:,:)
	real(8), allocatable, public ::  global_state(:,:,:,:)
	real(8), allocatable, public ::  global_state_t(:,:,:,:)
c     2nd Piola-Kirchoff stress at 't' (vectorized)
	real(8), allocatable, public ::  global_S(:,:,:)
	real(8), allocatable, public ::  global_S_t(:,:,:)
c   Damaged 2nd Piola-Kirchoff stress at 't' (vectorized)
	real(8), allocatable, public ::  global_S_damaged(:,:,:)
	real(8), allocatable, public ::  global_S_damaged_t(:,:,:)
c	Global orientation matrix
	real(8), allocatable, public ::  global_ori(:,:,:,:)
c	Global total slip per slip system (-) --- output
	real(8), allocatable, public ::	 global_gamma(:,:,:)
      real(8), allocatable, public ::	 global_gamma_t(:,:,:)
c	Global total overall slip (-) --- output
	real(8), allocatable, public ::	 global_gamma_sum(:,:)
      real(8), allocatable, public ::	 global_gamma_sum_t(:,:)
c	Global total sliprate per slip system (-) --- output
	real(8), allocatable, public ::	 global_gammadot(:,:,:)
	real(8), allocatable, public ::	 global_gammadot_t(:,:,:)
c	Global old jacobian matrix
	real(8), allocatable, public ::	 global_jacob_t(:,:,:,:)
c	Global current jacobian matrix
	real(8), allocatable, public ::	 global_jacob(:,:,:,:)
c	Global current stress vector - Cauchy
      real(8), allocatable, public ::	 global_sigma(:,:,:)
      real(8), allocatable, public ::	 global_sigma_t(:,:,:)
c   Global current damaged stress vector - Cauchy
      real(8), allocatable, public ::    global_sigma_damaged(:,:,:)
      real(8), allocatable, public ::    global_sigma_damaged_t(:,:,:)
c     Global IP coordinates
      real(8), allocatable, public ::	 global_coords(:,:,:)

c     Damage variable
      real(8), allocatable, public :: global_damage(:,:)
c     tensile part of the Helmholtz free energy
c	  which is degraded by damage
	  real(8), allocatable, public :: global_F_pos(:,:)
c     compressive part of the Helmholtz free energy
c	  which is degraded by damage
      real(8), allocatable, public :: global_F_neg(:,:)
c     tensile part of the second Piola-Kirchoff stress
c     this is needed for calculating _dstress_dc
      real(8), allocatable, public :: global_pk2_pos(:,:,:,:)
c     plastic work in finite strain formulation
      real(8), allocatable, public :: global_Wp(:,:)
      real(8), allocatable, public :: global_Wp_t(:,:)
c     Creep damage Parameters
      real(8), allocatable, public :: global_f_ep_c(:,:) 
      real(8), allocatable, public :: global_f_ep_c_t(:,:)
      integer, public :: creepphasefieldflag

c	Constants used everywhere
c     ______________________________________________________
c	Constant pi
	real(8), parameter, public :: pi = 3.141592654d+0
c	Taylor factor for a isotropic polycrytal aggregte
	real(8), parameter, public :: TF = 3.1d+0
c     Universal gas contant [J/mol/K]
      real(8), parameter, public :: Rgas = 8.31d+0
c     Boltzman contant [m2 kg s-2 K-1 ]
      real(8), parameter, public :: KB = 1.380649d-23
c     Number of slip systems
      integer , allocatable, public :: numslip(:)
c     Maximum number of slip systems (for array allocation)
      integer, public :: maxnumslip
c	Schmid tensor
      real(8)	, allocatable, public :: Schmid(:,:,:,:)
c	Climb tensor
      real(8)	, allocatable, public :: Climb(:,:,:,:)
c	Schmid tensor - transpose
      real(8)	, allocatable, public :: SchmidT(:,:,:,:)
c	Vectorized Schmid tensor
      real(8)	, allocatable, public :: Schmid_vec(:,:,:)
c	Vectorized Climb tensor
      real(8)	, allocatable, public :: Climb_vec(:,:,:)
c     Backstress tensor
      real(8)	, allocatable, public :: BS_dyad(:,:,:,:,:)
c	Dyadic product of Schmid tensor with Schmid tensor
      real(8)	, allocatable, public :: SchxSch(:,:,:,:,:,:)
c	Normalized slip directions
      real(8)	, allocatable, public :: b_slip(:,:,:)
c	Normalized slip plane normals
      real(8)	, allocatable, public :: n_slip(:,:,:)
c	Vector product of Burgers vector with slip plane normals
      real(8)	, allocatable, public :: l_slip(:,:,:)
c
C	CRYSTAL ORIENTATION ANGLES
c     ______________________________________________________
c
c	Euler angles (deg)
	real(8)	, allocatable, public :: Euler(:,:)
c
C	CONSTITUTIVE PARAMETERS
c     ______________________________________________________
c
c
c     Total number of materials
      integer, public :: totphase
c
c     Flag to indicate the existing material types
c     FCC/BCC/HCP (1/2/3)
c     Size: totphase
      integer	, allocatable, public ::	phases(:)
c     Index of the phase defined by the user
c     Size: total number of elements
      integer	, allocatable, public ::	phaseind(:)
c
c
c
c
c
c	Total number of state variables
      integer, public :: numstvar
c
c	Residual deformation defined or not
      integer, public :: resdef
c	Total amount of time for residual deformation
      integer, public :: tres
c
c	Thermally-coupled or mechanical problem
      integer, public :: thermo
c
c     Flag indicating coupling with fracture is active
      integer, public :: phasefielddamage
c
c     Isochoric free energy contribution to damage
c     1 = full contribution
c     0 = no contribution
      real(8), public :: cpl_contribution_to_damage
c
c     Hard coded plastic work energy contribution to damage
c     1 = full contribution
c     0 = no contribution
      real(8), public :: plastic_work_contribution_to_damage
c
c	Temperature (K) - initial temperature
      real(8), public :: temp0
c
c	Temperature dependent properties
c     "0": no temperature dependence
c     "1": linear temperature depedence
      integer, public :: tempdep
c
c
c     Flag to indicate material type
      integer, allocatable, public ::  phaseID(:)
c
c     Flag to indicate grain number
      integer, allocatable, public ::  grainID(:)
c
c	Elasticity constants (MPa)
      real(8), allocatable, public ::  elas_param(:,:)

c
c	Global variables that store the single crystal elasticity matrices
c	Elasticity tensor at the reference frame
      real(8), allocatable, public :: elas3333(:,:,:,:,:)

c	Work-conjugate elasticity matrix to stress at the reference frame
      real(8), allocatable, public :: elas66(:,:,:)
c
c	Work-conjugate elasticity matrix for isotropic phase
      real(8), allocatable, public :: elas66_iso(:,:,:)
c
c	modelno: 1/2/3/4/5/6
      integer, public :: modelno
c
c	modelno: 0/1/2/3/4
      integer, public :: creepno
c
c	interno: 0/1/2/3/4
      integer, public :: interno
c
C	LENGTHSCALE PARAMETERS
c     ______________________________________________________
c
c	Strain hardening interaction: 0 (OFF) /1 (ON)
      integer, public :: GSeffect

c	Strain rate parameters (1:param)
      real(8), allocatable, public :: sliprate_param(:,:)
c
c	Creep parameters (1:param)
      real(8), allocatable, public :: creep_param(:,:)

c	Strain hardening parameters (1:param)
      real(8), allocatable, public :: sliphard_param(:,:)
c
c	Strain hardening parameters (1:param)
      real(8), allocatable, public :: slipint_param(:,:)
c
c	Hardening interaction matrix
      real(8), allocatable, public :: intmat(:,:,:)
c	Irradiated Hardening matrices for gnds and DLs
      real(8), allocatable, public :: intmat1(:,:,:)
      real(8), allocatable, public :: intmat2(:,:,:)

c	Length scale parameters (1:3)
      real(8), allocatable, public :: grainsize_param(:,:)
      integer, allocatable, public :: grainsize_init(:,:)
c     grainsize_param(1): "k", linear hardening coefficient
c     grainsize_param(2): "c", hardening exponent for grain boundary mismatch
c     grainsize_param(3): "SF", size conversion factor from mesh units to millimeter units
c
c     Global variables related with length scale calculation
c      real(8), allocatable, public :: oriensh(:,:)
      real(8), allocatable, public :: grainori(:,:)
      real(8), allocatable, public :: nodex(:,:)
      real(8), allocatable, public :: nodey(:,:)
      real(8), allocatable, public :: nodez(:,:)
      real(8), allocatable, public :: boundgrain(:,:)
      real(8), allocatable, public :: elcent(:,:)
      integer, public :: nodeout

c     Size to indicate the grain morphology
      real(8), allocatable, public :: grainmorph(:,:,:)

c     ______________________________________________________


C     STRAIN GRADIENT PARAMETERS
c     ______________________________________________________

c	Strain gradient: 0 (OFF) / 1 (slip gradients)
      integer, public :: GNDeffect

c     Global initialization flag for IP coordinates
      integer, allocatable, public ::	 coords_init(:,:)
c
c     Initialiation of elemental GND calculations
      integer GND_init
      data    GND_init     /0/
cc     GND time cutback
c      integer GND_cutback
c	data    GND_cutback     /0/
c     Integration point coordinates in the iso-parametric space
c     Parameters are consistent with ABAQUS - g, h, r
      real(8), allocatable, public :: IPghr(:,:)
c     Integration weights
      real(8), allocatable, public :: wtghr(:)
c     Inverse of Nmat: IP to node mapping
      real(8), allocatable, public :: invNmat(:,:)
c     Interpolation function derivatives
      real(8), allocatable, public :: dNmat(:,:,:)
c     Gradient mapping
      real(8), allocatable, public :: Gmat(:,:,:)
c
c     Gradient map for elements
      real(8), allocatable, public :: gradIP2IP(:,:,:,:)
c     Gradient map for elements
      real(8), allocatable, public :: drhoGND(:,:,:,:)
c     ______________________________________________________


c     Equivalent isotropic Young's modulus
      real(8), allocatable, public :: E(:)

c     Equivalent isotropic Poisson's ratio
      real(8), allocatable, public :: nu(:)

c     Equivalent isotropic Shear modulus
      real(8), allocatable, public :: G(:)

c     Equivalent isotropic Bulk modulus
      real(8), allocatable, public :: kappa(:)

c

c --- ED HORTON EDIT
c	OUTPUT VARIABLE CHECKS
c     1: misorientation angle
c     2: cumulative slip
c     3: average of state variables over slip systems
c     4: slip rates per slip system
c     5: state variables per slip system
      integer, public :: output_vars(5)
c --- ED HORTON EDIT END


c
C	NUMERICAL CONSTANTS
c     ______________________________________________________

c	Inner loop exit tolerance (ABSOLUTE)
      real(8), public :: innertol
c	Outer loop exit tolerance (RELATIVE)
      real(8), public :: outertol
c     Outer absolute tolerance (ABSOLUTE)
      real(8), allocatable, public :: outerabstol(:,:)
c	Maximum number of iterations allowed for inner loop
      integer, public :: innoitmax
c	Maximum number of iterations allowed for outer loop
      integer, public :: ounoitmax
c	Constant for the jacobian calculation
      real(8), public :: deps
c	Method of jacobian calculation (1: perturbation, 2: analytical)
      integer, public :: mtdjaco
c	Number of increments skipped for the jacobian calculation
      integer, public :: njaco

c     Critical threshold value for stress-update algorithm
      real(8), allocatable, public :: dS_cr(:)
c
c     Factor for critical threshold stress calculation
      real(8), public :: dSratio_cr

c     Specified amount of slip for time-stepping algorithm
      real(8), public :: dgamma_s

c     Upper bound for time-stepping
      real(8), public :: ratio_ub


c     Lower bound for time-stepping
      real(8), public :: ratio_lb


c     Forward fraction of time
      real(8), public :: tstep_forw


c     Backward fraction of time
      real(8), public :: tstep_back

c	Small real number
      real(8), parameter, public :: smallnum = 1.0d-20

c	Large real number
      real(8), parameter, public :: largenum = 1.0d+20
c     ______________________________________________________
c
c
C	GLOBAL FLAGS
c     ______________________________________________________
c
c	Increment number that is stored
      integer, allocatable, public :: inc_old(:,:)
c
c     time step - used as a flag for state update
      real(8), allocatable, public :: t_old(:,:)
c
c
c     ______________________________________________________
c
c
c     Global variables used within the global subroutines
c     ______________________________________________________
c	Order of 4th rank tensor components (not necessarily symmetric)
c	11-22-33-12-21-13-31-23-32
      integer, public :: order9x9(9,2)
      data order9x9(1,:)		/1, 1/
      data order9x9(2,:)		/1, 2/
      data order9x9(3,:)		/1, 3/
      data order9x9(4,:)		/2, 1/
      data order9x9(5,:)		/2, 2/
      data order9x9(6,:)		/2, 3/
      data order9x9(7,:)		/3, 1/
      data order9x9(8,:)		/3, 2/
      data order9x9(9,:)		/3, 3/
c
c    Order for 4th rank tensor components (symmetric)
c	11-22-33-12-13-23
      integer, public :: order6x6(6,2)
      data order6x6(1,:)		/1, 1/
      data order6x6(2,:)		/2, 2/
      data order6x6(3,:)		/3, 3/
      data order6x6(4,:)		/1, 2/
      data order6x6(5,:)		/1, 3/
      data order6x6(6,:)		/2, 3/
c
c
c	Identity tensors used
      real(8), public :: I3(3,3)
      real(8), public :: I6(6,6)
      real(8), public :: I9(9,9)
      real(8), public :: eijk(3,3,3)
      real(8), public :: I3333(3,3,3,3)
c
c
c	BH stress states
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
      data BHstress(21,:) /-0.5,0.5,0,0.5,0.5,0/
      data BHstress(20,:) /0,-0.5,0.5,0,-0.5,-0.5/
      data BHstress(22,:) /-0.5,0.5,0,-0.5,0.5,0/
      data BHstress(24,:) /-0.5,0.5,0,-0.5,-0.5,0/
      data BHstress(23,:) /-0.5,0.5,0,0.5,-0.5,0/
      data BHstress(25,:) /0,0,0,0.5,0.5,-0.5/
      data BHstress(26,:) /0,0,0,0.5,-0.5,0.5/
      data BHstress(27,:) /0,0,0,-0.5,0.5,0.5/
      data BHstress(28,:) /0,0,0,0.5,0.5,0.5/
c
c
c	Symmetry operators for cubic crystals
c     These are necessary for misorientation calculation
c     4-3-2 symmetry
      real(8), public :: cubsym(24,3,3)
      data cubsym(1,1,:) /1, 0, 0/
      data cubsym(1,2,:) /0, 1, 0/
      data cubsym(1,3,:) /0, 0, 1/
c
      data cubsym(2,1,:) /0, 0, 1/
      data cubsym(2,2,:) /1, 0, 0/
      data cubsym(2,3,:) /0, 1, 0/
c
      data cubsym(3,1,:) /0, 1, 0/
      data cubsym(3,2,:) /0, 0, 1/
      data cubsym(3,3,:) /1, 0, 0/
c
      data cubsym(4,1,:) /0, -1, 0/
      data cubsym(4,2,:) /0, 0, 1/
      data cubsym(4,3,:) /-1, 0, 0/
c
      data cubsym(5,1,:) /0, -1, 0/
      data cubsym(5,2,:) /0, 0, -1/
      data cubsym(5,3,:) /1, 0, 0/
c
      data cubsym(6,1,:) /0, 1, 0/
      data cubsym(6,2,:) /0, 0, -1/
      data cubsym(6,3,:) /-1, 0, 0/
c
      data cubsym(7,1,:) /0, 0, -1/
      data cubsym(7,2,:) /1, 0, 0/
      data cubsym(7,3,:) /0, -1, 0/
c
      data cubsym(8,1,:) /0, 0, -1/
      data cubsym(8,2,:) /-1, 0, 0/
      data cubsym(8,3,:) /0, 1, 0/
c
      data cubsym(9,1,:) /0, 0, 1/
      data cubsym(9,2,:) /-1, 0, 0/
      data cubsym(9,3,:) /0, -1, 0/
c
      data cubsym(10,1,:) /-1, 0, 0/
      data cubsym(10,2,:) /0, 1, 0/
      data cubsym(10,3,:) /0, 0, -1/
c
      data cubsym(11,1,:) /-1, 0, 0/
      data cubsym(11,2,:) /0, -1, 0/
      data cubsym(11,3,:) /0, 0, 1/
c
      data cubsym(12,1,:) /1, 0, 0/
      data cubsym(12,2,:) /0, -1, 0/
      data cubsym(12,3,:) /0, 0, -1/
c
      data cubsym(13,1,:) /0, 0, -1/
      data cubsym(13,2,:) /0, -1, 0/
      data cubsym(13,3,:) /-1, 0, 0/
c
      data cubsym(14,1,:) /0, 0, 1/
      data cubsym(14,2,:) /0, -1, 0/
      data cubsym(14,3,:) /1, 0, 0/
c
      data cubsym(15,1,:) /0, 0, 1/
      data cubsym(15,2,:) /0, 1, 0/
      data cubsym(15,3,:) /-1, 0, 0/
c
      data cubsym(16,1,:) /0, 0, -1/
      data cubsym(16,2,:) /0, 1, 0/
      data cubsym(16,3,:) /1, 0, 0/
c
      data cubsym(17,1,:) /-1, 0, 0/
      data cubsym(17,2,:) /0, 0, -1/
      data cubsym(17,3,:) /0, -1, 0/
c
      data cubsym(18,1,:) /1, 0, 0/
      data cubsym(18,2,:) /0, 0, -1/
      data cubsym(18,3,:) /0, 1, 0/
c
      data cubsym(19,1,:) /1, 0, 0/
      data cubsym(19,2,:) /0, 0, 1/
      data cubsym(19,3,:) /0, -1, 0/
c
      data cubsym(20,1,:) /-1, 0, 0/
      data cubsym(20,2,:) /0, 0, 1/
      data cubsym(20,3,:) /0, 1, 0/
c
      data cubsym(21,1,:) /0, -1, 0/
      data cubsym(21,2,:) /-1, 0, 0/
      data cubsym(21,3,:) /0, 0, -1/
c
      data cubsym(22,1,:) /0, 1, 0/
      data cubsym(22,2,:) /1, 0, 0/
      data cubsym(22,3,:) /0, 0, -1/
c
      data cubsym(23,1,:) /0, 1, 0/
      data cubsym(23,2,:) /1, 0, 0/
      data cubsym(23,3,:) /0, 0, -1/
c
      data cubsym(24,1,:) /0, -1, 0/
      data cubsym(24,2,:) /1, 0, 0/
      data cubsym(24,3,:) /0, 0, 1/
c
c
c	Symmetry operators for tetragonal materials
c	4-2 symmetry
      real(8), public :: tetsym(8,3,3)
      data tetsym(1,1,:) /1, 0, 0/
      data tetsym(1,2,:) /0, 1, 0/
      data tetsym(1,3,:) /0, 0, 1/
c
      data tetsym(2,1,:) /-1, 0, 0/
      data tetsym(2,2,:) /0, 1, 0/
      data tetsym(2,3,:) /0, 0, -1/
c
      data tetsym(3,1,:) /1, 0, 0/
      data tetsym(3,2,:) /0, -1, 0/
      data tetsym(3,3,:) /0, 0, -1/
c
      data tetsym(4,1,:) /-1, 0, 0/
      data tetsym(4,2,:) /0, -1, 0/
      data tetsym(4,3,:) /0, 0, 1/
c
      data tetsym(5,1,:) /0, 1, 0/
      data tetsym(5,2,:) /-1, 0, 0/
      data tetsym(5,3,:) /0, 0, 1/
c
      data tetsym(6,1,:) /0, -1, 0/
      data tetsym(6,2,:) /1, 0, 0/
      data tetsym(6,3,:) /0, 0, 1/
c
      data tetsym(7,1,:) /0, 1, 0/
      data tetsym(7,2,:) /1, 0, 0/
      data tetsym(7,3,:) /0, 0, -1/
c
      data tetsym(8,1,:) /0, -1, 0/
      data tetsym(8,2,:) /-1, 0, 0/
      data tetsym(8,3,:) /0, 0, -1/
c
c	Symmetry operators for hexagonal materials
c	Constant for hexagonal materials
      real(8), public ::  hcpsym(12,3,3)
c	6-2 symmetry
      data hcpsym(1,1,:) /1, 0, 0/
      data hcpsym(1,2,:) /0, 1, 0/
      data hcpsym(1,3,:) /0, 0, 1/
      data hcpsym(2,1,:) /-0.5, 0.866025403d+0, 0/
      data hcpsym(2,2,:) /-0.866025403d+0, -0.5 , 0/
      data hcpsym(2,3,:) /0, 0, 1/
      data hcpsym(3,1,:) /-0.5, -0.866025403d+0, 0/
      data hcpsym(3,2,:) /0.866025403d+0, -0.5 , 0/
      data hcpsym(3,3,:) /0, 0, 1/
      data hcpsym(4,1,:) /0.5, 0.866025403d+0, 0/
      data hcpsym(4,2,:) /-0.866025403d+0, 0.5 , 0/
      data hcpsym(4,3,:) /0, 0, 1/
      data hcpsym(5,1,:) /-1, 0, 0/
      data hcpsym(5,2,:) /0, -1 , 0/
      data hcpsym(5,3,:) /0, 0, 1/
      data hcpsym(6,1,:) /0.5, -0.866025403d+0, 0/
      data hcpsym(6,2,:) /0.866025403d+0, 0.5 , 0/
      data hcpsym(6,3,:) /0, 0, 1/
      data hcpsym(7,1,:) /-0.5, -0.866025403d+0, 0/
      data hcpsym(7,2,:) /-0.866025403d+0, 0.5 , 0/
      data hcpsym(7,3,:) /0, 0, -1/
      data hcpsym(8,1,:) /1, 0, 0/
      data hcpsym(8,2,:) /0, -1 , 0/
      data hcpsym(8,3,:) /0, 0, -1/
      data hcpsym(9,1,:) /-0.5, 0.866025403d+0, 0/
      data hcpsym(9,2,:) /0.866025403d+0, 0.5 , 0/
      data hcpsym(9,3,:) /0, 0, -1/
      data hcpsym(10,1,:) /0.5, 0.866025403d+0, 0/
      data hcpsym(10,2,:) /0.866025403d+0, -0.5 , 0/
      data hcpsym(10,3,:) /0, 0, -1/
      data hcpsym(11,1,:) /-1, 0, 0/
      data hcpsym(11,2,:) /0, 1 , 0/
      data hcpsym(11,3,:) /0, 0, -1/
      data hcpsym(12,1,:) /0.5, -0.866025403d+0, 0/
      data hcpsym(12,2,:) /-0.866025403d+0, -0.5 , 0/
      data hcpsym(12,3,:) /0, 0, -1/
c
c
c     ______________________________________________________
c
c
c
c
c
      end module globalvars
