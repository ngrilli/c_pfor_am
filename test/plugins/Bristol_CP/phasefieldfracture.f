      ! Nicolò Grilli
      ! University of Bristol
      ! 26 Gennaio 2022

      ! Phase field fracture model
	  
      module phasefieldfracture
      implicit none
      contains	  
	  
      ! Decompose elastic Green-Lagrange strain 
	  ! into volumetric and non-volumetric
      ! and calculate elastic energy and stress  
	  ! based on
	  ! Nicolo Grilli and Marisol Koslowski
	  ! The effect of crystal anisotropy and plastic response
	  ! on the dynamic fracture of energetic materials
	  ! Journal of Applied Physics 126, 155101 (2019)
	  
      subroutine computeStrainVolumetric(ee,ce,fe,
     & pk2_vec_new,el_no,ip_no)
	 
      use globalvars, only : elas66
      use globalvars, only : global_F_pos
      use globalvars, only : global_F_neg
      use globalvars, only : global_damage
      use globalvars, only : global_pk2_pos
	  
      use globalsubs, only : convert3x3to6
	  use globalsubs, only : convert6to3x3
      use globalsubs, only : invert3x3
	  
	  ! Elastic Green-Lagrange strain
      real(8), intent(in) :: ee(3,3)
	  
	  ! Elastic right Cauchy-Green tensor
      real(8), intent(in) :: ce(3,3)
	  
	  ! Elastic deformation gradient
      real(8), intent(in) :: fe(3,3)

      ! Element and IP number
      integer, intent(in) :: el_no, ip_no

      ! damaged second Piola-Kirchhoff stress
      real(8), intent(out) :: pk2_vec_new(6)	

      ! undamaged second Piola-Kirchhoff stress
      real(8) :: S_vec(6)	

	  ! inverse of the elastic right Cauchy-Green tensor
      real(8) :: invce(3,3)	

      ! determinant of the elastic right Cauchy-Green tensor
      real(8) :: detce
      
	  ! tensile part of the Helmholtz free energy
	  ! which is degraded by damage
      real(8) :: F_pos
	  
	  ! compressive part of the Helmholtz free energy
	  ! which is not degraded by damage
      real(8) :: F_neg	

	  ! tensile part of the second Piola-Kirchhoff stress
	  ! which is degraded by damage	  
      real(8) :: pk2_pos_vec(6)
      real(8) :: pk2_pos_mat(3,3)
	  
	  ! compressive part of the second Piola-Kirchhoff stress
	  ! which is not degraded by damage	  
      real(8) :: pk2_neg_vec(6)
      real(8) :: pk2_neg_mat(3,3)
	  
      ! Elastic Green-Lagrange strain in vector form
	  real(8) :: ee_vec(6)
	  
      ! degradation function
	  real(8) :: D
	  
	  ! local variable for damage phase field
      real(8) :: c
	  
 	  ! reference bulk modulus
      real(8) :: Kb
	  
      ! indices of iterations	  
      integer :: i, j
	  
      ! Je is relative elastic volume change
	  ! Je23 is Je^{2/3}
      real(8) :: Je, Je23 
	  
      ! delta is the trace of volumetric part of elastic Green-Lagrange strain
      real(8) :: delta
	  
      ! Positive and negative volumetric free energies
      real(8) :: a_pos_vol, a_neg_vol
	  
      ! Isochoric free energy
      real(8) :: a_pos_cpl
	  
      ! Anisotropic elasticity (Luscher 2017)
      ! Kb = K in (Luscher 2017)
      ! Kb = (1/9) I : C : I
      Kb = 0.0
      do i = 1,3
        do j = 1,3
          Kb = Kb + elas66(i,j)
        end do
      end do
      Kb = Kb / 9.0
	  
      call M33DET(fe,Je)

      Je23 = Je**(2.0/3.0)

      delta = 1.5 * (Je23 - 1.0)
	  
      ! Calculate volumetric free energy
      ! Equations 15 and 16 in Grilli, Koslowski, 2019
  
      if (Je >= 1.0) then ! expansion
	  
        a_pos_vol = 0.5 * Kb * delta * delta
        a_neg_vol = 0.0
	  
      else ! compression
	  
        a_pos_vol = 0.0
        a_neg_vol = 0.5 * Kb * delta * delta
	
      end if  

      call convert3x3to6(ee,ee_vec)

	  ! Calculate second Piola-Kirchhoff stress
	  ! in vector form, undamaged
      S_vec = matmul(elas66,ee_vec)
	  
	  ! calculate undamaged elastic energy
      a_pos_cpl = 0.5 * dot_product(S_vec,ee_vec)
	  
      ! Calculate isochoric free energy
      ! Equation 17 in Grilli, Koslowski, 2019
      a_pos_cpl = a_pos_cpl - 0.5 * Kb * delta * delta
	  
      ! Calculate positive and negative parts of the Piola-Kirchhoff stress
      ! Equation 18 in Grilli, Koslowski, 2019
	  
      if (Je >= 1.0) then ! expansion
  
        pk2_pos_vec = S_vec
        pk2_neg_vec = 0.0

      else ! compression
  
        call invert3x3(ce,invce,detce)
	
        pk2_neg_mat = Je23 * Kb * delta * invce
		
        call convert3x3to6(pk2_neg_mat,pk2_neg_vec)
	
        pk2_pos_vec = (-1.0) * pk2_neg_vec	
        pk2_pos_vec = pk2_pos_vec + S_vec

      end if

	  ! assign local damage phase field variable
      c = global_damage(el_no,ip_no)
	  
	  ! calculate degradation function
      D = (1.0-c)*(1.0-c)

	  ! calculate degraded stress	  
      ! positive part of the stress is degraded by D
      pk2_vec_new = D * pk2_pos_vec + pk2_neg_vec      
	  
      call convert6to3x3(pk2_pos_vec,pk2_pos_mat)
	  
      ! Positive and negative parts of the free energy
      ! Equations 13 and 14 in Grilli, Koslowski, 2019
      F_pos = a_pos_vol + a_pos_cpl
      F_neg = a_neg_vol
	  
	  ! assign to global variables
      global_F_pos(el_no,ip_no) = F_pos
      global_F_neg(el_no,ip_no) = F_neg
	  global_pk2_pos(el_no,ip_no,1:3,1:3) = pk2_pos_mat(1:3,1:3)
	  
	  return
      end subroutine computeStrainVolumetric
	  
	  
	  ! pass damage phase field to the UMAT through a state variable
	  
      subroutine moose_interface_input(elem,ip,STATEV,NSTATV)
	  
      use globalvars, only : global_damage

      ! element and integration point numbers
	  ! using moose convention, starting from zero
      integer, intent(in) :: elem, ip
	  
	  ! number of state variables
      integer, intent(in) :: NSTATV
	  
	  ! state variables
      real(8), intent(in) :: STATEV(NSTATV) 
	  
	  ! state variable 1 is by convention
	  ! the damage phase field
      global_damage(elem,ip) = STATEV(1)
	  
	  return
      end subroutine moose_interface_input
	  
	  
	  ! pass quantities needed for the damage phase field model 
	  ! from the UMAT using state variables
	  
      subroutine moose_interface_output(elem,ip,STATEV,NSTATV)
	  
      use globalvars, only : global_F_pos
      use globalvars, only : global_F_neg
      use globalvars, only : global_pk2_pos
	  
      ! element and integration point numbers
	  ! using moose convention, starting from zero
      integer, intent(in) :: elem, ip

	  ! number of state variables
      integer, intent(in) :: NSTATV
	  
	  ! state variables
      real(8), intent(out) :: STATEV(NSTATV) 

      ! indices of iterations
      integer :: i
      
      ! 2 and 3 are components of the Helmholtz free energy
	  ! by convention
      STATEV(2) = global_F_pos(elem,ip)
      STATEV(3) = global_F_neg(elem,ip)

      ! 4 to 9 are the components of the positive
	  ! part of the second Piola-Kirchhoff stress
	  ! which is a symmetric matrix
	  ! moose needs this for _dstress_dc
	  ! it follows Abaqus components order convention
	  ! for 6-vectors
      do i = 1,3
        STATEV(3+i) = global_pk2_pos(elem,ip,i,i)
      end do
      STATEV(7) = global_pk2_pos(elem,ip,1,2)
      STATEV(8) = global_pk2_pos(elem,ip,1,3)
      STATEV(9) = global_pk2_pos(elem,ip,2,3)	  

	  return
      end subroutine moose_interface_output
	  
	  
	  ! calculate determinant of 3x3 matrix
      ! A = input 3x3 matrix
      ! detA = output determinant	  
	  
      subroutine M33DET(A,detA)	  
	  
      real(8), intent(in) :: A(3,3)

      real(8), intent(out) :: detA
	  
      detA = 0.0
	  
      detA = detA + A(1,1)*A(2,2)*A(3,3)
      detA = detA - A(1,1)*A(2,3)*A(3,2)
      detA = detA - A(1,2)*A(2,1)*A(3,3)
      detA = detA + A(1,2)*A(2,3)*A(3,1)
      detA = detA + A(1,3)*A(2,1)*A(3,2)	  
      detA = detA - A(1,3)*A(2,2)*A(3,1)

      return
      end subroutine M33DET			

      end module phasefieldfracture
