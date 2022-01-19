      ! Nicolò Grilli
      ! University of Bristol
      ! 19 Gennaio 2022

      ! Phase field fracture model
	  
      module phasefieldfracture
      implicit none
      contains	  
	  
      ! Decompose elastic Green-Lagrange strain 
	  ! into volumetric and non-volumetric
      ! and calculate elastic energy and stress  
      subroutine computeStrainVolumetric(ee,ce,
     & pk2_vec_new,el_no,ip_no)
	 
      use globalvars, only : elas66
      use globalvars, only : global_F_pos
      use globalvars, only : global_F_neg
      use globalvars, only : global_damage
      use globalvars, only : global_pk2_pos
	  
      use globalsubs, only : convert3x3to6
	  use globalsubs, only : convert6to3x3
	  
	  ! Elastic Green-Lagrange strain
      real(8), intent(in) :: ee(3,3)
	  
	  ! Elastic right Cauchy-Green tensor
      real(8), intent(in) :: ce(3,3)

      ! Element and IP number
      integer, intent(in) :: el_no, ip_no

      ! damaged second Piola-Kirchhoff stress
      real(8), intent(out) :: pk2_vec_new(6)	

      ! undamaged second Piola-Kirchhoff stress
      real(8) :: S_vec(6)	  
      
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
	  
      ! Elastic Green-Lagrange strain in vector form
	  real(8) :: ee_vec(6)
	  
      ! degradation function
	  real(8) :: D
	  
	  ! local variable for damage phase field
      real(8) :: c
	  
      call convert3x3to6(ee,ee_vec)
	  
	  ! Calculate second Piola-Kirchhoff stress
	  ! in vector form, undamaged
      S_vec = matmul(elas66,ee_vec)
	  
	  ! assign local damage phase field variable
      c = global_damage(el_no,ip_no)
	  
	  ! calculate degradation function
      D = (1.0-c)*(1.0-c)
	  
      ! this is temporary
	  ! now I am degrading all the stress
      pk2_neg_vec = 0.0
	  pk2_pos_vec = S_vec
	  
	  ! 
      call convert6to3x3(pk2_pos_vec,pk2_pos_mat)
	  
	  ! calculate degraded stress
      pk2_vec_new = D * pk2_pos_vec + pk2_neg_vec
	  
	  ! calculate degraded elastic energy
      ! this is temporary
	  ! now I am degrading all the elastic energy
      F_pos = 0.5 * dot_product(pk2_vec_new,ee_vec)

      ! this is temporary
	  ! now I am degrading all the elastic energy
      F_neg = 0.0
	  
	  ! assign to global variables
      global_F_pos(el_no,ip_no) = F_pos
      global_F_neg(el_no,ip_no) = F_neg
	  global_pk2_pos(el_no,ip_no,1:3,1:3) = pk2_pos_mat(1:3,1:3)
	  
	  return
      end subroutine computeStrainVolumetric
	  
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
      integer :: i, j
      
      ! 2 and 3 are components of the Helmholtz free energy
	  ! by convention
      STATEV(2) = global_F_pos(elem,ip)
      STATEV(3) = global_F_pos(elem,ip)

      ! 4 to 12 are the components of the positive
	  ! part of the second Piola-Kirchhoff stress
	  ! moose needs this for _dstress_dc
      do i = 1,3
        do j = 1,3
          STATEV(3+3*(i-1)+j) = global_pk2_pos(elem,ip,i,j)
        end do
      end do	  
	  
	  return
      end subroutine moose_interface_output

      end module phasefieldfracture
