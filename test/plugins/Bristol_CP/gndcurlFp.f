      ! Nicolò Grilli
      ! University of Bristol
      ! 30 Dicembre 2021

      ! Compute Curl of a 2nd order tensor
	  
      module gndcurlFp
      implicit none
      contains	  
	  
	  ! calculate the curl of the plastic deformation gradient Fp
      subroutine CalculateCurlFp(curlFp,NOEL)
	  
	  ! import coordinates and plastic deformation gradient
	  ! coordinates and Fp refer to the Gauss points
      use globalvars, only: global_coords, global_Fp,
     & smallnum
      use globalsubs, only: invertnxn, invert3x3
	  
	  ! element number
	  ! integration point is not needed because
	  ! this subroutine calculates curlFp
	  ! for all integration points of this element
      integer, intent(in) :: NOEL
	  
	  ! number of nodes
      integer, parameter :: nnodes = 8
	  
	  ! number of Gauss points
      integer, parameter :: ngp = 8
	  
	  ! curl of the plastic deformation gradient
      real(8), intent(out) :: curlFp(ngp,3,3)

      integer :: i, j
	  
	  ! temporary indices for Gauss points and nodes
      integer :: gp, node

	  ! Nodal coordinates in the reference element
      real(8) :: noderef(nnodes,3)
	  
	  ! Gauss point coordinates in the reference element
      real(8) :: gpref(ngp,3)
	  
	  ! temporary variables for the coordinates
	  ! of the node under consideration 
      ! in the reference element
	  ! when cycling over nodes
      real(8) :: ri, si, ti
	  
	  ! temporary variables for the factors
	  ! of the shape functions
      real(8) :: gr, gs, gt
	  
	  ! temporary variables for the derivatives
	  ! of the factors constituting
	  ! derivatives of the shape functions
      real(8) :: dgr, dgs, dgt
	  
	  ! shape functions calculated at the Gauss points
	  ! the first index indicates the Gauss point at
	  ! which the shape function is calculated
	  ! the second index indicate the shape function index
      ! corresponding to a specific node	  
	  ! they refer to the reference element
	  ! so their value is 1 on a single node of the reference element
	  real(8) :: shape_functions(ngp,nnodes)
	  
	  ! corresponding derivatives
	  ! with respect to x, y, z of the reference element
	  real(8) :: dshape_functions_dr(ngp,nnodes,3)
	  
	  ! this is the corresponding temporary variable
	  ! calculated at each Gauss point in the loop
      real(8) :: temp_dshape_functions_dr(nnodes,3)
	  
	  ! inverse of shape function matrix
	  ! note: a vector representing a quantity at the nodes
	  ! and multiplied by the shape function matrix
	  ! gives the same quantity at the Gauss points
	  ! therfore, the inverse, multiplied by a vector
	  ! representing a quantity at the Gauss points
	  ! gives the same quantity at the nodes
      real(8) :: inv_shape_functions(nnodes,ngp)
	  
	  ! derivatives of the shape functions 
	  ! in the global coordinates, with respect to X,Y,Z
	  ! this is a temporary variable
	  ! calculated at each Gauss point in the loop
      real(8) :: dshape_dbigR(nnodes,3)
	  
	  ! global system coordinates R=(X,Y,Z) of the Gauss points
	  ! of this element
      real(8) :: GP_coords(ngp,3)
	  
	  ! global system coordinates R=(X,Y,Z) of the nodes
	  ! of this element
      real(8) :: nodal_coords(nnodes,3)
	  
	  ! Jacobian of the global system coordinate function
	  ! Imagine expressing the real coordinates R=(X,Y,Z)
	  ! as a function of the coordinates in the reference
	  ! element (x,y,z), interpolating using the shape functions.
	  ! Then J_coords is a 3x3 matrix expressing 
	  ! the derivative of R=(X,Y,Z) with respect to r=(x,y,z)
	  ! this Jacobian is a temporary variable
	  ! calculated at each Gauss point in the loop
      real(8) :: J_coords(3,3)
	  
	  ! the inverse: derivatives of the coordinates
	  ! in the reference element r=(x,y,z)
	  ! with respect to the global system coordinates R=(X,Y,Z)
      real(8) :: inv_J_coords(3,3)
	  
	  ! check if matrix is invertible
c      logical :: invert_flag
      real(8) det
	  
	  ! plastic deformation gradient
      ! at the Gauss points read row by row
	  ! order: 11, 12, 13, 21, 22, 23, 31, 32, 33
	  ! first index is the Gauss point
	  ! this is a temporary variable
	  ! calculated at each Gauss point in the loop
      real(8) :: GP_Fp(ngp,9)
	  
	  ! plastic deformation gradient at the nodes
	  ! same ordering as GP_Fp
      real(8) :: nodal_Fp(ngp,9)
	  
	  ! temporary 3x3 matrix to store the following derivatives
	  ! of Fp with respect to the global coordinates R=(X,Y,Z)
	  ! in different parts of the code
	  ! (dFp_xx , dFp_xy, dFp_xz) / (dX, dY, dZ)
	  ! (dFp_yx , dFp_yy, dFp_yz) / (dX, dY, dZ)
	  ! (dFp_zx , dFp_zy, dFp_zz) / (dX, dY, dZ)
      real(8) :: dFp_dR_row(3,3)

      ! initialise node coordinates in the reference element
	  ! this works only for hexahedral elements
      noderef(1,:) = (/-1.,1.,1./)
      noderef(2,:) = (/-1.,-1.,1./)
      noderef(3,:) = (/-1.,1.,-1./)
      noderef(4,:) = (/-1.,-1.,-1./)
      noderef(5,:) = (/1.,1.,1./)
      noderef(6,:) = (/1.,-1.,1./)
      noderef(7,:) = (/1.,1.,-1./)
      noderef(8,:) = (/1.,-1.,-1./)

      gpref = 0.577350269189626 * noderef
	  
      ! assign real coordinates of the Gauss points
	  ! of this element
      do gp = 1,ngp
        do i = 1,3
          GP_coords(gp,i) = global_coords(NOEL,gp,i)
        end do
	  end do
	  
      ! assign row by row Fp for every Gauss point
      do gp = 1,ngp
        do i = 1,3
          do j = 1,3
            GP_Fp(gp,3*(i-1)+j) = global_Fp(NOEL,gp,i,j)		
          end do
        end do
      end do
	  
      ! iterate over Gauss points
	  ! in this first iteration
	  ! shape functions and their derivatives are calculated
	  ! and stored in the reference element system
      do gp = 1,ngp
	  
		! calculate the value of the shape functions at gp
		! this is done by cycling over the shape functions
        ! one for each node
        do node = 1,nnodes
		
		  ! temporary variables for the coordinates
		  ! of the node under consideration 
		  ! in the reference element
          ri = noderef(node,1)
          si = noderef(node,2)
          ti = noderef(node,3)       

          ! gr * gs * gt is the shape function for node "node"
		  ! calculated at the coordinates of the
		  ! Gauss point gp
		  ! everything in the reference element
          gr = 0.5*(1.0 + ri*gpref(i,1))
          dgr = 0.5*ri
      
          gs = 0.5*(1.0 + si*gpref(i,2))
          dgs = 0.5*si       
   
          gt = 0.5*(1.0 + ti*gpref(i,3))
          dgt = 0.5*ti

          shape_functions(gp,i) = gr*gs*gt
		  
          dshape_functions_dr(gp,i,1) = dgr*gs*gt
          dshape_functions_dr(gp,i,2) = gr*dgs*gt
          dshape_functions_dr(gp,i,3) = gr*gs*dgt

        end do ! end iteration over nodes
		
      end do ! end first iteration over integration points
	  
	  ! the matrix inv_shape_functions, multiplied by a vector
	  ! representing a quantity at the Gauss points
	  ! gives the same quantity at the nodes
      call invertnxn(shape_functions,inv_shape_functions,nnodes)
	  
	  ! calculate coordinates of the nodes
	  ! in the global reference frame
      nodal_coords = matmul(inv_shape_functions,GP_coords)

      ! calculate plastic deformation gradient at the nodes
      nodal_Fp = matmul(inv_shape_functions,GP_Fp)	  
	  
      ! iterate over Gauss points
	  ! in this second iteration
	  ! calculations in the global reference frame
	  ! are carried out
      do gp = 1,ngp	

      ! assign dS/dr for this Gauss point
      temp_dshape_functions_dr(1:nnodes,1:3) = 
     + dshape_functions_dr(gp,1:nnodes,1:3)
	  
      ! calculate Jacobian dR/dr
      ! at this Gauss point
      J_coords = matmul(transpose(nodal_coords),
     + temp_dshape_functions_dr)	

      ! invert Jacobian to find dr/dR
	  ! at this Gauss point
      call invert3x3(J_coords, inv_J_coords, det)	 
	  
      if (det.lt.smallnum) then
	  
          
          ! this choice guarantees that output curlFp is zero
	    ! if inversion fails
          dshape_dbigR = 0.0
          write(6,*) "Inverse failure in CalculateCurlFp"
          write(6,*) "Setting zero curlFp"
 
	  
      else
		
         ! calculate the derivatives of the shape functions 
	    ! in the real coordinate system
	    ! calculated at the Gauss points
	    ! dS/dr * dr/dR = dS/dR
        dshape_dbigR = matmul(temp_dshape_functions_dr,inv_J_coords)	 
        
      end if	  
	  
      ! derivatives with respect to the global coordinates R=(x,y,z)
      ! of the first row of the plastic deformation gradient
	  ! dFp/dR = nodal_Fp * dS/dR
	  ! (dFp_xx , dFp_xy, dFp_xz) / (dX, dY, dZ)
      dFp_dR_row = matmul(transpose(nodal_Fp(1:ngp,1:3)),
     + dshape_dbigR)

      ! calculate and store 
	  ! curlFp at this Gauss point
      curlFp(gp,1,1) = dFp_dR_row(3,2) - dFp_dR_row(2,3)
      curlFp(gp,2,1) = dFp_dR_row(1,3) - dFp_dR_row(3,1)
      curlFp(gp,3,1) = dFp_dR_row(2,1) - dFp_dR_row(1,2)

	  ! (dFp_yx , dFp_yy, dFp_yz) / (dX, dY, dZ)
      dFp_dR_row = matmul(transpose(nodal_Fp(1:ngp,4:6)),
     + dshape_dbigR)

      ! calculate and store 
	  ! curlFp at this Gauss point
      curlFp(gp,1,2) = dFp_dR_row(3,2) - dFp_dR_row(2,3)
      curlFp(gp,2,2) = dFp_dR_row(1,3) - dFp_dR_row(3,1)
      curlFp(gp,3,2) = dFp_dR_row(2,1) - dFp_dR_row(1,2) 

	  ! (dFp_zx , dFp_zy, dFp_zz) / (dX, dY, dZ)
      dFp_dR_row = matmul(transpose(nodal_Fp(1:ngp,7:9)),
     + dshape_dbigR)

      ! calculate and store 
	  ! curlFp at this Gauss point
      curlFp(gp,1,3) = dFp_dR_row(3,2) - dFp_dR_row(2,3)
      curlFp(gp,2,3) = dFp_dR_row(1,3) - dFp_dR_row(3,1)
      curlFp(gp,3,3) = dFp_dR_row(2,1) - dFp_dR_row(1,2) 	 
	  
      end do ! end second iteration over integration points
     
      return
      end subroutine CalculateCurlFp

    

	  
	  end module gndcurlFp
