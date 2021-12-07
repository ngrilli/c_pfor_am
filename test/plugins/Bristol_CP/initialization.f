! Chris Allen
! Edward Horton
! Eralp Demir
! Aug. 12th, 2021 - 1st working version
!        
      module initialization
      implicit none
      
      contains
     
      
      
      subroutine initialize_all(str)
      use globalvars, only : GSeffect
      implicit none
      character(len=:), allocatable :: str
	
!
!     If it is the first call for the routine initialize all required variables   
!	This has to be done ONLY ONCE!!!

!     Initializations required for the subroutines in globalsubs.f 
      call read_inputs(str)
      write(6,*) 'read_inputs completed!'
	call allocate_arrays
      write(6,*) 'allocate_arrays completed!'
	call initialize_identity
      write(6,*) 'initialize_identity completed!'
	call initialize_vectors
      write(6,*) 'initialize_vectors completed!'
      call initialize_orientations
      write(6,*) 'initialize_orientations completed!'
	call initialize_elasticity
      write(6,*) 'initialize_elasticity completed!'
	call initialize_statevars
      write(6,*) 'initialize_statevars completed!'      
	call initialize_hardeningmatrix
      write(6,*) 'initialize_hardeningmatrix completed!'
      call initialize_jacobian
      write(6,*) 'initialize_jacobian completed!'

!     Read grain size files if length scale analysis is "ON".
      if (GSeffect.eq.1) then
          
          call read_grainsize(str)
           write(6,*) 'read_grainsize completed!'
          
      endif
      
      
      call initialize_outputfiles(str)
      write(6,*) 'initialize_files completed!'


      
	return
      end subroutine initialize_all
	







      subroutine initialize_jacobian
      use globalvars, only : numel, numip, global_jacob, global_jacob_t
      use globalvars, only : elas66_iso, elas3333, phaseID, global_ori
      use globalsubs, only : transform4, convert3x3x3x3to6x6
      implicit none
      
      integer i,j,k,l
      real(8) elas3333_(3,3,3,3), elas66_(6,6), ori(3,3)
      
      
      write(6,*) 'transformed elasticities for jacobian'
      
      do i=1,numel
          
          do j=1,numip
      
!             If the material is isotropic              
              if (phaseID(i).eq.0) then
              
                  
                  global_jacob(i,j,:,:)=elas66_iso
                  global_jacob_t(i,j,:,:)=elas66_iso
                  
                  
!                  write(6,*) 'el_no',i
!                  write(6,*) 'ip_no',j
!                  do k=1,6
!                      write(6,*) (elas66_iso(k,l),l=1,6)
!                  enddo
                  
                  
                  
                  
              else
                  
      
!                 Get the orientation
                  ori = global_ori(i,j,:,:)
                  
!                 Transform elasticity tensor    
                  call transform4(ori,elas3333,elas3333_)
      
!                 Convert the transformed 4th rank tensor to 6x6 matrix      
                  call convert3x3x3x3to6x6(elas3333_,elas66_)

                  
!                 DO NOT USE DOUBLE SHEAR TERMS - WRONG!
!                 Double the shear terms
!                  elas66_(4,4)=2.0*elas66_(4,4)
!                  elas66_(5,5)=2.0*elas66_(5,5)
!                  elas66_(6,6)=2.0*elas66_(6,6)
                  
                  
                  global_jacob(i,j,:,:)=elas66_
                  global_jacob_t(i,j,:,:)=elas66_
                  
!                  write(6,*) 'el_no',i
!                  write(6,*) 'ip_no',j
!                  do k=1,6
!                      write(6,*) (elas66_(k,l),l=1,6)
!                  enddo
                  
                  
                  
      
              endif
              
      
          enddo
          
          
      enddo
          
              
      return
      end subroutine initialize_jacobian




!	This subroutine assigns identity tensors
!	USES: I3(3,3),I6(6,6),I9(9,9),eijk(3,3,3)
	subroutine initialize_identity
	use globalvars, only : I3,I6,I9,eijk
	implicit none
	integer i
	I3=0.d0
      I6=0.d0
      I9=0.d0
!	Identity matrix (3x3)
	do i=1,3
         I3(i,i)=+1.d0
      enddo 
!	Identity matrix (6x6)
      do i=1,6
         I6(i,i)=+1.d0
      enddo 
!	Identity matrix (9x9)
      do i=1,9
         I9(i,i)=+1.d0
      enddo 
      eijk=0.d0
      eijk(1,2,3)=+1.d0
      eijk(2,3,1)=+1.d0
      eijk(3,1,2)=+1.d0
      eijk(3,2,1)=-1.d0							
      eijk(2,1,3)=-1.d0
      eijk(1,3,2)=-1.d0
	return
	end subroutine initialize_identity
!
!	This subroutine normalizes the slip vectors and forms the line directions: l = n x b
!	USES: n_slip(12,3), b_slip(12,3), eijk(3,3,3)
!	OUTPUTS: l_slip(12,3), Schmid(12,3,3), SchxSch(12,3,3,3,3)
	subroutine initialize_vectors
      use globalvars, only: n_slip, b_slip, l_slip, Schmid
      use globalvars, only: Climb, SchmidT
      use globalvars, only: SchxSch, Schmid_vec, Climb_vec 
      use globalvars, only: eijk, numslip, mattyp
	implicit none
	integer i,j,k,l,m



      allocate (n_slip(numslip,3))
      allocate (b_slip(numslip,3))
      allocate (l_slip(numslip,3))
      
      allocate (Schmid(numslip,3,3))
      allocate (Climb(numslip,3,3))
      allocate (SchmidT(numslip,3,3))
      allocate (SchxSch(numslip,3,3,3,3))
      allocate (Schmid_vec(numslip,6))
      allocate (Climb_vec(numslip,6))


!     For FCC materials
!     mattyp = 1
      if (mattyp.eq.1) then

          
!         Schmid-Boas Notation:
! -----------------------------
!          B2
!         -B4
!          B6
!         -C1
!          C3
!         -C6
!         -A2
!         -A3
!          A5
!          D1
!          D4
!         -D5
! -----------------------------       
          
!         Slip plane normals
	    n_slip(1,:) = (/1., 1., 1./)
	    n_slip(2,:) = (/1., 1., 1./)
	    n_slip(3,:) = (/1., 1., 1./)
	    n_slip(4,:) = (/-1., -1., 1./)
	    n_slip(5,:) = (/-1., -1., 1./)
	    n_slip(6,:) = (/-1., -1., 1./)
	    n_slip(7,:) = (/1., -1., -1./)
	    n_slip(8,:) = (/1., -1., -1./)
	    n_slip(9,:) = (/1., -1., -1./)
	    n_slip(10,:) = (/-1., 1., -1./)
	    n_slip(11,:) = (/-1., 1., -1./)
	    n_slip(12,:) = (/-1., 1., -1./)
!	
!         Slip directions
	    b_slip(1,:) = (/0.,  1., -1./)
	    b_slip(2,:) = (/-1.,  0.,  1./)
	    b_slip(3,:) = (/1., -1.,  0./)
	    b_slip(4,:) = (/0., -1., -1./)
	    b_slip(5,:) = (/1.,  0.,  1./)
	    b_slip(6,:) = (/-1.,  1.,  0./)
	    b_slip(7,:) = (/0., -1.,  1./)
	    b_slip(8,:) = (/-1.,  0., -1./)
	    b_slip(9,:) = (/1.,  1.,  0./)
	    b_slip(10,:) = (/0.,  1.,  1./)
	    b_slip(11,:) = (/1.,  0., -1./)
	    b_slip(12,:) = (/-1., -1.,  0./)



!         Display slip directions          
          do i=1,numslip
              write(6,*) 'system',i
              write(6,*) 'slip direction', (b_slip(i,j),j=1,3)
          enddo
              

!         Display slip plane normals         
          do i=1,numslip
              write(6,*) 'system',i
              write(6,*) 'slip plane normals', (n_slip(i,j),j=1,3)
          enddo          
              
          
          

	    n_slip=n_slip/dsqrt(3.d+0)
	    b_slip=b_slip/dsqrt(2.d+0)
          
        
          
         
          
          
          
!     For BCC materials
!     mattyp = 2              
      elseif (mattyp.eq.2) then
              
!         Schmid-Boas Notation: 1-12
! -----------------------------
!          A2
!          A3
!          A6
!          B2
!          B4
!          B5
!         -C1
!         -C3
!         -C5
!          D1
!          D4
!          D6
! ----------------------------- 


!         Schmid-Boas Notation: 13-24
! -----------------------------
!          A2' (A)
!         -A3' (T)
!         -A6' (T)
!          B2" (A)
!          B4' (A)
!          B5' (A)
!         -C1' (T)
!         -C3" (T)
!          C5" (A)
!         -D1" (T)
!          D4" (A)
!         -D6" (T)
! ----------------------------- 

          

!	    Slip plane normals and slip directions
	    n_slip(1,:) = (/0., -1., 1./)
	    n_slip(2,:) = (/1., 0., 1./)
	    n_slip(3,:) = (/1., 1., 0./)
	    n_slip(4,:) = (/0., -1., 1./)
	    n_slip(5,:) = (/-1., 0., 1./)
	    n_slip(6,:) = (/-1., 1., 0./)
!         corrected ref. Strainer et al., JMPS 50 (2002) 1511-1545          
	    n_slip(7,:) = (/0., 1., 1./)
	    n_slip(8,:) = (/1., 0., 1./)
	    n_slip(9,:) = (/-1., 1., 0./)
!         corrected ref. Strainer et al., JMPS 50 (2002) 1511-1545                   
	    n_slip(10,:) = (/0., 1., 1./)
	    n_slip(11,:) = (/-1., 0., 1./)
	    n_slip(12,:) = (/1., 1., 0./)
      
          if (numslip.gt.12) then
          
              n_slip(13,:) = (/2., 1., 1./)
              n_slip(14,:) = (/-1., -2., 1./)
              n_slip(15,:) = (/-1., 1., -2./)
              n_slip(16,:) = (/-2., 1., 1./)
              n_slip(17,:) = (/1., -2., 1./)
              n_slip(18,:) = (/1., 1., -2./)
              n_slip(19,:) = (/2., -1., 1./)
              n_slip(20,:) = (/-1., 2., 1./)
              n_slip(21,:) = (/-1., -1., -2./)
              n_slip(22,:) = (/-2., -1., 1./)
              n_slip(23,:) = (/1., 2., 1./)
	        n_slip(24,:) = (/1., -1., -2./)     
              
          endif
          
      
!	    Slip direction
	    b_slip(1,:) = (/-1.,  1., 1./)
	    b_slip(2,:) = (/-1.,  1., 1./)
	    b_slip(3,:) = (/-1.,  1., 1./)
	    b_slip(4,:) = (/1.,  1., 1./)
	    b_slip(5,:) = (/1.,  1., 1./)
	    b_slip(6,:) = (/1.,  1., 1./)
	    b_slip(7,:) = (/-1., -1., 1./)
	    b_slip(8,:) = (/-1., -1., 1./)
	    b_slip(9,:) = (/-1., -1., 1./)
	    b_slip(10,:) = (/1., -1.,  1./)
	    b_slip(11,:) = (/1., -1.,  1./)
	    b_slip(12,:) = (/1., -1.,  1./)
      
          
          if (numslip.gt.12) then
          
              b_slip(13,:) = (/-1.,  1., 1./)
              b_slip(14,:) = (/-1.,  1., 1./)
              b_slip(15,:) = (/-1.,  1., 1./)
              b_slip(16,:) = (/1.,  1., 1./)
              b_slip(17,:) = (/1.,  1., 1./)
              b_slip(18,:) = (/1.,  1., 1./)
              b_slip(19,:) = (/-1., -1., 1./)
              b_slip(20,:) = (/-1., -1., 1./)
              b_slip(21,:) = (/-1., -1., 1./)
              b_slip(22,:) = (/1., -1.,  1./)
              b_slip(23,:) = (/1., -1.,  1./)
	        b_slip(24,:) = (/1., -1.,  1./)     
              
          endif
          
!
          

!         Display slip directions          
          do i=1,numslip
              write(6,*) 'system',i
              write(6,*) 'slip direction', (b_slip(i,j),j=1,3)
          enddo
              

!         Display slip plane normals         
          do i=1,numslip
              write(6,*) 'system',i
              write(6,*) 'slip plane normals', (n_slip(i,j),j=1,3)
          enddo     




!         Normalize slip directions          
          b_slip = b_slip / dsqrt(3.0d+0)    
          
!         Normalize slip vectors
          n_slip(1:12,:) = n_slip(1:12,:) / dsqrt(2.0d+0) 
          
          if (numslip.gt.12) then
          
              n_slip(13:24,:) = n_slip(13:24,:) / dsqrt(6.0d+0) 
              
          endif
              
!     For HCP materials
!     mattyp = 3
      elseif (mattyp.eq.3) then    
              
!     NEEDS TO BE FILLED HERE!!!

      endif
      
      
      
!     Schmid tensors      
      do i=1,numslip
          do j=1,3
		    l_slip(i,j)=0.0
		    do k=1,3
			    Schmid(i,j,k)=b_slip(i,j)*n_slip(i,k)
                  Climb(i,j,k)=n_slip(i,j)*n_slip(i,k)
                  SchmidT(i,k,j)=b_slip(i,j)*n_slip(i,k)
			    do l=1,3
      l_slip(i,j)=l_slip(i,j)+(eijk(j,k,l)*b_slip(i,k)*n_slip(i,l))
			    enddo
		    enddo
          enddo
		do j=1,3
		    do k=1,3
			    do l=1,3
				    do m=1,3
					    SchxSch(i,j,k,l,m)=Schmid(i,j,k)*Schmid(i,l,m)
				    enddo
			    enddo
		    enddo
		enddo
		Schmid_vec(i,1)=Schmid(i,1,1)
	    Schmid_vec(i,2)=Schmid(i,2,2)
	    Schmid_vec(i,3)=Schmid(i,3,3)
	    Schmid_vec(i,4)=Schmid(i,1,2)+Schmid(i,2,1)
	    Schmid_vec(i,5)=Schmid(i,3,1)+Schmid(i,1,3)
	    Schmid_vec(i,6)=Schmid(i,2,3)+Schmid(i,3,2)
 
        Climb_vec(i,1)=Climb(i,1,1)
	    Climb_vec(i,2)=Climb(i,2,2)
	    Climb_vec(i,3)=Climb(i,3,3)
	    Climb_vec(i,4)=Climb(i,1,2)+Climb(i,2,1)
	    Climb_vec(i,5)=Climb(i,3,1)+Climb(i,1,3)
	    Climb_vec(i,6)=Climb(i,2,3)+Climb(i,3,2)
      enddo                   
          

	return
	end subroutine initialize_vectors
!
!	This subroutine forms the hardening matrix for the given interaction coefficient
!	USES: sliphard_law, sliphard_param, intmat,numslip
!	OUTPUTS: intmat
	  subroutine initialize_hardeningmatrix
      use globalvars, only: interno, slipint_param, intmat, numslip
      use globalvars, only: mattyp
	  implicit none
      integer i, j, k
      real(8) g0, g1, g2, g3, g4, g5
      
      allocate (intmat(numslip,numslip))

      
!     No latent hardeing (Identity matrix)
      if (interno.eq.0) then
          
          intmat = 0.0d+0
          do i =1,numslip
              intmat(i,i) = 1.
          enddo
      
!     Latent hardening      
      elseif (interno.eq.1) then
          
          
         intmat = slipint_param(1)
         do k = 1, int(numslip/3.)       
	        do i = 1, 3
                  do j = 1, 3
	                intmat(3*(k-1)+i, 3*(k-1)+j)=1.
                  enddo
              enddo
         enddo
         
!     Interaction matrix
      elseif (interno.eq.2) then
         
!         REF.:          
!         FCC - structure
          if (mattyp.eq.1) then
              
              g0 = slipint_param(1)
              g1 = slipint_param(2)
              g2 = slipint_param(3)
              g3 = slipint_param(4)
              g4 = slipint_param(5)
              g5 = slipint_param(6)
              
              
               
              intmat(1,:) = (/g0,g1,g1,g4,g5,g3,g2,g3,g3,g4,g3,g5/)
	        intmat(2,:) = (/g1,g0,g1,g5,g4,g3,g3,g4,g5,g3,g2,g3/)
	        intmat(3,:) = (/g1,g1,g0,g3,g3,g2,g3,g5,g4,g5,g3,g4/)
	        intmat(4,:) = (/g4,g5,g3,g0,g1,g1,g4,g3,g5,g2,g3,g3/)
	        intmat(5,:) = (/g5,g4,g3,g1,g0,g1,g3,g2,g3,g3,g4,g5/)
	        intmat(6,:) = (/g3,g3,g2,g1,g1,g0,g5,g3,g4,g3,g5,g4/)
	        intmat(7,:) = (/g2,g3,g3,g4,g3,g5,g0,g1,g1,g4,g5,g3/)
	        intmat(8,:) = (/g3,g4,g5,g3,g2,g3,g1,g0,g1,g5,g4,g3/)
	        intmat(9,:) = (/g3,g5,g4,g5,g3,g4,g1,g1,g0,g3,g3,g2/)
	        intmat(10,:) = (/g4,g3,g5,g2,g3,g3,g4,g5,g3,g0,g1,g1/)
	        intmat(11,:) = (/g3,g2,g3,g3,g4,g5,g5,g4,g3,g1,g0,g1/)
	        intmat(12,:) = (/g5,g3,g4,g3,g5,g4,g3,g3,g2,g1,g1,g0/)


!         BCC - structure
          elseif (mattyp.eq.2) then
          
!             REF.:               
!             Slip system set: 1-12
              if (numslip.eq.12) then

                intmat(1,:) = (/g0,g1,g1,g2,g2,g2,g3,g3,g3,g3,g3,g3/)
	            intmat(2,:) = (/g1,g0,g1,g3,g3,g3,g2,g2,g2,g3,g3,g3/)
	            intmat(3,:) = (/g1,g1,g0,g3,g3,g3,g3,g3,g3,g2,g2,g2/)
	            intmat(4,:) = (/g2,g2,g2,g0,g1,g1,g3,g3,g3,g3,g3,g3/)
	            intmat(5,:) = (/g3,g3,g3,g1,g0,g1,g3,g3,g3,g2,g2,g2/)
	            intmat(6,:) = (/g3,g3,g3,g1,g1,g0,g2,g2,g2,g3,g3,g3/)
	            intmat(7,:) = (/g4,g4,g4,g4,g4,g4,g0,g1,g1,g5,g5,g5/)
	            intmat(8,:) = (/g2,g2,g2,g3,g3,g3,g1,g0,g1,g3,g3,g3/)
	            intmat(9,:) = (/g3,g3,g3,g2,g2,g2,g1,g1,g0,g3,g3,g3/)
	            intmat(10,:) = (/g4,g4,g4,g4,g4,g4,g5,g5,g5,g0,g1,g1/)
	            intmat(11,:) = (/g3,g3,g3,g2,g2,g2,g3,g3,g3,g1,g0,g1/)
	            intmat(12,:) = (/g2,g2,g2,g3,g3,g3,g3,g3,g3,g1,g1,g0/)

!             REF.: Strainer et al., JMPS 50 (2002) 1511-1545
!             Slip system set: 1-24
              elseif (numslip.eq.24) then
!     
      intmat(1,:)=(/g0,g1,g1,g2,g2,g2,g3,g3,g3,g3,g3,g3,
     + g1,g1,g1,g2,g2,g2,g3,g3,g3,g3,g3,g3/)
!      
      intmat(2,:)=(/g1,g0,g1,g3,g3,g3,g2,g2,g2,g3,g3,g3,
     + g1,g1,g1,g3,g3,g3,g2,g2,g2,g3,g3,g3/)
      
      intmat(3,:)=(/g1,g1,g0,g3,g3,g3,g3,g3,g3,g2,g2,g2,
     + g1,g1,g1,g3,g3,g3,g3,g3,g3,g2,g2,g2/)

      intmat(4,:)=(/g2,g2,g2,g0,g1,g1,g3,g3,g3,g3,g3,g3,
     + g2,g2,g2,g1,g1,g1,g3,g3,g3,g3,g3,g3/)
     
      intmat(5,:)=(/g3,g3,g3,g1,g0,g1,g3,g3,g3,g2,g2,g2,
     + g3,g3,g3,g1,g1,g1,g3,g3,g3,g2,g2,g2/)
      
      intmat(6,:)=(/g3,g3,g3,g1,g1,g0,g2,g2,g2,g3,g3,g3,
     + g3,g3,g3,g1,g1,g1,g2,g2,g2,g3,g3,g3/)
      
      intmat(7,:)=(/g4,g4,g4,g4,g4,g4,g0,g1,g1,g5,g5,g5,
     + g4,g4,g4,g4,g4,g4,g1,g1,g1,g5,g5,g5/)
      
      intmat(8,:)=(/g2,g2,g2,g3,g3,g3,g1,g0,g1,g3,g3,g3,
     + g2,g2,g2,g3,g3,g3,g1,g1,g1,g3,g3,g3/)
      
      intmat(9,:)=(/g3,g3,g3,g2,g2,g2,g1,g1,g0,g3,g3,g3,
     + g3,g3,g3,g2,g2,g2,g1,g1,g1,g3,g3,g3/)
      
      intmat(10,:)=(/g4,g4,g4,g4,g4,g4,g5,g5,g5,g0,g1,g1,
     + g4,g4,g4,g4,g4,g4,g5,g5,g5,g1,g1,g1/)
      
      intmat(11,:)=(/g3,g3,g3,g2,g2,g2,g3,g3,g3,g1,g0,g1,
     + g3,g3,g3,g2,g2,g2,g3,g3,g3,g1,g1,g1/)
      
      intmat(12,:)=(/g2,g2,g2,g3,g3,g3,g3,g3,g3,g1,g1,g0,
     + g2,g2,g2,g3,g3,g3,g3,g3,g3,g1,g1,g0/)
      
      intmat(13,:)=(/g1,g1,g1,g5,g5,g5,g4,g4,g4,g4,g4,g4,
     + g0,g1,g1,g5,g5,g5,g4,g4,g4,g4,g4,g4/)
      
      intmat(14,:)=(/g1,g1,g1,g4,g4,g4,g5,g5,g5,g4,g4,g4,
     + g1,g0,g1,g4,g4,g4,g5,g5,g5,g4,g4,g4/)
      
      intmat(15,:)=(/g1,g1,g1,g4,g4,g4,g4,g4,g4,g5,g5,g5,
     + g1,g1,g0,g4,g4,g4,g4,g4,g4,g5,g5,g5/)
      
      intmat(16,:)=(/g5,g5,g5,g1,g1,g1,g4,g4,g4,g4,g4,g4,
     + g5,g5,g5,g0,g1,g1,g4,g4,g4,g4,g4,g4/)
      
      intmat(17,:)=(/g4,g4,g4,g1,g1,g1,g4,g4,g4,g5,g5,g5,
     + g4,g4,g4,g1,g0,g1,g4,g4,g4,g5,g5,g5/)
      
      intmat(18,:)=(/g4,g4,g4,g1,g1,g1,g5,g5,g5,g4,g4,g4,
     + g4,g4,g4,g1,g1,g0,g5,g5,g5,g4,g4,g4/)
      
      intmat(19,:)=(/g4,g4,g4,g4,g4,g4,g1,g1,g1,g5,g5,g5,
     + g4,g4,g4,g4,g4,g4,g0,g1,g1,g5,g5,g5/)
      
      intmat(20,:)=(/g5,g5,g5,g4,g4,g4,g1,g1,g1,g4,g4,g4,
     + g5,g5,g5,g4,g4,g4,g1,g0,g1,g4,g4,g4/)
      
      intmat(21,:)=(/g4,g4,g4,g5,g5,g5,g1,g1,g1,g4,g4,g4,
     + g4,g4,g4,g5,g5,g5,g1,g1,g0,g4,g4,g4/)
      
      intmat(22,:)=(/g4,g4,g4,g4,g4,g4,g5,g5,g5,g1,g1,g1,
     + g4,g4,g4,g4,g4,g4,g5,g5,g5,g0,g1,g1/)
      
      intmat(23,:)=(/g4,g4,g4,g5,g5,g5,g4,g4,g4,g1,g1,g1,
     + g4,g4,g4,g5,g5,g5,g4,g4,g4,g1,g0,g1/)
      
      intmat(24,:)=(/g5,g5,g5,g4,g4,g4,g4,g4,g4,g1,g1,g1,
     + g5,g5,g5,g4,g4,g4,g4,g4,g4,g1,g1,g0/)
                  
              endif
              
              
         
          
          endif
          
         
      endif
      
       write(6,*) 'intmat'
      do i=1,numslip 
          write(6,*) (intmat(i,j),j=1,numslip)
      enddo
      
          
	return
      end subroutine initialize_hardeningmatrix
      
      
      

!	This subroutine forms the isotropic elasticity tensor
!	USES: C11, C12, C44, I3(3,3)
!	INPUTS: Crystal orientation ori(3,3)
!	OUTPUTS: zeta3333(3,3,3,3), zeta66(6,6)
	subroutine initialize_elasticity
      use globalvars, only: elas3333, elas66, elas66_iso
	  use globalvars, only: elas_param, mattyp
      use globalvars, only: nu, E, G
	implicit none
	integer i, j
	real(8) C11, C12, C44, con




!     IMPORTANT NOTE: THERE MUST BE A MULTIPLIER OF TWO (2.0) IN FRONT OF SHEAR MODULI
!     THIS IS VERY IMPORTANT TO GET CONVERGENCE


!	Isotropic elasticity
!	C44=(C11-C12)/2.0
      
!     FCC material
      if (mattyp.eq.1) then
      
          C11 = elas_param(1)
          C12 = elas_param(2)
          C44 = elas_param(3)



          elas3333=0.0
!         Normal terms      
          elas3333(1,1,1,1)=C11
          elas3333(2,2,2,2)=C11
          elas3333(3,3,3,3)=C11
      
!         Shear terms      
          elas3333(1,2,1,2)=C44
          elas3333(2,1,1,2)=C44
          elas3333(2,1,2,1)=C44
          elas3333(1,2,2,1)=C44
      
          elas3333(2,3,2,3)=C44
          elas3333(3,2,2,3)=C44
          elas3333(3,2,3,2)=C44
          elas3333(2,3,3,2)=C44
      
          elas3333(1,3,1,3)=C44
          elas3333(3,1,1,3)=C44
          elas3333(3,1,3,1)=C44
          elas3333(1,3,3,1)=C44
      
!         Transverse terms
          elas3333(1,1,2,2)=C12
          elas3333(2,2,1,1)=C12
      
          elas3333(1,1,3,3)=C12
          elas3333(3,3,1,1)=C12
      
          elas3333(3,3,2,2)=C12
          elas3333(2,2,3,3)=C12
 
      
      

      
      
!	    Form 6x6 elasticity matrix
	    elas66=0.0d+0
	    do i=1,3
		    do j=1,3
			    if (i.eq.j) then
				    elas66(i,j)=C11
				    elas66(i+3,j+3)=2.*C44
			    else
				    elas66(i,j)=C12
			    endif
              enddo 
          enddo
      


      
          write(6,*) 'SC elasticity'
          do i=1,6
                write(6,*) (elas66(i,j),j=1,6)
          enddo           

      
      
!	    Shear modulus
	    G=dsqrt((C11-C12)*C44)
          
!         Poisson's ratio
!          nu = C11/2./G - 1.
          nu = 0.3
 
!         Homogenized Youngs modulus
          E = 2.*G*(1.+nu)          
          
          
          
          write(6,*) 'Youngs Modulus - iso', E
          write(6,*) 'Poissons Ratio - iso', nu
          write(6,*) 'Shear Modulus - iso', G
          
          
          
          con = E/(1.0+nu)/(1.0-2.0*nu)
      
      
      
      
          elas66_iso = 0.0
      
          elas66_iso(1,1) = 1.0-nu
          elas66_iso(2,2) = 1.0-nu
          elas66_iso(3,3) = 1.0-nu
          elas66_iso(1,2) = nu
          elas66_iso(1,3) = nu
          elas66_iso(2,1) = nu
          elas66_iso(2,3) = nu
          elas66_iso(3,1) = nu
          elas66_iso(3,2) = nu
          elas66_iso(4,4) = (1.0-2.0*nu)/2.0
          elas66_iso(5,5) = (1.0-2.0*nu)/2.0
          elas66_iso(6,6) = (1.0-2.0*nu)/2.0
      
      
      
      
!      write(6,*) 'E', E
!      write(6,*) 'nu', nu
!      write(6,*) 'con', con
      
      
      
          elas66_iso = con * elas66_iso
          
          
          write(6,*) 'isotropic elasticity'
          do i=1,6
                write(6,*) (elas66_iso(i,j),j=1,6)
          enddo          
          
          
          
!     BCC material          
      elseif (mattyp.eq.2) then
          
          
!	    write(*,*) 'Elasticity matrix'
!	    write(*,*) zeta66
      
          C11 = elas_param(1)
          C12 = elas_param(2)
          C44 = elas_param(3)



          elas3333=0.
!         Normal terms      
          elas3333(1,1,1,1)=C11
          elas3333(2,2,2,2)=C11
          elas3333(3,3,3,3)=C11
      
!         Shear terms      
          elas3333(1,2,1,2)=C44
          elas3333(2,1,1,2)=C44
          elas3333(2,1,2,1)=C44
          elas3333(1,2,2,1)=C44
      
          elas3333(2,3,2,3)=C44
          elas3333(3,2,2,3)=C44
          elas3333(3,2,3,2)=C44
          elas3333(2,3,3,2)=C44
      
          elas3333(1,3,1,3)=C44
          elas3333(3,1,1,3)=C44
          elas3333(3,1,3,1)=C44
          elas3333(1,3,3,1)=C44
      
!         Transverse terms
          elas3333(1,1,2,2)=C12
          elas3333(2,2,1,1)=C12
      
          elas3333(1,1,3,3)=C12
          elas3333(3,3,1,1)=C12
      
          elas3333(3,3,2,2)=C12
          elas3333(2,2,3,3)=C12
 
      
      
!	Forming 4th order elasticity tensor
!	do i=1,3
!		do j=1,3
!			do k=1,3
!				do l=1,3
!					dummy=0.0
!					do r=1,3
!						dummy=dummy+(I3(i,r)*I3(j,r)*I3(k,r)*I3(l,r))
!					enddo
!					zeta3333(i,j,k,l)=(C12*I3(i,j)*I3(k,l))+
!     &				(C44*((I3(i,k)*I3(j,l))+(I3(i,l)*I3(j,k))))+
!     &				(dummy*(C11-C12-2.0*C44))		
!				enddo
!			enddo
!		enddo
!	enddo
      
      
      
      
!	    Form 6x6 elasticity matrix
	    elas66=0.d+0
	    do i=1,3
		    do j=1,3
			    if (i.eq.j) then
				    elas66(i,j)=C11
				    elas66(i+3,j+3)=2.0*C44
			    else
				    elas66(i,j)=C12
			    endif
              enddo 
          enddo
      


!	    Shear modulus
	    G=C44
          
!         Poisson's ratio
          nu = C11/2./G - 1.
          
          
          
 
!         Homogenized Youngs modulus
          E = 2.*G*(1.+nu)   

          
          con = E/(1.0+nu)/(1.0-2.0*nu)
      
      
      
      
          elas66_iso = 0.0
      
          elas66_iso(1,1) = 1.0-nu
          elas66_iso(2,2) = 1.0-nu
          elas66_iso(3,3) = 1.0-nu
          elas66_iso(1,2) = nu
          elas66_iso(1,3) = nu
          elas66_iso(2,1) = nu
          elas66_iso(2,3) = nu
          elas66_iso(3,1) = nu
          elas66_iso(3,2) = nu
          elas66_iso(4,4) = (1.0-2.0*nu)/2.0
          elas66_iso(5,5) = (1.0-2.0*nu)/2.0
          elas66_iso(6,6) = (1.0-2.0*nu)/2.0
      
      
      
      
!      write(6,*) 'E', E
!      write(6,*) 'nu', nu
!      write(6,*) 'con', con
      
      
      
          elas66_iso = con * elas66_iso
                    

!     HCP material
      elseif (mattyp.eq.3) then


      endif
          
      

	return
      end subroutine initialize_elasticity

!	This subroutine assigns orientations
!	USES: I3(3,3),I6(6,6),I9(9,9),eijk(3,3,3)
	subroutine initialize_orientations
      use globalvars, only : Euler, phaseID, global_Fp_t
      use globalvars, only : global_Fe_t, numel
      use globalvars, only : numip, I3, global_Fp, global_Fe, global_ori
      use globalsubs, only: ang2ori
	implicit none
	integer el, ip, i, j
      real(8) ori(3,3)

      
      
!     For each element
      do el = 1, numel
                      
 
!             If material is prescribed to be isotropic          
              if (phaseID(el).eq.0) then
              
                  ori = I3
                      
              else
                      
                  
              
                  call ang2ori(Euler(el,1:3),ori)
                      
                      
	                
                      
              endif          
          
          
          
              do ip = 1, numip
          
          

 
                  
!                 Initial deformation gradient                    
                  global_ori(el,ip,:,:) = ori
                  global_Fp(el,ip,:,:) = transpose(ori)
                  global_Fp_t(el,ip,:,:)= transpose(ori)
                  global_Fe(el,ip,:,:) = ori
                  global_Fe_t(el,ip,:,:) = ori
      
      
      
      
              enddo
              
              
              
!              write(6,*) 'el_no',el
!              write(6,*) 'Euler: ', Euler(el,1:3)
!              write(6,*) 'orientation matrix'
!              do i=1,3
!                  write(6,*) (ori(i,j),j=1,3)
!              enddo
              
              
      enddo
      
              
      
      
      return
	end subroutine initialize_orientations



!	This subroutine is written to initialize crystal orientations
!	USES: Euler(:,:), ngrain
	subroutine read_inputs(str)
      use globalvars, only: numel, numip, mattyp, elas_param, numslip
      use globalvars, only: modelno, interno
      use globalvars, only: sliprate_param, dSratio_cr, dgamma_s
      use globalvars, only: sliphard_param, slipint_param
	  use globalvars, only: innoitmax, ounoitmax, thermo
      use globalvars, only: innertol, outertol, njaco, mtdjaco
	  use globalvars, only: deps, temp0, tempdep
      use globalvars, only: Euler, phaseID, dS_cr
	  use globalvars, only: ratio_lb, ratio_ub, numstvar
      use globalvars, only: GSeffect, grainsize_param
	  use globalvars, only: grainID, tstep_forw, tstep_back
      use globalvars, only: numgrain, nodeout, grainori, output_vars
      implicit none	
      integer i, j, dummy, iele, ind
      real(8) dum, dum1, dum2, dum3, phi1, PHI, phi2
      character(len=:), allocatable   :: str
      character(len=17) :: param
      character(len=50) :: line1
      character(len=50) :: line2
      character(len=50) :: line3
      character*8 :: ch


      
      
      
      
      

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(100,file=str // '/inputs.dat',action='read',status='old')
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
!     Read input data for analysis
      
      write(6,*) 'inputs.dat'
      
!     number of elements      
      read(100,*) dum
      numel = int(dum)
      write(6,*) 'numel: ', numel
      
      
!     number of integration points   
      read(100,*) dum
      numip = int(dum)
      write(6,*) 'numip: ', numip
     

!     mechanical (0) or thermomechanical problem (1)
      read(100,*) thermo
      write(6,*) 'thermo: ', thermo
      
      
      
!     temperature dependent properties
!     0: NO / 1: YES
      read(100,*) tempdep
      write(6,*) 'tempdep: ', tempdep

      
      
!     initial temperature [K]
      read(100,*) temp0
      write(6,*) 'temp0: ', temp0      
      
      
!     material type
      read(100,*) dum
      mattyp = int(dum)
      write(6,*) 'mattyp: ', mattyp
      
      
      
      

      
      

!     Sizing of elasticity parameters
      if (mattyp.eq.1) then
          
           allocate (elas_param(3))
           
      elseif (mattyp.eq.2) then
          
           allocate (elas_param(3))   
           
      elseif (mattyp.eq.3) then
          
           allocate (elas_param(5))
          
      endif
      
          

!     Number of slip systems
      read(100,*) dum
      numslip = int(dum)
      write(6,*) 'numslip: ', numslip
             

          

!     Constititive law
!     Only one type of law can be chosen in the model
      read(100,*) dum
      modelno = int(dum)
      write(6,*) 'modelno: ', modelno
         
          
!     slip interaction law
!     Only one type of law can be chosen in the model
      read(100,*) dum
      interno = int(dum)
      write(6,*) 'interno: ', interno
      
      
!     grain size effect
!     Only one type of law can be chosen in the model
      read(100,*) dum
      GSeffect = int(dum)
      write(6,*) 'GSeffect: ', GSeffect    
      
      
!     This part should be expanded to cover different material types      
!     if model = 1 (slip or creep)
      if (modelno.eq.1) then
          allocate (sliprate_param(2))
          
          
                    
          
!     model = 2 (slip and creep law with backstress)          
      elseif (modelno.eq.2) then
          allocate (sliprate_param(4))    
          
     
!     model = 3 (P.R. Dawson's model)                    
      elseif (modelno.eq.3) then
          allocate (sliprate_param(6))           
          
          
          
!     model = 4 (Thermal activation: slip or creep law)                    
      elseif (modelno.eq.4) then
          allocate (sliprate_param(6))           

          
      endif    
   
      
      
      
      
      
!     This part should be expanded to cover different material types      
!     if slip hardening law = 1
!     REF: Kalidindi, S.R., Bronkhorst, C.A. and Anand, L., 1992. 
!     Journal of the Mechanics and Physics of Solids, 40(3), pp.537-569.
      if (modelno.eq.1) then
          
!         Number of state variables          
          numstvar = 1
          
          allocate (sliphard_param(4))
          
          
          
!     if slip hardening law = 2 (Isotropic + Kinematic hardening with Armstrong-Fredrick softening model)
!     REF: Agius, D., Al Mamun, A., Simpson, C.A., Truman, C., Wang, Y., Mostafavi, M. and Knowles, D., 2020. 
!     Computational Materials Science, 183, p.109823.
      elseif (modelno.eq.2) then
          
!         Number of state variables
          numstvar = 2
          
          allocate (sliphard_param(8))
          
          
             
          
          
      endif                   
      
      
      
      
      
         
!     if slip interaction law = 1 (latent hardening)
!     This part should be expanded to cover different material types
      if (interno.eq.1) then
          allocate (slipint_param(1))
          
!     Interaction matrix          
      elseif (interno.eq.2) then    
          allocate (slipint_param(6))
      endif              
          
          
!     Maximum number of iterations at the inner loop
      read(100,*) dum
      innoitmax = int(dum)
      write(6,*) 'innoitmax: ', innoitmax

            
      
      
!     Maximum number of iterations at the outer loop
      read(100,*) dum
      ounoitmax = int(dum)
      write(6,*) 'ounoitmax: ', ounoitmax
       
!     Convergence tolerance for the inner loop
      read(100,*) innertol
      write(6,*) 'innertol: ', innertol

      
!     Convergence tolerance for the outer loop
      read(100,*) outertol            
      write(6,*) 'outertol: ', outertol

!     Ratio of critical resolved shear stress for correction
      read(100,*) dSratio_cr
      write(6,*) 'dSratio_cr: ', dSratio_cr
      
      
      
!     Specified amount of allowed total slip increment for time stepping algorithm
      read(100,*) dgamma_s
      write(6,*) 'dgamma_s: ', dgamma_s          
!     Specified amount of slip increment for time stepping algorithm
!      dgamma_s = 0.01      
      
      
!     Allowed upper bound for slip ratio
      read(100,*) ratio_ub
      write(6,*) 'ratio_ub: ', ratio_ub
      

!     Allowed lower bound for slip ratio
      read(100,*) ratio_lb
      write(6,*) 'ratio_lb: ', ratio_lb      
      
      

!     Fraction of forward time stepping
      read(100,*)  tstep_forw
      write(6,*) 'tstep_forw: ', tstep_forw
      

!     Allowed lower bound for slip ratio
      read(100,*) tstep_back
      write(6,*) 'tstep_back: ', tstep_back            
      
      
           
      
      

!     Frequency of Jacobian calculation 
!     "1": Every time step
      read(100,*) dum
      njaco = int(dum)
      write(6,*) 'njaco: ', njaco

      
      
!     Method of Material Tangent calculation
!     "1": Perturbation
!     "2": Analytical
      read(100,*) dum
      mtdjaco = int(dum)
      write(6,*) 'mtdjaco: ', mtdjaco


!     Strain increment for jacobian calculation
!     If the method is perturbation 
      if (mtdjaco.eq.1) then
          read(100,*) deps
          write(6,*) 'deps: ', deps
      endif
       
          
          
          
          
          

              
!     read elasticity
!     FCC type material
      if (mattyp.eq.1) then
!         C11            
          read(100,*) elas_param(1)
!         C12                  
          read(100,*) elas_param(2)
!         C44
          read(100,*) elas_param(3)
!     BCC type material                  
      elseif (mattyp.eq.2) then       
!         C11            
          read(100,*) elas_param(1)
!         C12                  
          read(100,*) elas_param(2)
!         C44
          read(100,*) elas_param(3)
!     HCP type material
      elseif (mattyp.eq.3) then 
!         Needs to be filled for HCP
      endif
              
       
                  
                  
!     read slip rate parameters
!     Power Law
      if (modelno.eq.1) then
!         Reference slip rate - gammadot0    
          read(100,*) sliprate_param(1)
!         Rate sensitivity exponent - m             
          read(100,*) sliprate_param(2)
!         Factor to set threshold for slip activation - thres
!          read(100,*) sliprate_param(3)
                
!         Other strain rate laws are to be added here!    

!     Power Law + Creep 
      elseif (modelno.eq.2) then
!         Reference slip rate - gammadot0(s) 
          read(100,*) sliprate_param(1)
!         Rate sensitivity exponent - m(s)      
          read(100,*) sliprate_param(2)
!         Reference slip rate - gammadot0(c) 
          read(100,*) sliprate_param(3)
!         Rate sensitivity exponent - m(c)        
          read(100,*) sliprate_param(4)          
!         Factor to set threshold for slip activation - thres
!          read(100,*) sliprate_param(5)


                 
      endif
              
              
              
           
          
!     read strain hardening parameters
!     Voce type hardening law
      if (modelno.eq.1) then
!         initial slip resistance - tauc0                
          read(100,*) sliphard_param(1)
!         hardening rate - h0
          read(100,*) sliphard_param(2)
!         saturation slip resistance - ss
          read(100,*) sliphard_param(3)
!         strain hardening exponent - a
          read(100,*) sliphard_param(4)
           
                  


!     Isotropic + Kinematic hardening with Armstrong-Fredrick softening model
      elseif (modelno.eq.2) then
          
!         isotropic hardening constants          
!         initial slip resistance - tauc0
          read(100,*) sliphard_param(1)
!         hardening rate - h0
          read(100,*) sliphard_param(2)
!         strain hardening exponent - m
          read(100,*) sliphard_param(3)
!      

!         
!         Armstrong-Fredrick softening
!         activation energy for creep - Q
          read(100,*) sliphard_param(4)
!         slip resistance exponent - d
          read(100,*) sliphard_param(5)

!         fitting parameter - A
          read(100,*) sliphard_param(6)
          
!         kinematic hardening constants
!         hardening term - h
          read(100,*) sliphard_param(7)
!         softening term - hD
          read(100,*) sliphard_param(8)          
                  
      endif     
          
          

      
!     read hardening interaction parameters   
!     latent hardening
      if (interno.eq.1) then
!         latent hardening coefficient    
          read(100,*) slipint_param(1)
!
!     interaction matrix coeffcients
      elseif (interno.eq.2) then
!         interaction hardening coefficients
!         g0: self interaction
          read(100,*) slipint_param(1)     
!         g1: co-planar interaction
          read(100,*) slipint_param(2)               
!         g2: cross-slip interaction
          read(100,*) slipint_param(3)    
!         g3: glissile dislocation interaction
          read(100,*) slipint_param(4)            
!         g4: Hirth lock interaction
          read(100,*) slipint_param(5)   
!         g5: Hirth lock interaction
          read(100,*) slipint_param(6)              
      endif
     
      
      
      
      
!     read length scale parameters   
      if (GSeffect.eq.1) then
      
          
!         linear hardening coefficient - "k" [MPa mm^0.5]
          read(100,*) grainsize_param(1)
          
!         grain boundary mismatch exponent - "c" [-]
          read(100,*) grainsize_param(2)          
          
          
      endif
      
      
      

      
      
! ----ED HORTON EDIT ---
!     Outputs to extract: 
!     1: misorientation angle
!     2: cumulative slip
!     3: average of state variables over slip systems
!     4: slip rates per slip system
!     5: state variables per slip system

!     Misorientations
      read(100,*) dum
      output_vars(1) = int(dum)
      write(6,*) 'output_vars(1)', output_vars(1)

!     cumulative slip
      read(100,*) dum
      output_vars(2) = int(dum)
      write(6,*) 'output_vars(2)', output_vars(2)

!     average state variables over slip systems
      read(100,*) dum
      output_vars(3) = int(dum)
      write(6,*) 'output_vars(3)', output_vars(3)

!     slip rates per slip system
      read(100,*) dum
      output_vars(4) = int(dum)
      write(6,*) 'output_vars(4)', output_vars(4)

!     state variables per slip system
      read(100,*) dum
      output_vars(5) = int(dum)
      write(6,*) 'output_vars(5)', output_vars(5)

!---ED HORTON EDIT END ---      
      
      
      
      
      
      
      
      
      
      
      
!     close inputs.dat file     
      close(100)       
      
      
      
      
      
      
 

!     write elasticity
!     FCC type material                  
      if (mattyp.eq.1) then
!         C11            
          write(6,*) 'elas_param(1): ', elas_param(1)
!         C12                  
          write(6,*) 'elas_param(2): ', elas_param(2)
!         C44
          write(6,*) 'elas_param(3): ', elas_param(3)
!     BCC type material               
      elseif (mattyp.eq.2) then       
!         C11            
          write(6,*) 'elas_param(1): ', elas_param(1)
!         C12                  
          write(6,*) 'elas_param(2): ', elas_param(2)
!         C44
          write(6,*) 'elas_param(3): ', elas_param(3)
!     HCP type material
      elseif (mattyp.eq.3) then 
!         Needs to be filled for HCP
      endif
          
          
!     write slip rate parameters
!     Power Law
      if (modelno.eq.1) then
!         Reference slip rate                
          write(6,*) 'sliprate_param(1): ',(sliprate_param(1))
!         Rate sensitivity exponent              
          write(6,*) 'sliprate_param(2): ',(sliprate_param(2))
!         Factor to set threshold for slip activation    
!          write(6,*) 'sliprate_param(3): ',
!     &                (sliprate_param(3))
!         Other slip rate laws are to be added here!                  
                 

      elseif (modelno.eq.2) then
!         Reference slip rate                
          write(6,*) 'sliprate_param(1): ',(sliprate_param(1))
!         Rate sensitivity exponent for slip            
          write(6,*) 'sliprate_param(2): ',(sliprate_param(2))
!         Reference creep rate                
          write(6,*) 'sliprate_param(3): ',(sliprate_param(3))
!         Rate sensitivity exponent for creep
          write(6,*) 'sliprate_param(4): ',(sliprate_param(4))          
!         Factor to set threshold for slip activation    
!          write(6,*) 'sliprate_param(5): ',
!     &                (sliprate_param(5))
!         Other slip rate laws are to be added here!  

      endif          
          
          
          
!     write slip hard parameters
!     Voce type hardening law
      if (modelno.eq.1) then
!         initial slip resistance                
          write(6,*) 'sliphard_param(1): ',sliphard_param(1)
!         hardening rate             
          write(6,*) 'sliphard_param(2): ',sliphard_param(2)
!         saturation slip resistance    
          write(6,*) 'sliphard_param(3): ',sliphard_param(3)
!         slip hardening exponent    
          write(6,*) 'sliphard_param(4): ',sliphard_param(4)
      
                  
!         Other slip rate laws are to be added here!  


      elseif (modelno.eq.2) then
!         initial slip resistance - tauc0
          write(6,*) 'sliphard_param(1): ',sliphard_param(1)
!         hardening rate - h0
          write(6,*) 'sliphard_param(2): ',sliphard_param(2)
!         slip hardening exponent - m
          write(6,*) 'sliphard_param(3): ',sliphard_param(3)
!         activation energy for creep - Q
          write(6,*) 'sliphard_param(4): ',sliphard_param(4)             
!         softening exponent for creep - d          
          write(6,*) 'sliphard_param(5): ',sliphard_param(5)
!         fitting parameter for creep - AF
          write(6,*) 'sliphard_param(6): ',sliphard_param(6)
!         backstress coefficient - h 
          write(6,*) 'sliphard_param(7): ',sliphard_param(7)
!         backstress coefficient - hD    
          write(6,*) 'sliphard_param(8): ',sliphard_param(8)                     
!         Other slip rate laws are to be added here!  

                  
      endif           
          

      
      
  
      
            
      
      
      
      
      
      
      
!     Slip interactions      
      if (interno.eq.0) then
          
          write(6,*) 'no slip interactions'
          
      elseif (interno.eq.1) then
          
      
          write(6,*) 'latent hardening matrix'
      
!         latent hardening coefficient    
          write(6,*) 'slipint_param(1): ',slipint_param(1)                 
          
          
          
      elseif (interno.eq.2) then
          
          
          write(6,*) 'interaction matrix'
          
          
!         latent hardening coefficient    
          write(6,*) 'slipint_param(1): ',slipint_param(1)    
          write(6,*) 'slipint_param(2): ',slipint_param(2)  
          write(6,*) 'slipint_param(3): ',slipint_param(3)  
          write(6,*) 'slipint_param(4): ',slipint_param(4)  
          write(6,*) 'slipint_param(5): ',slipint_param(5)    
          
      endif
      
 
      
!     write length scale parameters   
      if (GSeffect.eq.1) then
      
          
!         linear hardening coefficient - "k" [MPa mm^0.5]
          write(6,*) 'grainsize_param(1)', grainsize_param(1)
          
!         grain boundary mismatch exponent - "c" [-]
          write(6,*) 'grainsize_param(2)', grainsize_param(2)          
          
          
      endif    
      
          

!     Threshold value for stress update algorithm
      dS_cr = dSratio_cr * sliphard_param(1)
      
      
      
      

      
    
!     Length scale model parameters
!     Read the value for "nodeout" from param_array.inc file
      if (GSeffect.eq.1) then
          
!         Read param_array.inc file    
      open(150,file=str // '/param_array.inc',
     + action='read',status='old')

          read(150,*) param, line1
          read(150,*) param, line2
          read(150,*) param, line3
          
          close(150)         
          
          
          


          ind=index(line3, ')')



!          write(6,*) 'ind: ', ind
          
!          write(6,*) 'line1: ', line1
          
!          write(6,*) 'line2: ', line2
          
!          write(6,*) 'line3: ', line3
         
          ch = line3(10:ind-1)
          
!          write(6,*) 'ch: ', ch
          
          
          
          read(ch,*) nodeout  
        
          
          

          write(6,*) 'nodeout', nodeout     

          

          
          
      endif      
      
      
      
      
      
      
      
      
      
      

      
      
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(200,file=str // '/materials.dat',action='read',status='old')
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
!     Read Euler angles, grain ID,  material ID
      allocate (Euler(numel,3))
      allocate (phaseID(numel))
      allocate (grainID(numel))
      

      
      do i=1,numel        
!         dummy first read - Euler angles
          read(200,*) dum1, phi1, PHI, phi2, dum2, dum3
          iele = int(dum1)
          Euler(iele,1) = phi1
          Euler(iele,2) = PHI
          Euler(iele,3) = phi2
          grainID(iele) = int(dum2)
          phaseID(iele) = int(dum3)          
!            call flush(6)
      enddo
      
      close(200) 
      
      
      
      
!     Process Euler angles to find the length scale parameters
!     "numgrain": former "totalfeat"
!     "grainori": former "oriensh"
      
      numgrain = maxval(grainID(:))
      
!     Allocate the array that stores the Euler angles of the grains      
      allocate (grainori(numgrain,3))

!     Assign the Euler angles of the grains
!     Angles are in degrees
      do i=1,numel
          
          grainori(grainID(i),:)=Euler(i,:)
          
      enddo
      
      
      
!     Number of grains in the mesh 
      write(6,*) 'numgrain', numgrain
      
      
!     Write Euler angles of the grains      
      write(6,*) 'grain no', 'Euler angles (deg)' 
      do i=1,numgrain
      
      write(6,*) i, ' ', grainori(i,:)
          
      enddo
      
      
      
      
      
!      write(6,*) 'materials.dat'
!     Output grain orientations and weights
!      do i=1,numel
!          write(6,*) 'Element no.: ', i
!          write(6,*)'Euler angles (deg): ', (Euler(i,j),j=1,3)
!          write(6,*)'Grain ID: ', (grainID(i))
!          write(6,*)'Phase ID: ', (phaseID(i))
!          call flush(6)
!      enddo     
      
      
      
	return 
	end subroutine read_inputs
      











      subroutine read_grainsize(str)
      use globalvars, only : nodex, nodey, nodez, boundgrain
      use globalvars, only : elcent, numel, numgrain, grainori, nodeout
      implicit none


      character(len=:), allocatable :: str
      integer*8 :: nr, nc, i
      real(8) :: numrowval,numcolval


      real(8), allocatable :: rowdata(:)

      
      
!     allocate arrays
      allocate (boundgrain(numgrain,nodeout))
      allocate (elcent(numel,3))
      allocate (nodex(numgrain,nodeout))
      allocate (nodey(numgrain,nodeout))
      allocate (nodez(numgrain,nodeout))
      
      



      
      
      ! Read in bound feature array
      open(300,file=str // '/boundfeat.bin',
     + form='unformatted',status='old',access='stream')
      
      read(300) numrowval
      nr=int(numrowval)
     
      read(300) numcolval
      nc=int(numcolval)
     
      allocate (rowdata(nc))
 
      do i=1,numgrain
        	read(300) rowdata
        	boundgrain(i,:)=rowdata(:)
      end do
	
      deallocate(rowdata)

      
      
      
      
      
      !Read in element centroid coordinates
      open(400,file=str // '/el_centroid.bin',
     + form='unformatted',status='old',access='stream')
                
                
      read(400) numrowval
      nr=int(numrowval)
     
      read(400) numcolval
      nc=int(numcolval)
     
      allocate (rowdata(nc))
     
      do i=1,numel
          read(400) rowdata
          elcent(i,:)=rowdata(:)
      end do
     
      
      deallocate(rowdata)

      
      
      

      ! reading in the x, y, and z coordinates for nodes on the outside
      ! of the grains
      open(500, file=str // '/xvalues.bin',form='unformatted',
     + status='old',access='stream')
                
      read(500) numrowval
      nr=int(numrowval)

      read(500) numcolval
      nc=int(numcolval)
      
      
  
          
          
      open(600,file=str // '/yvalues.bin',form='unformatted',
     + status='old',access='stream')

      read(600) numrowval
      nr=int(numrowval)
     
      read(600) numcolval
      nc=int(numcolval)          
          
      
      
      
      open(700,file=str // '/zvalues.bin',form='unformatted',
     + status='old',access='stream')

      read(700) numrowval
      nr=int(numrowval)
     
      read(700) numcolval
      nc=int(numcolval)
     
      allocate (rowdata(nc))    


      
      
      
      
      
      
      do i=1,nr
        	read(500) rowdata
        	nodex(i,:)=rowdata(:)
	
        
        	read(600) rowdata
        	nodey(i,:)=rowdata(:)
        
        	read(700) rowdata
        	nodez(i,:)=rowdata(:)
       
      end do
      
      deallocate(rowdata)



	close(300)
	close(400)
	close(500)
	close(600)
	close(700)

      
      
      
      
!      write(6,*) 'intotalfeat', intotalfeat
      
 
!     read Euler angles      
!      include 'orien.inc'
    

      
	!reshape array
!      allocate(oriensh(intotalfeat,3))
      
!	oriensh=reshape(orient, (/intotalfeat,3/), order=(/2,1/))     
      
!      write(6,*) 'oriensh'
!      do i=1,intotalfeat
!          write(6,*) oriensh(i,1:3)    
!      enddo
      
!      write(6,*) 'nodex'
!      do i=1,intotalfeat
!          write(6,*) nodex(i,:)    
!      enddo
      
!      write(6,*) 'nodey'
!      do i=1,intotalfeat
!          write(6,*) nodey(i,:)    
!      enddo
      
      
!      write(6,*) 'nodez'
!      do i=1,intotalfeat
!          write(6,*) nodez(i,:)    
!      enddo
      

      
!      write(6,*) 'elcent'
!      do i=1,totalels
!          write(6,*) elcent(i,:)    
!      enddo
      
      
      
	return 
      end subroutine read_grainsize

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
!
!	This subroutine allocates the variables that have to be stored
!	INPUTS: Element number; el, integration point number; ip
	subroutine allocate_arrays
      use globalvars, only: numel, numip, numslip, global_Fp
      use globalvars, only: global_state,global_Fp_t,sliphard_param
      use globalvars, only: global_Fe,global_Fe_t
      use globalvars, only: global_state_t,global_ori,global_jacob_t
	  use globalvars, only: global_gamma
      use globalvars, only: global_jacob,global_sigma,global_S
	  use globalvars, only: global_S_t,global_gammadot
      use globalvars, only: global_gamma_t,global_gamma_sum
	  use globalvars, only: global_gamma_sum_t,numstvar
      use globalvars, only: global_sigma_t, grainsize_init
	  use globalvars, only: global_state0
	  implicit none
      integer i
!
	allocate (global_Fp(numel,numip,3,3))
	global_Fp=0.0d+0
	allocate (global_Fe(numel,numip,3,3))
	global_Fe=0.0d+0
      do i=1,3
          global_Fp(:,:,i,i)=1.0d+0
          global_Fe(:,:,i,i)=1.0d+0
      enddo
      
      allocate (global_state0(numel,numip,numslip,numstvar))
      allocate (global_state(numel,numip,numslip,numstvar))
      allocate (global_state_t(numel,numip,numslip,numstvar))     



	
	
      
      allocate (global_S(numel,numip,6))
      global_S = 0.0d+0
      allocate (global_S_t(numel,numip,6))
      global_S_t = 0.0d+0
	allocate (global_jacob_t(numel,numip,6,6))
	global_jacob_t=0.0d+0
      allocate (global_jacob(numel,numip,6,6))
	global_jacob=0.0d+0
      allocate (global_sigma(numel,numip,6))
	global_sigma=0.0d+0
      allocate (global_sigma_t(numel,numip,6))
	global_sigma_t=0.0d+0
	allocate (global_ori(numel,numip,3,3))
	global_ori=0.0d+0
      allocate(global_gammadot(numel,numip,numslip))
	global_gammadot=0.0d+0

      
      allocate(grainsize_init(numel,numip))
      grainsize_init=0d+0
      
      
      allocate (global_Fp_t(numel,numip,3,3))
      global_Fp_t=0.0d+0
	allocate (global_Fe_t(numel,numip,3,3))
      global_Fe_t=0.0d+0
      do i=1,3
          global_Fp_t(:,:,i,i)=1.0d+0
          global_Fe_t(:,:,i,i)=1.0d+0
          global_ori(:,:,i,i)=1.0d+0
      enddo
      
      

	      
      allocate(global_gamma(numel,numip,numslip))
      allocate(global_gamma_t(numel,numip,numslip))
	global_gamma=0.0d+0
      global_gamma_t=0.0d+0
      

      allocate(global_gamma_sum(numel,numip))
      allocate(global_gamma_sum_t(numel,numip))
	  global_gamma_sum=0.0d+0
      global_gamma_sum_t=0.0d+0
      
!      allocate(global_R(numel,numip))
!	global_R=1.0d+0      
      
	return
	end subroutine allocate_arrays
      
      
      
!	This initializes the state variables
	subroutine initialize_statevars
      use globalvars, only: modelno, sliphard_param
      use globalvars, only: global_state_t, global_state
      use globalvars, only: global_state0
      implicit none
      integer i





      if (modelno.eq.1) then
	    
          global_state0=sliphard_param(1)
          global_state=sliphard_param(1)
          global_state_t=sliphard_param(1)
          
      elseif (modelno.eq.2) then
          
!         Strain/slip hardening term
          global_state0(:,:,:,1)=sliphard_param(1)
          global_state(:,:,:,1)=sliphard_param(1)
          global_state_t(:,:,:,1)=sliphard_param(1)
          
!         Kinematic hardening term
          global_state0(:,:,:,2)=0.0d+0
          global_state(:,:,:,2)=0.0d+0
          global_state_t(:,:,:,2)=0.0d+0
          
          
          
      endif
      
      
      

      
	return
      end subroutine initialize_statevars   
      
      
      
      
      
      
      

      
!	This initializes the state variables related with the grain size
	subroutine initialize_grainsize(el_no,ip_no,gr_no,coords)
!--------------------------------------------------------------------
!       Include file to read in an array containing
!       all information on grain orientations.
      use lengthscale, only: grainsize
      use globalvars, only: grainsize_param, numslip, b_slip, n_slip
      use globalvars, only: global_state, global_state_t
	  use globalvars, only: sliphard_param, global_state0



	implicit none
      
      
      
	integer el_no, ip_no, gr_no, i
      real(8) coords(3)
      real(8) k, c, kval, L, R, val
      
!      real(8) nodex(totalfeat,nodeout), nodey(totalfeat,nodeout),    
!     &nodez(totalfeat,nodeout), elcent(totalels,3), oriensh(totalfeat,3)



	real(8) lm(numslip), ldistance(numslip), rdistance(numslip)


!	common nodex,nodey,nodez,boundgrain,elcent,intotalfeat

!     Length scale parameters
      k = grainsize_param(1)
      c = grainsize_param(2)

      
      
!      write(6,*) '-------------'
!      write(6,*) 'el_no',el_no
!      write(6,*) 'ip_no',ip_no
!      write(6,*) 'gr_no',gr_no
      
      
	    
!     Run the length scale code for analysis         
      call grainsize(coords,ldistance,rdistance,lm,
     + b_slip,n_slip,numslip,gr_no,el_no)


          
!     Assign the strain/slip hardening term
!      write(6,*) 'global_state'
      do i=1,numslip
          
          

!         L-distance (L in ref.)
          L = ldistance(i)*1.d-3
          
!         R-distance (X in ref.)
          R = rdistance(i)*1.d-3          
          
!         Strength coefficient
          kval = k * (1.-lm(i))**(c)
              
          
!         Value of length scale effect
      val = kval/dsqrt(L)*((R+0.5*L)/
     + dsqrt((R+0.5*L)**2.0 - (0.5*L)**2.0)-1.0)
          

!         val is infinity or NaN
          if ((val.gt.1.0d+4).or.(isnan(val))) then
!              
              global_state0(el_no,ip_no,i,1)=sliphard_param(1)
!          
              global_state(el_no,ip_no,i,1)=sliphard_param(1)
!              
              global_state_t(el_no,ip_no,i,1)=sliphard_param(1)
!          
          else
!         Constants are for MPa.mm units so the distance


!             Micrometers are converted to mm
              global_state0(el_no,ip_no,i,1)=sliphard_param(1) + val

          
          
              global_state(el_no,ip_no,i,1)=sliphard_param(1) + val
      
      
              global_state_t(el_no,ip_no,i,1)=sliphard_param(1) + val
      

          endif
          
          
!          write(6,*) global_state(el_no,ip_no,i,1)
               
          
          
      enddo
      
          
!      write(6,*) 'ldistance'
!      write(6,*) ldistance
!      write(6,*) 'rdistance'
!      write(6,*) rdistance      
!      write(6,*) 'lm'
!      write(6,*) lm
!      
!      write(6,*) '-------------'
	return
	end subroutine initialize_grainsize        





      
!	This subroutine initializes the output .txt files
	subroutine initialize_outputfiles(str)
      use globalvars, only: output_vars, numslip, numstvar
	implicit none
	character(len=:), allocatable :: str
      integer count, iss, isv
      
      
!c     Write the captions of the text file for average stresses
!      open(800,file=str //'/avCauchyStress.txt',action='write',
!     &status='replace')
!c      
!      write(800,*) 'INC    ', 'TIME    ', 'S11     ', 
!     &'S22     ', 'S33     ', 'S12     ', 'S13     ', 
!     &'S23     ', 'SeVM'
!c      
!c
!      close(800)     

!     Write the legend of the output variables (UVARM)
!     Outputs to extract: 
!     1: misorientation angle
!     2: cumulative slip
!     3: average of state variables over slip systems
!     4: slip rates per slip system
!     5: state variables per slip system      
      count = 0d+0
      open(850,file=str //'/UVARM_legend.txt',
     + action='write',status='replace')
      
      if (output_vars(1) .eq. 1d+0) then

          count = count + 1d+0
          
      write(850,'(A6,I2,A28)') 'UVARM-', count, 
     + ': misorientation angle [deg]'
          
      endif
      
      
      if (output_vars(2) .eq. 1d+0) then

          count = count + 1d+0
          
          write(850,'(A6,I2,A12)') 'UVARM-', count, ': total slip'
          
      endif
      
      
      if (output_vars(3) .eq. 1d+0) then

          do isv =1,numstvar
          
              count = count + 1d+0
          
      write(850,'(A6,I2,A25,I2)') 'UVARM-', count,
     + ': average state variable-', isv
              
          enddo    
          
      endif
      
      
       if (output_vars(4) .eq. 1d+0) then

          do iss =1,numslip
          
              count = count + 1d+0
          
      write(850,'(A6,I2,A27,I2)') 'UVARM-',count, 
     + ': slip rate of slip system-', iss
              
          enddo    
          
       endif     
       
       
       
      if (output_vars(5) .eq. 1d+0) then

          do isv =1,numstvar
          
              do iss=1,numslip
              
                  count = count + 1d+0
          
      write(850,'(A6,I2,A17,I1,A20,I2)') 'UVARM-', count, 
     + ': state variable-', isv,' of the slip system-', iss
                  
              enddo
              
              
          enddo    
          
      endif
       
      
      
      
      close(850)
      
      

      
      
      
	return
	end subroutine initialize_outputfiles
      
      end module initialization