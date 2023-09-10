c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c Nicolò Grilli
c
c Aug. 12th, 2021 - 1st working version
c
      module initialization
      implicit none

      contains



      subroutine initialize_all(str)
      use globalvars, only : GSeffect, GNDeffect, resdef
      implicit none
      character(len=:), allocatable :: str

c
c     If it is the first call for the routine initialize all required variables
c	This has to be done ONLY ONCE!!!

c     Initializations required for the subroutines in globalsubs.f
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
      call initialize_SCelasticity
      write(6,*) 'initialize_SCelasticity completed!'
      call initialize_J2elasticity
      write(6,*) 'initialize_isoelasticity completed!'
      call initialize_statevars
      write(6,*) 'initialize_statevars completed!'
      call initialize_hardeningmatrix
      write(6,*) 'initialize_hardeningmatrix completed!'
      call initialize_relativetolerances
      write(6,*) 'initialize_relativetolerances completed!'
      call initialize_jacobian
      write(6,*) 'initialize_jacobian completed!'

c     Read grain size files if length scale analysis is "ON".
      if ((GSeffect.eq.1d+0).or.(GSeffect.eq.2d+0)) then

          call read_grainsize(str)
          write(6,*) 'read_grainsize completed!'

      endif

c     Read grain size files if length scale analysis is "ON".
      if (GSeffect.gt.2d+0) then

          call initialize_grainshape(str)
          write(6,*) 'read_grainshape completed!'

      endif


c     Initialize residual distortions
      if (resdef.eq.1d+0) then
          call initialize_residualdeformation(str)
          write(6,*) 'initialize_residualdeformation completed!'
      endif


      call initialize_outputfiles(str)
      write(6,*) 'initialize_files completed!'


c     Initialize GND related constants and mappings
      if ((GNDeffect.eq.1d+0).or.(GNDeffect.eq.2d+0)) then

          call initialize_gndslipgrad0
          write(6,*) 'initialize_gndslipgrad0 completed!'


      endif


      return
      end subroutine initialize_all
c







      subroutine initialize_jacobian
      use globalvars, only : numel, numip, global_jacob, global_ori,
     & global_jacob_t, phaseind, elas66_iso, elas3333, phaseID
      use globalsubs, only : transform4, convert3x3x3x3to6x6
      implicit none

      integer i,j,k,l, ph_no
      real(8) elas3333_(3,3,3,3), elas66_(6,6), ori(3,3)


      write(6,*) 'transformed elasticities for jacobian'

      do i=1,numel


c         Phase index
        ph_no = phaseind(i)

        do j=1,numip

c             If the material is isotropic
          if (phaseID(i).eq.0) then


            global_jacob(i,j,:,:)=elas66_iso(ph_no,:,:)
            global_jacob_t(i,j,:,:)=elas66_iso(ph_no,:,:)


c                  write(6,*) 'el_no',i
c                  write(6,*) 'ip_no',j
c                  do k=1,6
c                      write(6,*) (elas66_iso(k,l),l=1,6)
c                  enddo



c             If crystalline material (cubic)
          else


c                 Get the orientation
            ori = global_ori(i,j,:,:)

c                 Transform elasticity tensor
c                 Transpose of the orientation matrix: crystal to sample reference transformation is used!!!
            call transform4(ori,elas3333,elas3333_)

c                 Convert the transformed 4th rank tensor to 6x6 matrix
            call convert3x3x3x3to6x6(elas3333_,elas66_)


c                 DO NOT USE DOUBLE SHEAR TERMS - WRONG!
c                 Double the shear terms
!                  elas66_(4,4)=2.0*elas66_(4,4)
!                  elas66_(5,5)=2.0*elas66_(5,5)
!                  elas66_(6,6)=2.0*elas66_(6,6)


            global_jacob(i,j,:,:)=elas66_
            global_jacob_t(i,j,:,:)=elas66_

c                  write(6,*) 'el_no',i
c                  write(6,*) 'ip_no',j
c                  do k=1,6
c                      write(6,*) (elas66_(k,l),l=1,6)
c                  enddo




          endif


        enddo


      enddo


      return
      end subroutine initialize_jacobian




c	This subroutine assigns identity tensors
c	USES: I3(3,3),I6(6,6),I9(9,9),eijk(3,3,3)
	subroutine initialize_identity
	use globalvars, only : I3,I6,I9,eijk,I3333
	implicit none
	integer i,j,k,l
      I3333=0.0d+0
	I3=0.d+0
      I6=0.d+0
      I9=0.d+0
c	Identity matrix (3x3)
	do i=1,3
         I3(i,i)=+1.d0
      enddo
c     4th order identity tensor


      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      I3333(i,j,k,l)=I3(i,k)*I3(j,l)
                  enddo
              enddo
          enddo
      enddo
c
c
c
c	Identity matrix (6x6)
      do i=1,6
         I6(i,i)=+1.d0
      enddo
c	Identity matrix (9x9)
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
c
c	This subroutine normalizes the slip vectors and forms the line directions: l = n x b
c	USES: n_slip(12,3), b_slip(12,3), eijk(3,3,3)
c	OUTPUTS: l_slip(12,3), Schmid(12,3,3), SchxSch(12,3,3,3,3)
      subroutine initialize_vectors

      use globalvars, only: n_slip, b_slip, l_slip,
     & Schmid, Climb, SchmidT, SchxSch, Schmid_vec,
     & Climb_vec, eijk, numslip, phases, BS_dyad,
     & totphase, maxnumslip

      implicit none
      integer i,j,k,l,m,iph
c


      allocate (n_slip(totphase,maxnumslip,3))
      n_slip=0.0d+0
      allocate (b_slip(totphase,maxnumslip,3))
      b_slip=0.0d+0
      allocate (l_slip(totphase,maxnumslip,3))
      l_slip=0.0d+0

      allocate (Schmid(totphase,maxnumslip,3,3))
      Schmid=0.0d+0
      allocate (Climb(totphase,maxnumslip,3,3))
      Climb=0.0d+0
      allocate (SchmidT(totphase,maxnumslip,3,3))
      SchmidT=0.0d+0
      allocate (SchxSch(totphase,maxnumslip,3,3,3,3))
      SchxSch=0.0d+0
      allocate (Schmid_vec(totphase,maxnumslip,6))
      Schmid_vec=0.0d+0
      allocate (Climb_vec(totphase,maxnumslip,6))
      Climb_vec=0.0d+0
      allocate (BS_dyad(totphase,maxnumslip,6,3,3))
      BS_dyad=0.0d+0



      do iph=1,totphase

          write(6,*) 'phase no: ', iph
          write(6,*) 'phase ID: ', phases(iph)

c         For FCC materials
c         mattyp = 1
          if (phases(iph).eq.1d+0) then


c             Schmid-Boas Notation:
c -----------------------------
c              B2
c             -B4
c              B6
c             -C1
c              C3
c             -C6
c             -A2
c             -A3
c              A5
c              D1
c              D4
c             -D5
c -----------------------------

c             Slip plane normals
	        n_slip(iph,1,:) = (/1., 1., 1./)
	        n_slip(iph,2,:) = (/1., 1., 1./)
	        n_slip(iph,3,:) = (/1., 1., 1./)
	        n_slip(iph,4,:) = (/-1., -1., 1./)
	        n_slip(iph,5,:) = (/-1., -1., 1./)
	        n_slip(iph,6,:) = (/-1., -1., 1./)
	        n_slip(iph,7,:) = (/1., -1., -1./)
	        n_slip(iph,8,:) = (/1., -1., -1./)
	        n_slip(iph,9,:) = (/1., -1., -1./)
	        n_slip(iph,10,:) = (/-1., 1., -1./)
	        n_slip(iph,11,:) = (/-1., 1., -1./)
	        n_slip(iph,12,:) = (/-1., 1., -1./)
c
c             Slip directions
	        b_slip(iph,1,:) = (/0.,  1., -1./)
	        b_slip(iph,2,:) = (/-1.,  0.,  1./)
	        b_slip(iph,3,:) = (/1., -1.,  0./)
	        b_slip(iph,4,:) = (/0., -1., -1./)
	        b_slip(iph,5,:) = (/1.,  0.,  1./)
	        b_slip(iph,6,:) = (/-1.,  1.,  0./)
	        b_slip(iph,7,:) = (/0., -1.,  1./)
	        b_slip(iph,8,:) = (/-1.,  0., -1./)
	        b_slip(iph,9,:) = (/1.,  1.,  0./)
	        b_slip(iph,10,:) = (/0.,  1.,  1./)
	        b_slip(iph,11,:) = (/1.,  0., -1./)
	        b_slip(iph,12,:) = (/-1., -1.,  0./)



c             Display slip directions
              do i=1,numslip(iph)
                  write(6,*) 'system',i
                  write(6,*) 'slip direction',
     &            (b_slip(iph,i,j),j=1,3)
              enddo


c             Display slip plane normals
              do i=1,numslip(iph)
                  write(6,*) 'system',i
                  write(6,*) 'slip plane normals',
     &            (n_slip(iph,i,j),j=1,3)
              enddo




	        n_slip(iph,:,:)=n_slip(iph,:,:)/dsqrt(3.d+0)
	        b_slip(iph,:,:)=b_slip(iph,:,:)/dsqrt(2.d+0)







c     For BCC materials
c     mattyp = 2
      elseif (phases(iph).eq.2d+0) then

c             Schmid-Boas Notation: 1-12
c -----------------------------
c              A2
c              A3
c              A6
c              B2
c              B4
c              B5
c             -C1
c             -C3
c             -C5
c              D1
c              D4
c              D6
c -----------------------------


c             Schmid-Boas Notation: 13-24
c -----------------------------
c              A2' (A)
c             -A3' (T)
c             -A6' (T)
c              B2" (A)
c              B4' (A)
c              B5' (A)
c             -C1' (T)
c             -C3" (T)
c              C5" (A)
c             -D1" (T)
c              D4" (A)
c             -D6" (T)
c -----------------------------



c	        Slip plane normals and slip directions  (12 systems)
	        n_slip(iph,1,:) = (/0., -1., 1./)
	        n_slip(iph,2,:) = (/1., 0., 1./)
	        n_slip(iph,3,:) = (/1., 1., 0./)
	        n_slip(iph,4,:) = (/0., -1., 1./)
	        n_slip(iph,5,:) = (/-1., 0., 1./)
	        n_slip(iph,6,:) = (/-1., 1., 0./)
c             corrected ref. Strainer et al., JMPS 50 (2002) 1511-1545
	        n_slip(iph,7,:) = (/0., 1., 1./)
	        n_slip(iph,8,:) = (/1., 0., 1./)
	        n_slip(iph,9,:) = (/-1., 1., 0./)
c             corrected ref. Strainer et al., JMPS 50 (2002) 1511-1545
	        n_slip(iph,10,:) = (/0., 1., 1./)
	        n_slip(iph,11,:) = (/-1., 0., 1./)
	        n_slip(iph,12,:) = (/1., 1., 0./)

c             24 slip systems
              if (numslip(iph).gt.12) then
                  n_slip(iph,13,:) = (/2., 1., 1./)
                  n_slip(iph,14,:) = (/-1., -2., 1./)
                  n_slip(iph,15,:) = (/-1., 1., -2./)
                  n_slip(iph,16,:) = (/-2., 1., 1./)
                  n_slip(iph,17,:) = (/1., -2., 1./)
                  n_slip(iph,18,:) = (/1., 1., -2./)
                  n_slip(iph,19,:) = (/2., -1., 1./)
                  n_slip(iph,20,:) = (/-1., 2., 1./)
                  n_slip(iph,21,:) = (/-1., -1., -2./)
                  n_slip(iph,22,:) = (/-2., -1., 1./)
                  n_slip(iph,23,:) = (/1., 2., 1./)
	            n_slip(iph,24,:) = (/1., -1., -2./)
              endif


c	        Slip direction (12 systems)
	        b_slip(iph,1,:) = (/-1.,  1., 1./)
	        b_slip(iph,2,:) = (/-1.,  1., 1./)
	        b_slip(iph,3,:) = (/-1.,  1., 1./)
	        b_slip(iph,4,:) = (/1.,  1., 1./)
	        b_slip(iph,5,:) = (/1.,  1., 1./)
	        b_slip(iph,6,:) = (/1.,  1., 1./)
	        b_slip(iph,7,:) = (/-1., -1., 1./)
	        b_slip(iph,8,:) = (/-1., -1., 1./)
	        b_slip(iph,9,:) = (/-1., -1., 1./)
	        b_slip(iph,10,:) = (/1., -1.,  1./)
	        b_slip(iph,11,:) = (/1., -1.,  1./)
	        b_slip(iph,12,:) = (/1., -1.,  1./)

c             24 slip systems
              if (numslip(iph).gt.12) then

                  b_slip(iph,13,:) = (/-1.,  1., 1./)
                  b_slip(iph,14,:) = (/-1.,  1., 1./)
                  b_slip(iph,15,:) = (/-1.,  1., 1./)
                  b_slip(iph,16,:) = (/1.,  1., 1./)
                  b_slip(iph,17,:) = (/1.,  1., 1./)
                  b_slip(iph,18,:) = (/1.,  1., 1./)
                  b_slip(iph,19,:) = (/-1., -1., 1./)
                  b_slip(iph,20,:) = (/-1., -1., 1./)
                  b_slip(iph,21,:) = (/-1., -1., 1./)
                  b_slip(iph,22,:) = (/1., -1.,  1./)
                  b_slip(iph,23,:) = (/1., -1.,  1./)
	            b_slip(iph,24,:) = (/1., -1.,  1./)

              endif

c


c             Display slip directions
              do i=1,numslip(iph)
                  write(6,*) 'system',i
                  write(6,*) 'slip direction',
     &            (b_slip(iph,i,j),j=1,3)
              enddo


c             Display slip plane normals
              do i=1,numslip(iph)
                  write(6,*) 'system',i
                  write(6,*) 'slip plane normals',
     &            (n_slip(iph,i,j),j=1,3)
              enddo




c             Normalize slip directions
              b_slip(iph,:,:)=b_slip(iph,:,:)/dsqrt(3.0d+0)

c             Normalize slip vectors
              n_slip(iph,1:12,:)=n_slip(iph,1:12,:)/dsqrt(2.0d+0)

              if (numslip(iph).gt.12) then

                  n_slip(iph,13:24,:)=n_slip(iph,13:24,:)/dsqrt(6.0d+0)

              endif

c         For HCP materials
c         phaseID = 3
          elseif (phases(iph).eq.3d+0) then

c             NEEDS TO BE FILLED HERE!!!



c         For BCT materials
c         phaseID = 4
          elseif (phases(iph).eq.4+0) then

c             NEEDS TO BE FILLED HERE!!!

          endif

      enddo



      do iph=1,totphase


c         Schmid tensors
          do i=1,numslip(iph)
              do j=1,3
		        l_slip(iph,i,j)=0.0
		        do k=1,3
			        Schmid(iph,i,j,k)=b_slip(iph,i,j)*n_slip(iph,i,k)
                      Climb(iph,i,j,k)=n_slip(iph,i,j)*n_slip(iph,i,k)
                      SchmidT(iph,i,k,j)=b_slip(iph,i,j)*n_slip(iph,i,k)
			        do l=1,3
				        l_slip(iph,i,j)=l_slip(iph,i,j)+
     &				    (eijk(j,k,l)*b_slip(iph,i,k)*n_slip(iph,i,l))
                      enddo
		        enddo
              enddo

              do j=1,3
                  do k=1,3
c                     Dyads for screw dislocations
                      BS_dyad(iph,i,1,k,j)= l_slip(iph,i,j)*
     &                b_slip(iph,i,k) + b_slip(iph,i,j)*l_slip(iph,i,k)
                      BS_dyad(iph,i,2,k,j)= n_slip(iph,i,j)*
     &                b_slip(iph,i,k) + b_slip(iph,i,j)*n_slip(iph,i,k)
c                     Dyads for edge dislocations
                      BS_dyad(iph,i,3,k,j)= b_slip(iph,i,j)*
     &                b_slip(iph,i,k)
                      BS_dyad(iph,i,4,k,j)= n_slip(iph,i,j)*
     &                n_slip(iph,i,k)
                      BS_dyad(iph,i,5,k,j)= l_slip(iph,i,j)*
     &                l_slip(iph,i,k)
                      BS_dyad(iph,i,6,k,j)= n_slip(iph,i,j)*
     &                b_slip(iph,i,k) + b_slip(iph,i,j)*n_slip(iph,i,k)
                  enddo
              enddo


		    do j=1,3
		        do k=1,3
			        do l=1,3
				        do m=1,3
					        SchxSch(iph,i,j,k,l,m)=
     &                        Schmid(iph,i,j,k)*Schmid(iph,i,l,m)
				        enddo
			        enddo
		        enddo
		    enddo
		    Schmid_vec(iph,i,1)=Schmid(iph,i,1,1)
	        Schmid_vec(iph,i,2)=Schmid(iph,i,2,2)
	        Schmid_vec(iph,i,3)=Schmid(iph,i,3,3)
	        Schmid_vec(iph,i,4)=Schmid(iph,i,1,2)+Schmid(iph,i,2,1)
	        Schmid_vec(iph,i,5)=Schmid(iph,i,3,1)+Schmid(iph,i,1,3)
	        Schmid_vec(iph,i,6)=Schmid(iph,i,2,3)+Schmid(iph,i,3,2)

              Climb_vec(iph,i,1)=Climb(iph,i,1,1)
	        Climb_vec(iph,i,2)=Climb(iph,i,2,2)
	        Climb_vec(iph,i,3)=Climb(iph,i,3,3)
	        Climb_vec(iph,i,4)=Climb(iph,i,1,2)+Climb(iph,i,2,1)
	        Climb_vec(iph,i,5)=Climb(iph,i,3,1)+Climb(iph,i,1,3)
	        Climb_vec(iph,i,6)=Climb(iph,i,2,3)+Climb(iph,i,3,2)


          enddo



      enddo




	return
	end subroutine initialize_vectors
c
c
c	This subroutine forms the hardening matrix for the given interaction coefficient
c	USES: sliphard_law, sliphard_param, intmat,numslip
c	OUTPUTS: intmat
      subroutine initialize_hardeningmatrix
      use globalvars, only: interno, slipint_param, intmat, numslip,
     & phases, totphase, intmat1, intmat2, totphase, maxnumslip
      implicit none
      integer i, j, k, iph
      real(8) g0, g1, g2, g3, g4, g5, g6
c Model with 6 interaction hardening matrix coefficients
      real(8) a, b, ai, bi, omega1, omega2, omegai1, omegai2

      allocate (intmat(totphase,maxnumslip,maxnumslip))
      allocate (intmat1(totphase,maxnumslip,maxnumslip))
      allocate (intmat2(totphase,maxnumslip,maxnumslip))


      do iph=1,totphase

          write(6,*) 'phase no: ', iph
          write(6,*) 'phase ID: ', phases(iph)


c         No latent hardeing (Identity matrix)
          if (interno.eq.0d+0) then

              intmat = 0.0d+0
              do i = 1,numslip(iph)
                  intmat(iph,i,i) = 1.0d+0
              enddo

              write(6,*) 'intmat'
              do i=1,numslip(iph)
                  write(6,*) (intmat(iph,i,j),j=1,numslip(iph))
              enddo

c         Latent hardening
          elseif (interno.eq.1d+0) then


              intmat = slipint_param(iph,1)
              do k = 1, int(numslip(iph)/3.)
	            do i = 1,3
                      do j = 1,3
	                    intmat(iph,3*(k-1)+i,3*(k-1)+j)=1.0d+0
                      enddo
                  enddo
              enddo

              write(6,*) 'intmat'
              do i=1,numslip(iph)
                  write(6,*) (intmat(iph,i,j),j=1,numslip(iph))
              enddo

c         Interaction matrix
          elseif (interno.eq.2d+0) then


              g0 = slipint_param(iph,1)
              g1 = slipint_param(iph,2)
              g2 = slipint_param(iph,3)
              g3 = slipint_param(iph,4)
              g4 = slipint_param(iph,5)
              g5 = slipint_param(iph,6)



c             REF.:
c             FCC - structure
              if (phases(iph).eq.1d+0) then



                  intmat(iph,1,:) = (/g0,g1,g1,g4,g5,g3,
     &            g2,g3,g3,g4,g3,g5/)
                  intmat(iph,2,:) = (/g1,g0,g1,g5,g4,g3,
     &            g3,g4,g5,g3,g2,g3/)
                  intmat(iph,3,:) = (/g1,g1,g0,g3,g3,g2,
     &            g3,g5,g4,g5,g3,g4/)
                  intmat(iph,4,:) = (/g4,g5,g3,g0,g1,g1,
     &            g4,g3,g5,g2,g3,g3/)
                  intmat(iph,5,:) = (/g5,g4,g3,g1,g0,g1,
     &            g3,g2,g3,g3,g4,g5/)
                  intmat(iph,6,:) = (/g3,g3,g2,g1,g1,g0,
     &            g5,g3,g4,g3,g5,g4/)
                  intmat(iph,7,:) = (/g2,g3,g3,g4,g3,g5,
     &            g0,g1,g1,g4,g5,g3/)
                  intmat(iph,8,:) = (/g3,g4,g5,g3,g2,g3,
     &            g1,g0,g1,g5,g4,g3/)
                  intmat(iph,9,:) = (/g3,g5,g4,g5,g3,g4,
     &            g1,g1,g0,g3,g3,g2/)
                  intmat(iph,10,:) = (/g4,g3,g5,g2,g3,g3,
     &            g4,g5,g3,g0,g1,g1/)
                  intmat(iph,11,:) = (/g3,g2,g3,g3,g4,g5,
     &            g5,g4,g3,g1,g0,g1/)
                  intmat(iph,12,:) = (/g5,g3,g4,g3,g5,g4,
     &            g3,g3,g2,g1,g1,g0/)

                  write(6,*) 'intmat'
                  do i=1,numslip(iph)
                      write(6,*) (intmat(iph,i,j),j=1,numslip(iph))
                  enddo


c             BCC - structure
              elseif (phases(iph).eq.2d+0) then

c                 REF.:
c                 Slip system set: 1-12
                  if (numslip(iph).eq.12d+0) then

                      intmat(iph,1,:)=(/g0,g1,g1,g2,g2,g2,
     &                g3,g3,g3,g3,g3,g3/)
	                intmat(iph,2,:)=(/g1,g0,g1,g3,g3,g3,
     &                g2,g2,g2,g3,g3,g3/)
	                intmat(iph,3,:)=(/g1,g1,g0,g3,g3,g3,
     &                g3,g3,g3,g2,g2,g2/)
	                intmat(iph,4,:)=(/g2,g2,g2,g0,g1,g1,
     &                g3,g3,g3,g3,g3,g3/)
	                intmat(iph,5,:)=(/g3,g3,g3,g1,g0,g1,
     &                g3,g3,g3,g2,g2,g2/)
	                intmat(iph,6,:)=(/g3,g3,g3,g1,g1,g0,
     &                g2,g2,g2,g3,g3,g3/)
	                intmat(iph,7,:)=(/g4,g4,g4,g4,g4,g4,
     &                g0,g1,g1,g5,g5,g5/)
	                intmat(iph,8,:)=(/g2,g2,g2,g3,g3,g3,
     &                g1,g0,g1,g3,g3,g3/)
	                intmat(iph,9,:)=(/g3,g3,g3,g2,g2,g2,
     &                g1,g1,g0,g3,g3,g3/)
	                intmat(iph,10,:)=(/g4,g4,g4,g4,g4,g4,
     &                g5,g5,g5,g0,g1,g1/)
	                intmat(iph,11,:)=(/g3,g3,g3,g2,g2,g2,
     &                g3,g3,g3,g1,g0,g1/)
	                intmat(iph,12,:)=(/g2,g2,g2,g3,g3,g3,
     &                g3,g3,g3,g1,g1,g0/)

                      write(6,*) 'intmat'
                      do i=1,numslip(iph)
                          write(6,*) (intmat(iph,i,j),j=1,numslip(iph))
                      enddo


c                 REF.: Strainer et al., JMPS 50 (2002) 1511-1545
c                 Slip system set: 1-24
                  elseif (numslip(iph).eq.24d+0) then
c
                      intmat(iph,1,:)=(/ g0, g1, g1, g2, g2, g2, g3, g3,
     & g3, g3, g3, g3, g1, g1, g1, g2, g2, g2, g3, g3, g3, g3, g3, g3 /)
c
                      intmat(iph,2,:)=(/ g1, g0, g1, g3, g3, g3, g2, g2,
     & g2, g3, g3, g3, g1, g1, g1, g3, g3, g3, g2, g2, g2, g3, g3, g3 /)
c
                      intmat(iph,3,:)=(/ g1, g1, g0, g3, g3, g3, g3, g3,
     & g3, g2, g2, g2, g1, g1, g1, g3, g3, g3, g3, g3, g3, g2, g2, g2 /)
c
                      intmat(iph,4,:)=(/ g2, g2, g2, g0, g1, g1, g3, g3,
     & g3, g3, g3, g2, g2, g2, g1, g1, g1, g3, g3, g3, g3, g3, g3 /)
c
                      intmat(iph,5,:)=(/ g3, g3, g3, g1, g0, g1, g3, g3,
     & g3, g2, g2, g2, g3, g3, g3, g1, g1, g1, g3, g3, g3, g2, g2, g2 /)
c
                      intmat(iph,6,:)=(/ g3, g3, g3, g1, g1, g0, g2, g2,
     & g2, g3, g3, g3, g3, g3, g3, g1, g1, g1, g2, g2, g2, g3, g3, g3 /)
c
                      intmat(iph,7,:)=(/ g4, g4, g4, g4, g4, g4, g0, g1,
     & g1, g5, g5, g5, g4, g4, g4, g4, g4, g4, g1, g1, g1, g5, g5, g5 /)
c
                      intmat(iph,8,:)=(/ g2, g2, g2, g3, g3, g3, g1, g0,
     & g1, g3, g3, g3, g2, g2, g2, g3, g3, g3, g1, g1, g1, g3, g3, g3 /)
c
                      intmat(iph,9,:)=(/ g3, g3, g3, g2, g2, g2, g1, g1,
     & g0, g3, g3, g3, g3, g3, g3, g2, g2, g2, g1, g1, g1, g3, g3, g3 /)
c
                      intmat(iph,10,:)=(/ g4,	g4,	g4,	g4,	g4,	g4,	g5,	g5
     &,g5, g0, g1, g1, g4, g4, g4, g4, g4, g4, g5, g5, g5, g1, g1, g1 /)
c
                      intmat(iph,11,:)=(/ g3,	g3,	g3,	g2,	g2,	g2,	g3,	g3
     &,g3, g1 ,g0, g1, g3, g3, g3, g2, g2, g2, g3, g3, g3, g1, g1, g1 /)
c
                      intmat(iph,12,:)=(/ g2,	g2,	g2,	g3,	g3,	g3,	g3,	g3
     &,g3, g1 ,g1, g0, g2, g2, g2, g3, g3, g3, g3, g3, g3, g1, g1, g0 /)
c
                      intmat(iph,13,:)=(/ g1,	g1,	g1,	g5,	g5,	g5,	g4,	g4
     &,g4, g4, g4, g4, g0, g1, g1, g5, g5, g5, g4, g4, g4, g4, g4, g4 /)
c
                      intmat(iph,14,:)=(/ g1,	g1,	g1,	g4,	g4,	g4,	g5,	g5
     &,g5, g4, g4, g4, g1, g0, g1, g4, g4, g4, g5, g5, g5, g4, g4, g4 /)
c
                      intmat(iph,15,:)=(/ g1,	g1,	g1,	g4,	g4,	g4,	g4,	g4
     &,g4, g5, g5, g5, g1, g1, g0, g4, g4, g4, g4, g4, g4, g5, g5, g5 /)
c
                      intmat(iph,16,:)=(/ g5,	g5,	g5,	g1,	g1,	g1,	g4,	g4
     &,g4, g4, g4, g4, g5, g5, g5, g0, g1, g1, g4, g4, g4, g4, g4, g4 /)
c
                      intmat(iph,17,:)=(/ g4,	g4,	g4,	g1,	g1,	g1,	g4,	g4
     &,g4, g5, g5, g5, g4, g4, g4, g1, g0, g1, g4, g4, g4, g5, g5, g5 /)
c
                      intmat(iph,18,:)=(/ g4,	g4,	g4,	g1,	g1,	g1,	g5,	g5
     &,g5, g4, g4, g4, g4, g4, g4, g1, g1, g0, g5, g5, g5, g4, g4, g4 /)
c
                      intmat(iph,19,:)=(/ g4,	g4,	g4,	g4,	g4,	g4,	g1,	g1
     &,g1, g5, g5, g5, g4, g4, g4, g4, g4, g4, g0, g1, g1, g5, g5, g5 /)
c
                      intmat(iph,20,:)=(/ g5,	g5,	g5,	g4,	g4,	g4,	g1,	g1
     &,g1, g4, g4, g4, g5, g5, g5, g4, g4, g4, g1, g0, g1, g4, g4, g4 /)
c
                      intmat(iph,21,:)=(/ g4,	g4,	g4,	g5,	g5,	g5,	g1,	g1
     &,g1, g4, g4, g4, g4, g4, g4, g5, g5, g5, g1, g1, g0, g4, g4, g4 /)
c
                      intmat(iph,22,:)=(/ g4,	g4,	g4,	g4,	g4,	g4,	g5,	g5
     &,g5, g1, g1, g1, g4, g4, g4, g4, g4, g4, g5, g5, g5, g0, g1, g1 /)
c
                      intmat(iph,23,:)=(/ g4,	g4,	g4,	g5,	g5,	g5,	g4,	g4
     &,g4, g1, g1, g1, g4, g4, g4, g5, g5, g5, g4, g4, g4, g1, g0, g1 /)
c
                      intmat(iph,24,:)=(/ g5,	g5,	g5,	g4,	g4,	g4,	g4,	g4
     &,g4, g1, g1, g1, g5, g5, g5, g4, g4, g4, g4, g4, g4, g1, g1, g0 /)


                      write(6,*) 'intmat'
                      do i=1,numslip(iph)
                          write(6,*) (intmat(iph,i,j),j=1,numslip(iph))
                      enddo

                  endif




              endif





c         Interaction matrix - Code Aster

          elseif (interno.eq.3d+0) then


              g0 = slipint_param(iph,1)
              g1 = slipint_param(iph,2)
              g2 = slipint_param(iph,3)
              g3 = slipint_param(iph,4)
              g4 = slipint_param(iph,5)
              g5 = slipint_param(iph,6)
              g6 = slipint_param(iph,7)

c             REF.:
c             FCC - structure
c             Only defined for FCC tpye of materials
              if (phases(iph).eq.1d+0) then





                  intmat(iph,1,:) = (/g0,g1,g2,g2,g3,g4,
     &            g5,g6,g5,g6,g4,g3/)
	            intmat(iph,2,:) = (/g1,g0,g2,g2,g6,g5,
     &            g4,g3,g4,g3,g5,g6/)
	            intmat(iph,3,:) = (/g2,g2,g0,g1,g5,g6,
     &            g3,g4,g6,g5,g3,g4/)
	            intmat(iph,4,:) = (/g2,g2,g1,g0,g4,g3,
     &            g6,g5,g3,g4,g6,g5/)
	            intmat(iph,5,:) = (/g3,g4,g5,g6,g0,g1,
     &            g2,g2,g6,g5,g4,g3/)
	            intmat(iph,6,:) = (/g6,g5,g4,g3,g1,g0,
     &            g2,g2,g3,g4,g5,g6/)
	            intmat(iph,7,:) = (/g5,g6,g3,g4,g2,g2,
     &            g0,g1,g5,g6,g3,g4/)
	            intmat(iph,8,:) = (/g4,g3,g6,g5,g2,g2,
     &            g1,g0,g4,g3,g6,g5/)
	            intmat(iph,9,:) = (/g5,g6,g4,g3,g4,g3,
     &            g5,g6,g0,g1,g2,g2/)
	            intmat(iph,10,:) = (/g4,g3,g5,g6,g5,g6,
     &            g4,g3,g1,g0,g2,g2/)
	            intmat(iph,11,:) = (/g6,g5,g3,g4,g6,g5,
     &            g3,g4,g2,g2,g0,g1/)
	            intmat(iph,12,:) = (/g3,g4,g6,g5,g3,g4,
     &            g6,g5,g2,g2,g1,g0/)


                  write(6,*) 'intmat'
                  do i=1,numslip(iph)
                      write(6,*) (intmat(iph,i,j),j=1,numslip(iph))
                  enddo


              endif


c         Irradiation defects interaction matrix
          elseif (interno.eq.4d+0) then

              omega1 = slipint_param(iph,1)
              omega2 = slipint_param(iph,2)
              omegai1 = slipint_param(iph,3)
              omegai2 = slipint_param(iph,4)

              a = omega1
              b = 1.0d+0 + (omega1 - omega2)
              ai = omegai1
              bi = 1.0d+0 + (omegai1 - omegai2)

c             FCC - structure
c             Only defined for FCC tpye of materials
              if (phases(iph).eq.1d+0) then

c              do i = 1, numslip
c                  do j = 1, numslip
c                      if (i == j) then
c                          intmat1(i,j) = omega1
c                      else
c                          intmat1(i,j) = omega2
c                      endif
c                  enddo
c              enddo

c              do i = 1, numslip
c                  do j = 1, numslip
c                      if (i == j) then
c                          intmat2(i,j) = omegai1
c                      else
c                          intmat2(i,j) = omegai2
c                      endif
c                  enddo
c              enddo


                  intmat1(iph,1,:) = (/a,b,b,b,b,b,b,b,b,b,b,b/)
	            intmat1(iph,2,:) = (/b,a,b,b,b,b,b,b,b,b,b,b/)
	            intmat1(iph,3,:) = (/b,b,a,b,b,b,b,b,b,b,b,b/)
	            intmat1(iph,4,:) = (/b,b,b,a,b,b,b,b,b,b,b,b/)
	            intmat1(iph,5,:) = (/b,b,b,b,a,b,b,b,b,b,b,b/)
	            intmat1(iph,6,:) = (/b,b,b,b,b,a,b,b,b,b,b,b/)
	            intmat1(iph,7,:) = (/b,b,b,b,b,b,a,b,b,b,b,b/)
	            intmat1(iph,8,:) = (/b,b,b,b,b,b,b,a,b,b,b,b/)
	            intmat1(iph,9,:) = (/b,b,b,b,b,b,b,b,a,b,b,b/)
	            intmat1(iph,10,:) = (/b,b,b,b,b,b,b,b,b,a,b,b/)
	            intmat1(iph,11,:) = (/b,b,b,b,b,b,b,b,b,b,a,b/)
	            intmat1(iph,12,:) = (/b,b,b,b,b,b,b,b,b,b,b,a/)

                  intmat2(iph,1,:) = (/ai,bi,bi,bi,bi,bi,
     &            bi,bi,bi,bi,bi,bi/)
	            intmat2(iph,2,:) = (/bi,ai,bi,bi,bi,bi,
     &            bi,bi,bi,bi,bi,bi/)
	            intmat2(iph,3,:) = (/bi,bi,ai,bi,bi,bi,
     &            bi,bi,bi,bi,bi,bi/)
	            intmat2(iph,4,:) = (/bi,bi,bi,ai,bi,bi,
     &            bi,bi,bi,bi,bi,bi/)
	            intmat2(iph,5,:) = (/bi,bi,bi,bi,ai,bi,
     &             bi,bi,bi,bi,bi,bi/)
	            intmat2(iph,6,:) = (/bi,bi,bi,bi,bi,ai,
     &            bi,bi,bi,bi,bi,bi/)
	            intmat2(iph,7,:) = (/bi,bi,bi,bi,bi,bi,
     &            ai,bi,bi,bi,bi,bi/)
	            intmat2(iph,8,:) = (/bi,bi,bi,bi,bi,bi,
     &            bi,ai,bi,bi,bi,bi/)
	            intmat2(iph,9,:) = (/bi,bi,bi,bi,bi,bi,
     &            bi,bi,ai,bi,bi,bi/)
	            intmat2(iph,10,:) = (/bi,bi,bi,bi,bi,bi,
     &            bi,bi,bi,ai,bi,bi/)
	            intmat2(iph,11,:) = (/bi,bi,bi,bi,bi,bi,
     &            bi,bi,bi,bi,ai,bi/)
	            intmat2(iph,12,:) = (/bi,bi,bi,bi,bi,bi,
     &            bi,bi,bi,bi,bi,ai/)

                  write(6,*) 'intmat1'
                  do i=1,numslip(iph)
                      write(6,*) (intmat1(iph,i,j),j=1,numslip(iph))
                  enddo


                  write(6,*) 'intmat2'
                  do i=1,numslip(iph)
                      write(6,*) (intmat2(iph,i,j),j=1,numslip(iph))
                  enddo



              endif

          endif



      enddo
c     End of phase loop









      return
      end subroutine initialize_hardeningmatrix




c
c	This subroutine initializes the anisotropic elasticity tensor
c	USES: C11, C12, C44, I3(3,3)
c	INPUTS: Crystal orientation ori(3,3)
c	OUTPUTS: zeta3333(3,3,3,3), zeta66(6,6)
      subroutine initialize_SCelasticity
      use globalvars, only: elas3333, elas66,
     & elas_param, phases, totphase, phases
      use globalsubs, only: convert6x6to3x3x3x3
      implicit none
      integer i, j, iph
      real(8) C11, C12, C44
c
c     Allocate arrays
      allocate (elas3333(totphase,3,3,3,3))
      elas3333=0.0d+0
      allocate (elas66(totphase,6,6))
      elas66=0.d+0


c     IMPORTANT NOTE: THERE MUST BE A MULTIPLIER OF TWO (2.0) IN FRONT OF SHEAR MODULI
c     THIS IS VERY IMPORTANT TO GET CONVERGENCE

c     Material type
      do iph=1,totphase

          write(6,*) 'phase no: ', iph
          write(6,*) 'phase ID: ', phases(iph)

c         FCC material
          if (phases(iph).eq.1d+0) then



              C11 = elas_param(iph,1)
              C12 = elas_param(iph,2)
              C44 = elas_param(iph,3)





c	        Form 6x6 elasticity matrix
	        do i=1,3
		        do j=1,3
			        if (i.eq.j) then
				        elas66(iph,i,j)=C11
				        elas66(iph,i+3,j+3)=C44
			        else
				        elas66(iph,i,j)=C12
			        endif
                  enddo
              enddo



c             4th order elasticity
              call convert6x6to3x3x3x3(elas66(iph,:,:),
     &        elas3333(iph,:,:,:,:))


c             Shear corrections after conversion
              elas66(iph,4,4) = 2.0d+0 * elas66(iph,4,4)
              elas66(iph,5,5) = 2.0d+0 * elas66(iph,5,5)
              elas66(iph,6,6) = 2.0d+0 * elas66(iph,6,6)


              write(6,*) 'SC elasticity'
              do i=1,6
                write(6,*) (elas66(iph,i,j),j=1,6)
              enddo





c         BCC material
          elseif (phases(iph).eq.2d+0) then


c	        write(*,*) 'Elasticity matrix'
c	        write(*,*) zeta66

              C11 = elas_param(iph,1)
              C12 = elas_param(iph,2)
              C44 = elas_param(iph,3)





c	        Form 6x6 elasticity matrix
	        do i=1,3
		        do j=1,3
			        if (i.eq.j) then
				        elas66(iph,i,j)=C11
				        elas66(iph,i+3,j+3)=C44
			        else
				        elas66(iph,i,j)=C12
			        endif
                  enddo
              enddo


c             4th order elasticity
              call convert6x6to3x3x3x3(elas66(iph,:,:),
     &        elas3333(iph,:,:,:,:))


c             Shear corrections after conversion
              elas66(iph,4,4) = 2.0d+0 * elas66(iph,4,4)
              elas66(iph,5,5) = 2.0d+0 * elas66(iph,5,5)
              elas66(iph,6,6) = 2.0d+0 * elas66(iph,6,6)


              write(6,*) 'SC elasticity'
              do i=1,6
                write(6,*) (elas66(iph,i,j),j=1,6)
              enddo





c         HCP material
          elseif (phases(iph).eq.3d+0) then

c             FILL IN HERE!


c         BCT material
          elseif (phases(iph).eq.4d+0) then

c             FILL IN HERE!



          endif

c     end the phase loop
      enddo



      return
      end subroutine initialize_SCelasticity
cc
c
c
c
c
c	This subroutine initializes the isotropic elasticity tensor
c	USES: C11, C12, C44, I3(3,3)
c	INPUTS: Crystal orientation ori(3,3)
c	OUTPUTS: zeta3333(3,3,3,3), zeta66(6,6)
      subroutine initialize_J2elasticity
      use globalvars, only: elas66_iso, elas_param,
     & nu, E, G, kappa, elas3333, totphase
      use globalsubs, only: convert6x6to3x3x3x3, angax2ori,
     & transform4, convert3x3x3x3to6x6, invertnxn
      implicit none
      integer i, j, no, iph
      real(8) con, Str(6,6), S(6,6)
      real(8) elas3333tr(3,3,3,3), Ctr(6,6), C(6,6)
      real(8) u, v, w, ax(3), ang, R(3,3), RT(3,3)
c

c     Define number of random orientations
c     Assumed as 100
      no = 100d+0


c     Allocate arrays
      allocate (elas66_iso(totphase,6,6))
      allocate (E(totphase))
      allocate (nu(totphase))
      allocate (G(totphase))
      allocate (kappa(totphase))


      E=0.0d+0
      nu=0.0d+0
      G=0.0d+0


c     Repeat this procedure for all the materials
      do iph=1,totphase






c     IMPORTANT NOTE: THERE MUST BE A MULTIPLIER OF TWO (2.0) IN FRONT OF SHEAR MODULI
c     THIS IS VERY IMPORTANT TO GET CONVERGENCE








c         Randomly generate orientations (100)
          do i = 1, no


c             Randomly generate the axis
              call random_number(u)

              call random_number(v)

              call random_number(w)


              ax(1) = u / dsqrt(u**2.0 + v**2.0 + w**2.0)

              ax(2) = v / dsqrt(u**2.0 + v**2.0 + w**2.0)

              ax(3) = w / dsqrt(u**2.0 + v**2.0 + w**2.0)



c             Randomly generate the angle[0-62.8] degrees
              call random_number(ang)
              ang = 62.8d+0 * ang

c             Transformation matrix from sample to crystal
              call angax2ori(ang,ax,R)


c             Crystal to sample transformation
              RT = transpose(R)


c             Transform elasticity
              call transform4(RT,elas3333(iph,:,:,:,:),elas3333tr)


c             Convert to 6x6 matrix
              call convert3x3x3x3to6x6(elas3333tr,Ctr)

c             Invert to get compliance
              call invertnxn(Ctr,Str,6)

c             Young's modulus
              E(iph) = E(iph) + (Ctr(1,1)+Ctr(2,2)+Ctr(3,3))/3.0d+0

c             Shear modulus
              G(iph) = G(iph) + (Ctr(4,4)+Ctr(5,5)+Ctr(6,6))/3.0d+0

c             Poisson ratio
              nu(iph) = nu(iph) - (Str(1,2)/Str(1,1) + Str(1,3)/Str(1,1)
     &        + Str(2,1)/Str(2,2) + Str(2,3)/Str(2,2)
     &        + Str(3,1)/Str(3,3) + Str(3,2)/Str(3,3))/6.0d+0


          enddo


c         Obtain the average values
          E(iph) = E(iph) / no

          G(iph) = G(iph) / no

          nu(iph) = nu(iph) / no

c         Calculate Bulk modulus
          kappa(iph) = E(iph)/3.0d+0/(1.0d+0-2.0d+0*nu(iph))


          write(6,*) 'Material no: ', iph

          write(6,*) 'Youngs Modulus - iso', E(iph)
          write(6,*) 'Poissons Ratio - iso', nu(iph)
          write(6,*) 'Shear Modulus - iso', G(iph)
          write(6,*) 'Bulk Modulus - iso', kappa(iph)


          S=0.0d+0
          S(1,1) = 1.0d+0/E(iph)
          S(1,2) = -nu(iph)/E(iph)
          S(1,3) = -nu(iph)/E(iph)
          S(2,1) = -nu(iph)/E(iph)
          S(2,2) = 1.0d+0/E(iph)
          S(2,3) = -nu(iph)/E(iph)
          S(3,1) = -nu(iph)/E(iph)
          S(3,2) = -nu(iph)/E(iph)
          S(3,3) = 1.0d+0/E(iph)
          S(4,4) = 1.0d+0/G(iph)
          S(5,5) = 1.0d+0/G(iph)
          S(6,6) = 1.0d+0/G(iph)


          write(6,*) 'isotropic compliance'
          do i=1,6
              write(6,*) (S(i,j),j=1,6)
          enddo

c         Invert to get elasticity
          call invertnxn(S,C,6)

c         Correct the shear terms
          C(4,4) = 2.0d+0 * C(4,4)
          C(5,5) = 2.0d+0 * C(5,5)
          C(6,6) = 2.0d+0 * C(6,6)

          write(6,*) 'isotropic elasticity'
          do i=1,6
              write(6,*) (C(i,j),j=1,6)
          enddo

c         Assign the global variable
          elas66_iso(iph,:,:) = C




      enddo
c     End of material loop


      return
      end subroutine initialize_J2elasticity

cc
c
c
c
c
c
c	This subroutine assigns orientations
c	USES: I3(3,3),I6(6,6),I9(9,9),eijk(3,3,3)
      subroutine initialize_orientations
      use globalvars, only : Euler, phaseID, global_Fp_t,
     & global_Fe_t, numel,
     & numip, I3, global_Fp, global_Fe, global_ori
      use globalsubs, only: ang2ori
      implicit none
      integer el, ip, i, j
      real(8) ori(3,3)



c     For each element
      do el = 1, numel


c             If material is prescribed to be isotropic
              if (phaseID(el).eq.0d+0) then

                  ori = I3

              else



                  call ang2ori(Euler(el,1:3),ori)


                  do ip = 1, numip





c                     Initial deformation gradient
                      global_ori(el,ip,:,:) = ori
                      global_Fp(el,ip,:,:) = ori
                      global_Fp_t(el,ip,:,:)= ori
                      global_Fe(el,ip,:,:) = transpose(ori)
                      global_Fe_t(el,ip,:,:) = transpose(ori)
c
c


                  enddo



c                 write(6,*) 'el_no',el
c                 write(6,*) 'Euler: ', Euler(el,1:3)
c                 write(6,*) 'orientation matrix'
c                 do i=1,3
c                     write(6,*) (ori(i,j),j=1,3)
c                 enddo






              endif






      enddo
c
c
c
c
c
c
      return
      end subroutine initialize_orientations


c
c
c	This subroutine is written to initialize crystal orientations
c	USES: Euler(:,:), ngrain
      subroutine read_inputs(str)
      use globalvars, only: numel, numip, elas_param, numslip,
     &modelno, creepno, interno, sliprate_param, dSratio_cr, dgamma_s,
     &sliphard_param, slipint_param, innoitmax, ounoitmax, thermo,
     &innertol, outertol, njaco, mtdjaco, deps, temp0, tempdep, nnpe,
     &Euler, phaseID, dS_cr, ratio_lb, ratio_ub, numstvar, GNDeffect,
     &GSeffect, grainsize_param, grainID, tstep_forw, tstep_back,
     &numgrain, nodeout, grainori, output_vars, resdef, tres,
     &creep_param, phasefielddamage, totphase, phases, phaseind,
     &maxnumslip, cpl_contribution_to_damage,
     &plastic_work_contribution_to_damage, creepphasefieldflag

	implicit none
      integer i, j, dummy, iele, ind, sa
      real(8), allocatable   ::  da(:)
      real(8) dum, dum1, dum2, dum3, phi1, PHI, phi2
      character(len=:), allocatable   :: str
      character(len=17) :: param
      character(len=50) :: line1
      character(len=50) :: line2
      character(len=50) :: line3
      character*8 :: ch








ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(100,file=str // '/inputs.dat',action='read',status='old')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Read input data for analysis

      write(6,*) 'inputs.dat'

c     number of elements
      read(100,*) dum
      numel = int(dum)
      write(6,*) 'numel: ', numel


c     number of integration points
      read(100,*) dum
      numip = int(dum)
      write(6,*) 'numip: ', numip



c     Decide on the element type based on number of integration points
c     4-node linear tetrahedron - C3D4 (nnpe=4, numip=1)
      if (numip.eq.1d+0) then
          nnpe=4d+0
c     10-node quadratic tetrahedron - C3D10 (nnpe=10, numip=4)
      elseif (numip.eq.4d+0) then
          nnpe=10d+0
c     8-node linear brick - C3D8 (nnpe=8, numip=8)
      elseif (numip.eq.8d+0) then
          nnpe=8d+0
c     20-node quadratic brick - C3D20 (nnpe=20, numip=27)
      elseif (numip.eq.20d+0) then
          nnpe=20d+0
      endif

c     phase field damage model flag
c     0: NO / 1: YES
      read(100,*) dum
      phasefielddamage = int(dum)
      write(6,*) 'phasefielddamage: ', phasefielddamage

      if (phasefielddamage.eq.1d+0) then
        read(100,*) cpl_contribution_to_damage
        write(6,*) 'cpl_contribution_to_damage: ',
     & cpl_contribution_to_damage
        read(100,*) plastic_work_contribution_to_damage
        write(6,*) 'plastic_work_contribution_to_damage: ',
     & plastic_work_contribution_to_damage
      end if

c     creep damage model flag
c     0: NO / 1: YES
      read(100,*) dum
      creepphasefieldflag = int(dum)
      write(6,*) 'creep damage: ', creepphasefieldflag

c     residual stresses are defined: (0) no / (1) yes
      read(100,*) dum
      resdef = int(dum)
      write(6,*) 'resdef: ', resdef


c     total amount of time that residual deformation is applied
c     BC will apply after that time

      if (resdef.eq.1d+0) then
          read(100,*) tres
          write(6,*) 'tres: ', tres
      else
          tres=0.0d+0
      endif



c     mechanical (0) or thermomechanical problem (1)
      read(100,*) dum
      thermo = int(dum)
      write(6,*) 'thermo: ', thermo



c     temperature dependent properties
c     0: NO / 1: YES
      read(100,*) dum
      tempdep = int(dum)
      write(6,*) 'tempdep: ', tempdep





c     initial temperature [K]
      read(100,*) temp0
      write(6,*) 'temp0: ', temp0



c     number of different materials in the mesh
      read(100,*) dum
      totphase = int(dum)
      write(6,*) 'totphase: ', totphase


c     allocate material array
      allocate (phases(totphase))

c     allocate dummy array
      allocate (da(totphase))

c     material type(s)
      read(100,*) (da(i), i=1,totphase)
      phases = int(da)

      write(6,*) 'phases: ', (phases(i),i=1,totphase)






c     Find the array size
c     Based on the phase
c     1:FCC / 2:BCC / 3:HCP / 4:BCT



c     Sizing of elasticity parameters
c     Based on maximum amount of constants (BCT: 6)
      allocate (elas_param(totphase,6))
      elas_param=0.0d+0



c     For each material
c     Allocate numslip array
      allocate (numslip(totphase))

c     Read the number of slip systems
c     Number of slip systems
      read(100,*) (da(i), i=1,totphase)
      numslip = int(da)
      write(6,*) 'numslip: ', (numslip(i),i=1,totphase)


      maxnumslip=maxval(numslip)



c     Constititive law
c     Only one type of law can be chosen in the model
      read(100,*) dum
      modelno = int(dum)
      write(6,*) 'modelno: ', modelno


c     Creep law no
c     Only one type of law can be chosen in the model
      read(100,*) dum
      creepno = int(dum)
      write(6,*) 'creepno: ', creepno


c     slip interaction law
c     Only one type of law can be chosen in the model
      read(100,*) dum
      interno = int(dum)
      write(6,*) 'interno: ', interno


c     grain size effect
c     1: Dylan's original model
c     2: Dylan's modified model
c     3: Grain shapes based on ellipsoid fittings
c     4: Average grain size (uniform for all slip systems)
c     5: Maximum grain size (uniform for all slip systems)
      read(100,*) dum
      GSeffect = int(dum)
      write(6,*) 'GSeffect: ', GSeffect


c     GND calculations
c     Only one type of law can be chosen in the model
c     0: none / 1: by curl of Fp / 2: by slip gradients
      read(100,*) dum
      GNDeffect = int(dum)
      write(6,*) 'GNDeffect: ', GNDeffect


c     This part should be expanded to cover different material types
c     if model = 1 (slip or creep)
      if (modelno.eq.1d+0) then
        allocate (sliprate_param(totphase,2))




c     model = 2 (slip and creep law with backstress)
      elseif (modelno.eq.2d+0) then
        allocate (sliprate_param(totphase,2))


c     model = 3 (Code Aster constitutive model)
      elseif (modelno.eq.3d+0) then
        allocate (sliprate_param(totphase,2))



c     model = 4 (Hugh and Eralp's slip gradient model)
      elseif (modelno.eq.4d+0) then
        allocate (sliprate_param(totphase,2))


c     model = 5 (slip or creep with backstress)
      elseif (modelno.eq.5d+0) then
        allocate (sliprate_param(totphase,2))

c     model = 6 Irradiation hardening addition to model 4
      elseif (modelno.eq.6d+0) then
        allocate (sliprate_param(totphase,3))





      endif


c     This part should be expanded to cover different material types
c     if creep = 1 (Dylan and Abdullah - creep model)
      if (creepno.eq.1d+0) then
        allocate (creep_param(totphase,2))




c     creep = 2 (power law creep, Takeuchi and Argon)
      elseif (creepno.eq.2d+0) then
        allocate (creep_param(totphase,4))


c     creep = 3 (sinh law)
      elseif (creepno.eq.3d+0) then
        allocate (creep_param(totphase,3))



c     creep = 4 (exponential law)
      elseif (creepno.eq.4d+0) then
        allocate (creep_param(totphase,3))


c     creep = 0 (nocreep)
      else

        allocate (creep_param(totphase,1))
        creep_param=0.0d+0

      endif



c     This part should be expanded to cover different material types
c     if slip hardening law = 1
c     REF: Kalidindi, S.R., Bronkhorst, C.A. and Anand, L., 1992.
c     Journal of the Mechanics and Physics of Solids, 40(3), pp.537-569.
      if (modelno.eq.1d+0) then

c         Number of state variables
        numstvar = 1d+0

        allocate (sliphard_param(totphase,4))



c     if slip hardening law = 2 (Isotropic + Kinematic hardening with Armstrong-Fredrick softening model)
c     REF: Agius, D., Al Mamun, A., Simpson, C.A., Truman, C., Wang, Y., Mostafavi, M. and Knowles, D., 2020.
c     Computational Materials Science, 183, p.109823.
      elseif (modelno.eq.2d+0) then

c         Number of state variables
        numstvar = 2d+0

        allocate (sliphard_param(totphase,8))

c     if slip hardening law = 3 (Isotropic + Kinematic hardening model)
c     REF: Helfer, T., Michel, B., Proix, J.M., Salvo, M., Sercombe, J. and Casella, M., 2015.
c     Introducing the open-source mfront code generator: Application to mechanical behaviours and
c     material knowledge management within the PLEIADES fuel element modelling platform.
c     Computers & Mathematics with Applications, 70(5), pp.994-1023.
      elseif (modelno.eq.3d+0) then

c         Number of state variables
        numstvar = 2d+0

        allocate (sliphard_param(totphase,5))



c     Slip gradient model
      elseif (modelno.eq.4d+0) then

c         Number of state variables (SSD/GNDe/GNDs/backstress)
        numstvar = 4d+0

        allocate (sliphard_param(totphase,9))





c     Slip hardening with slip distance - grain size
      elseif (modelno.eq.5d+0) then

c         Number of state variables
        numstvar = 2d+0

        allocate (sliphard_param(totphase,10))

c     Irradiation hardening addition to Slip gradient model
      elseif (modelno.eq.6d+0) then

c 1 - SSD density
c 2 - GND edge dislocation density
c 3 - GND screw dislocation density
c 4 - backstress
c 5 - frank partial dislocation loop density
        numstvar = 5d+0

        allocate (sliphard_param(totphase,12))






      endif


      write(6,*) 'numstvar: ', numstvar



c     if slip interaction law = 1 (latent hardening)
c     This part should be expanded to cover different material types
      if (interno.eq.1d+0) then
        allocate (slipint_param(totphase,1))

c     Interaction matrix - FCC/BCC upto 24 slip systems are available
      elseif (interno.eq.2d+0) then
        allocate (slipint_param(totphase,6))

c     Interaction matrix - (Code Aster) - only defined for FCC material
      elseif (interno.eq.3d+0) then
        allocate (slipint_param(totphase,7))

c     Interaction matrix - Irradiation damage
      elseif (interno.eq.4d+0) then
        allocate (slipint_param(totphase,4))

c     NO interaction matrix- still allocate this array
      else
        allocate (slipint_param(totphase,1))
        slipint_param = 0.0d+0
      endif


c     Maximum number of iterations at the inner loop
      read(100,*) dum
      innoitmax = int(dum)
      write(6,*) 'innoitmax: ', innoitmax




c     Maximum number of iterations at the outer loop
      read(100,*) dum
      ounoitmax = int(dum)
      write(6,*) 'ounoitmax: ', ounoitmax

c     Convergence tolerance for the inner loop
      read(100,*) innertol
      write(6,*) 'innertol: ', innertol


c     Convergence tolerance for the outer loop
      read(100,*) outertol
      write(6,*) 'outertol: ', outertol

c     Ratio of critical resolved shear stress for correction
      read(100,*) dSratio_cr
      write(6,*) 'dSratio_cr: ', dSratio_cr



c     Specified amount of allowed total slip increment for time stepping algorithm
      read(100,*) dgamma_s
      write(6,*) 'dgamma_s: ', dgamma_s
cc     Specified amount of slip increment for time stepping algorithm
c      dgamma_s = 0.01


c     Allowed upper bound for slip ratio
      read(100,*) ratio_ub
      write(6,*) 'ratio_ub: ', ratio_ub


c     Allowed lower bound for slip ratio
      read(100,*) ratio_lb
      write(6,*) 'ratio_lb: ', ratio_lb



c     Fraction of forward time stepping
      read(100,*)  tstep_forw
      write(6,*) 'tstep_forw: ', tstep_forw


c     Allowed lower bound for slip ratio
      read(100,*) tstep_back
      write(6,*) 'tstep_back: ', tstep_back






c     Frequency of Jacobian calculation
c     "1": Every time step
      read(100,*) dum
      njaco = int(dum)
      write(6,*) 'njaco: ', njaco



c     Method of Material Tangent calculation
c     "1": Perturbation
c     "2": Analytical
      read(100,*) dum
      mtdjaco = int(dum)
      write(6,*) 'mtdjaco: ', mtdjaco



c     If the method is perturbation
      if (mtdjaco.eq.1) then
c         Strain increment for jacobian calculation
        read(100,*) deps
        write(6,*) 'deps: ', deps
      endif




c     Elastic constants


c     C11 (FCC or BCC or HCP or BCT)
      read(100,*) (da(i), i=1,totphase)
      elas_param(1:totphase,1) = int(da)
c     C12 (FCC or BCC or HCP or BCT)
      read(100,*) (da(i), i=1,totphase)
      elas_param(1:totphase,2) = int(da)
c     C44 (FCC or BCC or HCP or BCT)
      read(100,*) (da(i), i=1,totphase)
      elas_param(1:totphase,3) = int(da)
c     C33 (HCP or BCT)
      read(100,*) (da(i), i=1,totphase)
      elas_param(1:totphase,4) = int(da)
c     C14 (for HCP) / C13(for BCT)
      read(100,*) (da(i), i=1,totphase)
      elas_param(1:totphase,5) = int(da)
c     C66 (BCT)
      read(100,*) (da(i), i=1,totphase)
      elas_param(1:totphase,6) = int(da)








c     read slip rate parameters
c     Power Law
      if (modelno.eq.1d+0) then

c         Reference slip rate - gammadot0
        read(100,*) (da(i), i=1,totphase)
        sliprate_param(1:totphase,1) = da
c         Rate sensitivity exponent - m
        read(100,*) (da(i), i=1,totphase)
        sliprate_param(1:totphase,2) = da







c     Power Law
      elseif (modelno.eq.2d+0) then

c         Reference slip rate - gammadot0
        read(100,*) (da(i), i=1,totphase)
        sliprate_param(1:totphase,1) = da
c         Rate sensitivity exponent - m
        read(100,*) (da(i), i=1,totphase)
        sliprate_param(1:totphase,2) = da



c     Power Law
      elseif (modelno.eq.3d+0) then

c         Reference slip rate - gammadot0
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,1) = da
c         Rate sensitivity exponent - m
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,2) = da



c     Power Law
      elseif (modelno.eq.4d+0) then


c         Reference slip rate - gammadot0
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,1) = da
c         Rate sensitivity exponent - m
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,2) = da



c     Power Law
      elseif (modelno.eq.5d+0) then


c         Reference slip rate - gammadot0
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,1) = da
c         Rate sensitivity exponent - m
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,2) = da



c     Power Law with irradiation hardening
      elseif (modelno.eq.6d+0) then



c         Reference slip rate - gammadot0
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,1) = da
c         Rate sensitivity exponent - m
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,2) = da

c         statistical coefficient - lambda
      read(100,*) (da(i), i=1,totphase)
      sliprate_param(1:totphase,2) = da



      endif




c     Read creep parameters

c     Abdullah and Dylan
      if (creepno.eq.1d+0) then


c         Reference creep rate - gammadot0
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,1) = da
c         Creep rate sensitivity exponent - m
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,2) = da




c     power law creep (Takeuchi and Argon)
      elseif (creepno.eq.2d+0) then

c         Reference creep rate - A
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,1) = da
c         Burgers vector - b
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,2) = da
c         Diffusion rate - D
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,3) = da
c         Creep power exponent - n
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,4) = da



c     sinh law
      elseif (creepno.eq.3d+0) then

c         Reference creep rate - CD
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,1) = da
c         Creep exponent - n
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,2) = da
c         Activation volume - V
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,3) = da




c     exponential law
      elseif (creepno.eq.4d+0) then

c         Coefficient - B
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,1) = da
c         Creep activation energy - Qc
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,2) = da
c         Creep activation volume - V
      read(100,*) (da(i), i=1,totphase)
      creep_param(1:totphase,3) = da




      endif










c     read strain hardening parameters
c     Voce type hardening law
      if (modelno.eq.1d+0) then


c         initial slip resistance - tauc0
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,1) = da
c         hardening rate - h0
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,2) = da
c         saturation slip resistance - ss
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,3) = da
c         strain hardening exponent - n
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,4) = da





c     Isotropic + Kinematic hardening with Armstrong-Fredrick softening model
      elseif (modelno.eq.2d+0) then

c         isotropic hardening constants
c         initial slip resistance - tauc0
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,1) = da
c         hardening rate - h0
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,2) = da
c         strain hardening exponent - m
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,3) = da

c         Armstrong-Fredrick softening
c         activation energy for creep - Q
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,4) = da
c         slip resistance exponent - d
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,5) = da
c         fitting parameter - A
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,6) = da

c         kinematic hardening constants
c         hardening term - h
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,7) = da
c         softening term - hD
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,8) = da




c     Isotropic + Kinematic hardening (Code Aster)
      elseif (modelno.eq.3d+0) then


c         isotropic hardening constants
c         initial slip resistance - tauc0
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,1) = da
c         hardening rate coefficient - b
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,2) = da
c         strain hardening factor - Q
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,3) = da

c         kinematic hardening constants
c         hardening term - C
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,4) = da
c         softening term - D
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,5) = da




      elseif (modelno.eq.4d+0) then

c         Initial critical strength (MPa)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,1) = da
c         Initial SSD density (1/mm^2)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,2) = da
c         Burgers vector (mm)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,3) = da
c         Geometric factor - alpha for tauc
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,4) = da
c         Mean-free-path hardening - K
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,5) = da
c         Capture radii - yc (mm)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,6) = da
c         Effective radius - R for backstress
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,7) = da
c         Size factor - mesh size (micron to mm)
c         i.e. Dream3D output is in micrometers so, the conversion to mm is 1000
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,8) = da
c         Hardening scale factor
c         Scaling factor for strain hardening
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,9) = da


c     Voce type hardening law + slip distance based hardening
      elseif (modelno.eq.5d+0) then

c         Initial critical strength (MPa)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,1) = da
c         Initial SSD density (1/mm^2)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,2) = da
c         Burgers vector (mm)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,3) = da
c         Mean-free-path hardening - K
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,4) = da
c         Capture radii - yc (mm)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,5) = da
c         Geometric factor - alpha for tauc
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,6) = da
c         kinematic hardening constants
c         hardening term - h
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,7) = da
c         softening term - hD
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,8) = da
c         Hardening scale factor
c         Scaling factor for size-depedent strain hardening
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,9) = da
c         i.e. Dream3D output is in micrometers so, the conversion to mm is 1/1000
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,10) = da


c     Irradiation model
      elseif (modelno.eq.6d+0) then

c         Initial critical strength (MPa)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,1) = da
c         Initial SSD density (1/mm^2)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,2) = da
c         Burgers vector (mm)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,3) = da
c         Geometric factor - alpha for tauc
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,4) = da
c         Mean-free-path hardening - K
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,5) = da
c         Capture radii - yc (mm)
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,6) = da
c         Effective radius - R for backstress
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,7) = da
c         Size factor - mesh size (micron to mm)
c         i.e. Dream3D output is in micrometers so, the conversion to mm is 1000
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,8) = da
c         Hardening scale factor
c         Scaling factor for strain hardening
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,9) = da
c         Annihilation area for frank loop dislocations - AL
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,10) = da
c         Saturation density for frank loop dislocations - rhoSATURATION
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,11) = da
c         Initial frank loop dislocation density
      read(100,*) (da(i), i=1,totphase)
      sliphard_param(1:totphase,12) = da




      endif




c     read hardening interaction parameters
c     latent hardening
      if (interno.eq.1d+0) then

c         latent hardening coefficient
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,1) = da

c
c     interaction matrix coeffcients
      elseif (interno.eq.2d+0) then

c         interaction hardening coefficients
c         g0: self interaction
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,1) = da
c         g1: co-planar interaction
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,2) = da
c         g2: cross-slip interaction
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,3) = da
c         g3: glissile dislocation interaction
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,4) = da
c         g4: Hirth lock interaction
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,5) = da
c         g5: Hirth lock interaction
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,6) = da




c     interaction matrix coeffcients (Code Aster - 7 constants)
      elseif (interno.eq.3d+0) then

c         interaction hardening coefficients
c         g0
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,1) = da
c         g1
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,2) = da
c         g2
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,3) = da
c         g3
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,4) = da
c         g4
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,5) = da
c         g5
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,6) = da
c         g6
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,7) = da


c     interaction matrix coeffcients for irradiation model
      elseif (interno.eq.4d+0) then

c         interaction hardening coefficients
c         omega1 - self hardening coefficient
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,1) = da
c         omega2 - latent hardening coefficient
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,2) = da
c         omegai1 - irradiation self hardening coefficient
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,3) = da
c         omegai2 - irradiation latent hardening coefficient
      read(100,*) (da(i), i=1,totphase)
      slipint_param(1:totphase,4) = da



      endif



c     read length scale parameters
      if ((GSeffect.eq.1d+0).or.(GSeffect.eq.2d+0)) then


      allocate (grainsize_param(totphase,3))



c         linear hardening coefficient - "k" [MPa mm^0.5]
      read(100,*) (da(i), i=1,totphase)
      grainsize_param(1:totphase,1) = da

c         grain boundary mismatch exponent - "C" [-]
      read(100,*) (da(i), i=1,totphase)
      grainsize_param(1:totphase,2) = da

c         i.e. Dream3D output is in micrometers so, the conversion to mm is 1/1000
      read(100,*) (da(i), i=1,totphase)
      grainsize_param(1:totphase,3) = da





      endif






c ----ED HORTON EDIT ---
c     Outputs to extract:
c     1: misorientation angle
c     2: cumulative slip
c     3: average of state variables over slip systems
c     4: slip rates per slip system
c     5: state variables per slip system

c     Misorientations
      read(100,*) dum
      output_vars(1) = int(dum)
      write(6,*) 'output_vars(1)', output_vars(1)

c     cumulative slip
      read(100,*) dum
      output_vars(2) = int(dum)
      write(6,*) 'output_vars(2)', output_vars(2)

c     average state variables over slip systems
      read(100,*) dum
      output_vars(3) = int(dum)
      write(6,*) 'output_vars(3)', output_vars(3)

c     slip rates per slip system
      read(100,*) dum
      output_vars(4) = int(dum)
      write(6,*) 'output_vars(4)', output_vars(4)

c     state variables per slip system
      read(100,*) dum
      output_vars(5) = int(dum)
      write(6,*) 'output_vars(5)', output_vars(5)

c---ED HORTON EDIT END ---











c     close inputs.dat file
      close(100)






      do i=1,totphase

          write(6,*) 'phase no: ', i
          write(6,*) 'phase ID: ', phases(i)

c         write elastic constants

c         C11
          write(6,*) 'elas_param(1): ', elas_param(i,1)
c         C12
          write(6,*) 'elas_param(2): ', elas_param(i,2)
c         C44
          write(6,*) 'elas_param(3): ', elas_param(i,3)
c         C33
          write(6,*) 'elas_param(4): ', elas_param(i,4)
c         C14-C13
          write(6,*) 'elas_param(5): ', elas_param(i,5)
c         C66
          write(6,*) 'elas_param(6): ', elas_param(i,6)

      enddo



c     write slip rate parameters
c     Power Law
      if (modelno.eq.1d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference slip rate
              write(6,*) 'sliprate_param(1): ',
     &                (sliprate_param(i,1))
c             Rate sensitivity exponent
              write(6,*) 'sliprate_param(2): ',
     &                (sliprate_param(i,2))

          enddo


      elseif (modelno.eq.2d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference slip rate
              write(6,*) 'sliprate_param(1): ',
     &                (sliprate_param(i,1))
c             Rate sensitivity exponent for slip
              write(6,*) 'sliprate_param(2): ',
     &                (sliprate_param(i,2))


          enddo


      elseif (modelno.eq.3d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference slip rate
              write(6,*) 'sliprate_param(1): ',
     &                (sliprate_param(i,1))
c             Rate sensitivity exponent for slip
              write(6,*) 'sliprate_param(2): ',
     &                (sliprate_param(i,2))



          enddo


      elseif (modelno.eq.4d+0) then



          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference slip rate
              write(6,*) 'sliprate_param(1): ',
     &                (sliprate_param(i,1))
c             Rate sensitivity exponent for slip
              write(6,*) 'sliprate_param(2): ',
     &                (sliprate_param(i,2))

          enddo


      elseif (modelno.eq.5d+0) then



          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference slip rate
              write(6,*) 'sliprate_param(1): ',
     &                (sliprate_param(i,1))
c             Rate sensitivity exponent for slip
              write(6,*) 'sliprate_param(2): ',
     &                (sliprate_param(i,2))


          enddo


c Irradiation hardening addition to model 4
      elseif (modelno.eq.6d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference slip rate - gdot0
              write(6,*) 'sliprate_param(1): ',
     &                (sliprate_param(i,1))
c             Rate sensitivity exponent for slip - m
              write(6,*) 'sliprate_param(2): ',
     &                (sliprate_param(i,2))
c             Statistical coefficient - lambda
              write(6,*) 'sliprate_param(3): ',
     &                (sliprate_param(i,3))


          enddo

c         Other slip rate laws are to be added here!

      endif




c     write slip rate parameters
c     Power Law
      if (creepno.eq.1d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference slip rate
              write(6,*) 'creep_param(1): ',
     &                (creep_param(i,1))
c             Rate sensitivity exponent
              write(6,*) 'creep_param(2): ',
     &                (creep_param(i,2))


          enddo

      elseif (creepno.eq.2d+0) then

          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference creep rate - A
              write(6,*) 'creep_param(1): ',
     &                (creep_param(i,1))
c             Burgers vector - b
              write(6,*) 'creep_param(2): ',
     &                (creep_param(i,2))
c             Diffusion rate - D
              write(6,*) 'creep_param(3): ',
     &                (creep_param(i,3))
c             Power law exponent - n
              write(6,*) 'creep_param(4): ',
     &                (creep_param(i,4))


          enddo


      elseif (creepno.eq.3d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference creep rate - CD
              write(6,*) 'creep_param(1): ',
     &                (creep_param(i,1))
c             Exponent for creep sinh function - n
              write(6,*) 'creep_param(2): ',
     &                (creep_param(i,2))
c             Activation volume - V
              write(6,*) 'creep_param(3): ',
     &                (creep_param(i,3))

          enddo

      elseif (creepno.eq.4d+0) then



          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Reference creep rate - B
              write(6,*) 'creep_param(1): ',
     &                (creep_param(i,1))
c             Activation energy for creep - Qc
              write(6,*) 'creep_param(2): ',
     &                (creep_param(i,2))
c             Activation volume - V
              write(6,*) 'creep_param(3): ',
     &                (creep_param(i,3))


          enddo


c         Other creep laws are to be added here!

      endif




c     write slip hard parameters
c     Voce type hardening law
      if (modelno.eq.1d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             initial slip resistance
              write(6,*) 'sliphard_param(1): ',
     &                sliphard_param(i,1)
c             hardening rate
              write(6,*) 'sliphard_param(2): ',
     &                sliphard_param(i,2)
c             saturation slip resistance
              write(6,*) 'sliphard_param(3): ',
     &                sliphard_param(i,3)
c             slip hardening exponent
               write(6,*) 'sliphard_param(4): ',
     &                sliphard_param(i,4)


          enddo

c         Other slip rate laws are to be added here!


      elseif (modelno.eq.2d+0) then



          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)
c             initial slip resistance - tauc0
              write(6,*) 'sliphard_param(1): ',
     &                sliphard_param(i,1)
c             hardening rate - h0
              write(6,*) 'sliphard_param(2): ',
     &                sliphard_param(i,2)
c             slip hardening exponent - m
              write(6,*) 'sliphard_param(3): ',
     &                sliphard_param(i,3)
c             activation energy for creep - Q
              write(6,*) 'sliphard_param(4): ',
     &                sliphard_param(i,4)
c             softening exponent for creep - d
              write(6,*) 'sliphard_param(5): ',
     &                sliphard_param(i,5)
c             fitting parameter for creep - AF
              write(6,*) 'sliphard_param(6): ',
     &                sliphard_param(i,6)
c             backstress coefficient - h
              write(6,*) 'sliphard_param(7): ',
     &                sliphard_param(i,7)
c             backstress coefficient - hD
              write(6,*) 'sliphard_param(8): ',
     &                sliphard_param(i,8)

          enddo

      elseif (modelno.eq.3d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             initial slip resistance - tauc0
              write(6,*) 'sliphard_param(1): ',
     &                sliphard_param(i,1)
c             hardening rate coefficient - b
              write(6,*) 'sliphard_param(2): ',
     &                sliphard_param(i,2)
c             strain hardening factor - Q
              write(6,*) 'sliphard_param(3): ',
     &                sliphard_param(i,3)
c             backstress coefficient - C
              write(6,*) 'sliphard_param(4): ',
     &                sliphard_param(i,4)
c             backstress coefficient - D
              write(6,*) 'sliphard_param(5): ',
     &                sliphard_param(i,5)

          enddo

      elseif (modelno.eq.4d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Initial critical resolved shear strength
              write(6,*) 'sliphard_param(1): ',
     &                sliphard_param(i,1)
c             Initial density of statistically stored dislocations
              write(6,*) 'sliphard_param(2): ',
     &                sliphard_param(i,2)
c             Burgers vector
              write(6,*) 'sliphard_param(3): ',
     &                sliphard_param(i,3)
c             Geometric factor for tauc
              write(6,*) 'sliprhard_param(4): ',
     &                (sliphard_param(i,4))
c             Mean-free-path hardening constant
              write(6,*) 'sliphard_param(5): ',
     &                sliphard_param(i,5)
c             Capture radii
              write(6,*) 'sliphard_param(6): ',
     &                sliphard_param(i,6)
c             Effective radius for backstress
              write(6,*) 'sliphard_param(7): ',
     &                sliphard_param(i,7)
c             Size factor between the unit of calculations (MPa-mm) and mesh size (i.e. Dream3D; micrometers)
              write(6,*) 'sliphard_param(8): ',
     &                sliphard_param(i,8)
c             Scaling factor for hardening
              write(6,*) 'sliphard_param(9): ',
     &                sliphard_param(i,9)


          enddo

      elseif (modelno.eq.5d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             initial slip resistance
              write(6,*) 'sliphard_param(1): ',
     &                sliphard_param(i,1)
c             initial dislocation density
              write(6,*) 'sliphard_param(2): ',
     &                sliphard_param(i,2)
c             Burgers vector
              write(6,*) 'sliphard_param(3): ',
     &                sliphard_param(i,3)
c             Mean-free-path hardening coeff.
              write(6,*) 'sliphard_param(4): ',
     &                sliphard_param(i,4)
c             Softening / Annihiliation coeff.
              write(6,*) 'sliphard_param(5): ',
     &                sliphard_param(i,5)
c             Geometric factor for tauc
              write(6,*) 'sliphard_param(6): ',
     &                sliphard_param(i,6)
c             backstress coefficient - h
              write(6,*) 'sliphard_param(7): ',
     &                sliphard_param(i,7)
c             backstress coefficient - hD
              write(6,*) 'sliprhard_param(8): ',
     &                sliphard_param(i,8)
c             Scaling factor for hardening
              write(6,*) 'sliphard_param(9): ',
     &                sliphard_param(i,9)
c             Size conversion
              write(6,*) 'sliphard_param(10): ',
     &                sliphard_param(i,10)

          enddo

c irradiation hardening addition to model 4

      elseif (modelno.eq.6d+0) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             Initial critical resolved shear strength
              write(6,*) 'sliphard_param(1): ',
     &                sliphard_param(i,1)
c             Initial density of statistically stored dislocations
              write(6,*) 'sliphard_param(2): ',
     &                sliphard_param(i,2)
c             Burgers vector
              write(6,*) 'sliphard_param(3): ',
     &                sliphard_param(i,3)
c             Geometric factor for tauc
              write(6,*) 'sliprhard_param(4): ',
     &                (sliphard_param(i,4))
c             Mean-free-path hardening constant
              write(6,*) 'sliphard_param(5): ',
     &                sliphard_param(i,5)
c             Capture radii
              write(6,*) 'sliphard_param(6): ',
     &                sliphard_param(i,6)
c             Effective radius for backstress
              write(6,*) 'sliphard_param(7): ',
     &                sliphard_param(i,7)
c             Size factor between the unit of calculations (MPa-mm) and mesh size (i.e. Dream3D; micrometers)
              write(6,*) 'sliphard_param(8): ',
     &                sliphard_param(i,8)
c             Scaling factor for hardening
              write(6,*) 'sliphard_param(9): ',
     &                sliphard_param(i,9)
c             Frank loop dislocation annihilation area
              write(6,*) 'sliphard_param(10): ',
     &                sliphard_param(i,10)
c             saturation dislocation density for frank loop dislocations
              write(6,*) 'sliphard_param(11): ',
     &                sliphard_param(i,11)
c             Initial frank loop dislocation density
              write(6,*) 'sliphard_param(12): ',
     &                sliphard_param(i,12)

          enddo

c         Other slip rate laws are to be added here!


      endif














c     Slip interactions
      if (interno.eq.0d+0) then

          write(6,*) 'no slip interactions'

      elseif (interno.eq.1d+0) then

          write(6,*) 'latent hardening matrix'


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)



c             latent hardening coefficient
              write(6,*) 'slipint_param(1): ',
     &                slipint_param(i,1)
          enddo




      elseif (interno.eq.2d+0) then


          write(6,*) 'interaction matrix'


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             latent hardening coefficient
              write(6,*) 'slipint_param(1): ',
     &                slipint_param(i,1)
              write(6,*) 'slipint_param(2): ',
     &                slipint_param(i,2)
              write(6,*) 'slipint_param(3): ',
     &                slipint_param(i,3)
              write(6,*) 'slipint_param(4): ',
     &                slipint_param(i,4)
              write(6,*) 'slipint_param(5): ',
     &                slipint_param(i,5)
              write(6,*) 'slipint_param(6): ',
     &                slipint_param(i,6)


          enddo


      elseif (interno.eq.3d+0) then


          write(6,*) 'interaction matrix'


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             latent hardening coefficient
              write(6,*) 'slipint_param(1): ',
     &                slipint_param(i,1)
              write(6,*) 'slipint_param(2): ',
     &                slipint_param(i,2)
              write(6,*) 'slipint_param(3): ',
     &                slipint_param(i,3)
              write(6,*) 'slipint_param(4): ',
     &                slipint_param(i,4)
              write(6,*) 'slipint_param(5): ',
     &                slipint_param(i,5)
              write(6,*) 'slipint_param(6): ',
     &                slipint_param(i,6)
              write(6,*) 'slipint_param(7): ',
     &                slipint_param(i,7)


          enddo


c     Hardening coefficient for irradiation model
      elseif (interno.eq.4d+0) then

          write(6,*) 'interaction matrix'


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)

c             self and latent hardening coefficients
              write(6,*) 'slipint_param(1): ',
     &                slipint_param(i,1)
              write(6,*) 'slipint_param(2): ',
     &                slipint_param(i,2)
              write(6,*) 'slipint_param(3): ',
     &                slipint_param(i,3)
              write(6,*) 'slipint_param(4): ',
     &                slipint_param(i,4)


          enddo

      endif




c     write length scale parameters
      if ((GSeffect.eq.1d+0).or.(GSeffect.eq.2d+0)) then


          do i=1,totphase

              write(6,*) 'phase no: ', i
              write(6,*) 'phase ID: ', phases(i)


c             linear hardening coefficient - "k" [MPa mm^0.5]
              write(6,*) 'grainsize_param(1)', grainsize_param(i,1)

c             grain boundary mismatch exponent - "c" [-]
              write(6,*) 'grainsize_param(2)', grainsize_param(i,2)

c             size correction factor micrometer to mm - "SF" [-]
              write(6,*) 'grainsize_param(3)', grainsize_param(i,3)


          enddo


      endif



c     Threshold value for stress update algorithm
c     Allocate first
      allocate (dS_cr(totphase))
      do i=1,totphase
          dS_cr(i) = dSratio_cr * sliphard_param(i,1)
      enddo








c     Length scale model parameters
c     Read the value for "nodeout" from param_array.inc file
      if ((GSeffect.eq.1d+0).or.(GSeffect.eq.2d+0)) then

c         Read param_array.inc file
          open(150,file=str // '/param_array.inc',action='read',
     &status='old')

          read(150,*) param, line1
          read(150,*) param, line2
          read(150,*) param, line3

          close(150)





          ind=index(line3, ')')



c          write(6,*) 'ind: ', ind

c          write(6,*) 'line1: ', line1

c          write(6,*) 'line2: ', line2

c          write(6,*) 'line3: ', line3

          ch = line3(10:ind-1)

c          write(6,*) 'ch: ', ch



          read(ch,*) nodeout




          write(6,*) 'nodeout', nodeout





      endif














ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(200,file=str // '/materials.dat',action='read',
     &status='old')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Read Euler angles, grain ID,  material ID
      allocate (Euler(numel,3))
      allocate (phaseID(numel))
      allocate (grainID(numel))



      do i=1,numel
c         dummy first read - Euler angles
          read(200,*) dum1, phi1, PHI, phi2, dum2, dum3
          iele = int(dum1)
          Euler(iele,1) = phi1
          Euler(iele,2) = PHI
          Euler(iele,3) = phi2
          grainID(iele) = int(dum2)
          phaseID(iele) = int(dum3)
c            call flush(6)
      enddo

      close(200)




c     Process Euler angles to find the length scale parameters
c     "numgrain": former "totalfeat"
c     "grainori": former "oriensh"

      numgrain = maxval(grainID(:))

c     Allocate the array that stores the Euler angles of the grains
      allocate (grainori(numgrain,3))

c     Assign the Euler angles of the grains
c     Angles are in degrees
      do i=1,numel

          grainori(grainID(i),:)=Euler(i,:)

      enddo



c     Number of grains in the mesh
      write(6,*) 'numgrain', numgrain


c     Write Euler angles of the grains
      write(6,*) 'grain no', 'Euler angles (deg)'
      do i=1,numgrain

          write(6,*) i, ' ', grainori(i,:)



      enddo



c     The entered phase indices may not be 1 and 2 always, and need not to be in order
c     So, every element needs to have a phase index based on the assignment defined by the user
      allocate (phaseind(numel))
      phaseind=0d+0

      do i=1,numel

c         Loop through the entered phases and find the index
          do j=1,totphase

              if (phaseID(i).eq.phases(j)) then

                  phaseind(i) = j



              endif



          enddo


      enddo




c      write(6,*) 'materials.dat'
cc     Output grain orientations and weights
c      do i=1,numel
c          write(6,*) 'Element no.: ', i
c          write(6,*)'Euler angles (deg): ', (Euler(i,j),j=1,3)
c          write(6,*)'Grain ID: ', (grainID(i))
c          write(6,*)'Phase ID: ', (phaseID(i))
c          write(6,*)'Phase index: ', (phaseind(i))
c          call flush(6)
c      enddo



	return
	end subroutine read_inputs
c











      subroutine read_grainsize(str)
      use globalvars, only : nodex, nodey, nodez, boundgrain,
     &elcent, numel, numgrain, grainori, nodeout
      implicit none


      character(len=:), allocatable :: str
      integer*8 :: nr, nc, i
      real(8) :: numrowval,numcolval


      real(8), allocatable :: rowdata(:)



c     allocate arrays
      allocate (boundgrain(numgrain,nodeout))
      allocate (elcent(numel,3))
      allocate (nodex(numgrain,nodeout))
      allocate (nodey(numgrain,nodeout))
      allocate (nodez(numgrain,nodeout))







      ! Read in bound feature array
      open(300,file=str // '/boundfeat.bin',
     &form='unformatted',status='old',access='stream')

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
      open(400, file=str //  '/el_centroid.bin',
     &form='unformatted',status='old',access='stream')


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
      open(500, file=str //  '/xvalues.bin',
     &form='unformatted',status='old',access='stream')

      read(500) numrowval
      nr=int(numrowval)

      read(500) numcolval
      nc=int(numcolval)





      open(600, file=str //  '/yvalues.bin',
     &form='unformatted',status='old',access='stream')

      read(600) numrowval
      nr=int(numrowval)

      read(600) numcolval
      nc=int(numcolval)




      open(700, file=str //  '/zvalues.bin',
     &form='unformatted', status='old',access='stream')

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





c      write(6,*) 'intotalfeat', intotalfeat


cc     read Euler angles
c      include 'orien.inc'



	!reshape array
c      allocate(oriensh(intotalfeat,3))

c	oriensh=reshape(orient, (/intotalfeat,3/), order=(/2,1/))

c      write(6,*) 'oriensh'
c      do i=1,intotalfeat
c          write(6,*) oriensh(i,1:3)
c      enddo

c      write(6,*) 'nodex'
c      do i=1,intotalfeat
c          write(6,*) nodex(i,:)
c      enddo

c      write(6,*) 'nodey'
c      do i=1,intotalfeat
c          write(6,*) nodey(i,:)
c      enddo


c      write(6,*) 'nodez'
c      do i=1,intotalfeat
c          write(6,*) nodez(i,:)
c      enddo



c      write(6,*) 'elcent'
c      do i=1,totalels
c          write(6,*) elcent(i,:)
c      enddo



	return
      end subroutine read_grainsize



















c
c	This subroutine allocates the variables that have to be stored
c	INPUTS: Element number; el, integration point number; ip
	subroutine allocate_arrays

      use globalvars, only: numel,numip,numslip,
     &numgrain,numstvar,largenum,global_state,
     &global_Fp_t,sliphard_param,global_Fe,
     &global_Fe_t,global_state_t,global_ori,
     &global_jacob_t,global_gamma,gradIP2IP,
     &global_jacob,global_sigma,
     &global_sigma_damaged,global_sigma_damaged_t,
     &global_S,global_S_t,
     &global_S_damaged,global_S_damaged_t,
     &global_gammadot,coords_init,
     &global_gamma_t,global_gamma_sum,drhoGND,
     &global_gamma_sum_t,global_Fp,outerabstol,
     &global_sigma_t, grainsize_init,maxnumslip,
     &global_state0,global_Fr0,t_old,inc_old,
     &global_gammadot_t,global_coords,
     &grainmorph,totphase

      use globalvars, only: global_damage
      use globalvars, only: global_F_pos
      use globalvars, only: global_F_neg
      use globalvars, only: global_pk2_pos
      use globalvars, only: global_Wp
      use globalvars, only: global_Wp_t
      use globalvars, only: global_f_ep_c
      use globalvars, only: global_f_ep_c_t

	implicit none
	integer i
c



	allocate (t_old(numel,numip))
	t_old=0.0d+0

      allocate (inc_old(numel,numip))
	inc_old=0d+0

      allocate(gradIP2IP(numel,numip+1,3,numip))
      gradIP2IP=0.0

	allocate (global_Fp(numel,numip,3,3))
	global_Fp=0.0d+0
	allocate (global_Fe(numel,numip,3,3))
	global_Fe=0.0d+0
      allocate (global_Fr0(numel,numip,3,3))
	global_Fr0=0.0d+0
      do i=1,3
          global_Fp(:,:,i,i)=1.0d+0
          global_Fe(:,:,i,i)=1.0d+0
          global_Fr0(:,:,i,i)=1.0d+0
      enddo

      allocate (global_state0(numel,numip,maxnumslip,numstvar))
      global_state0=0.0d+0
      allocate (global_state(numel,numip,maxnumslip,numstvar))
      global_state=0.0d+0
      allocate (global_state_t(numel,numip,maxnumslip,numstvar))
      global_state_t=0.0d+0

      allocate (outerabstol(totphase,numstvar))



      allocate (global_S(numel,numip,6))
      global_S = 0.0d+0
      allocate (global_S_t(numel,numip,6))
      global_S_t = 0.0d+0
      allocate (global_S_damaged(numel,numip,6))
      global_S_damaged = 0.0d+0
      allocate (global_S_damaged_t(numel,numip,6))
	allocate (global_jacob_t(numel,numip,6,6))
	global_jacob_t=0.0d+0
      allocate (global_jacob(numel,numip,6,6))
	global_jacob=0.0d+0
      allocate (global_sigma(numel,numip,6))
	global_sigma=0.0d+0
      allocate (global_sigma_damaged(numel,numip,6))
      global_sigma_damaged = 0.0d+0
      allocate (global_sigma_damaged_t(numel,numip,6))
      global_sigma_damaged_t = 0.0d+0
      allocate (global_sigma_t(numel,numip,6))
	global_sigma_t=0.0d+0
	allocate (global_ori(numel,numip,3,3))
	global_ori=0.0d+0
      allocate(global_gammadot(numel,numip,maxnumslip))
	global_gammadot=0.0d+0
      allocate(global_gammadot_t(numel,numip,maxnumslip))
	global_gammadot_t=0.0d+0

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


c
c
      allocate(global_gamma(numel,numip,maxnumslip))
      allocate(global_gamma_t(numel,numip,maxnumslip))
	global_gamma=0.0d+0
      global_gamma_t=0.0d+0


      allocate(global_gamma_sum(numel,numip))
      allocate(global_gamma_sum_t(numel,numip))
	global_gamma_sum=0.0d+0
      global_gamma_sum_t=0.0d+0

      allocate(global_coords(numel,numip,3))
	global_coords=0.0d+0

      allocate(global_damage(numel,numip))
	  ! damage variables is initialised by MOOSE
	  ! this is just a local copy of the MOOSE variable
	  global_damage = 0.0d+0

      allocate(global_F_pos(numel,numip))
	  global_F_pos = 0.0d+0
	  
      allocate(global_Wp(numel,numip))
      global_Wp = 0.0d+0

      allocate(global_Wp_t(numel,numip))
      global_Wp_t = 0.0d+0

      allocate(global_F_neg(numel,numip))
      global_F_neg = 0.0d+0

      allocate(global_pk2_pos(numel,numip,3,3))
      global_pk2_pos = 0.0d+0

      allocate(coords_init(numel,numip))
      coords_init=0d+0
    
      allocate(global_f_ep_c(numel,numip))
      global_f_ep_c = 1d+0
      allocate(global_f_ep_c_t(numel,numip))
      global_f_ep_c_t = 1d+0

c     A large number is assigned since 1/X is used in the calculations
      allocate(grainmorph(numel,numip,maxnumslip))
	grainmorph=largenum


      allocate(drhoGND(numel,numip,maxnumslip,2))
	drhoGND=0.0d+0

	return
	end subroutine allocate_arrays



c	This initializes the state variables
	subroutine initialize_statevars
	use globalvars, only: modelno, sliphard_param,
     &global_state_t, global_state, global_state0,numel,
     &numip,numslip,numstvar,phaseind
	implicit none
	integer ie, ip, is, iv, iph
c



      do ie = 1,numel
          iph = phaseind(ie)
          do ip = 1,numip
              do is = 1,numslip(iph)



      if (modelno.eq.1d+0) then

          global_state0(ie,ip,is,1)=sliphard_param(iph,1)
          global_state(ie,ip,is,1)=sliphard_param(iph,1)
          global_state_t(ie,ip,is,1)=sliphard_param(iph,1)



      elseif (modelno.eq.2d+0) then

c         Strain/slip hardening term
          global_state0(ie,ip,is,1)=sliphard_param(iph,1)
          global_state(ie,ip,is,1)=sliphard_param(iph,1)
          global_state_t(ie,ip,is,1)=sliphard_param(iph,1)

c         Kinematic hardening term
          global_state0(ie,ip,is,2)=0.0d+0
          global_state(ie,ip,is,2)=0.0d+0
          global_state_t(ie,ip,is,2)=0.0d+0

       elseif (modelno.eq.3d+0) then

c         Strain/slip hardening term
          global_state0(ie,ip,is,1)=sliphard_param(iph,1)
          global_state(ie,ip,is,1)=sliphard_param(iph,1)
          global_state_t(ie,ip,is,1)=sliphard_param(iph,1)

c         Kinematic hardening term
          global_state(ie,ip,is,2)=0.0d+0
          global_state(ie,ip,is,2)=0.0d+0
          global_state_t(ie,ip,is,2)=0.0d+0


       elseif (modelno.eq.4d+0) then



c         Statistically stored initial dislocations (SSD)
          global_state0(ie,ip,is,1)=sliphard_param(iph,2)
          global_state(ie,ip,is,1)=sliphard_param(iph,2)
          global_state_t(ie,ip,is,1)=sliphard_param(iph,2)

c         All the state variables (GNDe/GNDs/backstress)
          global_state0(ie,ip,is,2)=0.0d+0
          global_state(ie,ip,is,2)=0.0d+0
          global_state_t(ie,ip,is,2)=0.0d+0

          global_state0(ie,ip,is,3)=0.0d+0
          global_state(ie,ip,is,3)=0.0d+0
          global_state_t(ie,ip,is,3)=0.0d+0

          global_state0(ie,ip,is,4)=0.0d+0
          global_state(ie,ip,is,4)=0.0d+0
          global_state_t(ie,ip,is,4)=0.0d+0


       elseif (modelno.eq.5d+0) then

c         Strain/slip hardening term
          global_state0(ie,ip,is,1)=sliphard_param(iph,2)
          global_state(ie,ip,is,1)=sliphard_param(iph,2)
          global_state_t(ie,ip,is,1)=sliphard_param(iph,2)

c         Kinematic hardening term
          global_state0(ie,ip,is,2)=0.0d+0
          global_state(ie,ip,is,2)=0.0d+0
          global_state_t(ie,ip,is,2)=0.0d+0


c      Initialize defects for irradiation model
       elseif (modelno.eq.6d+0) then


c         Irradiation induced frank dislocation loops
          global_state0(ie,ip,is,5)=sliphard_param(iph,12)
          global_state(ie,ip,is,5)=sliphard_param(iph,12)
          global_state_t(ie,ip,is,5)=sliphard_param(iph,12)

c         Statistically stored initial dislocations (SSD)
          global_state0(ie,ip,is,1)=sliphard_param(iph,2)
          global_state(ie,ip,is,1)=sliphard_param(iph,2)
          global_state_t(ie,ip,is,1)=sliphard_param(iph,2)



c         All the state variables (GNDe/GNDs/backstress)
          global_state0(ie,ip,is,2)=0.0d+0
          global_state(ie,ip,is,2)=0.0d+0
          global_state_t(ie,ip,is,2)=0.0d+0

          global_state0(ie,ip,is,3)=0.0d+0
          global_state(ie,ip,is,3)=0.0d+0
          global_state_t(ie,ip,is,3)=0.0d+0

          global_state0(ie,ip,is,4)=0.0d+0
          global_state(ie,ip,is,4)=0.0d+0
          global_state_t(ie,ip,is,4)=0.0d+0


      endif

              enddo
          enddo
      enddo




	return
      end subroutine initialize_statevars









c	This initializes the state variables related with the grain size
	subroutine initialize_grainsize(el_no,ip_no,gr_no,ph_no,coords)
!--------------------------------------------------------------------
!       Include file to read in an array containing
!       all information on grain orientations.
      use lengthscale, only: grainsize
      use globalvars, only: grainsize_param, numslip, b_slip, n_slip,
     &global_state, global_state_t, sliphard_param, global_state0,
     &GSeffect, largenum, totphase



	implicit none



	integer el_no, ip_no, gr_no, ph_no, i
      real(8) coords(3), nslip(numslip(ph_no),3)
      real(8) bslip(numslip(ph_no),3)
      real(8) k, c, SF, kval, L, R, val
      real(8) lm(numslip(ph_no)), ldistance(numslip(ph_no))
      real(8) rdistance(numslip(ph_no))

c      real(8) nodex(totalfeat,nodeout), nodey(totalfeat,nodeout),
c     &nodez(totalfeat,nodeout), elcent(totalels,3), oriensh(totalfeat,3)
c




c	common nodex,nodey,nodez,boundgrain,elcent,intotalfeat






c     Length scale parameters
      k = grainsize_param(ph_no,1)
      c = grainsize_param(ph_no,2)
      SF = grainsize_param(ph_no,3)


c      write(6,*) '-------------'
c      write(6,*) 'el_no',el_no
c      write(6,*) 'ip_no',ip_no
c      write(6,*) 'gr_no',gr_no


c     Slip directions
      bslip = b_slip(ph_no,:,:)
      nslip = n_slip(ph_no,:,:)




c     Run the length scale code for analysis
      call grainsize(coords,ldistance,rdistance,lm,
     &bslip,nslip,numslip(ph_no),gr_no,el_no)





c     Assign the strain/slip hardening term
c      write(6,*) 'global_state'
      do i=1,numslip(ph_no)



c         L-distance (L in ref.)
          L = ldistance(i)*SF

c         R-distance (X in ref.)
          R = rdistance(i)*SF

c         Strength coefficient
          kval = k * (1.0d+0-lm(i))**(c)




          if (GSeffect.eq.1d+0) then

c             Value of length scale effect (without stress distribution)
              val = kval/dsqrt(L)


          elseif (GSeffect.eq.2d+0) then

c             Value of length scale effect (full equation)
              val=kval/dsqrt(L)*((R+0.5d+0*L)/dsqrt((R+0.5d+0*L)**2.0d+0
     &- (0.5d+0*L)**2.0d+0)-1.0d+0)

          endif





c         val is infinity or NaN
          if ((val.gt.largenum).or.(isnan(val))) then
c
              global_state0(el_no,ip_no,i,1)=sliphard_param(ph_no,1)
c
              global_state(el_no,ip_no,i,1)=sliphard_param(ph_no,1)
c
              global_state_t(el_no,ip_no,i,1)=sliphard_param(ph_no,1)
c
          else
c         Constants are for MPa.mm units so the distance


c             Micrometers are converted to mm
              global_state0(el_no,ip_no,i,1)=sliphard_param(ph_no,1)
     &        +val

              global_state(el_no,ip_no,i,1)=sliphard_param(ph_no,1)
     &        +val

              global_state_t(el_no,ip_no,i,1)=sliphard_param(ph_no,1)
     &        +val


          endif


c          write(6,*) global_state(el_no,ip_no,i,1)



      enddo


c      write(6,*) 'ldistance'
c      write(6,*) ldistance
c      write(6,*) 'rdistance'
c      write(6,*) rdistance
c      write(6,*) 'lm'
c      write(6,*) lm
c
c      write(6,*) '-------------'
	return
	end subroutine initialize_grainsize
c
c
c
c
c
c	This subroutine reads the text file for grainshape analysis and
c     initializes the slip distances
c	USES:
	subroutine initialize_grainshape(str)
	use globalvars, only : numel, numip, numgrain, numslip, b_slip,
     & n_slip, grainori, grainmorph, sliphard_param, grainID,
     & smallnum, largenum, GSeffect, phaseind
      use globalsubs, only: ang2ori, cross, normvec
	implicit none
      character(len=:), allocatable   :: str
      integer ig, is, ie, ip, el_no, ph
      real(8) D1(numgrain,3), D2(numgrain,3), D3(numgrain,3)
      real(8) d, dx, dy, dz, d1vec(3), d2vec(3), d3vec(3), ang(3)
      real(8) R(3,3), SlipDist
      real(8) d1v(3), d2v(3), d3v(3)
      real(8) d1len, d2len, d3len d1u(3), d2u(3), d3u(3)
      real(8) Dave(numgrain),Dmax(numgrain), conv

      Dave=0.0d+0
      Dmax=0.0d+0





      if (GSeffect.eq.3d+0) then

c     Read the file
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          open(1000,file=str // '/diameter_xyz.txt',action='read',
     &    status='old')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c         Read input data for analysis

          write(6,*) 'grain size from ellipsoids'

c         Read for all grains
          do ig=1,numgrain
              write(6,*) 'grain no: ', ig
c
              read(1000,*) D1(ig,1:3)
              write(6,*) '1st principal axis: ', D1(ig,1:3)
c
              read(1000,*) D2(ig,1:3)
              write(6,*) '2nd principal axis: ', D2(ig,1:3)
c
              read(1000,*) D3(ig,1:3)
              write(6,*) '3rd principal axis: ', D3(ig,1:3)

          enddo

          close(1000)












c         Assign the slip distance to the elements
          do ie=1,numel
              ig = grainID(ie)

c             Find phasind of the grain
              ph = phaseind(ie)


c             Conversion factor for the phase
              conv = sliphard_param(ph,10)



c             1st axis
              d1vec = D1(ig,:)

c             2nd axis
              d2vec = D2(ig,:)

c             3rd axis
              d3vec = D3(ig,:)



c             Transform from sample to crystal
              ang=grainori(ig,:)

              call ang2ori(ang,R)

              d1vec = matmul(R,d1vec)

              d2vec = matmul(R,d2vec)

              d3vec = matmul(R,d3vec)



              do is=1,numslip(ph)




c                 Cross product with the slip plane normal
                  call cross(d1vec, n_slip(ph,is,:), d1v)

                  call cross(d2vec, n_slip(ph,is,:), d2v)

                  call cross(d3vec, n_slip(ph,is,:), d3v)

c                 Magnitude
                  dx = dot_product(d1v, d1v)

                  dy = dot_product(d2v, d2v)

                  dz = dot_product(d3v, d3v)


                  d = dsqrt(dx + dy + dz)

c                 if grain size is equal to zero
                  if (d.lt.smallnum) then
                      d=largenum
                  endif


c                 This scheme only works for model-5
c                 Size factor - conversion to mm
                  SlipDist = d * conv

                  write(6,*) 'Slip system: ', is
                  write(6,*) 'd (size parameter): ', SlipDist


                  do ip=1,numip


                      grainmorph(ie,ip,is) = SlipDist

                  enddo
              enddo
          enddo
















      elseif (GSeffect.eq.4d+0) then


c     Read the file
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          open(1100,file=str // '/diameter_ave.txt',action='read',
     &    status='old')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c         Read input data for analysis

          write(6,*) 'grain size from average diameters'

c         Read for all grains
          do ig=1,numgrain
              write(6,*) 'grain no: ', ig
c
              read(1100,*) Dave(ig)
              write(6,*) 'Av. grain size: ', Dave(ig)

          enddo

          close(1100)



c         Calculate projections along slip directions
          write(6,*) 'Grain shape parameter:'
          do ie=1,numel
              ig = grainID(ie)


c             Find phasind of the grain
              ph = phaseind(ie)


c             Conversion factor for the phase
              conv = sliphard_param(ph,10)


c             Size parameter
              d = Dave(ig)

c             if grain size is equal to zero
              if (d.lt.smallnum) then
                  d=largenum
              endif


              do ip=1,numip
                  do is=1,numslip(ph)





c                     This scheme only works for model-5
c                     Size factor - conversion to mm
                      SlipDist = d * conv

                      grainmorph(ie,ip,is) = SlipDist

                      write(6,*) 'Slip system: ', is
                      write(6,*) 'd (size parameter): ', SlipDist


                  enddo



              enddo



          enddo





      elseif (GSeffect.eq.5d+0) then


c     Read the file
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          open(1200,file=str // '/diameter_max.txt',action='read',
     &    status='old')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c         Read input data for analysis

          write(6,*) 'grain size from average diameters'

c         Read for all grains
          do ig=1,numgrain
              write(6,*) 'grain no: ', ig
c
              read(1200,*) Dmax(ig)
              write(6,*) 'Max. grain size: ', Dmax(ig)

          enddo

          close(1200)



c         Calculate projections along slip directions
          write(6,*) 'Grain shape parameter:'
          do ie=1,numel
              ig = grainID(ie)


c             Find phasind of the grain
              ph = phaseind(ie)


c             Conversion factor for the phase
              conv = sliphard_param(ph,10)








              d = Dmax(ig)

c             if grain size is equal to zero
              if (d.lt.smallnum) then
                  d=largenum
              endif

c             For each slip system - the same value is used!!!
              do ip=1,numip
                  do is=1,numslip(ph)


c                     This scheme only works for model-5
c                     Size factor - conversion to mm
                      SlipDist = d * conv


                      grainmorph(ie,ip,is) = SlipDist


                      write(6,*) 'Slip system: ', is
                      write(6,*) 'd (size parameter): ', SlipDist







                  enddo



              enddo




          enddo



      endif














	return
      end subroutine initialize_grainshape
c
c
c
c
c

c	This subroutine initializes the output .txt files
	subroutine initialize_outputfiles(str)
      use globalvars, only: output_vars, maxnumslip, numstvar
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

c     Write the legend of the output variables (UVARM)
c     Outputs to extract:
c     1: misorientation angle
c     2: cumulative slip
c     3: average of state variables over slip systems
c     4: slip rates per slip system
c     5: state variables per slip system
      count = 0d+0
      open(850,file=str //'/UVARM_legend.txt',action='write',
     &status='replace')

      if (output_vars(1) .eq. 1d+0) then

          count = count + 1d+0

          write(850,'(A6,I2,A28)') 'UVARM-', count,
     &': misorientation angle [deg]'

      endif


      if (output_vars(2) .eq. 1d+0) then

          count = count + 1d+0

          write(850,'(A6,I2,A12)') 'UVARM-', count, ': total slip'

      endif


      if (output_vars(3) .eq. 1d+0) then

          do isv =1,numstvar

              count = count + 1d+0

              write(850,'(A6,I2,A25,I2)') 'UVARM-', count,
     &': average state variable-', isv

          enddo

      endif


       if (output_vars(4) .eq. 1d+0) then

          do iss =1,maxnumslip

              count = count + 1d+0

              write(850,'(A6,I2,A27,I2)') 'UVARM-',count,
     &': slip rate of slip system-', iss

          enddo

       endif



      if (output_vars(5) .eq. 1d+0) then

          do isv =1,numstvar

              do iss=1,maxnumslip

                  count = count + 1d+0

                  write(850,'(A6,I2,A17,I1,A20,I2)') 'UVARM-', count,
     &': state variable-', isv,' of the slip system-', iss

              enddo


          enddo

      endif




      close(850)






	return
      end subroutine initialize_outputfiles









c	This subroutine assigns residual deformation
c	USES: I3(3,3),I6(6,6),I9(9,9),eijk(3,3,3)
	subroutine initialize_residualdeformation(str)
	use globalvars, only : resdef, global_Fe_t, global_Fe, global_Fr0,
     &numel, numip
      use globalsubs, only: convert9to3x3, invert3x3
	implicit none
	integer i, iele, iquad
      real(8) Fe0_vec(9), Fe0(3,3), invFe0(3,3)
	  real(8) det, dum1, dum2, R(3,3)
      character(len=:), allocatable   :: str

      if (resdef.eq.1d+0) then

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          open(900,file=str // '/resdef.dat',action='read',status='old')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c          write(6,*) 'residual deformation'



c         The size of the file shall be the same as size of elements*numip
          do i=1,numel*numip

c             dummy first read - Euler angles
              read(900,*) dum1, dum2, Fe0_vec(:)
              iele = int(dum1)
              iquad = int(dum2)

c             write(6,*) iele, iquad, Fe0_vec(:)

c             calculate 3x3 deformation gradient
              call convert9to3x3(Fe0_vec,Fe0)

cc             store the initially assigned orientatios
c              R = global_Fe_t(iele,iquad,:,:)

cc             assign the initial elastic deformation togther with rotation
cc             PRE-MULTIPLY ROTATIONS
c              global_Fe_t(iele,iquad,:,:) = matmul(R,Fe0)
c              global_Fe(iele,iquad,:,:) = matmul(R,Fe0)

c             assign the inverse of this to the residual deformation
              call invert3x3(Fe0,invFe0,det)

              global_Fr0(iele,iquad,:,:)=invFe0




c            call flush(6)
          enddo

          close(900)



      endif





	return
      end subroutine initialize_residualdeformation






c	This subroutine calculates the interpolation mappings necessaru for gnd analysis by slip gradients
c	USES:
	subroutine initialize_gndslipgrad0
	use globalvars, only : nnpe, numip
      use gndslipgrad, only : shafunmap
	implicit none

c     Initialize mappings at the isoparametric space for the element type under concern
      call shafunmap(nnpe,numip)

      write(6,*) 'Shape functions are initialized!'

	return
      end subroutine initialize_gndslipgrad0







c	This subroutine calculates the interpolation mappings necessary for gnd analysis by slip gradients
c	USES:
	subroutine initialize_gndslipgradel
      use globalvars, only : numel, numip,
     &nnpe, invNmat, dNmat, sliphard_param, phaseind,
     &Gmat, global_coords, global_ori, gradIP2IP
      use globalsubs, only: invert3x3
	implicit none
	integer iele, iqpt, islp, i, j, iph
      real(8) IPel(numip,3), Nodeel(nnpe,3), JT(3,3), invJT(3,3), detJ
      real(8) grad(3,numip), dN(3,nnpe), Gmati(3,numip), ori(3,3)
      real(8) c2s(3,3), SF, IPsub(numip,3)






c     Elemental gradient calculation
      do iele=1,numel


          iph=phaseind(iele)

c         Size factor of gradient (mesh size, i.e. micrometer to mm)
          SF = sliphard_param(iph,8)

c          write(6,*) 'element no: ', iele


c         IP coordinates of element NOEL
          IPel = global_coords(iele,:,:)*SF






c          write(6,*) 'Overall IP coords.'
c          do i=1,numip
c              write(6,*) (IPel(i,j), j=1,3)
c          enddo





          Nodeel = matmul(invNmat,IPel)

c          write(6,*) 'node coords.'
c          do i=1,nnpe
c              write(6,*) (Nodeel(i,j), j=1,3)
c          enddo



c         Calculation at each integration point
          do iqpt=1,numip


c             Crystal orientation matrix
              ori = global_ori(iele,iqpt,:,:)


c             Calculate jacobian transpose
c              Gmati(1:3,1:numip)=Gmat(iqpt,:,:)
              dN = dNmat(iqpt,1:3,1:nnpe)
              JT = matmul(dN,Nodeel)




c             Calculate the inverse of the jacobian
              call invert3x3(JT,invJT,detJ)



c             Gradient at the global coordinate reference/physical space
c             Maps from the nodal values
c              dN(1:3,1:nnpe) = dNmat(iqpt,1:3,1:nnpe)
              Gmati = Gmat(iqpt,:,:)
              grad = matmul(invJT, Gmati)



cc             Maps from the IP values
c              grad = matmul(grad, invNmat)



c             Transform the gradient into the crystal reference
              grad = matmul(ori, grad)






c             Save the value of the gradient
              gradIP2IP(iele,iqpt,:,:) = grad



          enddo


c         Gradient at the element center
c         using NODE values

c         Calculate jacobian transpose
c          Gmati(1:3,1:numip)=Gmat(numip+1,:,:)
          dN = dNmat(numip+1,1:3,1:nnpe)


c         Switch the IP coordinates to match the node ordering in shape functions
          IPsub=IPel
          IPsub(3,:)=IPel(4,:)
          IPsub(4,:)=IPel(3,:)
          IPsub(7,:)=IPel(8,:)
          IPsub(8,:)=IPel(7,:)


          JT = matmul(dN,IPsub)



c         Calculate the inverse of the jacobian
          call invert3x3(JT,invJT,detJ)


c         Gradient at the global coordinate reference/physical space
          grad = matmul(invJT, dN)

!c         Gradient at the global coordinate reference/physical space
!c         Maps from the nodal values
!c          dN(1:3,1:nnpe) = dNmat(numip+1,1:3,1:nnpe)
!          Gmati = Gmat(numip+1,:,:)
!          grad = matmul(invJT, Gmati)
!
!
!cc         Maps from the IP values
!c          grad = matmul(grad, invNmat)
!


c         Transform the gradient into the crystal reference
          grad = matmul(ori, grad)


c         Store
          gradIP2IP(iele,numip+1,:,:) = grad







      enddo

      write(6,*) 'Gradient operators are initialized!'

	return
      end subroutine initialize_gndslipgradel
c
c
c
c
c
c
c	This subroutine initializes the relative tolerances
	subroutine initialize_relativetolerances
      use globalvars, only : outertol, outerabstol,
     & numstvar, sliphard_param, G, modelno, phases, totphase
	implicit none
	integer i, j, iph
      real(8) alpha, b


      do iph=1,totphase


          write(6,*) 'phase no: ', iph
          write(6,*) 'phase ID: ', phases(iph)


c         The relative tolerance specs depend on the model type
c         tauc
          if (modelno.eq.1d+0) then

              outerabstol(iph,1) = outertol * sliphard_param(iph,1)

c         tauc
          elseif (modelno.eq.2d+0) then

              outerabstol(iph,:) = outertol * sliphard_param(iph,1)

c         tauc
          elseif (modelno.eq.3d+0) then

              outerabstol(iph,:) = outertol * sliphard_param(iph,1)

c         rho0==>(1:3), tauc==>(4)
          elseif (modelno.eq.4d+0) then

              outerabstol(iph,1:3) = outertol * sliphard_param(iph,2)

c             Geometric factor for slip resistance calculation
              alpha = sliphard_param(iph,4)

c             Burgers vector for slip resistance calculation
              b = sliphard_param(iph,3)

              outerabstol(iph,4) = outertol * alpha * G(iph) * b *
     &        dsqrt(sliphard_param(iph,2))

c         rho0==>(1), tauc==>(2)
          elseif (modelno.eq.5d+0) then

              outerabstol(iph,1) = outertol * sliphard_param(iph,2)

c             Geometric factor for slip resistance calculation
              alpha = sliphard_param(iph,4)

c             Burgers vector for slip resistance calculation
              b = sliphard_param(iph,3)

              outerabstol(iph,2) = outertol * alpha * G(iph) * b *
     &        dsqrt(sliphard_param(iph,2))

c         output for irradiation model
c         rho0==>(1:3 & 5), tauc==>(4)
          elseif (modelno.eq.6d+0) then

              outerabstol(iph,1:3) = outertol * sliphard_param(iph,2)


c             Added by Eralp ---- PLS. CHECK!

c             Geometric factor for slip resistance calculation
              alpha = sliphard_param(iph,4)

c             Burgers vector for slip resistance calculation
              b = sliphard_param(iph,3)

              outerabstol(iph,4) = outertol * alpha * G(iph) * b *
     &        dsqrt(sliphard_param(iph,2))




              outerabstol(iph,5) = outertol * sliphard_param(iph,2)





      endif



          write(6,*) 'Absolute values of tolerances'
          do i=1,numstvar
              write(6,*) 'tol. for state var. no.', i, ': ',
     &        outerabstol(iph,i)
          enddo



      enddo










      return
      end subroutine initialize_relativetolerances




      end module initialization
