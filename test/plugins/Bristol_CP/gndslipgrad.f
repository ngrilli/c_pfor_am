c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c
c
      module gndslipgrad
      implicit none
      contains





      subroutine calculategnds(dt)
      use globalvars, only : numel, numip, numslip, gradIP2IP,
     & global_state, global_state_t, sliphard_param, global_gammadot,
     & b_slip, l_slip, drhoGND, phaseind
      implicit none
      integer iele, iqpt, islp, iq, nss, ph_no
      real(8) dt, b, SF, C, gammadot(numip), grad(3,numip)
      real(8) gradgammadot(3), drhoGNDe, drhoGNDs

c     Needed for backstress increment
      drhoGND=0.0d+0

c     Burger's vector
      b = sliphard_param(ph_no,3)



c     Scaling factor for GNDs
      C = sliphard_param(ph_no,9)


c     Loop over the elements
      do iele=1,numel


c         Phase index (not phase id)
          ph_no = phaseind(iele)


c         Number of slip systems
          nss=numslip(ph_no)



          do iqpt=1,numip

c             Gradient operator
              grad = gradIP2IP(iele,iqpt,1:3,1:numip)





c             write(6,*) 'grad'
c             write(6,*) grad

c             Loop over the slip systems
              do islp=1,nss



c                 Vectorize element slip rates
                  gammadot=0.0
                  do iq=1,numip

                      gammadot(iq) = global_gammadot(iele,iq,islp)

                  enddo



c                 Slip gradient
                  gradgammadot = matmul(grad,gammadot)*C/b



c                 Edge dislocation
                  drhoGNDe=-dot_product(b_slip(ph_no,islp,:),
     &            gradgammadot)*dt

                  drhoGND(iele,iqpt,islp,1)=drhoGNDe




c                 Screw dislocation
                  drhoGNDs=dot_product(l_slip(ph_no,islp,:),
     &            gradgammadot)*dt

                  drhoGND(iele,iqpt,islp,2)=drhoGNDs




c                 Store and update the state variables
c                 Edge dislocaitons
                  global_state(iele,iqpt,islp,2)=
     & global_state_t(iele,iqpt,islp,2) + drhoGNDe


c                 Store and update the state variables
c                 Screw dislocaitons
                  global_state(iele,iqpt,islp,3)=
     & global_state_t(iele,iqpt,islp,3) + drhoGNDs


              enddo

          enddo


c



      enddo



      return
      end subroutine calculategnds





      subroutine calculatebackstress1
      use globalvars, only : global_state, numel, numip, numslip,
     & gradIP2IP, G, nu, sliphard_param, b_slip, l_slip, global_state_t,
     & drhoGND, phaseind
      implicit none
      integer iele, iqpt, islp, nss, ph_no
      real(8) rhoGNDe(numip), rhoGNDs(numip), grad(3,numip)
      real(8) nabla_rhoGNDe(3), nabla_rhoGNDs(3), b, R
      real(8) rhoGNDe_x, rhoGNDs_x, cst1, cst2, dX



c     Burger's vector
      b = sliphard_param(ph_no,3)

c     Effective radius
      R = sliphard_param(ph_no,7)






c     Loop over the elements
      do iele=1,numel

c         Phase index (not phase id)
          ph_no = phaseind(iele)

c         Coefficients
          cst1=G(ph_no)*b*(R**2.0)/8.0/(1.0-nu(ph_no))

          cst2=G(ph_no)*b*(R**2.0)/4.0


c         Number of slip systems
          nss=numslip(ph_no)

c         Gradient operator at the element center
          grad = gradIP2IP(iele,numip+1,1:3,1:numip)

c         Loop over the slip systems
          do islp=1,nss


c             Vectorize element GNDs
              do iqpt=1,numip

c                 Edge dislocations
c                  rhoGNDe(iqpt) = global_state(iele,iqpt,islp,2)
                  rhoGNDe(iqpt) = drhoGND(iele,iqpt,islp,1)

c                 Screw dislocations
c                  rhoGNDs(iqpt) = global_state(iele,iqpt,islp,3)
                  rhoGNDs(iqpt) = drhoGND(iele,iqpt,islp,2)

              enddo




c             gradient of edge dislocatinos
              nabla_rhoGNDe = matmul(grad,rhoGNDe)

c             gradient along slip direction
              rhoGNDe_x=dot_product(b_slip(ph_no,islp,:),nabla_rhoGNDe)

c             gradient of screw dislocatinos
              nabla_rhoGNDs = matmul(grad,rhoGNDs)

c             gradient along line direction
              rhoGNDs_x=dot_product(l_slip(ph_no,islp,:),nabla_rhoGNDs)





c             Calculate and assign the same backstress value to all of the ips of an element
              do iqpt=1,numip

c                 Backstress
                  dX = rhoGNDe_x*cst1 - rhoGNDs_x*cst2

c                 Store the backstress as another state variable
                  global_state(iele,iqpt,islp,4) = dX +
     &            global_state_t(iele,iqpt,islp,4)




              enddo






          enddo


      enddo






      return
      end subroutine calculatebackstress1








      subroutine calculatebackstress2
      use globalvars, only : global_state, numel, numip, numslip,
     & gradIP2IP, BS_dyad, Schmid_vec, sliphard_param, G, nu,
     & b_slip, l_slip, n_slip, global_state_t, drhoGND, phaseind
      use globalsubs, only : convert3x3to6
      implicit none
      integer iele, iqpt, islp, i, ph_no, nss
      real(8) rhoGNDe(numip), rhoGNDs(numip), BS(3,3), BS_vec(6)
      real(8) rhoGNDe_x, rhoGNDs_x, rhoGNDe_y, rhoGNDs_y, b, R
      real(8) nabla_rhoGNDe(3), nabla_rhoGNDs(3), dX
      real(8) cst1, cst2, grad(3,numip)



c     Burger's vector
      b = sliphard_param(ph_no,3)

c     Effective radius
      R = sliphard_param(ph_no,7)






c     Loop over the elements
      do iele=1,numel

c         Phase index (not phase id)
          ph_no = phaseind(iele)


c         Coefficients
          cst1=G(ph_no)*b*(R**2.0)/8.0/(1.0-nu(ph_no))

          cst2=G(ph_no)*b*(R**2.0)/4.0

c         Number of slip systems
          nss=numslip(ph_no)

c         Gradient operator at the element center
          grad = gradIP2IP(iele,numip+1,1:3,1:numip)


c         Sum the backstress tensor over the slip systems
          BS = 0.0d+0
c         Loop over the slip systems
          do islp=1,nss


c             Vectorize element GNDs
              do iqpt=1,numip

c                 Edge dislocations
c                  rhoGNDe(iqpt) = global_state(iele,iqpt,islp,2)
                  rhoGNDe(iqpt) = drhoGND(iele,iqpt,islp,1)

c                 Screw dislocations
c                  rhoGNDs(iqpt) = global_state(iele,iqpt,islp,3)
                  rhoGNDs(iqpt) = drhoGND(iele,iqpt,islp,2)

              enddo




c             gradient of edge dislocatinos
              nabla_rhoGNDe = matmul(grad,rhoGNDe)

c             gradient along slip direction
              rhoGNDe_x=dot_product(b_slip(ph_no,islp,:),nabla_rhoGNDe)

c             gradient along normal direction
              rhoGNDe_y=dot_product(n_slip(ph_no,islp,:),nabla_rhoGNDe)



c             gradient of screw dislocatinos
              nabla_rhoGNDs = matmul(grad,rhoGNDs)

c             gradient along line direction
              rhoGNDs_x=dot_product(l_slip(ph_no,islp,:),nabla_rhoGNDs)

c             gradient along normal direction
              rhoGNDs_y=dot_product(n_slip(ph_no,islp,:),nabla_rhoGNDs)




c             Calculate backstress TENSOR

c             Screw dislocation - xz stress component
              BS = BS + BS_dyad(ph_no,islp,1,:,:)*rhoGNDs_y*cst2


c             Screw dislocation - yz stress component
              BS = BS - BS_dyad(ph_no,islp,2,:,:)*rhoGNDs_x*cst2


c             Edge dislocation - xx stress component
              BS = BS - 3.0*BS_dyad(ph_no,islp,3,:,:)*rhoGNDe_y*cst1


c             Edge dislocation - yy stress component
              BS = BS - BS_dyad(ph_no,islp,4,:,:)*rhoGNDe_y*cst1


c             Edge dislocation - zz stress component
              BS = BS - 4.0*BS_dyad(ph_no,islp,5,:,:)*rhoGNDe_y*cst1


c             Edge dislocation - xy stress component
              BS = BS + BS_dyad(ph_no,islp,6,:,:)*rhoGNDe_x*cst1






          enddo


c         Vectorize 3x3 backstress tensor
          call convert3x3to6(BS,BS_vec)







c         Project it to each slip system again using Schmid tensor
          do islp=1,nss

              dX = 0.0d+0
              do i=1,6
                  dX = dX + (BS_vec(i)*Schmid_vec(ph_no,islp,i))
              enddo


c             Vectorize element GNDs
              do iqpt=1,numip

c                 Store the backstress as another state variable
                  global_state(iele,iqpt,islp,4) = dX +
     &            global_state_t(iele,iqpt,islp,4)






              enddo





          enddo





      enddo








      return
      end subroutine calculatebackstress2










c     Calculates the shape function mapping for gradients
      subroutine shafunmap(nnpe,numip)
      use globalvars, only: IPghr, wtghr, invNmat, dNmat, Gmat
      use globalsubs, only: invertnxn
      implicit none
      integer nnpe, i, j, numip, iqpt
      real(8) Nmat(numip,nnpe)
      real(8) g, h, r, N(nnpe), dN(3,nnpe), dNc(3,nnpe)

c     nnpe: number of nodes per element

c     allocate arrays
      allocate (IPghr(numip,3))
      IPghr=0.0d+0
      allocate (wtghr(numip))
      wtghr=0.0d+0
      allocate(invNmat(nnpe,numip))
      invNmat=0.0d+0
c     array size has "numip+1" because of the gradient at the center of the element
      allocate(dNmat(numip+1,3,nnpe))
      dNmat=0.0d+0
      allocate(Gmat(numip+1,3,numip))
      Gmat=0.0d+0

c     4-node linear tetrahedron - C3D4 (nnpe=4, numip=1)

c     10-node quadratic tetrahedron - C3D10 (nnpe=10, numip=4)

c     8-node linear brick - C3D8 (nnpe=8, numip=8)

c     20-node quadratic brick - C3D20 (nnpe=20, numip=27)




c     Solid Brick element - C3D8
      if (nnpe.eq.8d+0) then



c         Quad point locations

          IPghr(1,1:3) = (/-1.0, -1.0, -1.0 /)

          IPghr(2,1:3) = (/ 1.0, -1.0, -1.0 /)

          IPghr(3,1:3) = (/-1.0,  1.0, -1.0 /)

          IPghr(4,1:3) = (/ 1.0,  1.0, -1.0 /)

          IPghr(5,1:3) = (/-1.0, -1.0,  1.0 /)

          IPghr(6,1:3) = (/ 1.0, -1.0,  1.0 /)

          IPghr(7,1:3) = (/-1.0,  1.0,  1.0 /)

          IPghr(8,1:3) = (/ 1.0,  1.0,  1.0 /)


          IPghr = IPghr/dsqrt(3.0d+0)



          write(6,*) 'IP coords for C3D8'
          do i=1,numip
              write(6,*) (IPghr(i,j), j=1,3)
          enddo


c         Integration weights
          wtghr = 1.0



          do iqpt=1,numip

              g = IPghr(iqpt,1)
              h = IPghr(iqpt,2)
              r = IPghr(iqpt,3)


              call shapefunctions(nnpe,g,h,r,N)



              Nmat(iqpt,:) = N


          enddo


          write(6,*) 'Nmat matrix'
          do i=1,numip
              write(6,*) (Nmat(i,j), j=1,nnpe)
          enddo


          call invertnxn(Nmat,invNmat,numip)



          write(6,*) 'invNmat matrix'
          do i=1,nnpe
              write(6,*) (invNmat(i,j), j=1,numip)
          enddo



          do iqpt=1,numip

              g = IPghr(iqpt,1)
              h = IPghr(iqpt,2)
              r = IPghr(iqpt,3)



              call shapefunctionderivatives(nnpe,g,h,r,dN)

              dNmat(iqpt,:,:) = dN

              write(6,*) 'IP no: ', iqpt
              write(6,*) 'dN'
              do i=1,3
                  write(6,*) (dN(i,j), j=1,nnpe)
              enddo

              Gmat(iqpt,:,:) = matmul(dN,invNmat)




          enddo




c         gradient at the element center
          call shapefunctionderivatives(nnpe,0.0d+0,0.0d+0,0.0d+0,dNc)


          dNmat(numip+1,:,:) = dNc

          Gmat(numip+1,:,:) = matmul(dNc,invNmat)



      endif









      end subroutine shafunmap













      subroutine shapefunctions(nnpe,g,h,r,N)
      implicit none
      integer nnpe
      real(8) N(nnpe), g, h, r


c     Solid Brick element - C3D8
      if (nnpe.eq.8) then

          N(1) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0-h) * (1.0d+0-r)

          N(2) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0-h) * (1.0d+0-r)

          N(3) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0+h) * (1.0d+0-r)

          N(4) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0+h) * (1.0d+0-r)




          N(5) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0-h) * (1.0d+0+r)

          N(6) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0-h) * (1.0d+0+r)

          N(7) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0+h) * (1.0d+0+r)

          N(8) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0+h) * (1.0d+0+r)




      endif



      end subroutine shapefunctions






      subroutine shapefunctionderivatives(nnpe,g,h,r,dN)
      implicit none

      integer nnpe
      real(8) dN(3,nnpe), g, h, r


c     Solid Brick element - C3D8
      if (nnpe.eq.8) then


c         dN_dg

          dN(1,1) = 1.0d+0/8.0d+0* (-1.0d+0) * (1.0d+0-h) * (1.0d+0-r)

          dN(1,2) = 1.0d+0/8.0d+0* (1.0d+0) * (1.0d+0-h) * (1.0d+0-r)

          dN(1,3) = 1.0d+0/8.0d+0* (1.0d+0) * (1.0d+0+h) * (1.0d+0-r)

          dN(1,4) = 1.0d+0/8.0d+0* (-1.0d+0) * (1.0d+0+h) * (1.0d+0-r)

          dN(1,5) = 1.0d+0/8.0d+0* (-1.0d+0) * (1.0d+0-h) * (1.0d+0+r)

          dN(1,6) = 1.0d+0/8.0d+0* (1.0d+0) * (1.0d+0-h) * (1.0d+0+r)

          dN(1,7) = 1.0d+0/8.0d+0* (1.0d+0) * (1.0d+0+h) * (1.0d+0+r)

          dN(1,8) = 1.0d+0/8.0d+0* (-1.0d+0) * (1.0d+0+h) * (1.0d+0+r)


c         dN_dh

          dN(2,1) = 1.0d+0/8.0d+0* (1.0d+0-g) * (-1.0d+0) * (1.0d+0-r)

          dN(2,2) = 1.0d+0/8.0d+0* (1.0d+0+g) * (-1.0d+0) * (1.0d+0-r)

          dN(2,3) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0) * (1.0d+0-r)

          dN(2,4) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0) * (1.0d+0-r)


          dN(2,5) = 1.0d+0/8.0d+0* (1.0d+0-g) * (-1.0d+0) * (1.0d+0+r)

          dN(2,6) = 1.0d+0/8.0d+0* (1.0d+0+g) * (-1.0d+0) * (1.0d+0+r)

          dN(2,7) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0) * (1.0d+0+r)

          dN(2,8) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0) * (1.0d+0+r)



c         dN_dr

          dN(3,1) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0-h) * (-1.0d+0)

          dN(3,2) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0-h) * (-1.0d+0)

          dN(3,3) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0+h) * (-1.0d+0)

          dN(3,4) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0+h) * (-1.0d+0)


          dN(3,5) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0-h) * (1.0d+0)

          dN(3,6) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0-h) * (1.0d+0)

          dN(3,7) = 1.0d+0/8.0d+0* (1.0d+0+g) * (1.0d+0+h) * (1.0d+0)

          dN(3,8) = 1.0d+0/8.0d+0* (1.0d+0-g) * (1.0d+0+h) * (1.0d+0)


      endif


      end subroutine shapefunctionderivatives












      end module gndslipgrad
