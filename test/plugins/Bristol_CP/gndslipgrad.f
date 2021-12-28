      module gndslipgrad
      implicit none
      contains
      
      
      
      
      
      subroutine calculategnds(dt)
      use globalvars, only : numel, numip, 
     + numslip, gradIP2IP,
     + global_state, global_state_t, sliphard_param, 
     + global_gammadot_t
      implicit none
      integer iele, iqpt, islp
      real(8) dt, b, SF, gdot(numip), grad(3,numip), rhodotGND(3)

      
c     Burger's vector
      b = sliphard_param(3)
      
c     Size factor (mesh size, i.e. micrometer to mm)
      SF = sliphard_param(8)
      
      
c     Loop over the elements
      do iele=1,numel
          
          
c         Loop over the slip systems
          do islp=1,numslip
          
              
                         
c             Vectorize element slip rates              
              gdot=0.0
              do iqpt=1,numip
                  
                  gdot(iqpt) = global_gammadot_t(iele,iqpt,islp)
              
              enddo
              
              do iqpt=1,numip
                  
c                 Gradient operator
                  grad = gradIP2IP(iele,iqpt,islp,1:3,1:numip)
              
c                 Slip gradient
                  rhodotGND = matmul(grad,gdot)/SF/b

              
c                 Store and update the state variables
c                 Edge dislocaitons
                  global_state(iele,iqpt,islp,2)=
     & global_state_t(iele,iqpt,islp,2) - rhodotGND(1)*dt
              
c                 Store and update the state variables
c                 Screw dislocaitons
                  global_state(iele,iqpt,islp,3)=
     & global_state_t(iele,iqpt,islp,3) + rhodotGND(3)*dt
                  
                  
              enddo       
              
          enddo
          
 
c                   
          
          
          
      enddo
      
      
      
      return
      end subroutine calculategnds
      
      
      
      
      
      
      subroutine calculatebackstress
      use globalvars, only : global_state, numel, numip, numslip,
     & gradIP2IP, G, nu, sliphard_param
      implicit none
      integer iele, iqpt, islp
      real(8) rhoGNDe(numip), rhoGNDs(numip), grad(3,numip)
      real(8) nabla_rhoGNDe(3), nabla_rhoGNDs(3), X, b, R
      
      
      
c     Burger's vector
      b = sliphard_param(3)      
      
c     Effective radius
      R = sliphard_param(7)
      
c     Loop over the elements
      do iele=1,numel
          
          
c         Loop over the slip systems
          do islp=1,numslip

              
              
c             Vectorize element GNDs              
              do iqpt=1,numip
                  
c                 Edge dislocations                  
                  rhoGNDe(iqpt) = global_state(iele,iqpt,islp,2)
                  
c                 Screw dislocations                  
                  rhoGNDs(iqpt) = global_state(iele,iqpt,islp,3)
              
              enddo              
              
              
              
c             Gradient operator
              grad = gradIP2IP(iele,numip+1,islp,1:3,1:numip)
              
              
              
c             gradient of edge dislocatinos              
              nabla_rhoGNDe = matmul(grad,rhoGNDe)              
              
c             gradient of screw dislocatinos              
              nabla_rhoGNDs = matmul(grad,rhoGNDs)                 
      
c             Calculate and assign the same backstress value to all of the ips of an element             
              do iqpt=1,numip      
                  
c                 Backstres              
                  X = nabla_rhoGNDe(1)*G*b*(R**2.)/8./(1.-nu)
     &              - nabla_rhoGNDs(3)*G*b*(R**2.)/4.
                  
c                 Store the backstress as another state variable
                  global_state(iele,iqpt,islp,4) = X
                  
              enddo
              
                  
          enddo
      
          
      enddo
      

                  
              
              
      
      return
      end subroutine calculatebackstress
      
      
      
      
      
      
      
      
      
      
      
      
      
      
c     Calculates the shape function mapping for gradients 
      subroutine shafunmap(nnpe,numip)
      use globalvars, only: IPghr, wtghr, invNmat, dNmat, Gmat
      use globalsubs, only: invertnxn
      implicit none
      integer nnpe, i, j, numip
      real(8) Nmat(numip,nnpe)
      real(8) g, h, r, N(nnpe), dN(3,nnpe), dNc(3,nnpe)
      
c     nnpe: number of nodes per element
      
c     allocate arrays
      allocate (IPghr(numip,3))
      IPghr=0.0
      allocate (wtghr(numip))
      wtghr=0.0
      allocate(invNmat(nnpe,numip))
      invNmat=0.0
c     array size has "numip+1" because of the gradient at the center of the element
      allocate(dNmat(numip+1,3,nnpe))
      dNmat=0.0
      allocate(Gmat(numip+1,3,numip))
      Gmat=0.0

      
c     4-node linear tetrahedron - C3D4 (nnpe=4, numip=1)

c     10-node quadratic tetrahedron - C3D10 (nnpe=10, numip=4)

c     8-node linear brick - C3D8 (nnpe=8, numip=8)
      
c     20-node quadratic brick - C3D20 (nnpe=20, numip=27)


      
      
c     Solid Brick element - C3D8      
      if (nnpe.eq.8) then
          

      
c         Quad point locations

          IPghr(1,1:3) = (/-1.0, -1.0, -1.0 /)
          
          IPghr(2,1:3) = (/ 1.0, -1.0, -1.0 /)
          
          IPghr(3,1:3) = (/-1.0,  1.0, -1.0 /)
          
          IPghr(4,1:3) = (/ 1.0,  1.0, -1.0 /)
          
          IPghr(5,1:3) = (/-1.0, -1.0,  1.0 /)
          
          IPghr(6,1:3) = (/ 1.0, -1.0,  1.0 /)
          
          IPghr(7,1:3) = (/-1.0,  1.0,  1.0 /)
          
          IPghr(8,1:3) = (/ 1.0,  1.0,  1.0 /)
              
          
          IPghr = 1.0d+0/dsqrt(3.0d+0) * IPghr

          
          
          write(6,*) 'IP coords for C3D8'
          do i=1,numip
              write(6,*) (IPghr(i,j), j=1,3)
          enddo 
          
          
c         Integration weights
          wtghr = 1.0
          
          
          
          do i=1,numip
          
              g = IPghr(i,1)
              h = IPghr(i,2)
              r = IPghr(i,3)
              
              
              call shapefunctions(nnpe,g,h,r,N)
              

              
              Nmat(i,:) = N
              
              
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
          
          
    
          do i=1,numip
          
              g = IPghr(i,1)
              h = IPghr(i,2)
              r = IPghr(i,3)
          
          
          
              call shapefunctionderivatives(nnpe,g,h,r,dN)
              
              dNmat(i,:,:) = dN
              
              
           
              Gmat(i,:,:) = matmul(dN,invNmat)
              
              
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
      
          N(1) = 1./8.* (1.-g) * (1.-h) * (1.-r)

          N(2) = 1./8.* (1.+g) * (1.-h) * (1.-r)

          N(3) = 1./8.* (1.+g) * (1.+h) * (1.-r)

          N(4) = 1./8.* (1.-g) * (1.+h) * (1.-r)




          N(5) = 1./8.* (1.-g) * (1.-h) * (1.+r)

          N(6) = 1./8.* (1.+g) * (1.-h) * (1.+r)

          N(7) = 1./8.* (1.+g) * (1.+h) * (1.+r)

          N(8) = 1./8.* (1.-g) * (1.+h) * (1.+r)
      
          

      
      endif
      
      
      
      end subroutine shapefunctions
      

      
      
      
      
      subroutine shapefunctionderivatives(nnpe,g,h,r,dN)
      implicit none
      
      integer nnpe
      real(8) dN(3,nnpe), g, h, r
      
      
c     Solid Brick element - C3D8      
      if (nnpe.eq.8) then
      
          
c         dN_dg          
          
          dN(1,1) = 1./8.* (-1.) * (1.-h) * (1.-r)

          dN(1,2) = 1./8.* (1.) * (1.-h) * (1.-r)

          dN(1,3) = 1./8.* (1.) * (1.+h) * (1.-r)

          dN(1,4) = 1./8.* (-1.) * (1.+h) * (1.-r)

          dN(1,5) = 1./8.* (-1.) * (1.-h) * (1.+r)

          dN(1,6) = 1./8.* (1.) * (1.-h) * (1.+r) 
          
          dN(1,7) = 1./8.* (1.) * (1.+h) * (1.+r)

          dN(1,8) = 1./8.* (-1.) * (1.+h) * (1.+r)
      
          
c         dN_dh       
          
          dN(2,1) = 1./8.* (1.-g) * (-1.) * (1.-r)

          dN(2,2) = 1./8.* (1.+g) * (-1.) * (1.-r)

          dN(2,3) = 1./8.* (1.+g) * (1.) * (1.-r)

          dN(2,4) = 1./8.* (1.-g) * (1.) * (1.-r)


          dN(2,5) = 1./8.* (1.-g) * (-1.) * (1.+r)

          dN(2,6) = 1./8.* (1.+g) * (-1.) * (1.+r)

          dN(2,7) = 1./8.* (1.+g) * (1.) * (1.+r)

          dN(2,8) = 1./8.* (1.-g) * (1.) * (1.+r)
          
          
          
c         dN_dr     
          
          dN(3,1) = 1./8.* (1.-g) * (1.-h) * (-1.) ;

          dN(3,2) = 1./8.* (1.+g) * (1.-h) * (-1.) ;

          dN(3,3) = 1./8.* (1.+g) * (1.+h) * (-1.) ;

          dN(3,4) = 1./8.* (1.-g) * (1.+h) * (-1.) ;


          dN(3,5) = 1./8.* (1.-g) * (1.-h) * (1.) ;

          dN(3,6) = 1./8.* (1.+g) * (1.-h) * (1.) ;

          dN(3,7) = 1./8.* (1.+g) * (1.+h) * (1.) ;

          dN(3,8) = 1./8.* (1.-g) * (1.+h) * (1.) ;         
          
      
      endif
      
      
      end subroutine shapefunctionderivatives
      
      
      
      
      
      
      
      
      
      
      
      
      end module gndslipgrad