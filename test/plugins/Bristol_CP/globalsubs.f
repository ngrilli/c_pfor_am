c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c
c Aug. 12th, 2021 - 1st working version
c
c This block stores the global variables & symmetry operators 
c used in misorientation calculation
c
	module globalsubs
      implicit none
	contains
c
c     This subroutine call all the other initialization routines
c	and assigns initial orientations to every integration point
c           
c    
c
c
c
c
c
c	This subroutine calculates the misorientation
c	INPUTS:
c	g1:		1st orientation 
c	g2:		2nd orientation
c	typ:	Type of crystal symmetry; 1:cubic, 2:tetragonal, 3:hexagonal
c	OUTPUTS:
c	uvw:	Misorientation axis
c	ang:	Misorientation angle (in degrees)
c	dg:		net rotation
c	USES:
c	Vars:	cubsym,tetsym,hcpsym
c	Consts:	pi
	subroutine misorientation(g1,g2,typ,uvw,ang,dg)
      use globalvars, only: cubsym,tetsym,hcpsym,pi,I3,smallnum
	implicit none
	real(8) g1(3,3),g2(3,3),uvw(3),ang,dg(3,3)
	integer typ,nosym,s1,s2
	real(8) g11(3,3),g22(3,3),dg1(3,3),dg2(3,3),sq,x1,x2,x3
	real(8) trace,temp_ang
	real(8) sym(24,3,3)
      real(8) sym1(3,3),sym2(3,3) 



	if (typ.eq.1d+0) then
		nosym=24d+0
		sym=cubsym
	elseif (typ.eq.2d+0) then
		nosym=8d+0
		sym(1:nosym,:,:)=tetsym
	elseif (typ.eq.3d+0) then
		nosym=12d+0
		sym(1:nosym,:,:)=hcpsym
	endif



      
c	Initially set the angle to be maximum (62.4 degrees indeed)
	ang=pi/2.0d+0

c	1: Crystal Symmetry to 1st orientation
	do s1=1,nosym
          sym1=sym(s1,:,:)
          !do i=1,3
          !    write(1,*) (s(i,j), j=1,3)
          !enddo
		g11=matmul(sym1, g1)
c	2: Crystal Symmetry to 2nd orientation
		do s2=1,nosym
              sym2 = sym(s2,:,:)
			g22=matmul(sym2,g2)
c	3: Net rotations (g2*inv(g1)=g2*transpose(g1))
			dg1=matmul(g22,transpose(g11))										
			dg2=matmul(g11,transpose(g22))		

     
c	5. CHECK IF THE MISORIENTATION IS IN FUNDAMENTAL ZONE
c	5.a. For 1st result: (dg1)
      sq=dsqrt((dg1(2,3)-dg1(3,2))**2.0d+0 + 
     &(dg1(1,3)-dg1(3,1))**2.0d+0 + (dg1(2,1)-dg1(1,2))**2.0d+0)
c      write(6,*) 'sq', sq
      if (dabs(sq).gt.smallnum) then
		x1=(dg1(2,3)-dg1(3,2))/sq
		x2=(dg1(3,1)-dg1(1,3))/sq
		x3=(dg1(1,2)-dg1(2,1))/sq
		if (x1.ge.0.0d+0) then
			if (x2.ge.0.0d+0) then
				if (x3.ge.0.0d+0) then
					if (x1.le.x2) then
						if (x2.le.x3) then
							trace=dg1(1,1)+dg1(2,2)+dg1(3,3)
							temp_ang=dacos((trace-1.0d+0)/2.0d+0)
							if (dabs(temp_ang).lt.ang) then
								ang=temp_ang
								dg=dg1
								uvw(1)=x1
                                  uvw(2)=x2
                                  uvw(3)=x3
							endif
						endif
					endif
				endif
			endif
            endif		
              
      endif
          

 
c	5.a. For the 2nd result: (dg2)
     	sq=dsqrt((dg2(2,3)-dg2(3,2))**2.0d+0 + 
     &(dg2(1,3)-dg2(3,1))**2.0d+0 + (dg2(2,1)-dg2(1,2))**2.0d+0)
c      write(6,*) 'sq', sq
      if (dabs(sq).gt.smallnum) then
   	    x1=(dg2(2,3)-dg2(3,2))/sq
	    x2=(dg2(3,1)-dg2(1,3))/sq
	    x3=(dg2(1,2)-dg2(2,1))/sq
	    if (x1.ge.0.0d+0) then
		    if (x2.ge.0.0d+0) then
			    if (x3.ge.0.0d+0) then
				    if (x1.le.x2) then
					    if (x2.le.x3) then
						    trace=dg2(1,1)+dg2(2,2)+dg2(3,3)
						    temp_ang=dacos((trace-1.0d+0)/2.0d+0)
						    if (dabs(temp_ang).lt.ang) then
							    ang=temp_ang
							    dg=dg2
							    uvw(1)=x1
                                  uvw(2)=x2
                                  uvw(3)=x3
						    endif
					    endif
				    endif
			    endif
		    endif
            endif       
      endif
              
 		enddo
      enddo
	
          
    

c     Check for zero misorientation      
      if (dabs(ang-pi/2.0).gt.smallnum) then
          ang = 0.0d+0
          uvw=0.0d+0
c         Assign a unit vector
          uvw(1)=1.0d+0
          dg=I3
      endif
       
c     Convert the misorientation angle into degrees      
	ang=ang*180.0d+0/pi
      
      

      
	return
	end subroutine misorientation
c
c	This subroutine calculates the Von-Mises stress
c	INPUTS: Stress						---	sigma(3,3)
c	OUTPUT: Von-Mises stress (scalar)	--- vms
	subroutine vonmises_stress(sigma,vms)
	use globalvars, only: I3
      implicit none
	real(8) sigma(3,3),vms,dev(3,3),hyd
	integer i,j
	real(8) sum
      
      hyd=0.0d+0
      do i=1,3
          hyd=hyd + sigma(i,i)
      enddo
      
      dev = sigma - hyd/3.0d+0*I3
      
c
c	Take the inner product
	sum=0.0d+0
	do i=1,3
		do j=1,3
			sum=sum+(3.0d+0*dev(i,j)*dev(i,j)/2.0d+0)
		enddo
	enddo
c	Find the Von-Mises stress
	vms=dsqrt(sum)
	return
	end subroutine vonmises_stress
c
c	This subroutine calculates the Von-Mises strain
c	INPUTS: Strain						---	eps(3,3)
c	OUTPUT: Von-Mises strain (scalar)	--- evm
	subroutine vonmises_strain(eps,evm)	
	implicit none
	real(8) eps(3,3),evm
	integer i,j
	real(8) sum
c

      


c	Take the inner product
	sum=0.0d+0
	do i=1,3
		do j=1,3
			sum=sum+(2.0d+0*eps(i,j)*eps(i,j)/3.0d+0)
		enddo
	enddo
c	Find the Von-Mises strain
	evm=dsqrt(sum)
	return 
	end subroutine vonmises_strain
c
c	This subroutine calculates the Taylor factor
c	INPUTS:
c	F: Deformation gradient at a certain reference 
c	g: Transformation matrix from the reference frame at which F is defined to the crystal frame
c	OUTPUT: 
c	M: Taylor factor
c	USES:
c	BHstress: Bishop Hill stress states
	subroutine taylorfactor(F,g,M)
	use globalvars, only: BHstress
	implicit none
	real(8) F(3,3),g(3,3),M
	integer i,j,k,l
	real(8) Fc(3,3),sum,evm,st(3,3),st_vec(6),w,max_w
c
c	Find the Von-Mises strain
	call vonmises_strain(F,evm)
	if (evm.eq.0d+0) then
		M=3.0d+0
		return
	endif
c	Transform the deformation gradient to crystal frame	
	do i=1,3
		do j=1,3
			sum=0
			do k=1,3
				do l=1,3
					sum=sum+g(i,k)*g(j,l)*F(k,l)
				enddo
			enddo
			Fc(i,j)=sum
		enddo
	enddo
c	Strain part of the transformed deformation
	st=(Fc+transpose(Fc))/2.0
c	Vectorizing stress (Convention used is different than in MARC)
c	Check Read's "Deformation Geometry" for details
	st_vec(1)=st(2,2)
	st_vec(2)=-st(1,1)
	st_vec(3)=st(3,3)
	st_vec(4)=st(1,2)+st(2,1)
	st_vec(5)=st(3,1)+st(1,3)
      st_vec(6)=st(2,3)+st(3,2)
c	Calculating work 
	max_w=0.0d+0
	do i=1,28
		w=0.0
		do j=1,6
			w=w+BHstress(i,j)*st_vec(j)
		enddo
		if (abs(w)>max_w) then
			max_w=abs(w)
		endif
	enddo
	M=sqrt(6.0d+0)*max_w/evm;
	return
	end subroutine taylorfactor
c
c
c	This subroutine calculates the Macauley function
c	INPUTS: scalar              		--- A
c	OUTPUT: Macauley function value		--- B
	subroutine Macauley(A,B)
	implicit none
	real(8) A, B
c
	if (A.lt.0.0d+0) then
          B = 0.0d+0
c     A >= 0          
      else
          B = A
      endif
      
	return
	end subroutine Macauley
c
c
c
c	This subroutine converts a symmetric 2nd order tensor to a vector
c	INPUTS: Matrix			--- s(3,3)
c	OUTPUT: Vectorized form	--- v(6)	
	subroutine convert3x3to6(s,v)
      use globalvars, only : order6x6
	implicit none
	real(8) s(3,3),v(6)
      integer i
c	
	
      do i=1,6

          v(i) = (s(order6x6(i,1), order6x6(i,2)) +
     &            s(order6x6(i,2), order6x6(i,1)) )/2.0d+0
          
      enddo
      
      
      
      
	return
	end subroutine convert3x3to6
c
c
c
c	This subroutine converts a vector to a SYMMETRIC 2nd order tensor
c	INPUTS: Vectorized form		--- v(6)
c	OUTPUT: Matrix				--- s(3,3)
	subroutine convert6to3x3(v,s)
      use globalvars, only : order6x6
	implicit none
	real(8) s(3,3),v(6)
      integer i
c

      do i=1,6

          s(order6x6(i,1),order6x6(i,2)) = v(i)
          
          s(order6x6(i,2),order6x6(i,1)) = v(i)
          
      enddo


	return
	end subroutine convert6to3x3
c
c	This subroutine is written to convert 4th order elasticity tensor to 6x6 matrix
c	INPUTS : 4th order tensor		--- C3333(3,3,3,3)
c	OUTPUTS: Matrix					---	C66(6,6)
	subroutine convert3x3x3x3to6x6(C3333,C66)
      use globalvars, only : order6x6
	implicit none
	real(8) C66(6,6), C3333(3,3,3,3)
	integer i,j

      
      C66 = 0.0d+0
      do i=1,6
          do j=1,6
              
              C66(i,j)=C3333(order6x6(i,1), order6x6(i,2), 
     &order6x6(j,1), order6x6(j,2))
              
          enddo
      enddo
c     
c
c
	return
	end subroutine convert3x3x3x3to6x6
c
c
c	This subroutine is written to convert 6x6 matrix to 4th order elasticity tensor
c	INPUTS : Matrix					---	C66(6,6)
c	OUTPUTS: 4th order tensor		--- C3333(3,3,3,3)
	subroutine convert6x6to3x3x3x3(C66,C3333)
      use globalvars, only : order6x6
	implicit none
	real(8) C66(6,6), C3333(3,3,3,3)
	integer i,j
c      
c      
c      
      C3333=0.0d+0
      do i=1,6
          do j=1,6
              
              C3333(  order6x6(i,1), order6x6(i,2), 
     &                order6x6(j,1), order6x6(j,2)  ) = C66(i,j)
              
              
              C3333(  order6x6(i,2), order6x6(i,1), 
     &                order6x6(j,1), order6x6(j,2)  ) = C66(i,j)
              
              
              C3333(  order6x6(i,1), order6x6(i,2), 
     &                order6x6(j,2), order6x6(j,1)  ) = C66(i,j)
              
              
              C3333(  order6x6(i,2), order6x6(i,1), 
     &                order6x6(j,2), order6x6(j,1)  ) = C66(i,j)     
              
              
          enddo
      enddo
      
	return
	end subroutine convert6x6to3x3x3x3
c     
c
c	This subroutine converts a 4th order tensor to 9x9 matrix
c	INPUTS: 4th order tensor	--- r(3,3,3,3)
c	OUTPUT: Matrix				--- s(9,9)
c	USES  :	Order of indices	---	order9x9(9,2)
	subroutine convert3x3x3x3to9x9(r,s)
	use globalvars, only : order9x9
	implicit none
	real(8) r(3,3,3,3), s(9,9)
	integer i,j,a,b,c,d
c	
	do i=1,9
		a=order9x9(i,1)
		b=order9x9(i,2)
		do j=1,9
			c=order9x9(j,1)
			d=order9x9(j,2)
			s(i,j)=r(a,b,c,d)
		enddo
	enddo
	return
      end subroutine convert3x3x3x3to9x9
c      
c	This subroutine converts a 9x9 matrix to 4th order tensor
c	INPUTS: 4th order tensor	--- s(9,9)
c	OUTPUT: Matrix				--- r(3,3,3,3)
c	USES  :	Order of indices	---	order9x9(9,2)
	subroutine convert9x9to3x3x3x3(s,r)
	use globalvars, only : order9x9
	implicit none
	real(8) r(3,3,3,3), s(9,9)
	integer i,j,a,b,c,d
c	
	do i=1,9
		a=order9x9(i,1)
		b=order9x9(i,2)
		do j=1,9
			c=order9x9(j,1)
			d=order9x9(j,2)
			r(a,b,c,d) = s(i,j)
		enddo
	enddo
	return
	end subroutine convert9x9to3x3x3x3      
c
c	This subroutine converts a vectorized 2nd order tensor to a matrix 
c	INPUTS: Vectorized form		--- s(9)
c	OUTPUT: Matrix				--- r(3,3)
c	USES  :	Order of indices	---	order9x9(9,2)
	subroutine convert9to3x3(s,r)
	use globalvars, only : order9x9
	implicit none
	real(8) r(3,3), s(9)
	integer i,a,b
c	
	do i=1,9
		a=order9x9(i,1)
		b=order9x9(i,2)
		r(a,b)=s(i)
	enddo
	return
	end subroutine convert9to3x3
c
c	This subroutine converts a matrix to a vectorized form
c	INPUTS: Matrix				--- r(3,3)
c	OUTPUT: Vectorized form		--- s(9)
c	USES  :	Order of indices	---	order9x9(9,2)
	subroutine convert3x3to9(r,s)
	use globalvars, only : order9x9
	implicit none
	real(8) r(3,3), s(9)
	integer i,a,b
c	
	do i=1,9
		a=order9x9(i,1)
		b=order9x9(i,2)
		s(i)=r(a,b)
	enddo
	return
	end subroutine convert3x3to9
c
c
c	This subroutine transforms a 2nd order tensor using a given transformation matrix
c	INPUTS : Transformation matrix, 2nd order tensor	--- g, A33
c	OUTPUTS: Transformed 2nd order tensor				--- A33_tr
	subroutine transform2(g,A33,A33_tr)
	implicit none
	real(8) g(3,3),A33_tr(3,3),A33(3,3)
	real(8) sum
	integer i,j,m,n
c
	do i=1,3
		do j=1,3
			sum=0.0
			do m=1,3
				do n=1,3
					sum=sum+
     &                     (g(i,m)*g(j,n)*A33(m,n))
				enddo	
			enddo
			A33_tr(i,j)=sum
		enddo
	enddo
	return
	end subroutine transform2
c
c	This subroutine transforms a 4th order tensor using a given transformation matrix
c	INPUTS : Transformation matrix, 4th order tensor	--- g, A3333
c	OUTPUTS: Transformed 4th order tensor				--- A3333_tr
	subroutine transform4(g,A3333,A3333_tr)
	implicit none
	real(8) g(3,3),A3333(3,3,3,3),A3333_tr(3,3,3,3)
	real(8) sum
	integer i,j,k,l,m,n,o,p
c		

	do i=1,3
		do j=1,3
			do k=1,3
				do l=1,3
					sum=0.0
					do m=1,3
						do n=1,3
							do o=1,3
								do p=1,3
									sum=sum+
     &                    (g(i,m)*g(j,n)*g(k,o)*g(l,p)*A3333(m,n,o,p))
								enddo
							enddo
						enddo	
					enddo
					A3333_tr(i,j,k,l)=sum
				enddo
			enddo
		enddo
      enddo
      
    
      
      
	return
	end subroutine transform4
c
c	This subroutine calculates the square norm of a matrix
c	INPUTS: Matrix, dimension(n)		--- A(nxn),dim
c	OUTPUT: Square norm				--- norm
	subroutine normmat(A,dim,norm)
	implicit none
	integer dim
	real(8) A(dim,dim),norm
	integer i,j
c
	norm=0.0d+0
	do i=1,dim
		do j=1,dim
			norm=norm+(A(i,j)*A(i,j))
		enddo
	enddo
	norm=dsqrt(norm)
	return
	end subroutine normmat
c
c	This subroutine calculates the trace of a matrix
c	INPUTS: Matrix, dimension(n)		--- A(nxn),dim
c	OUTPUT: Trace of a matrix			--- value
	subroutine trace(A,dim,value)
	implicit none
	integer dim
	real(8) A(dim,dim),value
	integer i
c
	value=0.0d+0
	do i=1,dim
		value=value+A(i,i)
	enddo
	return
	end subroutine trace
c
c	This subroutine calculates the square norm of a column vector
c	INPUTS: Vector, dimension(n)		--- A(nx1),dim
c	OUTPUT: Square norm				--- norm
	subroutine normvec(A,dim,norm)
	implicit none
	integer i,dim
	real(8) A(dim),norm
c
	norm=0.0d+0
	do i=1,dim
		norm=norm+(A(i)*A(i))
	enddo
	norm=dsqrt(norm)
	return
	end subroutine normvec
c
c	This subroutine multiplies a matrix with a matrix
c	INPUTS: Matrix (nxn), vector (nx1), dimension (n)	--- A(nxn),B(n,n),dim
c	OUTPUT: Result										--- C(nxn)
	subroutine matmulmat(A,B,dim,C)
	implicit none
	integer i,j,k,dim
	real(8) A(dim,dim),B(dim,dim),C(dim,dim)
c
	C=0.0d+0
	do i=1,dim
		do j=1,dim
			do k=1,dim
				C(i,j)=C(i,j)+(A(i,k)*B(k,j))
			enddo
		enddo
	enddo
	return
      end subroutine matmulmat
c
c
c	This subroutine computes the cross product of two vectors
c	INPUTS: vector (3), vector (3)
c	OUTPUT: vector(3)
	subroutine cross(a,b,c)
      use globalvars, only: eijk
	implicit none
	integer i,j,k
	real(8) a(3),b(3),c(3)
c
	c=0.0d+0
	do i=1,3
		do j=1,3
			do k=1,3
                  c(i) = c(i) + (eijk(i,j,k)*a(j)*b(k))
              enddo
		enddo
	enddo
	return
	end subroutine cross     
c
c	This subroutine multiplies a matrix with a vector
c	INPUTS: Matrix (nxn), vector (nx1), dimension (n)	--- A(nxn),b(n,1),dim
c	OUTPUT: Result										---	c(n,1)
	subroutine matmulvec(A,b,dim,c)
	implicit none
	integer i,j,k,dim
	real(8) A(dim,dim),b(dim),c(dim)
c
	c=0.0d+0
	do i=1,dim
		do j=1,dim
			c(i)=c(i)+(A(i,j)*b(j))
		enddo
	enddo
	return
	end subroutine matmulvec
c
c	This subroutine is written to check the determinant
c	INPUTS: 3x3 matrix	--- a(3,3)
c	OUTPUT: Determinant	--- det
	subroutine determinant(a,det)
	implicit none
	real(8) a(3,3),det
	real(8) v1,v2,v3
c	
	v1=a(1,1)*((a(2,2)*a(3,3))-(a(2,3)*a(3,2)))
	v2=-a(1,2)*((a(2,1)*a(3,3))-(a(2,3)*a(3,1)))
	v3=a(1,3)*((a(2,1)*a(3,2))-(a(2,2)*a(3,1)))
	det=v1+v2+v3
	return
	end subroutine determinant
c
c	This subroutine inverts a 3x3 matrix
c	INPUT:	Matrix								---	A(3,3)
c	OUTPUT:	Invereted matrix, determinant		---	invA(3,3),det
	subroutine invert3x3(A,invA,det)
      use globalvars, only: smallnum
	implicit none
	integer i,j
	real(8) A(3,3),invA(3,3),det
c
c	First calculate the determinant
	call determinant(A,det)
c	If the determinant is greater than certain value
	if (det.lt.smallnum) then
		invA=0.0d+0
	else
		invA(1,1)=((A(2,2)*A(3,3))-(A(2,3)*A(3,2)))/det
		invA(2,1)=-((A(2,1)*A(3,3))-(A(2,3)*A(3,1)))/det
		invA(3,1)=((A(2,1)*A(3,2))-(A(2,2)*A(3,1)))/det
		invA(1,2)=-((A(1,2)*A(3,3))-(A(1,3)*A(3,2)))/det
		invA(2,2)=((A(1,1)*A(3,3))-(A(1,3)*A(3,1)))/det
		invA(3,2)=-((A(1,1)*A(3,2))-(A(1,2)*A(3,1)))/det
		invA(1,3)=((A(1,2)*A(2,3))-(A(1,3)*A(2,2)))/det
		invA(2,3)=-((A(1,1)*A(2,3))-(A(2,1)*A(1,3)))/det
		invA(3,3)=((A(1,1)*A(2,2))-(A(1,2)*A(2,1)))/det
	endif
	return
	end subroutine invert3x3
c
c	This subroutine orientation matrix from Euler angles
c	INPUT:	Angles(deg)			---	ang(3)
c	OUTPUT:	Orientation matrix	---	R(3,3)
c     USES:     Number pi           --- pi
	subroutine ang2ori(ang,R)
	use globalvars, only : pi
	implicit none
	real(8) ang(3),R(3,3)
	real(8) phi1,phi2,PHI
c
	phi1=ang(1)*pi/180.0d+0
	phi2=ang(3)*pi/180.0d+0
	PHI=ang(2)*pi/180.0d+0
      R(1,1)=(dcos(phi1)*dcos(phi2))-(dsin(phi1)*dsin(phi2)*dcos(PHI))
      R(2,1)=-(dcos(phi1)*dsin(phi2))-(dsin(phi1)*dcos(phi2)*dcos(PHI))
      R(3,1)=dsin(phi1)*dsin(PHI)
      R(1,2)=(dsin(phi1)*dcos(phi2))+(dcos(phi1)*dsin(phi2)*dcos(PHI))
      R(2,2)=-(dsin(phi1)*dsin(phi2))+(dcos(phi1)*dcos(phi2)*dcos(PHI))
      R(3,2)=-dcos(phi1)*dsin(PHI)
      R(1,3)=dsin(phi2)*dsin(PHI)
      R(2,3)=dcos(phi2)*dsin(PHI)
      R(3,3)=dcos(PHI)	
	return
	end subroutine ang2ori
c
c	This subroutine orientation matrix from Euler angles
c	INPUT:	Orientation matrix	---	R(3)
c	OUTPUT:	Angles(deg)			---	ang(3)
c     USES:     Number pi           --- pi
	subroutine ori2ang(R,ang)
	use globalvars, only : pi
	implicit none
	real(8) ang(3),R(3,3)
	real(8) phi1,phi2,PHI
c
	if (R(3,3).eq.1) then
		PHI=0.0d+0
		phi1=atan2(R(1,2),R(1,1))
		phi2=0.0d+0
	else
		PHI=dacos(R(3,3))
		phi1=atan2(R(3,1)/dsin(PHI),-R(3,2)/dsin(PHI))
		phi2=atan2(R(1,3)/dsin(PHI),R(2,3)/dsin(PHI))
	endif
	ang=ang*180.0d+0/pi
	return
	end subroutine ori2ang
c
c	This subroutine calculates angle and axis pair from orientation matrix
c	INPUT:	Orientation matrix	--- R(3,3)
c	OUTPUT:	Angle(deg), axis	---	ang, ax(3)
c     USES:     Number pi           --- pi
	subroutine ori2angax(R,ang,ax)
	use globalvars, only : pi
	implicit none
	real(8) R(3,3),ang,ax(3)
	real(8) trace,norm,t,mag
c
	norm=dsqrt((R(2,3)-R(3,2))**2.0d+0 + (R(1,3)-R(3,1))**2.0d+0 +
     &(R(1,2)-R(2,1))**2.0d+0)
	trace=R(1,1)+R(2,2)+R(3,3)
	ang=2.0d+0*acos(dsqrt(1.0d+0+trace)/2.0d+0)
	t=dtan(ang/2.0d+0)
	ax(1)=t*(R(2,3)-R(3,2))/norm
	ax(2)=t*(R(3,1)-R(1,3))/norm
	ax(3)=t*(R(1,2)-R(2,1))/norm	
c	Normalize the axis vector
	mag=dsqrt(ax(1)**2.0d+0 + ax(2)**2.0d+0 + ax(3)**2.0d+0)
	ax=ax/mag		
	ang=ang*180.0d+0/pi
	return
	end subroutine ori2angax
c
c	This subroutine finds the symmetric part of a matirx
c	INPUT: 2nd order tensor	--- F(3,3)
c	OUPUT:	Symmetric part	--- S(3,3)
	subroutine sym(F,S)
	implicit none
	real(8) F(3,3),S(3,3)
	integer i,j
c
	do i=1,3
		do j=1,3
			S(i,j)=0.5d+0*(F(i,j)+F(j,i))
		enddo
	enddo
	return
	end subroutine sym
c
c	This subroutine finds the skew-symmetric part of a matirx
c	INPUT: 2nd order tensor	--- F(3,3)
c	OUPUT:	Symmetric part	--- S(3,3)
	subroutine skew(F,S)
	implicit none
	real(8) F(3,3),S(3,3)
	integer i,j
c
	do i=1,3
		do j=1,3
			S(i,j)=0.5d+0*(F(i,j)-F(j,i))
		enddo
	enddo
	return
	end subroutine skew
c
c
c
c    
c      
c      
c********1*********2*********3*********4*********5*********6*********7**
c
c	MATPP3(G,R,S)
c
c	positive polar decomposition of 3x3 matrix  G = R * S
c	forcing positive orthonormal rotation matrix
c
c	INPUTS
c	G = 3x3 general matrix
c
c	OUTPUTS
c	R = 3x3 positive orthonormal matrix
c	S = 3x3 symmetric matrix
c
c	PRECISION:	single
c	COMMONS:	none
c	CALLS:		none
c	FUNCTIONS:	ABS, SQRT
c	REFERENCE:	Veldpaus, F.E., H.J. Woltring, and L.J.M.G. Dortmans,
c				A Least-Squares Algorithm for the Equiform
c				Transformation from Spatial Marker Coordinates, 
c				J. Biomechanics, 21(1):45-54 (1988).
c	DATE:		10/8/92 - HJSIII
c
c
      SUBROUTINE polar(G,R,S)
c
      use globalvars, only: I3, smallnum
      implicit none
c	declarations
      REAL(8) G(3,3),R(3,3),S(3,3)
      REAL(8) COG(3,3),P(3,3),ADP(3,3),PBI(3,3)
      real(8) EPS, G1SQ, G1, G2SQ, G2, G3, H1, H2, X, Y
      real(8) DEN, RES1, RES2, DX, DY, BETA1, BETA2, DETPBI
      integer I, j, k, flag1, flag2
c
      flag1=0d+0
      flag2=0d+0
c     If G is identity
      do j=1,3
          do k=1,3
              if (dabs(G(j,k)-I3(j,k)).lt.smallnum) then
                  flag1 = 1d+0
              endif
          enddo
      enddo
      
      
      if (flag1.eq.0d+0) then




c	    constants
          EPS=1.0d-5
c
c	    cofactors and determinant of g
          COG(1,1)=G(2,2)*G(3,3)-G(2,3)*G(3,2)
          COG(2,1)=G(1,3)*G(3,2)-G(1,2)*G(3,3)
          COG(3,1)=G(1,2)*G(2,3)-G(1,3)*G(2,2)
          COG(1,2)=G(2,3)*G(3,1)-G(2,1)*G(3,3)
          COG(2,2)=G(1,1)*G(3,3)-G(1,3)*G(3,1)
          COG(3,2)=G(1,3)*G(2,1)-G(1,1)*G(2,3)
          COG(1,3)=G(2,1)*G(3,2)-G(2,2)*G(3,1)
          COG(2,3)=G(1,2)*G(3,1)-G(1,1)*G(3,2)
          COG(3,3)=G(1,1)*G(2,2)-G(1,2)*G(2,1)
          G3=G(1,1)*COG(1,1)+G(2,1)*COG(2,1)+G(3,1)*COG(3,1)
c
c	    P = trans(G) * G = S * S
          DO I=1,3
              P(I,1)=G(1,I)*G(1,1)+G(2,I)*G(2,1)+G(3,I)*G(3,1)
              P(I,2)=G(1,I)*G(1,2)+G(2,I)*G(2,2)+G(3,I)*G(3,2)
              P(I,3)=G(1,I)*G(1,3)+G(2,I)*G(2,3)+G(3,I)*G(3,3)
          ENDDO
      
      
c         If P is identity
          do j=1,3
              do k=1,3
                  if (dabs(P(j,k)-I3(j,k)).lt.smallnum) then
                      flag2 = 1
                  endif
              enddo
          enddo

c         Continue if non identity  
          if (flag2.eq.0d+0) then
      
c
c	        adjoint of P
              ADP(1,1)=P(2,2)*P(3,3)-P(2,3)*P(3,2)
              ADP(2,2)=P(1,1)*P(3,3)-P(1,3)*P(3,1)
              ADP(3,3)=P(1,1)*P(2,2)-P(1,2)*P(2,1)
c
c	        G invariants
              G1SQ=P(1,1)+P(2,2)+P(3,3)
              G1=SQRT(G1SQ)
              G2SQ=ADP(1,1)+ADP(2,2)+ADP(3,3)
              G2=SQRT(G2SQ)
c
c	        initialize iteration
              H1=G2/G1SQ
              H2=G3*G1/G2SQ
              X=1.0d+0
              Y=1.0d+0
c
c	        iteration loop
              DO WHILE (ABS(DX/X).GT.EPS.OR.ABS(DY/Y).GT.EPS)
                  DEN=2.0d+0*(X*Y-H1*H2)
                  RES1=1.0d+0 - X*X + 2.0d+0*H1*Y
                  RES2=1.0d+0 - Y*Y + 2.0d+0*H2*X
                  DX=(Y*RES1+H1*RES2)/DEN
                  DY=(H2*RES1+X*RES2)/DEN
                  X=X+DX
                  Y=Y+DY
              ENDDO
c
c	        BETA invariants
              BETA1=X*G1
              BETA2=Y*G2
c
c	        invert ( trans(G) * G + BETA2 * identity )
              P(1,1)=P(1,1)+BETA2
              P(2,2)=P(2,2)+BETA2
              P(3,3)=P(3,3)+BETA2
              PBI(1,1)=P(2,2)*P(3,3)-P(2,3)*P(3,2)
              PBI(1,2)=P(1,3)*P(3,2)-P(1,2)*P(3,3)
              PBI(1,3)=P(1,2)*P(2,3)-P(1,3)*P(2,2)
              PBI(2,1)=P(2,3)*P(3,1)-P(2,1)*P(3,3)
              PBI(2,2)=P(1,1)*P(3,3)-P(1,3)*P(3,1)
              PBI(2,3)=P(1,3)*P(2,1)-P(1,1)*P(2,3)
              PBI(3,1)=P(2,1)*P(3,2)-P(2,2)*P(3,1)
              PBI(3,2)=P(1,2)*P(3,1)-P(1,1)*P(3,2)
              PBI(3,3)=P(1,1)*P(2,2)-P(1,2)*P(2,1)
              DETPBI=P(1,1)*PBI(1,1)+P(2,1)*PBI(1,2)+P(3,1)*PBI(1,3)
c
c	        R = (cofac(G)+BETA1*G) * inv(trans(G)*G+BETA2*identity)
              DO I=1,3
                  R(I,1)=((COG(I,1)+BETA1*G(I,1))*PBI(1,1)
     &       +(COG(I,2)+BETA1*G(I,2))*PBI(2,1)
     &       +(COG(I,3)+BETA1*G(I,3))*PBI(3,1))/DETPBI
                  R(I,2)=((COG(I,1)+BETA1*G(I,1))*PBI(1,2)
     &       +(COG(I,2)+BETA1*G(I,2))*PBI(2,2)
     &       +(COG(I,3)+BETA1*G(I,3))*PBI(3,2))/DETPBI
                  R(I,3)=((COG(I,1)+BETA1*G(I,1))*PBI(1,3)
     &       +(COG(I,2)+BETA1*G(I,2))*PBI(2,3)
     &       +(COG(I,3)+BETA1*G(I,3))*PBI(3,3))/DETPBI
              ENDDO
c
c	        S = trans(R) * G
              DO I=1,3
                  S(I,1)=R(1,I)*G(1,1)+R(2,I)*G(2,1)+R(3,I)*G(3,1)
                  S(I,2)=R(1,I)*G(1,2)+R(2,I)*G(2,2)+R(3,I)*G(3,2)
                  S(I,3)=R(1,I)*G(1,3)+R(2,I)*G(2,3)+R(3,I)*G(3,3)
              ENDDO
              
          else
              
              S = I3
              R = G
              
          endif
          
      else
          
          S = I3
          R = I3
          
      endif
      
              
c
c	done
      RETURN
      END SUBROUTINE polar
c
c********1*********2*********3*********4*********5*********6*********7**      
c      
c      
c
c	This subroutine calculates the L-U decomp.
c	INPUT: Matrix, dimension					--- F(n,n), n
c	OUTPUT: Lower tri., Upper tri., determinant --- L(n,n), U(n,n), det
	subroutine LUdecomp(A,n,L,U,det)
	implicit none
	integer n,i,j,k
	real(8) A(n,n),L(n,n),U(n,n)
	real(8) factor(n,n), det
c
	U=A
	factor=0.0d+0
c	Triangularization
	do k=1,n-1
		do i=k+1,n
			factor(i,k)=U(i,k)/U(k,k)
			do j=1,n
				U(i,j)=U(i,j)-factor(i,k)*U(k,j)
			enddo	
		enddo
	enddo
c	Lower triangular matrix
	L=factor
	do i=1,n
		L(i,i)=1.0d+0
	enddo
c	Multiplication of the elements at the diagonal is the determinant
	det=1;
	do k=1,n
		det=det*L(k,k)*U(k,k)
	enddo
	return
	end subroutine LUdecomp
c
c
      subroutine invertnxn(a,c,n)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !===========================================================
      implicit none 
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k

      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0d+0
      U=0.0d+0
      b=0.0d+0

      ! step 1: forward elimination
      do k=1, n-1
          do i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                  a(i,j) = a(i,j)-coeff*a(k,j)
              end do
          end do
      end do

      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
      L(i,i) = 1.0d+0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
          do i=1,j
              U(i,j) = a(i,j)
          end do
      end do

      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
          b(k)=1.0d+0
          d(1) = b(1)
          ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
              d(i)=b(i)
              do j=1,i-1
                  d(i) = d(i) - L(i,j)*d(j)
              end do
          end do
          ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
              x(i) = d(i)
              do j=n,i+1,-1
                  x(i)=x(i)-U(i,j)*x(j)
              end do
              x(i) = x(i)/u(i,i)
          end do
          ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
              c(i,k) = x(i)
          end do
          b(k)=0.0d+0
      end do
      return
      end subroutine invertnxn
      
      
      
      

      

      end module globalsubs
