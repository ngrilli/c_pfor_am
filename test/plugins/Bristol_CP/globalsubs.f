! Chris Allen
! Edward Horton
! Eralp Demir
! Aug. 12th, 2021 - 1st working version
!
! This block stores the global variables & symmetry operators 
! used in misorientation calculation
!
	module globalsubs
      implicit none
	contains
!
!     This subroutine call all the other initialization routines
!	and assigns initial orientations to every integration point
           
    
!
!	This subroutine calculates the misorientation
!	INPUTS:
!	g1:		1st orientation 
!	g2:		2nd orientation
!	typ:	Type of crystal symmetry; 1:cubic, 2:tetragonal, 3:hexagonal
!	OUTPUTS:
!	uvw:	Misorientation axis
!	ang:	Misorientation angle (in degrees)
!	dg:		net rotation
!	USES:
!	Vars:	cubsym,tetsym,hcpsym
!	Consts:	pi
	subroutine misorientation(g1,g2,typ,uvw,ang,dg)
      use globalvars, only: cubsym,tetsym,hcpsym,pi,I3
	implicit none
	real(8) g1(3,3),g2(3,3),uvw(3),ang,dg(3,3)
	integer typ,nosym,s1,s2
	real(8) g11(3,3),g22(3,3),dg1(3,3),dg2(3,3),sq,x1,x2,x3
	real(8) trace,temp_ang
	real(8) sym(24,3,3)
      real(8) sym1(3,3),sym2(3,3) 



	if (typ.eq.1) then
		nosym=24
		sym=cubsym
	elseif (typ.eq.2) then
		nosym=8
		sym(1:nosym,:,:)=tetsym
	elseif (typ.eq.3) then
		nosym=12
		sym(1:nosym,:,:)=hcpsym
	endif



      
!	Initially set the angle to be maximum (62.4 degrees indeed)
	ang=pi/2.0

!	1: Crystal Symmetry to 1st orientation
	do s1=1,nosym
          sym1=sym(s1,:,:)
          !do i=1,3
          !    write(1,*) (s(i,j), j=1,3)
          !enddo
		g11=matmul(sym1, g1)
!	2: Crystal Symmetry to 2nd orientation
		do s2=1,nosym
              sym2 = sym(s2,:,:)
			g22=matmul(sym2,g2)
!	3: Net rotations (g2*inv(g1)=g2*transpose(g1))
      dg1=matmul(g22,transpose(g11))										
      dg2=matmul(g11,transpose(g22))		

     
!	5. CHECK IF THE MISORIENTATION IS IN FUNDAMENTAL ZONE
!	5.a. For 1st result: (dg1)
      sq=dsqrt((dg1(2,3)-dg1(3,2))**2.0+(dg1(1,3)-
     + dg1(3,1))**2.0+(dg1(2,1)-dg1(1,2))**2.0)
!      write(6,*) 'sq', sq
      if (dabs(sq).gt.1.0d-10) then
		x1=(dg1(2,3)-dg1(3,2))/sq
		x2=(dg1(3,1)-dg1(1,3))/sq
		x3=(dg1(1,2)-dg1(2,1))/sq
		if (x1.ge.0.0) then
			if (x2.ge.0.0) then
				if (x3.ge.0.0) then
					if (x1.le.x2) then
						if (x2.le.x3) then
							trace=dg1(1,1)+dg1(2,2)+dg1(3,3)
							temp_ang=dacos((trace-1.0)/2.0)
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
          

 
!	5.a. For the 2nd result: (dg2)
      sq=dsqrt((dg2(2,3)-dg2(3,2))**2.0+
     + (dg2(1,3)-dg2(3,1))**2.0+(dg2(2,1)-dg2(1,2))**2.0)
!      write(6,*) 'sq', sq
      if (dabs(sq).gt.1.0d-10) then
   	    x1=(dg2(2,3)-dg2(3,2))/sq
	    x2=(dg2(3,1)-dg2(1,3))/sq
	    x3=(dg2(1,2)-dg2(2,1))/sq
	    if (x1.ge.0.0) then
		    if (x2.ge.0.0) then
			    if (x3.ge.0.0) then
				    if (x1.le.x2) then
					    if (x2.le.x3) then
						    trace=dg2(1,1)+dg2(2,2)+dg2(3,3)
						    temp_ang=dacos((trace-1.0)/2.0)
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
	
          
    

!     Check for zero misorientation      
      if (ang.gt.0.9*pi/2.0) then
          ang = 0.0d+0
          uvw=0.0d+0
!         Assign a unit vector
          uvw(1)=1.0
          dg=I3
      endif
       
!     Convert the misorientation angle into degrees      
	ang=ang*180.0/pi
      
      

      
	return
	end subroutine misorientation
!
!	This subroutine calculates the Von-Mises stress
!	INPUTS: Stress						---	sigma(3,3)
!	OUTPUT: Von-Mises stress (scalar)	--- vms
	subroutine vonmises_stress(sigma,vms)
	use globalvars, only: I3
      implicit none
	real(8) sigma(3,3),vms,dev(3,3),hyd
	integer i,j
	real(8) sum
      
      hyd=0.
      do i=1,3
          hyd=hyd +sigma(i,i)
      enddo
      
      dev = sigma -hyd/3.0*I3
      
!
!	Take the inner product
	sum=0
	do i=1,3
		do j=1,3
			sum=sum+(3.0*dev(i,j)*dev(i,j)/2.0)
		enddo
	enddo
!	Find the Von-Mises stress
	vms=dsqrt(sum)
	return
	end subroutine vonmises_stress
!
!	This subroutine calculates the Von-Mises strain
!	INPUTS: Strain						---	eps(3,3)
!	OUTPUT: Von-Mises strain (scalar)	--- evm
	subroutine vonmises_strain(eps,evm)	
	implicit none
	real(8) eps(3,3),evm
	integer i,j
	real(8) sum
!

      


!	Take the inner product
	sum=0.0
	do i=1,3
		do j=1,3
			sum=sum+(2.0*eps(i,j)*eps(i,j)/3.0)
		enddo
	enddo
!	Find the Von-Mises strain
	evm=dsqrt(sum)
	return 
	end subroutine vonmises_strain
!
!	This subroutine calculates the Taylor factor
!	INPUTS:
!	F: Deformation gradient at a certain reference 
!	g: Transformation matrix from the reference frame at which F is defined to the crystal frame
!	OUTPUT: 
!	M: Taylor factor
!	USES:
!	BHstress: Bishop Hill stress states
	subroutine taylorfactor(F,g,M)
	use globalvars, only: BHstress
	implicit none
	real(8) F(3,3),g(3,3),M
	integer i,j,k,l
	real(8) Fc(3,3),sum,evm,st(3,3),st_vec(6),w,max_w
!
!	Find the Von-Mises strain
	call vonmises_strain(F,evm)
	if (evm.eq.0) then
		M=3.0
		return
	endif
!	Transform the deformation gradient to crystal frame	
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
!	Strain part of the transformed deformation
	st=(Fc+transpose(Fc))/2.0
!	Vectorizing stress (Convention used is different than in MARC)
!	Check Read's "Deformation Geometry" for details
	st_vec(1)=st(2,2)
	st_vec(2)=-st(1,1)
	st_vec(3)=st(3,3)
	st_vec(4)=st(1,2)+st(2,1)
	st_vec(5)=st(3,1)+st(1,3)
      st_vec(6)=st(2,3)+st(3,2)
!	Calculating work 
	max_w=0.0
	do i=1,28
		w=0.0
		do j=1,6
			w=w+BHstress(i,j)*st_vec(j)
		enddo
		if (abs(w)>max_w) then
			max_w=abs(w)
		endif
	enddo
	M=sqrt(6.0)*max_w/evm;
	return
	end subroutine taylorfactor
!
!	This subroutine converts a symmetric 2nd order tensor to a vector
!	INPUTS: Matrix			--- s(3,3)
!	OUTPUT: Vectorized form	--- v(6)	
	subroutine convert3x3to6(s,v)
	implicit none
	real(8) s(3,3),v(6)
!	
	v(1)=s(1,1)
	v(2)=s(2,2)
	v(3)=s(3,3)
	v(4)=(s(1,2)+s(2,1))/2.0
	v(5)=(s(3,1)+s(1,3))/2.0
      v(6)=(s(2,3)+s(3,2))/2.0
	return
	end subroutine convert3x3to6
!
!	This subroutine converts a vector to a SYMMETRIC 2nd order tensor
!	INPUTS: Vectorized form		--- v(6)
!	OUTPUT: Matrix				--- s(3,3)
	subroutine convert6to3x3(v,s)
	implicit none
	real(8) s(3,3),v(6)
!
	s(1,1)=v(1)
	s(2,2)=v(2)
	s(3,3)=v(3)
	s(1,2)=v(4)
	s(2,1)=v(4)
	s(1,3)=v(5)
	s(3,1)=v(5)      
	s(2,3)=v(6)
	s(3,2)=v(6)
	return
	end subroutine convert6to3x3
!
!	This subroutine is written to convert 4th order elasticity tensor to 6x6 matrix
!	INPUTS : 4th order tensor		--- C4(3,3,3,3)
!	OUTPUTS: Matrix					---	C_6x6(6,6)
	subroutine convert3x3x3x3to6x6(C4,C_6x6)
	implicit none
	real(8) C_6x6(6,6), C4(3,3,3,3)
	integer i,j,k,l


      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      C4(i,j,k,l) = (C4(i,j,k,l) + C4(j,i,k,l))/2.0
                      C4(i,j,k,l) = (C4(i,j,k,l) + C4(j,i,l,k))/2.0
                  end do
              end do
          end do
      end do

      
      
!     1=11 / 2=22 / 3=33 / 4=12 / 5=13 / 6=23
      

!     Stiffness Matrix(6x6)

      C_6x6 = 0.0
      C_6x6(1,1) = C4(1,1,1,1)
      C_6x6(1,2) = C4(1,1,2,2) 
      C_6x6(1,3) = C4(1,1,3,3)
      C_6x6(1,4) = C4(1,1,1,2) 
      C_6x6(1,5) = C4(1,1,1,3) 
      C_6x6(1,6) = C4(1,1,2,3)
      
      C_6x6(2,1) = C4(2,2,1,1)
      C_6x6(2,2) = C4(2,2,2,2) 
      C_6x6(2,3) = C4(2,2,3,3) 
      C_6x6(2,4) = C4(2,2,1,2)
      C_6x6(2,5) = C4(2,2,1,3)
      C_6x6(2,6) = C4(2,2,2,3)
      
      C_6x6(3,1) = C4(3,3,1,1)
      C_6x6(3,2) = C4(3,3,2,2)
      C_6x6(3,3) = C4(3,3,3,3)
      C_6x6(3,4) = C4(3,3,1,2)
      C_6x6(3,5) = C4(3,3,1,3)
      C_6x6(3,6) = C4(3,3,2,3)
      
      C_6x6(4,1) = C4(1,2,1,1)
      C_6x6(4,2) = C4(1,2,2,2)
      C_6x6(4,3) = C4(1,2,3,3)
      C_6x6(4,4) = C4(1,2,1,2)
      C_6x6(4,5) = C4(1,2,1,3)
      C_6x6(4,6) = C4(1,2,2,3)
      
      C_6x6(5,1) = C4(1,3,1,1)
      C_6x6(5,2) = C4(1,3,2,2)
      C_6x6(5,3) = C4(1,3,3,3)
      C_6x6(5,4) = C4(1,3,1,2)
      C_6x6(5,5) = C4(1,3,1,3)
      C_6x6(5,6) = C4(1,3,2,3)
      
      C_6x6(6,1) = C4(2,3,1,1)
      C_6x6(6,2) = C4(2,3,2,2)
      C_6x6(6,3) = C4(2,3,3,3)
      C_6x6(6,4) = C4(2,3,1,2)
      C_6x6(6,5) = C4(2,3,1,3)
      C_6x6(6,6) = C4(2,3,2,3)
      



	return
	end subroutine convert3x3x3x3to6x6
!
!
!     
!
!	This subroutine converts a 4th order tensor to 9x9 matrix
!	INPUTS: 4th order tensor	--- r(3,3,3,3)
!	OUTPUT: Matrix				--- s(9,9)
!	USES  :	Order of indices	---	order(9,2)
	subroutine convert3x3x3x3to9x9(r,s)
	use globalvars, only : order
	implicit none
	real(8) r(3,3,3,3), s(9,9)
	integer i,j,a,b,c,d
!	
	do i=1,9
		a=order(i,1)
		b=order(i,2)
		do j=1,9
			c=order(j,1)
			d=order(j,2)
			s(i,j)=r(a,b,c,d)
		enddo
	enddo
	return
	end subroutine convert3x3x3x3to9x9
!
!	This subroutine converts a vectorized 2nd order tensor to a matrix 
!	INPUTS: Vectorized form		--- s(9)
!	OUTPUT: Matrix				--- r(3,3)
!	USES  :	Order of indices	---	order(9,2)
	subroutine convert9to3x3(s,r)
	use globalvars, only : order
	implicit none
	real(8) r(3,3), s(9)
	integer i,a,b
!	
	do i=1,9
		a=order(i,1)
		b=order(i,2)
		r(a,b)=s(i)
	enddo
	return
	end subroutine convert9to3x3
!
!	This subroutine converts a matrix to a vectorized form
!	INPUTS: Matrix				--- r(3,3)
!	OUTPUT: Vectorized form		--- s(9)
!	USES  :	Order of indices	---	order(9,2)
	subroutine convert3x3to9(r,s)
	use globalvars, only : order
	implicit none
	real(8) r(3,3), s(9)
	integer i,a,b
!	
	do i=1,9
		a=order(i,1)
		b=order(i,2)
		s(i)=r(a,b)
	enddo
	return
	end subroutine convert3x3to9
!
!
!	This subroutine transforms a 2nd order tensor using a given transformation matrix
!	INPUTS : Transformation matrix, 2nd order tensor	--- g, A33
!	OUTPUTS: Transformed 2nd order tensor				--- A33_tr
	subroutine transform2(g,A33,A33_tr)
	implicit none
	real(8) g(3,3),A33_tr(3,3),A33(3,3)
	real(8) sum
	integer i,j,m,n
!
	do i=1,3
		do j=1,3
			sum=0.0
			do m=1,3
				do n=1,3
					sum=sum+(g(i,m)*g(j,n)*A33(m,n))
				enddo	
			enddo
			A33_tr(i,j)=sum
		enddo
	enddo
	return
	end subroutine transform2
!
!	This subroutine transforms a 4th order tensor using a given transformation matrix
!	INPUTS : Transformation matrix, 4th order tensor	--- g, A3333
!	OUTPUTS: Transformed 4th order tensor				--- A3333_tr
	subroutine transform4(g,A3333,A3333_tr)
	implicit none
	real(8) g(3,3),A3333(3,3,3,3),A3333_tr(3,3,3,3)
	real(8) sum
	integer i,j,k,l,m,n,o,p
!		

	do i=1,3
		do j=1,3
			do k=1,3
				do l=1,3
					sum=0.0
					do m=1,3
						do n=1,3
							do o=1,3
								do p=1,3
									sum=sum+(g(i,m)*g(j,n)*g(k,o)*g(l,p)*A3333(m,n,o,p))
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
!
!	This subroutine calculates the square norm of a matrix
!	INPUTS: Matrix, dimension(n)		--- A(nxn),dim
!	OUTPUT: Square norm				--- norm
	subroutine normmat(A,dim,norm)
	implicit none
	integer dim
	real(8) A(dim,dim),norm
	integer i,j
!
	norm=0.0
	do i=1,dim
		do j=1,dim
			norm=norm+(A(i,j)*A(i,j))
		enddo
	enddo
	norm=dsqrt(norm)
	return
	end subroutine normmat
!
!	This subroutine calculates the trace of a matrix
!	INPUTS: Matrix, dimension(n)		--- A(nxn),dim
!	OUTPUT: Trace of a matrix			--- value
	subroutine trace(A,dim,value)
	implicit none
	integer dim
	real(8) A(dim,dim),value
	integer i
!
	value=0.0
	do i=1,dim
		value=value+A(i,i)
	enddo
	return
	end subroutine trace
!
!	This subroutine calculates the square norm of a column vector
!	INPUTS: Vector, dimension(n)		--- A(nx1),dim
!	OUTPUT: Square norm				--- norm
	subroutine normvec(A,dim,norm)
	implicit none
	integer i,dim
	real(8) A(dim),norm
!
	norm=0.0
	do i=1,dim
		norm=norm+(A(i)*A(i))
	enddo
	norm=dsqrt(norm)
	return
	end subroutine normvec
!
!	This subroutine multiplies a matrix with a matrix
!	INPUTS: Matrix (nxn), vector (nx1), dimension (n)	--- A(nxn),B(n,n),dim
!	OUTPUT: Result										--- C(nxn)
	subroutine matmulmat(A,B,dim,C)
	implicit none
	integer i,j,k,dim
	real(8) A(dim,dim),B(dim,dim),C(dim,dim)
!
	C=0.0
	do i=1,dim
		do j=1,dim
			do k=1,dim
				C(i,j)=C(i,j)+(A(i,k)*B(k,j))
			enddo
		enddo
	enddo
	return
	end subroutine matmulmat
!
!	This subroutine multiplies a matrix with a vector
!	INPUTS: Matrix (nxn), vector (nx1), dimension (n)	--- A(nxn),b(n,1),dim
!	OUTPUT: Result										---	c(n,1)
	subroutine matmulvec(A,b,dim,c)
	implicit none
	integer i,j,k,dim
	real(8) A(dim,dim),b(dim),c(dim)
!
	c=0.0
	do i=1,dim
		do j=1,dim
			c(i)=c(i)+(A(i,j)*b(j))
		enddo
	enddo
	return
	end subroutine matmulvec
!
!	This subroutine is written to check the determinant
!	INPUTS: 3x3 matrix	--- a(3,3)
!	OUTPUT: Determinant	--- det
	subroutine determinant(a,det)
	implicit none
	real(8) a(3,3),det
	real(8) v1,v2,v3
!	
	v1=a(1,1)*((a(2,2)*a(3,3))-(a(2,3)*a(3,2)))
	v2=-a(1,2)*((a(2,1)*a(3,3))-(a(2,3)*a(3,1)))
	v3=a(1,3)*((a(2,1)*a(3,2))-(a(2,2)*a(3,1)))
	det=v1+v2+v3
	return
	end subroutine determinant
!
!	This subroutine inverts a 3x3 matrix
!	INPUT:	Matrix								---	A(3,3)
!	OUTPUT:	Invereted matrix, determinant		---	invA(3,3),det
	subroutine invert3x3(A,invA,det)
	implicit none
	integer i,j
	real(8) A(3,3),invA(3,3),det
!
!	First calculate the determinant
	call determinant(A,det)
!	If the determinant is greater than certain value
	if (det.lt.1.d-10) then
		invA=0.0
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
!
!	This subroutine orientation matrix from Euler angles
!	INPUT:	Angles(deg)			---	ang(3)
!	OUTPUT:	Orientation matrix	---	R(3,3)
!     USES:     Number pi           --- pi
	subroutine ang2ori(ang,R)
	use globalvars, only : pi
	implicit none
	real(8) ang(3),R(3,3)
	real(8) phi1,phi2,PHI
!
	phi1=ang(1)*pi/180.0
	phi2=ang(3)*pi/180.0
	PHI=ang(2)*pi/180.0
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
!
!	This subroutine orientation matrix from Euler angles
!	INPUT:	Orientation matrix	---	R(3)
!	OUTPUT:	Angles(deg)			---	ang(3)
!     USES:     Number pi           --- pi
	subroutine ori2ang(R,ang)
	use globalvars, only : pi
	implicit none
	real(8) ang(3),R(3,3)
	real(8) phi1,phi2,PHI
!
	if (R(3,3).eq.1) then
		PHI=0.0
		phi1=atan2(R(1,2),R(1,1))
		phi2=0.0
	else
		PHI=dacos(R(3,3))
		phi1=atan2(R(3,1)/dsin(PHI),-R(3,2)/dsin(PHI))
		phi2=atan2(R(1,3)/dsin(PHI),R(2,3)/dsin(PHI))
	endif
	ang=ang*180.0/pi
	return
	end subroutine ori2ang
!
!	This subroutine calculates angle and axis pair from orientation matrix
!	INPUT:	Orientation matrix	--- R(3,3)
!	OUTPUT:	Angle(deg), axis	---	ang, ax(3)
!     USES:     Number pi           --- pi
	subroutine ori2angax(R,ang,ax)
	use globalvars, only : pi
	implicit none
	real(8) R(3,3),ang,ax(3)
	real(8) trace,norm,t,mag
!
      norm=dsqrt((R(2,3)-R(3,2))**2.0+(R(1,3)-
     + R(3,1))**2.0+(R(1,2)-R(2,1))**2.0)
	 
	trace=R(1,1)+R(2,2)+R(3,3)
	ang=2.0*acos(dsqrt(1+trace)/2.0)
	t=dtan(ang/2.0)
	ax(1)=t*(R(2,3)-R(3,2))/norm
	ax(2)=t*(R(3,1)-R(1,3))/norm
	ax(3)=t*(R(1,2)-R(2,1))/norm	
!	Normalize the axis vector
	mag=dsqrt(ax(1)**2.0+ax(2)**2.0+ax(3)**2.0)
	ax=ax/mag		
	ang=ang*180.0/pi
	return
	end subroutine ori2angax
!
!	This subroutine finds the symmetric part of a matirx
!	INPUT: 2nd order tensor	--- F(3,3)
!	OUPUT:	Symmetric part	--- S(3,3)
	subroutine sym(F,S)
	implicit none
	real(8) F(3,3),S(3,3)
	integer i,j
!
	do i=1,3
		do j=1,3
			S(i,j)=0.5*(F(i,j)+F(j,i))
		enddo
	enddo
	return
	end subroutine sym
!
!	This subroutine finds the skew-symmetric part of a matirx
!	INPUT: 2nd order tensor	--- F(3,3)
!	OUPUT:	Symmetric part	--- S(3,3)
	subroutine skew(F,S)
	implicit none
	real(8) F(3,3),S(3,3)
	integer i,j
!
	do i=1,3
		do j=1,3
			S(i,j)=0.5*(F(i,j)-F(j,i))
		enddo
	enddo
	return
	end subroutine skew
!
!	Polar decomposition of deformation
!	The method mentioned in 'Mechanics of deformable solids' by
!	Isaam Dogri was implemented
!	INPUT:	2nd order tensor	--- F(3,3)
!	USES :	Indentity matrix	--- I3(3,3)
!	OUPUT:	Rotation	--- R(3,3)
!			Stretch		--- U(3,3)
	subroutine polar(F,R,U,sing)
	use globalvars, only : I3
	implicit none
	real(8) F(3,3),R(3,3),U(3,3)
	real(8) RC(3,3),inv1,inv2,inv3,dummy,mat(3,3)
	real(8) p,q,p3,Cf,f3,lambda1,lambda2,lambda3
	real(8) ii1,ii2,ii3,check,invU(3,3)
	integer i,j,sing
!	Assign zero to nozero
	sing=0
!	Find the right Cauchy
	RC=matmul(transpose(F),F)
!	Find the principal invariants
!	First invariant is the trace along the diagonals
	call trace(RC,3,inv1)
!	Assign the multiplier sign in front of the equation
	inv1=-inv1
!	Second invariant
	mat=matmul(RC,RC)
	call trace(mat,3,dummy)
	inv2=((inv1**2.0)-dummy)/2.0
!	Third invariant
	call determinant(RC,inv3)
!	Assign the sing as well
	inv3=-inv3
!	The constants
	p=(-(inv1**2.0)/3.0)+inv2
	q=(-(inv2+(p*2.0))*inv1/9.0)+inv3;
	check=(4.0*(p**3.0))+(27.0*(q**2.0));
!	If the check is greater than zero then there exists only one
!	real eigenvalues and two complex eigenvalues
	if (check.gt.0) then
		sing=1
		return
	else
		if (p.lt.0) then
			P3=dsqrt(-p/3.0)
			Cf=3.0*q/(2.0*p*P3)
			if (Cf.gt.1) then
			  Cf=1.0
			elseif (Cf.lt.-1) then
			  Cf=-1.0
			endif
			f3=dacos(Cf)/3.0
      lambda1=dsqrt((-inv1/3.0)+(2.0*P3*dcos(f3)))
      lambda2=dsqrt((-inv1/3.0)-(P3*(dcos(f3)+
     + (dsin(f3)*dsqrt(3.0d+0)))))
      lambda3=dsqrt((-inv1/3.0)-(P3*(dcos(f3)-
     + (dsin(f3)*dsqrt(3.0d+0)))))
		else
			sing=1
			return
		endif
		ii1=lambda1+lambda2+lambda3
		ii2=(lambda1*lambda2)+(lambda3*lambda2)+(lambda1*lambda3)
		ii3=lambda1*lambda2*lambda3
!		Another check for the singularity
		check=(ii1*ii2)-ii3
		if (check.eq.0.0) then
			sing=1
			return
		else
			U=((-mat)+(((ii1**2.0)-ii2)*RC)+(ii1*ii3*I3))/check
			invU=(RC-(ii1*U)+(ii2*I3))/ii3
			R=matmul(F,invU)
		endif
	endif
	return
	end subroutine polar
!
!	This subroutine calculates the L-U decomp.
!	INPUT: Matrix, dimension					--- F(n,n), n
!	OUTPUT: Lower tri., Upper tri., determinant --- L(n,n), U(n,n), det
	subroutine LUdecomp(A,n,L,U,det)
	implicit none
	integer n,i,j,k
	real(8) A(n,n),L(n,n),U(n,n)
	real(8) factor(n,n), det
!
	U=A
	factor=0.0
!	Triangularization
	do k=1,n-1
		do i=k+1,n
			factor(i,k)=U(i,k)/U(k,k)
			do j=1,n
				U(i,j)=U(i,j)-factor(i,k)*U(k,j)
			enddo	
		enddo
	enddo
!	Lower triangular matrix
	L=factor
	do i=1,n
		L(i,i)=1.0
	enddo
!	Multiplication of the elements at the diagonal is the determinant
	det=1;
	do k=1,n
		det=det*L(k,k)*U(k,k)
	enddo
	return
	end subroutine LUdecomp
!	
!	This subroutine inverts a nxn matrix
!	INPUT:	Matrix, size of the matrix		---	A(n,n),n
!	OUTPUT:	Invereted matrix, determinant	---	invA(n,n),det
	subroutine invertnxn(A,n,invA,det)
	implicit none
	integer n,i,j,m
	real(8) A(n,n),invA(n,n),det,sum
	real(8) L(n,n),U(n,n),B(n,n),C(n,n),Ct(n,n),Lt(n,n)
!	First decompose into lower and upper triangles	
	call LUdecomp(A,n,L,U,det)
!	If the determinant is greater than zero
	if (det.gt.1.0d-10) then
!		Inversion algorithm (For the UPPER triangular matrix)
		do i=1,n
			do j=1,n
				if (i.gt.j) then
					B(i,j)=0.0
				elseif (i==j) then
					B(i,j)=1/U(i,j)			
				else
					sum=0.0
					do m=1,j-1
						sum=sum-(B(i,m)*U(m,j))
					enddo
					B(i,j)=sum/U(j,j)
				endif
			enddo
		enddo
!		First transpose (For the LOWER triangular matrix)
		Lt=transpose(L)
		do i=1,n
			do j=1,n
				if (i.gt.j) then
					Ct(i,j)=0
				elseif (i==j) then
					Ct(i,j)=1/Lt(i,j)			
				else
					sum=0.0
					do m=1,j-1
						sum=sum-(Ct(i,m)*Lt(m,j))
					enddo
					Ct(i,j)=sum/Lt(j,j)
				endif
			enddo
		enddo
!		Transpose the result back finally
		C=transpose(Ct)
!		Calculate the overall inverse
		invA=matmul(B,C)
	else
		invA=0.0
      endif
      
	return
	end subroutine invertnxn
!
      subroutine inverse(a,c,n)
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
      L=0.0
      U=0.0
      b=0.0

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
      L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
          do i=1,j
              U(i,j) = a(i,j)
          end do
      end do

      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
          b(k)=1.0
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
          b(k)=0.0
      end do
      end subroutine inverse
      
      
      
      

      

      end module globalsubs
