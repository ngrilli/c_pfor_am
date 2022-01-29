c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c
c Aug. 12th, 2021 - 1st working version
c
c
c This subroutine is taken from Dr. Dylan Agius (GitHub)
	module lengthscale
      implicit none
      contains
      
c     Includes the subroutines that are used to compute the length-scale effects      
      
   	!----------------------------------------------------------------------------
    	subroutine grainsize(coords,ldistance,rdistance,lm,  
     &		 slpdir,slpplane,nslptl,feature,noel)



	!subroutine containing the length scale calculation
	!INPUTS :
	!coords             - voxel coordinates
	!nodeout            - columns in nodex, nodey, nodez,boundfeat
	!numgrain          - total number of grains in the RVE
	!nodex,nodey,nodez  - x, y, and z coordinates for each boundary node
	!slpdir1             - slip directions in local (crystal) coordinate system
	!slpnor1             - slip plane normals in local (crystal) coordinate system
	!elcent             - array containing the centroids of all voxels
	!feature            - grain number
	!boundgrain         - boundary grain ids
	!grainori            - orientations in euler angles in each grain
	!noel               - voxel number

	!OUTPUTS:
	!rdistance          - distance from the voxel to grain boundary
	!ldistance          - distance along slip length
	!lm                 - Luster-Morris parameter

      use globalvars, only : pi, nodex, nodey,
     & nodez, boundgrain, elcent, grainori, nodeout
      use globalsubs, only : ang2ori
	implicit none

		
      
	   	
      integer:: i, arraysize,tester,grainb,    
     &arraysizeb,valb,NSLPTL,j,MININDEXSING,CHARCTERLENGT
      integer :: val,noel
      integer :: minidexinsing

      real*8, allocatable ::  mindist(:),vect(:,:),     		 
     &sortmindist(:),sortedbgrains(:),uniqu(:),uniindex(:),	      
     &uniqgraindist(:),zrot(:,:),slpdirrotate(:,:),                
     &grainboundindex(:) ,featureboundsnodes(:,:),		                
     &bxvalue(:),byvalue(:),bzvalue(:),btotal(:,:),bangle(:,:),      
     &normarray(:),normslp(:),minangleval(:),minangleval180(:),   
     &minangleactual(:),rnodes(:,:),vectr(:,:),			 
     &nodesnot(:,:),lxtotal(:,:,:),lxnormin(:,:),lxnorm(:),        
     &vectrnorm(:),lxangle(:,:),minxval0(:,:),minxval180(:,:),     
     &minxval0act(:,:),minxval180act(:,:),minxangle(:),            
     &minxangleact(:),lxnodes(:,:),lxnodesindex(:),                
     &slpdirrotatenorm(:),veca(:,:),vecb(:,:),checkangle(:),       
     &vanorm(:),vbnorm(:),slpplanrotate(:,:),         
     &slpdirrotatea(:,:),slpplanrotatea(:,:),slpplanrotatenorm(:), 
     &slpdirrotateanorm(:),slpplanrotateanorm(:),cosphi(:),cosk(:)
	  
	integer, allocatable :: minindex(:),minang0(:,:),		
     &minang180(:,:),boundnodegrainb(:),nodesnotindex(:), 
     &minangleactloc(:)
      
      integer FEATURE
	  
	real*8 :: coords(3),vectotal(3),                             
     &euler1,totalrot(3,3),groundbundid,charctlenght,           
     &twovoxdist(3),centrodes(3),totalrota(3,3),                   
     &sortmindistsing,uniqgraindistsing,testera

      real*8 :: lm(NSLPTL),ldistance(NSLPTL),rdistance(NSLPTL)


	  
	real*8 :: slpdir(NSLPTL,3), slpplane(NSLPTL,3)
!     &SLPDIR(NSLPTL,3),  slpplane(NSLPTL,3)


	LOGICAL, allocatable :: mk(:)
			
	logical :: locgrain= .False.


      real(8) Bunge(3)
      
      
      
	   
	!adjust the slip plane normal and slip direction to removing 
	 !any trailing zeros

	
	!allocate(slpdir(nslptl,3))
	!allocate(slpplane(nslptl,3))
	!	
	!do i=1,nslptl
	!!do i=1,nd
	!	do j=1,3
	!		slpdir(i,j)=slpdir1(i,j)
	!		slpplane(i,j)=slpnor1(i,j)
	!	end do
	!end do
		
	 
	! check to remove values of interest ** this for loop can be 
	! replaced with 'findloc'**

	

	checkerloop: DO i=1,nodeout
		
		if(int(nodey(feature,i)) < 0)then
			
			arraysize=i-1
			
			exit checkerloop
		else
			
			arraysize=i
				
		end if
	end do checkerloop


	 		
		 
 
	!find the closest boundary.  This will form 'Grain B' from which 
	!the information on orientation is extracted to rotate the slip
	!systems

	!create a vector from the current voxel position to boundary nodes
		allocate (vect(3,arraysize))
	
	
	vect(1,:)=nodex(feature,1:arraysize)-coords(1)
	
	vect(2,:)=nodey(feature,1:arraysize)-coords(2)
	
	vect(3,:)=nodez(feature,1:arraysize)-coords(3)
	
	vectotal=(/vect(1,:),vect(2,:),vect(3,:)/)
	

	allocate (mindist(arraysize))
	allocate (sortmindist(arraysize))
	allocate (mk(arraysize))
	allocate (minindex(arraysize))
	allocate (sortedbgrains(arraysize))
		  
	! calculate the distance of the current location to the grain boundary and
	! multiplying to bring distances to microns
      mindist=(((vect(1,:)**2.D0) + (vect(2,:)**2.D0)+           
     1              (vect(3,:)**2.D0))**0.5D0)

					   
	! sort distance by the minimum size and and keep indices of ordered minimum

	minindexsing=minloc(mindist,1)
	sortmindistsing=minval(mindist)

	
		 

	   ! mk= .TRUE.
	   ! DO i = 1, arraysize
	   !     sortmindist(i) = MINVAL(mindist,mk)
	   !     minindex(i)=MINLOC(mindist,1,mk)
	   !     mk(MINLOC(mindist,mk)) = .FALSE.
	   ! END DO
		
	! using the ordered indices, find the grains that these belong to
	!sortedbgrains=boundgrain(feature,minindex)

		

	!find the unique grains in this list
	! call remove_dups(sortedbgrains,arraysize,uniqu,uniindex)

		
		
	!using the index of unique grains, the corresponding distance value can be extracted
	allocate (uniqgraindist(size(uniindex)))
	uniqgraindistsing=int(boundgrain(feature,minindexsing))
		
	!loop over all indentified touching grains.  This is done to ensure that the 
	!voxel isn't situation in the middle which may result in an identified grain 
	!being not the most suitable.

	!***********leave this for later.  Let's just see how it works using the closest
	!grain boundary based on the euclidean distance.


	allocate( slpdirrotate(NSLPTL,3))
	allocate( slpplanrotate(NSLPTL,3))
	   

      grainb=uniqgraindistsing



		 
		!rotate the slip systems using the angles from each identified grain
c	call eulercosmatrix(grainori(grainb,:),totalrot)
c     This has changed to the ang2ori(ang,R)
c     The result is transposed since crystal to sample transformation is required
      Bunge = grainori(grainb,:)
      call ang2ori(Bunge,totalrot)
      totalrot=transpose(totalrot)
      
	do j =1, NSLPTL
		slpdirrotate(j,:)= matmul(totalrot, slpdir(j,:))
		slpplanrotate(j,:)=matmul(totalrot, slpplane(j,:))
	end do
		
	allocate( slpdirrotatea(NSLPTL,3))
	allocate( slpplanrotatea(NSLPTL,3))
	!rotate slip systems for the current grain
c	call eulercosmatrix(oriensh(feature,:),totalrota)
c     This has changed to the ang2ori(ang,R)
c     TRANSPOSE IS REVERT BACK TO THE ORIGINAL VERSION
c     The result is transposed since crystal to sample transformation is required
      Bunge = grainori(int(feature),:)
      call ang2ori(Bunge,totalrota)
c      totalrota=transpose(totalrota)
		
	do j =1, NSLPTL
		slpdirrotatea(j,:)= matmul(totalrota, slpdir(j,:))
		slpplanrotatea(j,:)=matmul(totalrota, slpplane(j,:))
	end do
		
	! end do idgrains    
		
	!using the identified closest grain, construct vectors extending from the
	!current voxel position to the grain boundary.

	!firstly extract only the nodes at the shared boundary between grain A and 
	!grain B.
	   
		
	!firtly find the size of the array at the indentified grain B
	! check to remove values of interest **this for loop can be replaced with 
	! 'findloc' if it is available**
	
c     THIS LINE IS CONVERTED TO ORIGINAL VERSION - IT DOES NOT EXIST
c      testera=minloc(nodey(int(grainb),:),dim=1)

	checkerloopb: DO i=1,nodeout
c         THIS LINE IS CONVERTED TO ORIGINAL VERSION          
c		if(nodex(int(grainb),i) < 0)then
          if(nodex(int(grainb),i) .eq. 9999999)then 
			arraysizeb=i-1
			exit checkerloopb
		else
			arraysizeb=i
				
		 end if
	end do checkerloopb
		
		
		

	!************************************************************        
	!find index of array locations which are true.
	!For loops are used here to make it more generically applicable
	!but you can also use 'all' first to create a boolean array
	! and the 'findloc' to find .TRUE. values.
		
	val=1
	do i=1,arraysizeb
		if(int(boundgrain(grainb,i)) .eq. feature)then
			val=val+1
		end if
	end do
		
	allocate (grainboundindex(val-1))
	val=1
	do i=1,arraysizeb
		if(int(boundgrain(grainb,i)) .eq. feature)then
			grainboundindex(val)=i
			val=val+1
		end if
	end do
		
	   !**********************************************************
	  
        
	   !find the coordinates of the identified shared boundary nodes
	   allocate (featureboundsnodes(val-1,3))
	   featureboundsnodes=reshape((/nodex(grainb,grainboundindex), 
     1   	 nodey(grainb,grainboundindex),                       
     2    	nodez(grainb,grainboundindex)/),(/val-1,3/))
		
		 
	   !calculate the vector from the boundary nodes to the current
	   !voxel
	   allocate (bxvalue(val-1))
	   allocate (byvalue(val-1))
	   allocate (bzvalue(val-1))
	   bxvalue=featureboundsnodes(:,1)-coords(1)
	   byvalue=featureboundsnodes(:,2)-coords(2)
	   bzvalue=featureboundsnodes(:,3)-coords(3)
		
	   allocate (btotal(val-1,3)) 
	   btotal=reshape((/bxvalue,byvalue,bzvalue/),(/val-1,3/))  
		
	   
	   !find the closest angle between the vectors from the interface
	   !boundary to the voxel and the slip direction in grain B.
	   
	   allocate (bangle(size(slpdirrotate,1),val-1))
	   
	   !these subroutine calls can just be replaced with the 
	   !intrisic function norm2*********
	   allocate (normarray(val-1))
	   allocate (normslp(size(slpdirrotate,1)))
	   
	   normarray=norm2(btotal,dim=2)
	   !call enorm(btotal,val-1,normarray)
		
	   !call enorm(slpdirrotate,size(slpdirrotate,1),normslp)
	   normslp=norm2(slpdirrotate,dim=2)
	   
	   !****************************************************
	   !determine if the current voxel is on the boundary.
	   !If this is true, this voxel should form the grain boundary voxel
	   !To do this, determine if any euclidean distances are less than 
	   !the elements characteristic length
	   
	   !calculate characteristic length of voxel
	   
	   twovoxdist=(elcent(3,:)-elcent(2,:))/2.D0
	   charcterlengt=(twovoxdist(1)**2.D0 + twovoxdist(2)**2.D0 +  
     1   twovoxdist(3)**2.D0)**0.5D0
	   
			
	  !this can be done with findloc if available
	  do i=1,size(normarray,1)
		if(normarray(i) .lt. charcterlengt)then
			locgrain=.True.
		 end if
		 
	  end do
	   
	    
	   allocate (minang0(size(slpdirrotate,1),1))
	   allocate (minangleval(size(slpdirrotate,1)))
	   allocate (minang180(size(slpdirrotate,1),1))
	   allocate (minangleval180(size(slpdirrotate,1)))
	   allocate (minangleactual(size(slpdirrotate,1)))
	   allocate (minangleactloc(size(slpdirrotate,1)))
	   allocate (boundnodegrainb(size(slpdirrotate,1)))
	   allocate (rnodes(size(slpdirrotate,1),3))
	  



		
	   if (locgrain .eqv. .FALSE.) then
	  
	   !voxel not on the boundary
	 
		   do i=1,size(slpdirrotate,1)
		   !calculate angle between the slipdirection in grian B and the vector
		   !created between the current voxel and the boundary nodes.
				bangle(i,:)=dacos(matmul(btotal,slpdirrotate(i,:))/    
     1       			(normarray*normslp(i))) 
			
		   !find minimum angle and the location index
		   !of the array for each slip direction 
				minang0(i,:)=minloc(bangle(i,:))
				
				minangleval(i)=bangle(i,minang0(i,1))
				 
			!find the minimums at angles of 180 degrees
				minang180(i,:)=minloc(abs(bangle(i,:)-pi))
				minangleval180(i)=abs(bangle(i,minang180(i,1))-pi)
			!find the smaller of the 0 and 180 degrees to determine
			!the true mininum value
				if (minangleval(i) .lt. minangleval180(i)) then
					minangleactual(i)=minangleval(i)
					minangleactloc(i)=minang0(i,1)
				else
					minangleactual(i)=minangleval180(i)
					minangleactloc(i)=minang180(i,1)
				end if        
		   end do 
			
		   
		   !using the minimum angles, the nodes at the boundary  which connect to the voxel location
		   !can be determined
		   
	  
		   
		   boundnodegrainb=grainboundindex(minangleactloc)
		   rnodes=reshape((/nodex(grainb,int(boundnodegrainb)),        
     1       		nodey(grainb,int(boundnodegrainb)),                
     2      		 nodez(grainb,int(boundnodegrainb))/),              
     3      		(/size(slpdirrotate,1),3/)) 
					
		   ! the corresponding vector for each slip direction
		   allocate (vectr(size(slpdirrotate,1),3))
		  
		   
		   vectr=reshape((/nodex(grainb,int(boundnodegrainb))-coords(1),
     1              	nodey(grainb,int(boundnodegrainb))-coords(2),
     2              	nodez(grainb,int(boundnodegrainb))-coords(3) 
     3              	/),(/size(slpdirrotate,1),3/))
						   
		   ! calculate the distance from the voxel centroid to the boundary nodes
c             WHY THIS IS IN 2D?             
		   rdistance=norm2(vectr,dim=2)
		   !call enorm(vectr,size(slpdirrotate,1),rdistance)
		   
	   else
		
	   !since the voxel is on the boundary, the distance to the boundary is taken 
	   !as half the voxel characteristic length
c         THIS IS CHANGED BY ERALP
		   rdistance(1:size(slpdirrotate,1))=charcterlengt*0.5D0
c             rdistance=0.5D0
cc             Need to multiply with the real voxelsize
c		   rdistance=0.5D0*voxelsize
		

		   do i=1,3
			  rnodes(1:size(slpdirrotate,1),i)=elcent(int(noel),i)
		   end do
	   end if
		 
	   ! calculate the distance from the boundary to boundary in grain B along
	   ! the slip direction
	   
	   !****************************************************************
	   !This is to be done by excluding the nodes at the boundary.
	   !For loops are used here to make it more generically applicable
	   !but you can also use 'all' first to create a boolean array
	   ! and the 'findloc' to find .TRUE. values.
	   
		valb=1
		do i=1,arraysizeb
			if(int(boundgrain(grainb,i)) .ne. feature)then
				valb=valb+1
			end if
		end do
		
		allocate (nodesnotindex(valb-1))
		valb=1
		do i=1,arraysizeb
			if(int(boundgrain(grainb,i)) .ne. feature)then
				nodesnotindex(valb)=i
				valb=valb+1
			end if
		end do

		!**********************************************************
		!node coordinates not at the shared boundary.
		allocate (nodesnot(valb-1,3))
		
		nodesnot=reshape((/nodex(grainb,int(nodesnotindex)),      
     1            	nodey(grainb,int(nodesnotindex)),          
     2           	 nodez(grainb,int(nodesnotindex))           
     3            	/),(/valb-1,3/))
					   
		
		allocate (lxtotal(size(slpdirrotate,1),valb-1,3))
		allocate (lxnorm(valb-1))
		allocate (lxnormin(valb-1,3))
		allocate (lxangle(size(slpdirrotate,1),valb-1))
		allocate (minxval0(size(slpdirrotate,1),1))
		allocate (minxval180(size(slpdirrotate,1),1))
		allocate (minxval180act(size(slpdirrotate,1),1))
		allocate (minxval0act(size(slpdirrotate,1),1))
		allocate (minxangle(size(slpdirrotate,1)))
		allocate (minxangleact(size(slpdirrotate,1)))
	   
		!calculate vector from the boundary to the nodes in Grain B not 
		!at the boundary
		
		! check to see if the current voxel is at the boundary or not
		! If not, progress in the following
		allocate (vectrnorm(size(slpdirrotate,1)))
		allocate (slpdirrotatenorm(size(slpdirrotate,1)))
		if (locgrain .eqv. .FALSE.)then
			
			!call enorm(vectr,size(slpdirrotate,1),vectrnorm)
			vectrnorm=norm2(vectr,dim=2)
			do i=1,size(slpdirrotate,1)
				lxtotal(i,:,:)=reshape((/nodex(grainb,int(nodesnotindex))   
     1            			-rnodes(i,1),                                
     2            			nodey(grainb,int(nodesnotindex))-rnodes(i,2),
     3           			 nodez(grainb,int(nodesnotindex))-rnodes(i,3) 
     4            			/),(/valb-1,3/))
				lxnormin=lxtotal(i,:,:)
				!******can be replaced with norm2 function
				!call enorm(lxnormin,valb-1,lxnorm)
				lxnorm=norm2(lxnormin,dim=2)

				!find the angle between the vectors from the shared boundary to 
				!boundary nodes in grain b and the vector developed from the voxel
				!to the shared boundary
				lxangle(i,:)= dacos(matmul(lxtotal(i,:,:),vectr(i,:))/(lxnorm*  
     1         				vectrnorm(i)))
				minxval0(i,:)=minloc(lxangle(i,:))
				minxval0act(i,:)=lxangle(i,int(minxval0(i,:)))
				minxval180(i,:)=minloc(abs(lxangle(i,:)-pi))
				minxval180act(i,:)=abs(pi-lxangle(i,int(minxval180(i,:))))
				
				minxangle(i)=minxval0(i,1)
				minxangleact(i)=minxval0act(i,1)
				
			 end do
		 else
			
			!call enorm(slpdirrotate,size(slpdirrotate,1),slpdirrotatenorm)
			slpdirrotatenorm=norm2(slpdirrotate,dim=2)
			
			do i=1,size(slpdirrotate,1)
				lxtotal(i,:,:)=reshape((/nodex(grainb,int(nodesnotindex))   
     1         				 -rnodes(i,1),                                
     2          			 nodey(grainb,int(nodesnotindex))-rnodes(i,2),
     3         				 nodez(grainb,int(nodesnotindex))-rnodes(i,3) 
     4          			/),(/valb-1,3/))
				lxnormin=lxtotal(i,:,:)
				!******can be replaced with norm2 function
				!call enorm(lxnormin,valb-1,lxnorm)
				lxnorm=norm2(lxnormin,dim=2)
				!find the angle between the vectors from the shared boundary to 
				!boundary nodes in grain b and the vector developed from the voxel
				!to the shared boundary

      lxangle(i,:)= dacos(
     & matmul(lxtotal(i,:,:),slpdirrotate(i,:))/(lxnorm
     & *slpdirrotatenorm(i)))

				minxval0(i,:)=minloc(lxangle(i,:))
				minxval0act(i,:)=lxangle(i,int(minxval0(i,:)))
				minxval180(i,:)=minloc(abs(lxangle(i,:)-pi))
				minxval180act(i,:)=abs(pi-lxangle(i,int(minxval180(i,:))))
				
				!since the vectors in this case case be parallel but in opposite directions
				!an 'if' loop is used all sort through the 0deg and 180deg possible 
				!minimums
				if(minxval0act(i,1)<minxval180act(i,1))then   
					minxangle(i)=minxval0(i,1)
					minxangleact(i)=minxval0act(i,1)
				else 
					minxangle(i)=minxval180(i,1)
					minxangleact(i)=minxval180act(i,1)
				end if
				
			 end do
		 end if
		

		 	 
		 
		 allocate (lxnodesindex(size(slpdirrotate,1)))
		 allocate (lxnodes(size(slpdirrotate,1),3))
		 !using the location of the closest vector in grain B to the developed vector
		 !along grain B's slip systems in grain A, the vector and corresponding nodes
		 !can be determined.
		 lxnodesindex=nodesnotindex(minxangle)
		 lxnodes=reshape((/nodex(grainb,lxnodesindex),               
     1  		nodey(grainb,lxnodesindex),                             
     2  		nodez(grainb,lxnodesindex)/),(/size(slpdirrotate,1),3/))

		 
		 !need to check the angle formed by the lxnodes vector and the slip vector
		 !in grain A.  If the difference is large (greater than a specified minimum
		 !allowable angle), the slip distance in grian B should be ignored since it
		 !is unreliable
		 allocate (vecb(size(slpdirrotate,1),3))
		 allocate (veca(size(slpdirrotate,1),3))
		 if (locgrain .eqv. .FALSE.)then
			!vector from shared bounday nodes to 
			vecb=lxnodes-rnodes
			!vector from voxel to shared boundary
			do i=1,size(slpdirrotate,1)
				veca(i,:)=elcent(int(noel),:)-rnodes(i,:)
			end do
		  else
			!vector from shared bounday nodes to 
			vecb=lxnodes-rnodes
			!voxel at the bounary so use slip direction in grain B
			veca=slpdirrotate
		  end if
		  
		  !calculate the angle between vector a and b
		  allocate (vanorm(size(slpdirrotate,1)))
		  allocate (vbnorm(size(slpdirrotate,1)))
		  allocate (checkangle(size(slpdirrotate,1)))
		  !call enorm(veca,size(slpdirrotate,1),vanorm)
		  vanorm=norm2(veca,dim=2)
		  !call enorm(vecb,size(slpdirrotate,1),vbnorm)
		   vbnorm=norm2(vecb,dim=2)
		   ldistance=vbnorm
!        -------------------------------------------------------------------------
!                 This section is added to check the angle to see if it is
!                 within a tolerance. This can be updated further by the user
!        -------------------------------------------------------------------------
!		  do i=1,size(slpdirrotate,1)
!			checkangle(i)=dacos(dot_product(veca(i,:),vecb(i,:)) 
!     1      			 /(vanorm(i)*vbnorm(i)))
		  !check to see if the angle is below the recommended minimum
!			if(checkangle(i)<abs(pi-(40.D0*pi/180.D0)))then
			   !lxnodes(i,:)=rnodes(i,:)
		   !distance is taken only as half the characteristic voxel length
			   
			   !ldistance(i)=charcterlengt*0.5D0
!			 else
!			   ldistance(i)=vbnorm(i)
!			 end if
!		  end do
!-----------------------------------------------------------------------------------
		  
		  !calculate the Luster-Morris parameter
		  allocate (slpplanrotatenorm(size(slpplanrotate,1)))
		  allocate (slpdirrotateanorm(size(slpdir,1)))
		  allocate (slpplanrotateanorm(size(slpplane,1)))
		  
		  !call enorm(slpplanrotate,size(slpplanrotate,1),slpplanrotatenorm)
		  slpplanrotatenorm=norm2(slpplanrotate,dim=2)
		  !call enorm(slpdirrotatea,size(slpdir,1),slpdirrotateanorm)
		   slpdirrotateanorm=norm2(slpdirrotatea,dim=2)
		  !call enorm(slpplanrotatea,size(slpplane,1),slpplanrotateanorm)
		   slpplanrotateanorm=norm2(slpplanrotatea,dim=2)
		  
		  allocate (cosphi(size(slpplane,1)))
		  allocate (cosk(size(slpplane,1)))
	
		  do i=1,size(slpplane,1)
			cosphi(i)=dot_product(slpplanrotate(i,:),slpplanrotatea(i,:))/   
     1        			(slpplanrotatenorm(i)*slpplanrotateanorm(i))
		   
			cosk(i)=dot_product(slpdirrotatea(i,:),slpdirrotate(i,:))/       
     1       			 (slpdirrotateanorm(i)*normslp(i))
			lm(i)=abs(cosphi(i))*abs(cosk(i))
		  end do

c		 ldistance=ldistance*(10.D0**-6.D0)
c		 rdistance=rdistance*(10.D0**-6.D0)

           
c         coords, nodex, nodey, nodez are all in mm
c         Therefore noscaling is used.


		  !need to make sure to deallocate arrays
      deallocate (vect,mindist  ,sortmindist, mk, 
     & minindex,sortedbgrains,  
     & uniqgraindist,slpdirrotate,slpplanrotate,slpdirrotatea,   
     & slpplanrotatea,grainboundindex,featureboundsnodes,        
     & bxvalue,byvalue,bzvalue,btotal,bangle,normarray,normslp,  
     & minangleval,minang180,minangleval180,minang0,             
     & minangleactual,minangleactloc,boundnodegrainb,rnodes,     
     & nodesnotindex,nodesnot,lxtotal,lxnorm,lxangle,lxnormin,   
     & minxval0,minxval180,minxval180act,minxval0act,minxangle,  
     & minxangleact,vectrnorm,slpdirrotatenorm,lxnodesindex,     
     & lxnodes,vecb,veca,vanorm,vbnorm,checkangle,     
     & slpplanrotatenorm,slpdirrotateanorm,slpplanrotateanorm,   
     & cosphi,cosk)		  
	  
	
	return
	end subroutine grainsize
	

	!----------------------------------
	subroutine remove_dups(sortedbgrains,arraysize, uniqu,uniindex)
      implicit none
      integer :: k, i, j,arraysize
      real*8, allocatable :: uniqu(:),uniindex(:)
      real*8 :: sortedbgrains(arraysize),res(arraysize),
     &savindex(arraysize)
		 


		  k = 1
		  res(1) = sortedbgrains(1)
		  savindex(1)=1
		  outer: do i=2,size(sortedbgrains)
			 do j=1,k
				if (res(j) .eq. sortedbgrains(i)) then
				   ! Found a match so start looking again
				   cycle outer
				end if
			 end do
			 ! No match found so add it to the output
			 k = k + 1
			 res(k) = sortedbgrains(i)
			 savindex(k)=i
		  end do outer

		
		

		  allocate (uniqu(1:k))
		  allocate (uniindex(1:k))
		  
		  uniqu=res(1:k)
		  uniindex=savindex(1:k)
            return
            
	end subroutine remove_dups


      
      
      
      
      
      
      
      
      
      end module lengthscale