PROGRAM SLICING
USE COMB

IMPLICIT NONE
integer,parameter::N=13394,clus=11
real(kind=8)::x(N),y(N),z(N),SM(N),Metal(N),ur(N),gr(N),Pos_data(N,3),origin_shift(3),origin_shiftx,origin_shifty,origin_shiftz
real(kind=8)::Rot_mat(3,3),temar1(3,N),temar2(3,N),Rot_inv(3,3)
real(kind=8)::xclus(11),yclus(11), zclus(11), Mclus(11), Rclus(11),clus1(3),clus2(3),dist,cluster_vec(3),trans_array(N,3),P_rot(N,3)
integer::row,col,i,M,stat,pcl,dum1,dum2,tot_points,j,copy=0,total_filaments,k
real(kind=8)::normal(3),xfil1,yfil1,xfil2,yfil2,Length
real(kind=8),allocatable::Slice_array(:,:)  
character:: no*4,directory*20,filename*40,NDfile*60,SKLfile*60,FILAMENT_FILE*80,line*100,FILDAT*80,newfil*80
logical::dir_ex,file_ex
character:: filno*3
integer::id(N),id_p
real(kind=8)::r_prime(3),r_p(3),TM(3)




pcl=1

	do row=1,36
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REDING CLUSTER FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CREATING LENGTH.dat FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  write(no,'(I2)') row
		  
		  call StripSpaces(no)
		  
		  directory = trim('SLICED_')//trim(no)



		  FILDAT = directory//'/'//'FILAMENTS.dat'
		  call StripSpaces(FILDAT)

		  
		  open(unit=12,file=FILDAT,status='old',action='read',iostat=stat)
		  
		 

		  read(12,*,iostat=stat)total_filaments

		  close(12)

		 
		  open(unit=98,file=trim(directory)//'/'//'Lengths.dat',action='write',status='replace')


		  do i=1,total_filaments

			
		    
			write(filno,'(I3)')i
			call StripSpaces(filno)
			newfil = trim(directory)//'/'//'filament_'//trim(filno)//'.dat'
			open(unit=99,file=newfil,status='old',action='read',iostat=stat)
			
						

			Length=0.0
			
			do
			  if(stat/=0)then
			      exit
			      
			  else
			      read(99,*,iostat=stat)xfil1,yfil1
			      read(99,*,iostat=stat)xfil2,yfil2
			      
			      Length = Length + sqrt((xfil2-xfil1)**2.0 + (yfil2-yfil1)**2.0)
			      
			  end if
			  
			end do
			close(99)
			
			write(98,*)Length
			
			
		  
		  end do
		  
		  close(98)

  
	end do


		  
		  
		  
	

close(1)
close(2)
close(3)
END PROGRAM SLICING







