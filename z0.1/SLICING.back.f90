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

open(unit=1,file='Data_all_mass.dat')
read(1,*)
do row=1,N
    read(1,*)id(row),x(row),y(row),z(row),SM(row),Metal(row),ur(row),gr(row)
   
end do



Pos_data(:,1) = x(:)
Pos_data(:,2) = y(:)
Pos_data(:,3) = z(:)




open(unit=2,file='Data_fortran_Cluster.dat')


do row=1,clus
    read(2,*)xclus(row),yclus(row), zclus(row), Mclus(row), Rclus(row)
end do


call Combinationset(clus,2,xclus,yclus,zclus)

open(unit=3,file='Cluster_pairs.dat',status='old',iostat=stat)

pcl=1

	do row=1,17
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REDING CLUSTER FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CALLING DISPERSE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		  
		  call system('/home/ankit/Python_Environments/EAGLE/DisPerSE/bin/delaunay_2D '//filename//&
			      &'-btype '//'periodic '//'-outDir '//directory)
		  NDfile = trim(filename)//'.NDnet'
		  call StripSpaces(NDfile)
		  call system('/home/ankit/Python_Environments/EAGLE/DisPerSE/bin/mse '//NDfile//' -nsig'//' 3'//&
			      &' -upSkl'//' -forceLoops'//' -outDir '//directory)
		  SKLfile = trim(NDfile)//'_s3.up.NDskl'
		  call StripSpaces(SKLfile)
		  call system('/home/ankit/Python_Environments/EAGLE/DisPerSE/bin/skelconv '//SKLfile//' -breakdown'//&
		  &' -assemble'//' 60'//' -to'//' vtp'//' -outDir '//directory)
		  call system('/home/ankit/Python_Environments/EAGLE/DisPerSE/bin/skelconv '//SKLfile//&
		  &' -to'//' NDskl_ascii'//' -outDir '//directory)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CREATING DIFFEENT FILAMENT FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		  FILAMENT_FILE = filename//'.NDnet_s3.up.NDskl.a.NDskl'
		  call StripSpaces(FILAMENT_FILE)
		  FILDAT = directory//'/'//'FILAMENTS.dat'
		  call StripSpaces(FILDAT)
		  
		  
		  open(unit=11,file=FILAMENT_FILE,status='old',iostat=stat)
		  open(unit=12,file=FILDAT,status='replace',action='write')
		  
		  do
		    if(stat /= 0) then
		      exit
		    else
			  read(11,'(A)',iostat=stat)line
			  
			  !write(*,*)line
			  if(line=='[FILAMENTS]') then
			    write(*,*)line
			    copy=1
			  else if(line=='[CRITICAL POINTS DATA]') then
			    write(*,*)line
			    copy=0
			    exit
			  else
			    continue
			  
			  if(copy==1) then
			    write(12,*)line
			  end if
			  
			  end if
			  
		    end if
		      
		  end do
		  
		  
		  close(11)
		  close(12)
		  
		  
		  
		  
		  open(unit=12,file=FILDAT,status='old',action='read',iostat=stat)
		  
		  read(12,*,iostat=stat)total_filaments
		  copy=0
		  
		  i=1
		  do 
		    
		    if(stat/=0) then
			exit
		    
		    else
			read(12,'(A)',iostat=stat)line
			
			if(line(1:2)=='  ') then
			  
			      
			      write(99,*)line
			      k=k+1
			      
			else  
			      
			      write(filno,'(I3)')i
			      call StripSpaces(filno)
			      newfil = trim(directory)//'/'//'filament'//trim(filno)//'.dat'
			      open(unit=99,file=newfil,status='replace')
			      i=i+1
			  
			end if
		    
		    end if
		  
		  end do
		  
		  
		  close(12)
		  close(99)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CREATING LENGTH.dat FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  
		  
		  open(unit=98,file=trim(directory)//'/'//'Lengths.dat',action='write',status='replace')

		  do i=1,total_filaments
		    
			write(filno,'(I3)')i
			call StripSpaces(filno)
			newfil = trim(directory)//'/'//'filament'//trim(filno)//'.dat'
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
		 deallocate(Slice_array)
  
	end do


		  
		  
		  
	

close(1)
close(2)
close(3)
END PROGRAM SLICING







