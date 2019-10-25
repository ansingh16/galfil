PROGRAM MATCHING
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





open(unit=3,file='Cluster_pairs.dat',status='old',iostat=stat)

pcl=1

	do row=1,17
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REDING CLUSTER FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		  read(3,*)dum1,dum2,clus1(1),clus1(2),clus1(3),&
			    & clus2(1),clus2(2),clus2(3),dist
		  write(*,*)dum1,dum2
		  
		  origin_shift = (clus1 + clus2)/2.0
		  origin_shiftx = origin_shift(1);origin_shifty = origin_shift(2);origin_shiftz = origin_shift(3)
		  
		  
		  clus1 = clus1 - origin_shift
		  clus2 = clus2 - origin_shift
		  
		  cluster_vec = clus2 - clus1
		  
		  !print *,cluster_vec
		  
		  normal(1) = rand()
		  normal(2) = rand()
		  normal(3) = -(cluster_vec(1)*normal(1) + cluster_vec(2)*normal(2))/cluster_vec(3)
		  
		  Rot_mat = ROTATION(normal)

		  
		  
		  write(no,'(I2)') row

		  call StripSpaces(no)
		  
		  directory = trim('SLICED_')//trim(no)
		  filename = './'//trim(directory)//'/'//'Clipped.dat'
		  call StripSpaces(filename)
		  
		  
		  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TRANSFORMING BACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		  inquire (file=directory//'/'//'Groups.dat',exist=file_ex)
		  if(dir_ex) then
		      print*,"yes"
		  else
		      exit
		      
		  endif
		  

		  
		  call inverse(Rot_mat,Rot_inv)
		  
		  call dgemm('N','N',3,N,3,1.0d0,Rot_inv,3,temar2,3,0.0d0,temar1,3)
		  
		  
		  trans_array = transpose(temar1)
		  
		  trans_array(:,1) = trans_array(:,1) + origin_shiftx
		  trans_array(:,2) = trans_array(:,2) + origin_shifty
		  trans_array(:,3) = trans_array(:,3) + origin_shiftz
		  
		  
		  
		  
		  !call print_array(trans_array,10,3)
		  
		  print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		  !call print_array(Pos_data,10,3)
		  
		  pcl=pcl+1
		  
	  
	end do


		  
		  
		  
	
END PROGRAM MATCHING







