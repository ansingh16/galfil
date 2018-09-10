Program Match

USE COMB
USE csv_file

IMPLICIT NONE
integer,parameter::N=13394,clus=11
real(kind=8)::x(N),y(N),z(N),SM(N),Metal(N),ur(N),gr(N),Pos_data(N,3),origin_shift(3),origin_shiftx,origin_shifty,origin_shiftz
real(kind=8)::Rot_mat(3,3),temar1(3,N),temar2(3,N),Rot_inv(3,3)
real(kind=8)::xclus(11),yclus(11), zclus(11), Mclus(11), Rclus(11),clus1(3),clus2(3),dist,cluster_vec(3),trans_array(N,3),P_rot(N,3)
integer::row,col,i,M,stat,pcl,dum1,dum2,tot_points,j,copy=0,total_filaments,k,stat2
real(kind=8)::normal(3),xfil1,yfil1,xfil2,yfil2,Length
real(kind=8),allocatable::Slice_array(:,:)  
character:: no*4,directory*20,filename*40,NDfile*60,SKLfile*60,FILAMENT_FILE*80,line*100,FILDAT*80,newfil*80,groupfile*80
logical::dir_ex,file_ex
character:: filno*3,header*80
integer(kind=8)::id(N),id_p
real(kind=8)::r_prime(3),r_p(3)






open(unit=3,file='Cluster_pairs.dat',status='old',iostat=stat)

pcl=1

	do row=1,17



		
		write(no,'(I2)') row
		  
		  call StripSpaces(no)
		  
		  directory = trim('SLICED_')//trim(no)
		  
		  read(3,*)dum1,dum2,clus1(1),clus1(2),clus1(3),&
			    & clus2(1),clus2(2),clus2(3),dist
		  
		  
		  origin_shift = (clus1 + clus2)/2.0
		  origin_shiftx = origin_shift(1);origin_shifty = origin_shift(2);origin_shiftz = origin_shift(3)
		  
		  
		  
		  open(999,file=trim(directory)//'/ROTATION'//trim(no)//'.dat')

		  read(999,*)Rot_mat

		  close(999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TRANSFORMING BACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		  
		call inverse(Rot_mat,Rot_inv)

		groupfile = trim(directory)//trim('/Groups.dat')
		call StripSpaces(groupfile)
		
		open(unit=1000,file=groupfile,status='old',action='read',iostat=stat)

		open(unit=2000,file=trim(directory)//trim('/Groups_transbk.dat'),status='replace',action='write')


		 	read(1000,'(A)',iostat=stat2)header
                        write(2000,*)header
    

			do
			  if(stat2 /= 0)then
				
			      exit
			      
			  else
			      read(1000,*,iostat=stat2)id_p,r_prime(1),r_prime(2),r_prime(3)
			      r_p = matmul(Rot_inv,r_prime)
			      r_p(1) = r_p(1) + origin_shiftx
			      r_p(2) = r_p(2) + origin_shifty
			      r_p(3) = r_p(3) + origin_shiftz
			      !write(*,*)id_p,r_prime(1),r_prime(2),r_prime(3)
			      write(2000,'(I15)')id_p,',',r_p(1),',',r_p(2),',',r_p(3)
			  end if
			  
			end do
			



		 	open(unit=3000,file=trim(directory)//trim('/Filaments.dat'),status='old',iostat=stat,action='read')
			open(unit=4000,file=trim(directory)//trim('/Filaments_transbk.dat'),status='replace',action='write')

			read(3000,*,iostat=stat)
			do
			  if(stat/=0)then
			      exit
			      
			  else
			      read(3000,*,iostat=stat)id_p,r_prime(1),r_prime(2),r_prime(3)
			      r_p = matmul(Rot_inv,r_prime)	
			      r_p(1) = r_p(1) + origin_shiftx
			      r_p(2) = r_p(2) + origin_shifty
			      r_p(3) = r_p(3) + origin_shiftz
			      write(4000,*)id_p,r_p(1),r_p(2),r_p(3)
			  end if
			  
			end do
			

		

		  


			
			open(unit=5000,file=trim(directory)//trim('/Field.dat'),status='old',iostat=stat,action='read')
			open(unit=6000,file=trim(directory)//trim('/Field_transbk.dat'),status='replace',action='write')

			read(5000,*,iostat=stat)
			do
			  if(stat/=0)then
			      exit
			      
			  else
			      read(5000,*,iostat=stat)id_p,r_prime(1),r_prime(2),r_prime(3)
			      r_p = matmul(Rot_inv,r_prime)
			      r_p(1) = r_p(1) + origin_shiftx
			      r_p(2) = r_p(2) + origin_shifty
			      r_p(3) = r_p(3) + origin_shiftz
			      write(6000,*)id_p,r_p(1),r_p(2),r_p(3)
			  end if
			  
			end do
			

		  

		  close(1000)
	      	  close(2000)
		  close(3000)
	      	  close(4000)
		  close(5000)
	      	  close(6000)

		  
		  
		  !call dgemm('N','N',3,N,3,1.0d0,Rot_inv,3,temar2,3,0.0d0,temar1,3)
		  
		  
		  !trans_array = transpose(temar1)
		  
		  !trans_array(:,1) = trans_array(:,1) + origin_shiftx
		  !trans_array(:,2) = trans_array(:,2) + origin_shifty
		  !trans_array(:,3) = trans_array(:,3) + origin_shiftz
		  
		  
		  
		  
		  !call print_array(trans_array,10,3)
		  
		  print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		  !call print_array(Pos_data,10,3)
		  
		  pcl=pcl+1
		  
		 
	end do

end program match
