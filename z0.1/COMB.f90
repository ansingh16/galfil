MODULE COMB

  contains
  
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
      SUBROUTINE inverse(A,Ainv)
	  real(kind=8), dimension(:,:), intent(in) :: A
	  real(kind=8), dimension(size(A,1),size(A,2)) :: Ainv

	  real(kind=8), dimension(size(A,1)) :: work  ! work array for LAPACK
	  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
	  integer :: n, info

	  ! External procedures defined in LAPACK
	  external DGETRF
	  external DGETRI

	  ! Store A in Ainv to prevent it from being overwritten by LAPACK
	  Ainv = A
	  n = size(A,1)

	  ! DGETRF computes an LU factorization of a general M-by-N matrix A
	  ! using partial pivoting with row interchanges.
	  call DGETRF(n, n, Ainv, n, ipiv, info)

	  if (info /= 0) then
	    stop 'Matrix is numerically singular!'
	  end if

	  ! DGETRI computes the inverse of a matrix using the LU factorization
	  ! computed by DGETRF.
	  call DGETRI(n, Ainv, n, ipiv, work, n, info)

	  if (info /= 0) then
	    stop 'Matrix inversion failed!'
	  end if
    END SUBROUTINE inverse
  
  
  
    SUBROUTINE print_array(array,n,m) 
      
      implicit none
      real(kind=8), intent(in) :: array(n,m)
      integer, intent(in) :: n,m
      integer :: i
      do i = 1,n
      print*, array(i,:)
      end do

    END SUBROUTINE
  
    
    FUNCTION cross_vec(a, b) result(cross)
      double precision, DIMENSION(3) :: cross(3)
      double precision, DIMENSION(3), INTENT(IN) :: a(3), b(3)

      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
    END FUNCTION cross_vec
      
    
    FUNCTION norm_vec(V1,n1) result(norm)
	    double precision:: norm,V1(3)
	    integer:: n1
	    norm = sqrt ( sum ( V1(:n1)*V1(:n1) ))
    END FUNCTION
  
  
    FUNCTION get_R(rcos,rsin,u,v,w,uvw) result(MAT)
	 double precision::rcos,rsin,u,v,w,uvw(3),zero=0.0
	 double precision, DIMENSION(3, 3) :: MAT,Iden,TMP,TMP1
	 
	 Iden = reshape((/ 1.0,0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /), shape(MAT))
	 TMP = reshape((/ zero,w, -v, -w, zero, u, v, -u, zero /), shape(MAT))
	 
	 TMP1 = r_cos*Iden(:,:) + rsin*TMP(:,:)
	 MAT = TMP1(:,:) + 1.0
	
    END FUNCTION
  
    RECURSIVE SUBROUTINE combinations(n,r,k,column,X1,Y1,Z1,fp)
    implicit none
    integer::column(r),k
    integer,intent(in)::n,r
    integer::j,fp
    real(kind=8),intent(in)::X1(:),Y1(:),Z1(:)
    real(kind=8)::d
      column(k) = column(k-1)
      
      do while(column(k)<n-r+k)
	column(k)=column(k)+1
	
	if(k<r) then
	  call combinations(n,r,k+1,column,X1,Y1,Z1,fp)
	else 
	  do j=1,r
	      write(*,*)column(j)
	  end do
	  d = SQRT((X1(column(1))-X1(column(2)))**2.0 + (Y1(column(1))-Y1(column(2)))**2.0 + (Z1(column(1))-Z1(column(2)))**2.0)
	  if((5.0<d).and.(d<50.0)) then
	  !write(*,*)"Above works",d
	  write(fp,*)column(1),column(2),X1(column(1)),Y1(column(1)),Z1(column(1)),&
	  & X1(column(2)),Y1(column(2)),Z1(column(2)),d
	  !write(*,*)column(1),column(2),X1(column(1)),Y1(column(1)),Z1(column(1)),&
	  !& X1(column(2)),Y1(column(2)),Z1(column(2)),d

	  end if
	  
	  write(*,*)
	end if
	
      end do
      
      
    END SUBROUTINE
    
  
  
  
    SUBROUTINE Combinationset(n,r,X1,Y1,Z1)
    IMPLICIT NONE
    real(kind=8),intent(in)::X1(:),Y1(:),Z1(:)
    integer::n,r
    integer::column(r),k
    open(unit=3,file='Cluster_pairs.dat',status='replace',action='write')

      column(1)=0
      k=1
      
      call combinations(n,r,k,column,X1,Y1,Z1,3)

    close(3)
    END SUBROUTINE Combinationset

    
    
    FUNCTION ROTATION(i_v) result(Rot_mat)
    
    implicit none
      integer::spacedim=3
      double precision::i_v(3),unit(3)
      !integer::dim1, dim2, dim3
      !parameter::(dim1 = 2, dim2 = 2, dim3 = 2)
      !double precision::A(dim1, dim2), B(dim2, dim3), C(dim1, dim3)
      double precision::alpha, beta,uvw(3),r_cos,r_sin,mag,u,v,w
      double precision::Rot_mat(3,3)
      
      unit(1)=0.0
      unit(2)=0.0
      unit(3)=1.0
      
      
      alpha = 1.0d0
      beta = 0.0d0
      
      mag = norm_vec(i_v,spacedim)   ! L2 norm
      
      i_v = i_v/mag
      
      uvw = cross_vec(i_v,unit)
      
      r_cos = dot_product(i_v,unit)
      r_sin = norm_vec(uvw,spacedim)

  
      if(r_sin .ne. 0.0)then
	    uvw = uvw/r_sin
      end if
      
      u=uvw(1);v=uvw(2);w=uvw(3);
      

        Rot_mat(1,1) =      r_cos + u*u*(1-r_cos);
	Rot_mat(2,1) =  w * r_sin + v*u*(1-r_cos);
	Rot_mat(3,1) = -v * r_sin + w*u*(1-r_cos);
	Rot_mat(1,2) = -w * r_sin + u*v*(1-r_cos);
	Rot_mat(2,2) =      r_cos + v*v*(1-r_cos);
	Rot_mat(3,2) =  u * r_sin + w*v*(1-r_cos);
	Rot_mat(1,3) =  v * r_sin + u*w*(1-r_cos);
	Rot_mat(2,3) = -u * r_sin + v*w*(1-r_cos);
	Rot_mat(3,3) =      r_cos + w*w*(1-r_cos);
	





      !Rot_mat = get_R(r_cos,r_sin,u,v,w,uvw)
      
      !call dgemm('N','N',dim1,dim3,dim2,alpha,A,dim1,B,dim2,beta,C,dim1)
      
      
      
    END FUNCTION
    
    
    subroutine StripSpaces(string)
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

    end subroutine

    
    
END MODULE COMB
