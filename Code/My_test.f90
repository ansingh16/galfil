Program My_test


real(kind=8)::normal(3),Rot_mat(3,3)







end program My_test






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
      
      Rot_mat = get_R(r_cos,r_sin,u,v,w,uvw)
      
      !call dgemm('N','N',dim1,dim3,dim2,alpha,A,dim1,B,dim2,beta,C,dim1)
      
      
      
    END FUNCTION
