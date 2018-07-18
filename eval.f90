module global
    implicit none
    integer, parameter :: N = 3 !order of matrix
    real*8, dimension(N,N) :: A, B, C

end module global

program lapack 
    use global
    implicit none
    
    !A(1,1) = 1
    !A(1,2) = 0
    !A(1,3) = 0
    !A(2,1) = 0 
    !A(2,2) = 0 
    !A(2,3) = 1 
    !A(3,1) = 0 
    !A(3,2) = 1 
    !A(3,3) = 0 
   
    A(1,1) = 1
    A(1,2) = -3
    A(1,3) = 3
    A(2,1) = 3 
    A(2,2) = -5 
    A(2,3) = 3 
    A(3,1) = 6 
    A(3,2) = -6 
    A(3,3) = 4 

    print*, "A=" 
    print*,"|", A(1,1:3) ,"|"
    print*,"|", A(2,1:3) ,"|"
    print*,"|", A(3,1:3) ,"|"

    call eval(A)

    !A(1,1) = 1
    !A(1,2) = -3
    !A(1,3) = 3
    !A(2,1) = 3 
    !A(2,2) = -5 
    !A(2,3) = 3 
    !A(3,1) = 6 
    !A(3,2) = -6 
    !A(3,3) = 4 

    !call decomp(A)
    
end program lapack 

!computes eigenvalues of a real double matrix
!over writes intial matrix
subroutine eval(M)
    use global
    implicit none
    !--------------------lapack--------------------------- 
    !integer, parameter :: lwmax = 1000
    integer, parameter :: lda = N !ld is leading dimension
    integer, parameter :: ldvr = N
    integer, parameter :: ldvl = N
    real*8, dimension(ldvl,N) :: vl 
    real*8, dimension(ldvr,N) :: vr 
    real*8, dimension(N) :: wr 
    real*8, dimension(N) :: wi 
    integer :: lwork = 4*N
    real*8, dimension(N*4) :: work 
    character*1 :: jobvl = 'N' !compute left eigenvectors : N = no, V = yes
    character*1 :: jobvr = 'V' ! ""     right "" 
    integer :: info ! if = 0 then successful, if = - i, ith paramter had illegal value

    !--------------------eval--------------------------- 
    real*8, dimension(N,N) :: M ! input matrix 

    vl = 0 
    vr = 0 
    wr = 0
    wi = 0

    call dgeev(jobvl, jobvr, N, M, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

    if( info.gt.0 ) then
         write(*,*)'the algorithm failed to compute eigenvalues.'
          stop
    end if

    !print*, work(1)

    print*, ""
    print*, "----------------EValues--------------------"
    print*,"re eigenvalues: ",  wr(:)
    print*,"im eigenvalues: ", wi(:)
    print*, ""
    print*, "-----------------EVectors--------------------"
    print*, vr(1,:)
    print*, vr(2,:)
    print*, vr(3,:)

    print*,sqrt(vr(1,3)**2+vr(2,3)**2+vr(3,3)**2)
    !print*,"left eigenvector: ", vl(1,:)
    !print*,"left eigenvector: ", vl(2,:)
    !print*,"left eigenvector: ", vl(3,:)

End subroutine eval

subroutine decomp(M)
    use global
    implicit none

    integer :: i, j 
    real*8, dimension(N,N) :: M ! input matrix 
    real*8, dimension(N,N) :: M_t ! transpose of M.
    real*8, dimension(N,N) :: M_sym ! symmetric part of M.
    real*8, dimension(N,N) :: M_asym ! anti symmetric part of M

    do i = 1, N 
        do j = 1, N 
            M_t(i,j) = M(j,i)
        end do
    end do

    M_sym = 0.5_8 * M + M_t
    M_asym = 0.5_8 * M - M_t
    
    print*, ""
    print*, "-----------------Matrix decomposition--------------------"
    print*, "M"
    print*,"|", M(1,1:3) ,"|"
    print*,"|", M(2,1:3) ,"|"
    print*,"|", M(3,1:3) ,"|"

    print*, ""
    print*, "M_t"
    print*,"|", M_t(1,1:3) ,"|"
    print*,"|", M_t(2,1:3) ,"|"
    print*,"|", M_t(3,1:3) ,"|"

    print*, ""
    print*, "M_sym"
    print*,"|", M_sym(1,1:3) ,"|"
    print*,"|", M_sym(2,1:3) ,"|"
    print*,"|", M_sym(3,1:3) ,"|"

    M = M_sym + M_asym

    print*, ""
    print*, "M = M_sym + M_asym"
    print*,"|", M(1,1:3) ,"|"
    print*,"|", M(2,1:3) ,"|"
    print*,"|", M(3,1:3) ,"|"
End subroutine decomp 


