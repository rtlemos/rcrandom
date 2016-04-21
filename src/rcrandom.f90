! For standard compilation, use
! R CMD SHLIB mrg32k3a.f90
!
! For extra optimization, use
! gfortran -fpic -g -O3  -c  mrg32k3a.f90 -o mrg32k3a.o
! gfortran -shared -Wl,-Bsymbolic-functions -Wl,-z,relro -o mrg32k3a.so mrg32k3a.o -L/usr/lib/R/lib -lR

!******************************************************************
!* Primary subroutines ********************************************
!******************************************************************

subroutine U01(n, s, Cg, anti, uu)
    implicit none
    integer,intent(in)                  :: n !number of variates per substream
    integer,intent(in)                  :: s !number of substreams
    real(8),dimension(6*s),intent(inout):: Cg
    logical,intent(in)                  :: anti
    real(8),dimension(n*s),intent(out)  :: uu
    real(8),dimension(s,6)              :: Cp
    real(8),dimension(6*s,1)            :: Cpt
    real(8),dimension(s,n)              :: up
    real(8),dimension(n*s,1)            :: upt
    integer                             :: i
    real(8),dimension(s)                :: p1, p2, u, k
    real(8),parameter                   :: norm=0.0000000002328306549295727688d0
    real(8),parameter                   :: m1=4294967087.0d0
    real(8),parameter                   :: m2=4294944443.0d0
    real(8),parameter                   :: a12=1403580.0d0
    real(8),parameter                   :: a13n=810728.0d0
    real(8),parameter                   :: a21=527612.0
    real(8),parameter                   :: a23n=1370589.0

    Cp = transpose( reshape(Cg, (/6,s/) ))
    do i=1,n
        !Component 1
        p1 = a12 * Cp(:,2) - a13n * Cp(:,1)
        k = int( p1 / m1 )
        p1 = p1 - k * m1
        where( p1 < 0d0) p1 = p1 + m1
        Cp(:,1) = Cp(:,2)
        Cp(:,2) = Cp(:,3)
        Cp(:,3) = p1

        !Component 2
        p2 = a21 * Cp(:,6) - a23n * Cp(:,4)
        k = int( p2 / m2 )
        p2 = p2 - k * m2
        where(p2 < 0d0) p2 = p2 + m2
        Cp(:,4) = Cp(:,5)
        Cp(:,5) = Cp(:,6)
        Cp(:,6) = p2

        !Combination
        where( p1 > p2)
            u = (p1 - p2) * norm
        elsewhere
            u = (p1 - p2 + m1) * norm
        end where
        if( anti ) then
            up(:,i) = 1d0 - u
        else
            up(:,i) = u
        endif
    enddo
    Cpt = reshape(transpose(Cp), (/6*s,1/))
    upt = reshape(transpose(up), (/n*s,1/))
    Cg = Cpt(:,1)
    uu = upt(:,1)
end subroutine

subroutine U01d(n, s, Cg, anti, uu)
    implicit none
    integer,intent(in)                  :: n !number of variates per substream
    integer,intent(in)                  :: s !number of substreams
    real(8),dimension(6*s),intent(inout):: Cg
    logical,intent(in)                  :: anti
    real(8),dimension(n*s),intent(out)  :: uu
    real(8),dimension(n*s)              :: u, u1, u2
    integer                             :: i
    real(8),parameter                   :: fact=0.000000059604644775390625d0

    do i=1,n
        call U01(1, s, Cg, anti, u1(i))
        call U01(1, s, Cg, anti, u2(i))
    enddo
    if( .not. (anti) ) then
        u = u1 + u2 * fact
        where(u < 1d0)
            uu = u
        elsewhere
            uu = u - 1d0
        end where
    else !Don't forget that U01() returns 1 - u in the antithetic case
        u = u1 + (u2 - 1.0d0) * fact
        where(u < 0.0d0)
            uu = u + 1.0d0
        elsewhere
            uu = u
        end where
    endif
end subroutine

!******************************************************************
!* Secondary subroutines ******************************************
!******************************************************************

subroutine get_advanced_seed(s)
    implicit none
    real(8),dimension(6), intent(inout) :: s
    real(8),parameter                   :: m1=4294967087.0d0
    real(8),parameter                   :: m2=4294944443.0d0
    real(8),dimension(3,3),parameter    :: A1p127 =reshape(  (/ 2427906178.0d0,  226153695.0d0,  1988835001.0d0,   &
                                                                3580155704.0d0, 1230515664.0d0,   986791581.0d0,   &
                                                                 949770784.0d0, 3580155704.0d0,  1230515664.0d0/), &
                                                                 (/ 3,3 /) )
    real(8),dimension(3,3),parameter    :: A2p127 =reshape(  (/ 1464411153.0d0,   32183930.0d0,  2824425944.0d0,   &
                                                                 277697599.0d0, 1464411153.0d0,    32183930.0d0,   &
                                                                1610723613.0d0, 1022607788.0d0,  2093834863.0d0/), &
                                                                 (/ 3,3 /) )
    call MatVecModM( A1p127, s(1:3), s(1:3), m1 )
    call MatVecModM( A2p127, s(4:6), s(4:6), m2 )
end subroutine

subroutine get_advanced_subseed( subseed, n_jumps )
    implicit none
    real(8),dimension(6),intent(inout)  :: subseed
    integer,intent(in)                  :: n_jumps
    integer                             :: i
    real(8),parameter                   :: m1=4294967087.0d0
    real(8),parameter                   :: m2=4294944443.0d0
    real(8),dimension(3,3),parameter    :: A1p76 = reshape((/     82758667.0d0, 3672831523.0d0,  3672091415.0d0,   &
                                                                 1871391091.0d0,   69195019.0d0,  3528743235.0d0,   &
                                                                 4127413238.0d0, 1871391091.0d0,    69195019.0d0/), &
                                                                 (/ 3,3 /) )
    real(8),dimension(3,3),parameter    :: A2p76 = reshape((/    1511326704.0d0, 4292754251.0d0,  3859662829.0d0,   &
                                                                 3759209742.0d0, 1511326704.0d0,  4292754251.0d0,   &
                                                                 1610795712.0d0, 3889917532.0d0,  3708466080.0d0/), &
                                                                 (/ 3,3 /) )
    do i=1,n_jumps
        call MatVecModM (A1p76, subseed(1:3), subseed(1:3), m1)
        call MatVecModM (A2p76, subseed(4:6), subseed(4:6), m2)
    enddo
end subroutine

subroutine get_advanced_state(ee, cc, Cg)
    implicit none
    real(8),dimension(6), intent(inout) :: Cg
    integer,intent(in)                  :: cc, ee
    integer                             :: c,e
    real(8),dimension(3,3)              :: B1, B2, C1, C2
    real(8),parameter                   :: m1=4294967087.0d0
    real(8),parameter                   :: m2=4294944443.0d0
    real(8),dimension(3,3),parameter    :: A1p0 =  reshape(  (/       0.0d0,          0.0d0,     -810728.0d0,   &
                                                                      1.0d0,          0.0d0,     1403580.0d0,   &
                                                                      0.0d0,          1.0d0,           0.0d0/), &
                                                                      (/ 3,3 /) )
    real(8),dimension(3,3),parameter    :: A2p0 =  reshape(  (/       0.0d0,          0.0d0,    -1370589.0d0,   &
                                                                      1.0d0,          0.0d0,           0.0d0,   &
                                                                      0.0d0,          1.0d0,      527612.0d0/), &
                                                                      (/ 3,3 /) )
    real(8),dimension(3,3),parameter    :: InvA1 = reshape( (/184888585.0d0,           1.0d0,           0.0d0,   &
                                                                      0.0d0,           0.0d0,           1.0d0,   &
                                                             1945170933.0d0,           0.0d0,           0.0d0/), &
                                                             (/ 3,3 /) )
    real(8),dimension(3,3),parameter    :: InvA2 = reshape((/         0.0d0,          1.0d0,           0.0d0,   &
                                                              360363334.0d0,          0.0d0,           1.0d0,   &
                                                             4225571728.0d0,          0.0d0,           0.0d0/), &
                                                             (/ 3,3 /) )
    if( ee>0 ) then
        e=ee
        call MatTwoPowModM (A1p0, B1, m1, e)
        call MatTwoPowModM (A2p0, B2, m2, e)
    elseif( ee<0 ) then
        e=-ee
        call MatTwoPowModM (InvA1, B1, m1, e)
        call MatTwoPowModM (InvA2, B2, m2, e)
    endif

    if( cc >= 0 ) then
        c=cc
        call MatPowModM (A1p0, C1, m1, c)
        call MatPowModM (A2p0, C2, m2, c)
    else
        c=-cc
        call MatPowModM (InvA1, C1, m1, c)
        call MatPowModM (InvA2, C2, m2, c)
    endif

    if( ee .ne. 0 ) then
        call MatMatModM (B1, C1, C1, m1)
        call MatMatModM (B2, C2, C2, m2)
    endif

    call MatVecModM (C1, Cg(1:3), Cg(1:3), m1)
    call MatVecModM (C2, Cg(4:6), Cg(4:6), m2)
end subroutine

subroutine get_next_substream(Bg)
    real(8),dimension(6),intent(inout)  :: Bg
    real(8),parameter                   :: m1=4294967087.0d0
    real(8),parameter                   :: m2=4294944443.0d0
    real(8),dimension(3,3),parameter    :: A1p76 = reshape((/     82758667.0d0, 3672831523.0d0,  3672091415.0d0,   &
                                                                 1871391091.0d0,   69195019.0d0,  3528743235.0d0,   &
                                                                 4127413238.0d0, 1871391091.0d0,    69195019.0d0/), &
                                                                 (/ 3,3 /) )
    real(8),dimension(3,3),parameter    :: A2p76 = reshape((/    1511326704.0d0, 4292754251.0d0,  3859662829.0d0,   &
                                                                 3759209742.0d0, 1511326704.0d0,  4292754251.0d0,   &
                                                                 1610795712.0d0, 3889917532.0d0,  3708466080.0d0/), &
                                                                 (/ 3,3 /) )
    call MatVecModM(A1p76, Bg(1:3), Bg(1:3), m1)
    call MatVecModM(A2p76, Bg(4:6), Bg(4:6), m2)
end subroutine get_next_substream

! Check that the seeds are legitimate values. Returns 0 if legal seeds, -1 otherwise
subroutine get_seed_check(seed,ok)
    implicit none
    integer,dimension(6),intent(in)     :: seed
    integer,intent(out)                 :: ok
    real(8),parameter                   :: m1=4294967087.0d0
    real(8),parameter                   :: m2=4294944443.0d0
    integer                             :: i

    ok=0
    do i=1,3
        if(seed(i) >= m1) then
            !write(*,*) 'ERROR: seed > m1'
            ok=-1
        endif
    enddo
    do i=4,6
        if(seed(i) >= m2) then
            !write(*,*) 'ERROR: seed > m2'
            ok=-1
        endif
    enddo
    if(seed(1) == 0 .and. seed(2)==0 .and. seed(3)==0) then
        !write(*,*) 'ERROR: first 3 seeds = 0'
        ok=-1
    endif
    if(seed(4) == 0 .and. seed(5)==0 .and. seed(6)==0) then
        !write(*,*) 'ERROR: last 3 seeds = 0'
        ok=-1
    endif
end subroutine

!******************************************************************
!* Auxiliary subroutines ******************************************
!******************************************************************

! Compute (a*s + c) % m. m must be < 2^35.  Works also for s, c < 0
subroutine MultModM (aa, s, c, m, v)
    implicit none
    real(8),intent(in)      :: aa, s, c, m
    real(8),intent(out)     :: v
    real(8)                 :: a
    real(8),parameter       :: two17=131072.0d0
    real(8),parameter       :: two53=9007199254740992.0d0
    integer                 :: a1

    a=aa
    v=a*s+c
    if( (v>= two53) .or. (v <= -two53)) then
        a1 = int(a/two17)
        a = a - a1 * two17
        v = a1 * s
        a1 = int(v/m)
        v = v - a1 * m
        v = v * two17 + a * s + c
    endif
    a1 = int(v/m)
    v = v - a1 * m
    if( v < 0d0) then
        v = v + m
    endif
end subroutine

! Returns v = A*s % m.  Assumes that -m < s[i] < m.
! Works even if v = s.
subroutine MatVecModM(A,s,v,m)
    implicit none
    real(8),dimension(3,3),intent(in)   :: A
    real(8),dimension(3),intent(inout)  :: s
    real(8),dimension(3),intent(out)    :: v
    real(8),intent(in)                  :: m
    integer                             :: i
    real(8),dimension(3)                :: x

    do i=1,3
        call MultModM(A(i,1),s(1),0d0,m, x(i))
        call MultModM(A(i,2),s(2),x(i),m,x(i))
        call MultModM(A(i,3),s(3),x(i),m,x(i))
    enddo
    v(:)=x(:)
end subroutine

! Returns C = A*B % m. Work even if A = C or B = C or A = B = C.
subroutine MatMatModM(A,B,C,m)
    implicit none
    real(8),dimension(3,3),intent(inout)      :: A
    real(8),dimension(3,3),intent(inout)      :: B
    real(8),dimension(3,3),intent(out)        :: C
    real(8),intent(in)                        :: m
    real(8),dimension(3,3)                    :: W
    real(8),dimension(3)                      :: V
    integer                                   :: i

    do i=1,3
        V(:)=B(:,i)
        call MatVecModM(A, V, V, m)
        W(:,i)=V(:)
    enddo
    C(:,:)=W(:,:)
end subroutine

! Compute matrix B = (A^(2^e) % m);  works even if A = B
subroutine MatTwoPowModM(A, B, m, e)
    implicit none
    real(8),dimension(3,3),intent(in)       :: A
    real(8),dimension(3,3),intent(inout)    :: B
    real(8),intent(in)                      :: m
    integer,intent(in)                      :: e
    integer                                 :: i

    B=A
    do i=0, e-1
        call MatMatModM (B, B, B, m)
    enddo
end subroutine

! Compute matrix B = A^n % m ;  works even if A = B
subroutine MatPowModM(A, B, m, p)
    implicit none
    real(8),dimension(3,3),intent(in)          :: A
    real(8),dimension(3,3),intent(inout)       :: B
    real(8),dimension(3,3)                     :: W
    real(8),intent(in)                         :: m
    integer,intent(in)                         :: p
    integer                                    :: n
    real(8),dimension(3,3),parameter           :: id3 =reshape((/1.0d0, 0.0d0, 0.0d0,   &
                                                                 0.0d0, 1.0d0, 0.0d0,   &
                                                                 0.0d0, 0.0d0, 1.0d0/), &
                                                                 (/ 3,3 /))
    W=A
    B=id3
    n=p
    ! Compute B = A^n % m using the binary decomposition of n
    do while(n>0)
        if( mod(n,2) .ne. 0 ) call MatMatModM (W, B, B, m)
        call MatMatModM (W, W, W, m)
        n = n/2
    end do
end subroutine
