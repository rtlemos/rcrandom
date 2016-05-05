subroutine continuous_betas_linear(npt, y, b0, b1)
    implicit none
    integer, intent(in)                      :: npt
    real(8), intent(in), dimension(npt)      :: y
    real(8), intent(inout), dimension(npt-1) :: b0, b1

    b0 = y(1:(npt-1))
    b1 = y(2:npt)-y(1:(npt-1))
end subroutine continuous_betas_linear

subroutine continuous_betas_quadratic(npt, x, y, b0, b1, b2)
    implicit none
    integer,intent(in)                      :: npt
    real(8),dimension(npt),intent(in)       :: x, y
    real(8),dimension(npt-1),intent(inout)  :: b0, b1, b2
    integer                                 :: i, m

    m         = (npt+1)/2 !position of the mode
    b0(m)     = y(m)
    b1(m)     = 0.0d0
    b2(m)     = y(m+1) - y(m)
    do i = m+1, npt-1
        b0(i) = y(i)
        b1(i) = ( b1(i-1) + 2.0d0 * b2(i-1) ) * (x(i+1)-x(i)) / (x(i)-x(i-1))
        b2(i) = y(i+1) - y(i) - b1(i)
    enddo
    do i = m-1, 1, -1
        b0(i) = y(i)
        b1(i) = 2.0d0 * (y(i+1) - y(i)) - b1(i+1) * (x(i+1)-x(i)) / (x(i+2)-x(i+1))
        b2(i) = y(i+1) - y(i) - b1(i)
    enddo
end subroutine continuous_betas_quadratic

subroutine continuous_betas_quadratic_old(npt, x, y, b0, b1, b2)
    implicit none
    integer,intent(in)                      :: npt
    real(8),dimension(npt),intent(in)       :: x, y
    real(8),dimension(npt-1),intent(inout)  :: b0, b1, b2
    integer                                 :: i
    real(8)                                 :: d

    b0(:) = y(1:(npt - 1))
    d = (x(3) - x(1)) / (x(2) - x(1))
    b1(1) = d / (d - 1) * (y(2) - y(3) / d ** 2 - y(1) * (1 - 1 / d ** 2))
    do i = 2, npt - 1
      b1(i) = (b1(i - 1) + 2 * (y(i) - y(i - 1) - b1(i - 1))) * &
        (x(i + 1) - x(i)) / (x(i) - x(i - 1))
    enddo
    b2(1:(npt - 1)) = y(2:npt) - y(1:(npt - 1)) - b1(1:(npt - 1))
end subroutine continuous_betas_quadratic_old

subroutine continuous_integrals_b2equal0(npt, beta0, beta1, xlow, xhigh, integrals)
    implicit none
    integer,intent(in)                      :: npt
    real(8),dimension(npt-1),intent(in)     :: beta0, beta1, xlow, xhigh
    real(8),dimension(npt-1),intent(out)    :: integrals

    integrals = exp(beta0) * (exp(beta1 * xhigh) - exp(beta1 * xlow) ) / beta1
end subroutine continuous_integrals_b2equal0

subroutine continuous_integrals(npt, beta0, beta1, beta2, xlow, xhigh, integrals)
    implicit none
    integer,intent(in)                      :: npt
    real(8),dimension(npt-1),intent(in)     :: beta0, beta1, beta2, xlow, xhigh
    real(8),dimension(npt-1),intent(out)    :: integrals
    real(8),dimension(npt-1)                :: imL, imU, dlt, sgB
    real(8),parameter                       :: pi = 3.141592653589793d0
    integer                                 :: i
    real(8)                                 :: erU, erL, twosqpi

    twosqpi = 2.0/sqrt(pi)
    sgB = sign(1.0d0, beta2)
    imL = -sgB * 0.5d0 * (2.0d0 * beta2 * xlow  + beta1) / sqrt(abs(beta2))
    imU = -sgB * 0.5d0 * (2.0d0 * beta2 * xhigh + beta1) / sqrt(abs(beta2))
    do i = 1, npt-1
        if (abs(beta2(i)) < 1d-16) then
          dlt(i) = 0d0
        else if (abs(beta1(i)/beta2(i)) < 100d0) then
            if(beta2(i) < 0.0d0) then
                call erfS(imL(i),erL)
                call erfS(imU(i),erU)
            else
                call erfiS(imL(i),erL)
                call erfiS(imU(i),erU)
            endif
            dlt(i) = erU - erL
        else
            if(beta2(i) < 0.0d0) then
                dlt(i) = exp(-imL(i))*twosqpi*(imU(i)-imL(i))* &
                    (1d0 - (imU(i)-imL(i))*imL(i))
            else if(beta2(i) > 0.0d0) then
                dlt(i) = exp(imL(i))*twosqpi*(imU(i)-imL(i))* &
                    (1d0 + (imU(i)-imL(i))*imL(i))
            else
              dlt(i) = 0d0
            endif
        endif
    enddo
    where (dlt == 0)
      integrals = 0
    elsewhere
      integrals = - 0.5d0 * sqrt(pi/abs(beta2)) * &
                exp( (4.0d0 * beta0 * beta2 - beta1**2)/(4.0d0 * beta2)) * dlt
    end where
end subroutine continuous_integrals

subroutine continuous_pdf(n, z, npt, x, dx, beta0, beta1, beta2, normct, dolog, pdf)
    implicit none
    integer,intent(in)                  :: n, npt
    logical, intent(in)                 :: dolog
    real(8),dimension(n),intent(in)     :: z
    real(8),dimension(npt),intent(in)   :: x
    real(8),dimension(npt-1),intent(in) :: dx, beta0, beta1, beta2
    real(8),intent(in)                  :: normct
    real(8),dimension(n),intent(inout)  :: pdf
    integer,dimension(n)                :: pos
    integer,dimension(1)                :: pvec
    real(8),dimension(n)                :: zscaled
    integer                             :: i
    real(8)                             :: b0, b1, b2

    if(dolog) then
        pdf = -50
    else
        pdf = 0
    endif
    do i = 1, n
        pvec   = minloc( abs(z(i) - x) )
        pos(i) = pvec(1)
    enddo
    where( x(pos) > z .and. pos > 1)
        pos = pos - 1
    end where
    where(pos < npt)
        zscaled = (z - x(pos)) / dx(pos)
    end where
    do i = 1, npt-1
        b0 = beta0(i)
        b1 = beta1(i)
        b2 = beta2(i)
        if(dolog) then
            where( pos == i )
                pdf = b0 + b1 * zscaled + b2 * zscaled**2 - log(normct)
            end where
        else
            where( pos == i )
                pdf = exp(b0 + b1 * zscaled + b2 * zscaled**2) / normct
            end where
        endif
    enddo

end subroutine continuous_pdf

subroutine continuous_cdf_quadratic(n, z, npt, x, dx, beta0, beta1, beta2, pcdf, normct, p)
    implicit none
    integer,intent(in)                  :: n, npt
    real(8),dimension(n),intent(in)     :: z
    real(8),dimension(npt),intent(in)   :: x
    real(8),dimension(npt-1),intent(in) :: dx, beta0, beta1, beta2, pcdf
    real(8),intent(in)                  :: normct
    real(8),dimension(n),intent(inout)  :: p
    real(8),parameter                   :: spi2 = sqrt(3.14159265359d0) / 2.0d0
    real(8),parameter                   :: a = 0.147d0
    integer,dimension(n)                :: pos
    real(8),dimension(n)                :: tU, imU, erU
    complex(8)                          :: uL, iL
    complex(8),dimension(n)             :: cU, iU
    real(8),dimension(npt)              :: lagged_pcdf
    integer,dimension(1)                :: pvec
    integer                             :: i
    real(8)                             :: b0, b1, b2, rd, xi, dxi
    real(8)                             :: imL, tL, erL, sb2, signb2, aux

    do i = 1, n
        pvec   = minloc( abs(z(i) - x) )
        pos(i) = pvec(1)
    enddo
    where (x(pos) > z .and. pos > 1)
        pos = pos - 1
    end where
    lagged_pcdf(1)     = 0.0d0
    lagged_pcdf(2:npt) = pcdf( 1:(npt-1) )
    do i = 1, npt-1
      b0 = beta0(i)
      b1 = beta1(i)
      b2 = beta2(i)
      rd = lagged_pcdf(i)
      xi = x(i)
      dxi = dx(i)
      sb2 = sqrt(abs(b2))
      signb2 = sign(1.0d0, b2)
      aux = -spi2 * dxi * exp( ( 4.0d0 * b0 * b2 - b1**2 )/(4.0d0 * b2) ) / &
        (sb2 * normct)
      imL = -signb2 * 0.5d0 * b1 / sb2
      if (abs(b2) < 1d-15) then
        where(pos == i)
          p = 0.0d0
        end where
      elseif (abs(b1/b2) < 100.0d0) then
          if( b2 < 0.0d0) then
              !local version of erf(imL) -> erL
              tL  = 1.0d0 / (1.0d0 + 0.327591100d0 * abs(imL))
              erL = sign(1.0d0 - &
                        (0.254829592d0 * tL - 0.284496736d0 * tL**2 + &
                        1.421413741d0 * tL**3 - &
                        1.453152027d0 * tL**4 + &
                        1.061405429d0 * tL**5) * exp(-imL**2), imL)
              where(pos == i)
                  !local version of erf(imU) -> erU
                  imU = -signb2 * (b2 * (z - xi) / dxi + 0.5d0 * b1) / sb2
                  tU  = 1.0d0 / (1.0d0 + 0.327591100d0 * abs(imU))
                  erU = sign(1.0d0 - &
                          (0.254829592d0 * tU - 0.284496736d0 * tU**2 + 1.421413741d0 * tU**3 - &
                          1.453152027d0 * tU**4 + 1.061405429d0 * tU**5) * exp(-imU**2), imU)
                  p   = rd + aux * (erU - erL)
              end where
          else
              !local version of erfi(imL) -> erL
              uL  = cmplx(0.0d0,abs(imL),KIND=8)
              iL  = 1.0d0 / (1.0d0 + 0.327591100d0 * uL )
              erL = sign(aimag(1.0d0 - &
                        (0.254829592d0 * iL - 0.284496736d0 * iL**2 + 1.421413741d0 * iL**3 - &
                        1.453152027d0 * iL**4 + 1.061405429d0 * iL**5) * exp(-uL**2) ), imL)
              where(pos == i)
                  !local version of erfi(imU) -> erU
                  imU = -signb2 * (b2 * (z - xi) / dxi + 0.5d0 * b1 ) / sb2
                  cU  = cmplx(0.0d0,abs(imU), KIND = 8)
                  iU  = 1.0d0 / (1.0d0 + 0.327591100d0 * cU )
                  erU = sign(aimag(1.0d0 - &
                        (0.254829592d0 * iU - 0.284496736d0 * iU**2 + 1.421413741d0 * iU**3 - &
                        1.453152027d0 * iU**4 + 1.061405429d0 * iU**5) * exp(-cU**2) ), imU)
                  p   = rd + aux * (erU - erL)
              end where
          endif
      elseif (b2 < 0.0d0) then
          where(pos == i)
              imU = -sb2 * (0.5d0 * b1 + b2 * (z - xi) / dxi) / sb2
              p   = rd + aux * exp(-imL**2)*(imU-imL)*( 1.0d0 - (imU-imL)*imL ) / spi2
          end where
      elseif (b2 > 0.0d0) then
          where(pos == i)
              imU = -sb2 * (0.5d0 * b1 + b2 * (z - xi) / dxi) / sb2
              p   = rd + aux * exp( imL**2)*(imU-imL)*( 1.0d0 + (imU-imL)*imL ) / spi2
          end where
      endif
    enddo
    where (pos >= npt)
        p = 1.0d0
    end where
end subroutine continuous_cdf_quadratic

subroutine continuous_cdf_linear(n, z, npt, x, dx, beta0, beta1, pcdf, normct, p)
    implicit none
    integer, intent(in)                     :: n, npt
    real(8), intent(in)                     :: normct
    real(8), intent(in), dimension(n)       :: z
    real(8), intent(in), dimension(npt)     :: x, dx
    real(8), intent(in), dimension(npt-1)   :: beta0, beta1, pcdf
    real(8), intent(inout), dimension(n)    :: p
    integer                                 :: i, pos(n), pvec(1)
    real(8)                                 :: rd, b0, b1, di, xi

    do i = 1, n
        pvec   = minloc( abs(z(i) - x) )
        pos(i) = pvec(1)
    enddo
    where( x(pos) > z )
        pos = pos - 1
    end where

    rd = 0.0d0
    do i=1,npt-1
        b0 = beta0(i)
        b1 = beta1(i)
        di = dx(i)
        xi = x(i)
        where( pos == i )
            p = rd + (di / normct) * &
                (exp(b0) / b1 ) * ( -1.0d0 + exp(b1 * &
                (z-xi) / di))
        end where
        rd = pcdf(i) !for next round
    enddo
    where( pos >= npt)
        p = 1.0d0
    end where
end subroutine continuous_cdf_linear

subroutine continuous_invcdf_quadratic(n, p, npt, x, dx, beta0, beta1, beta2, pcdf, normct, z)
    implicit none
    integer,intent(in)                  :: n, npt
    real(8),dimension(n),intent(in)     :: p
    real(8),dimension(npt),intent(in)   :: x
    real(8),dimension(npt-1),intent(in) :: dx, beta0, beta1, beta2, pcdf
    real(8),intent(in)                  :: normct
    real(8),dimension(n),intent(inout)  :: z
    real(8),parameter                   :: twosqpi  = 2.0d0 / sqrt(3.141592653589793d0)
    real(8),parameter                   :: spi2     = sqrt(3.141592653589793d0) / 2.0d0
    real(8),parameter                   :: twopi    = 2.0d0 / 3.141592653589793d0
    real(8),parameter                   :: inverfa  = 0.145d0
    real(8),parameter                   :: inverfia = 0.1382d0
    integer,dimension(n)                :: pos
    real(8),dimension(n)                :: xU, tU, imU
    complex(8)                          :: uL, iL
    complex(8),dimension(n)             :: cU
    real(8),dimension(npt)              :: lagged_pcdf
    integer,dimension(1)                :: pvec
    integer                             :: i, j
    real(8)                             :: b0, b1, b2, rd, xi
    real(8)                             :: dxi, imL, tL, erL, sb2, aux, signb2
    real(8)                             :: a, b, c, r1, r2, z1, z2

    do i = 1, n
        pvec   = minloc( abs(p(i) - pcdf) )
        pos(i) = pvec(1)
    enddo
    where( pcdf(pos) < p .and. pos < npt)
        pos = pos + 1
    end where
    lagged_pcdf(1)     = 0.0d0
    lagged_pcdf(2:npt) = pcdf( 1:(npt-1) )
    do i = 1, npt-1
        b0  = beta0(i)
        b1  = beta1(i)
        b2  = beta2(i)
        rd  = lagged_pcdf(i)
        xi  = x(i)
        dxi = dx(i)
        sb2 = sqrt(abs(b2))
        signb2 = sign(1.0d0, b2)
        aux = - twosqpi * sb2 * normct / (dxi * exp( (4.0d0 * b0 * b2 - b1**2) / (4.0d0 * b2) ) )
        imL = - signb2 * 0.5d0 * b1 / sb2
        if( abs(b1/b2) < 100.0d0 ) then
            if( b2 < 0.0d0) then
                !local version of erf(imL) -> erL
                tL  = 1.0d0 / (1.0d0 + 0.327591100d0 * abs(imL))
                erL = sign(1.0d0 - &
                        (0.254829592d0 * tL - 0.284496736d0 * tL**2 + 1.421413741d0 * tL**3 - &
                        1.453152027d0 * tL**4 + 1.061405429d0 * tL**5) * exp(-imL**2), imL)
                where(pos == i)
                    !vectors
                    xU  = aux * (p - rd) + erL
                    !local version of inv_erf(xU) -> imU
                    tU  = log(1.0d0 - xU**2)
                    imU = sign( sqrt( sqrt( (twopi/inverfa + 0.5d0 * tU)**2 - tU / inverfa) - &
                                (twopi/inverfa + 0.5d0 * tU) ), xU)
                    z   = xi + dxi * (sb2 * imU - 0.5d0 * b1) / b2
                end where
            else
                !local version of erfi(imL) -> erL
                uL  = cmplx(0.0d0,abs(imL),KIND=8)
                iL  = 1.0d0 / (1.0d0 + 0.327591100d0 * uL )
                erL = sign(aimag(1.0d0 - (0.254829592d0 * iL - 0.284496736d0 * iL**2 + &
                            1.421413741d0 * iL**3 - &
                            1.453152027d0 * iL**4 + 1.061405429d0 * iL**5) * exp(-uL**2) ), imL)
                where(pos == i)
                    xU  = aux * (p - rd) + erL
                    !local version of inv_erfi(xU) -> imU
                    cU  = log( 1.0d0 - cmplx(0.0d0, abs(xU),KIND=8)**2 )
                    imU = sign(aimag(sqrt(sqrt( &
                            (twopi/inverfia + 0.5d0 * cU)**2 - cU / inverfia) - &
                            (twopi/inverfia + 0.5d0 * cU) )), xU)
                    z   = xi + dxi * (sb2 * imU - 0.5d0 * b1) / b2
                end where
            end if
        else  if( b2 < 0.0d0 ) then
            a = imL
            b = -1.0d0 - 2.0d0 * imL**2
            do j=1,n
                if(pos(j) == i) then
                    c = imL**3 + imL + (p(j)-rd) * aux * spi2 / exp(-imL**2)
                    r1 = (-b - sqrt(b**2-4.0d0*a*c))/(2.0d0*a)
                    z1 = -r1/sb2 - 0.5d0*b1/b2
                    r2 = (-b + sqrt(b**2-4.0d0*a*c))/(2.0d0*a)
                    z2 = -r2/sb2 - 0.5d0*b1/b2
                    if( abs(z1 - 0.5d0) <= 0.5d0 ) then
                        z(j) = xi + dxi * z1
                    else
                        z(j) = xi + dxi * z2
                    endif
                endif
            enddo
        else
            a = imL
            b = 1.0d0 - 2.0d0 * imL**2
            do j=1,n
                if(pos(j) == i) then
                    c = imL**3 - imL - (p(j)-rd) * aux * spi2 / exp( imL**2)
                    r1 = (-b - sqrt(b**2-4.0d0*a*c))/(2.0d0*a)
                    z1 = -r1/sb2 - 0.5d0*b1/b2
                    r2 = (-b + sqrt(b**2-4.0d0*a*c))/(2.0d0*a)
                    z2 = -r2/sb2 - 0.5d0*b1/b2
                    if( abs(z1 - 0.5d0) <= 0.5d0 ) then
                        z(j) = xi + dxi * z1
                    else
                        z(j) = xi + dxi * z2
                    endif
                endif
            enddo
        endif
    enddo
end subroutine continuous_invcdf_quadratic

subroutine continuous_invcdf_quadratic_bkp(n, p, npt, x, dx, beta0, beta1, beta2, pcdf, normct, z)
    implicit none
    integer,intent(in)                  :: n, npt
    real(8),dimension(n),intent(in)     :: p
    real(8),dimension(npt),intent(in)   :: x
    real(8),dimension(npt-1),intent(in) :: dx, beta0, beta1, beta2, pcdf
    real(8),intent(in)                  :: normct
    real(8),dimension(n),intent(out)    :: z
    real(8),parameter                   :: twosqpi  = 2.0d0 / sqrt(3.141592653589793d0)
    real(8),parameter                   :: twopi    = 2.0d0 / 3.141592653589793d0
    real(8),parameter                   :: inverfa  = 0.145d0
    real(8),parameter                   :: inverfia = 0.1382d0
    integer,dimension(n)                :: pos
    real(8),dimension(n)                :: xU, tU, imU
    complex(8)                          :: uL, iL
    complex(8),dimension(n)             :: cU
    real(8),dimension(npt-1)            :: lagged_pcdf
    integer,dimension(1)                :: pvec
    integer                             :: i
    real(8)                             :: b0, b1, b2, rd, xi
    real(8)                             :: dxi, imL, tL, erL, sb2, aux

    do i = 1, n
        pvec   = minloc( abs(p(i) - pcdf) )
        pos(i) = pvec(1)
    enddo
    where( pcdf(pos) < p .and. pos < npt)
        pos = pos + 1
    end where
    lagged_pcdf(1) = 0.0d0
    lagged_pcdf(2:npt) = pcdf( 1:(npt-1) )
    do i = 1, npt-1
        b0  = beta0(i)
        b1  = beta1(i)
        b2  = beta2(i)
        rd  = lagged_pcdf(i)
        xi  = x(i)
        dxi = dx(i)
        sb2 = sqrt(abs(b2))
        imL = 0.5d0 * b1 / sb2
        if( b2 < 0.0d0) then
            !scalars
            aux = -sb2 * twosqpi * exp( (b1**2 - 4.0d0 * b0 * b2) / (4.0d0 * b2) ) * normct / dxi
            !local version of erf(imL) -> erL
            tL  = 1.0d0 / (1.0d0 + 0.327591100d0 * abs(imL))
            erL = sign(1.0d0 - &
                  (0.254829592d0 * tL - 0.284496736d0 * tL**2 + 1.421413741d0 * tL**3 - &
                  1.453152027d0 * tL**4 + 1.061405429d0 * tL**5) * exp(-imL**2), imL)
            where(pos == i)
                !vectors
                xU  = aux * (p - rd) + erL
                !local version of inv_erf(xU) -> imU
                tU  = log(1.0d0 - xU**2)
                imU = sign( sqrt( sqrt( (twopi/inverfa + 0.5d0 * tU)**2 - tU / inverfa) - &
                      (twopi/inverfa + 0.5d0 * tU) ), xU)
                z   = xi + dxi * (2.0d0 * sb2 * imU - b1) / ( 2.0d0 * b2 )
            end where
        else
            !scalars
            aux = sb2 * twosqpi * exp( (b1**2 - 4.0d0 * b0 * b2) / (4.0d0 * b2) ) * normct / dxi
            !local version of erfi(imL) -> erL
            uL  = cmplx(0.0d0,abs(imL),KIND=8)
            iL  = 1.0d0 / (1.0d0 + 0.327591100d0 * uL )
            erL = sign(aimag(1.0d0 - (0.254829592d0 * iL - 0.284496736d0 * iL**2 + 1.421413741d0 * iL**3 - &
                  1.453152027d0 * iL**4 + 1.061405429d0 * iL**5) * exp(-uL**2) ), imL)
            where(pos == i)
                xU  = aux * (p - rd) + erL
                !local version of inv_erfi(xU) -> imU
                cU  = log( 1.0d0 - cmplx(0.0d0, abs(xU),KIND=8)**2 )
                imU = sign(aimag(sqrt( sqrt((twopi/inverfia + 0.5d0 * cU)**2 - cU / inverfia) - &
                      (twopi/inverfia + 0.5d0 * cU) )), xU)
                z   = xi + dxi * (2.0d0 * sb2 * imU - b1) / ( 2.0d0 * b2 )
            end where
        end if
    enddo
end subroutine continuous_invcdf_quadratic_bkp

subroutine continuous_invcdf_linear(n, p, npt, x, dx, beta0, beta1, pcdf, normct, z)
    implicit none
    integer, intent(in)                     :: n, npt
    real(8), intent(in)                     :: normct
    real(8), intent(in), dimension(n)       :: p
    real(8), intent(in), dimension(npt)     :: x
    real(8), intent(in), dimension(npt-1)   :: dx, beta0, beta1, pcdf
    real(8), intent(inout), dimension(n)    :: z
    integer                                 :: i, pos(n), pvec(1)
    real(8)                                 :: rd, b0, b1, di, xi

    do i = 1, n
        pvec   = minloc( abs(p(i) - pcdf) )
        pos(i) = pvec(1)
    enddo
    where( pcdf(pos) > p .and. pos > 1 )
        pos = pos - 1
    end where

    rd = 0.0d0
    do i=1,npt-1
        b0 = beta0(i)
        b1 = beta1(i)
        di = dx(i)
        xi = x(i)
        where( pos == i )
            z = xi + di * log( (p-rd) * normct * b1 / &
                (di * exp(b0)) + 1.0d0 ) / b1
        end where
        rd = pcdf(i) !for the next round
    enddo
end subroutine continuous_invcdf_linear

! Subroutine erfV: computes the error function of a real number
! Arguments:
!    n   -- input, dimension of vector x
!    x   -- input, vector containing the n values to be evaluated
!    res -- output, vector containing erf(x)
! This is an approximate function, provided by
! Winitzki, Sergei (6 February 2008). "A handy approximation for the error function and its inverse"
! Wikipedia reference: http://en.wikipedia.org/wiki/Error_function
! Example:
! r <- .Fortran('erfF',n=as.integer(1),x=as.double(1), res=as.double(0)); r$res
! should provide 0.8427007, i.e., Erf[1]
subroutine erfV(n,x,res)
    implicit none
    integer,intent(in)                  :: n
    real(8),dimension(n),intent(in)     :: x
    real(8),dimension(n),intent(out)    :: res
    real(8),dimension(n)                :: t

    t   = 1.0d0 / (1.0d0 + 0.327591100d0 * abs(x))
    res = sign(1.0d0 - (0.254829592d0 * t - 0.284496736d0 * t**2 + 1.421413741d0 * t**3 - &
                        1.453152027d0 * t**4 + 1.061405429d0 * t**5) * exp(-x**2), x)
end subroutine erfV

subroutine erfS(x,res)
    implicit none
    real(8),intent(in)     :: x
    real(8),intent(out)    :: res
    real(8)                :: t

    t   = 1.0d0 / (1.0d0 + 0.327591100d0 * abs(x))
    res = sign(1.0d0 - (0.254829592d0 * t - 0.284496736d0 * t**2 + 1.421413741d0 * t**3 - &
                        1.453152027d0 * t**4 + 1.061405429d0 * t**5) * exp(-x**2), x)
end subroutine erfS

! Subroutine erfiV: computes the imaginary error function of a complex number
! with real part equal to zero and imaginary part equal to x
! Arguments:
!    n   -- input, dimension of vector x
!    x   -- input, vector containing the n values to be evaluated
!    res -- output, vector containing Im( erfi(x) )
! This is an approximate function, provided by
! Winitzki, Sergei (6 February 2008). "A handy approximation for the error function and its inverse"
! Wikipedia reference: http://en.wikipedia.org/wiki/Error_function
! Example:
! ri <- .Fortran('erfiF',n=as.integer(1),x=as.double(1), res=as.double(0)); ri$res
! should provide 1.651042, i.e., Erfi[1]
subroutine erfiV(n, x, res)
    implicit none
    integer,intent(in)                  :: n
    real(8),dimension(n),intent(in)     :: x
    real(8),dimension(n),intent(out)    :: res
    complex(8),dimension(n)             :: t, u

    u   = cmplx(0.0d0,abs(x),KIND=8)
    t   = 1.0d0 / (1.0d0 + 0.327591100d0 * u )
    res = sign(aimag(1.0d0 - (0.254829592d0 * t - 0.284496736d0 * t**2 + 1.421413741d0 * t**3 - &
                   1.453152027d0 * t**4 + 1.061405429d0 * t**5) * exp(-u**2) ), x)
end subroutine erfiV

subroutine erfiS(x, res)
    implicit none
    real(8),intent(in)     :: x
    real(8),intent(out)    :: res
    complex(8)             :: t, u

    u   = cmplx(0.0d0,abs(x),KIND=8)
    t   = 1.0d0 / (1.0d0 + 0.327591100d0 * u )
    res = sign(aimag(1.0d0 - (0.254829592d0 * t - 0.284496736d0 * t**2 + 1.421413741d0 * t**3 - &
                   1.453152027d0 * t**4 + 1.061405429d0 * t**5) * exp(-u**2) ), x)
end subroutine erfiS

! Subroutine inv_erfV: computes the inverse of the error function of a real number
! Arguments:
!    n   -- input, dimension of vector x
!    x   -- input, vector containing the n values to be evaluated
!    res -- output, vector containing erf^{-1}(x)
! This is an approximate function, provided by
! Winitzki, Sergei (6 February 2008). "A handy approximation for the error function and its inverse"
! Wikipedia reference: http://en.wikipedia.org/wiki/Error_function
! Example:
! ir <- .Fortran('inv_erfF',n=as.integer(1),x=as.double(0.8427008), res=as.double(0)); ir$res
! should provide 1, i.e., invErf[0.8427]
subroutine inv_erfV(n, x, res)
    implicit none
    integer,intent(in)                  :: n
    real(8),dimension(n),intent(in)     :: x
    real(8),dimension(n),intent(out)    :: res
    real(8),dimension(n)                :: t, q
    real(8),parameter                   :: a = 0.145d0
    real(8),parameter                   :: twopia = 2.0d0/(3.141592653589793d0*0.145d0)

    t   = log(1.0d0 - x**2)
    q   = twopia + 0.5d0 * t
    res = sign( sqrt( sqrt(q**2 - t / a) - q ), x)
end subroutine inv_erfV

subroutine inv_erfS(x, res)
    implicit none
    real(8),intent(in)     :: x
    real(8),intent(out)    :: res
    real(8)                :: t, q
    real(8),parameter      :: a = 0.145d0
    real(8),parameter      :: twopia = 2.0d0/(3.141592653589793d0*0.145d0)

    t   = log(1.0d0 - x**2)
    q   = twopia + 0.5d0 * t
    res = sign( sqrt( sqrt(q**2 - t / a) - q ), x)
end subroutine inv_erfS

! Subroutine inv_erfiV: computes the inverse of the imaginary error function of
! a complex number with imaginary part equal to zero and real part equal to x
! Arguments:
!    n   -- input, dimension of vector x
!    x   -- input, vector containing the n values to be evaluated
!    res -- output, vector containing erfi^{-1}(x)
! This is an approximate function, provided by
! Winitzki, Sergei (6 February 2008). "A handy approximation for the error function and its inverse"
! Wikipedia reference: http://en.wikipedia.org/wiki/Error_function
! Example:
! iri <- .Fortran('inv_erfiF',n=as.integer(1),x=as.double(1.650426), res=as.double(0)); iri$res
! should provide 1, i.e., erfi^{-1}(1.650426)
subroutine inv_erfiV(n, x, res)
    implicit none
    integer,intent(in)                  :: n
    real(8),dimension(n),intent(in)     :: x
    real(8),dimension(n),intent(out)    :: res
    complex(8),dimension(n)             :: t, q, u
    real(8), parameter                  :: a = 0.1382d0
    real(8),parameter                   :: itwopia = 2.0d0/(3.141592653589793d0 * a)

    u   = cmplx(0.0d0, abs(x),KIND=8)
    t   = log(1.0d0 - u**2)
    q   = itwopia + 0.5d0 * t
    res = sign(aimag(sqrt( sqrt(q**2 - t / a) - q )), x)
end subroutine inv_erfiV

subroutine inv_erfiS(x, res)
    implicit none
    real(8),intent(in)     :: x
    real(8),intent(out)    :: res
    complex(8)             :: t, q, u
    real(8),parameter      :: a = 0.1382d0
    real(8),parameter      :: itwopia = 2.0d0/(3.141592653589793d0 * a)

    u   = cmplx(0.0d0, abs(x),KIND=8)
    t   = log(1.0d0 - u**2)
    q   = itwopia + 0.5d0 * t
    res = sign(aimag(sqrt( sqrt(q**2 - t / a) - q )), x)
end subroutine inv_erfiS

!Subroutine slow_erf: computes the error function of a real number, using a slow but accurate algorithm
! Arguments:
!    n   -- input, dimension of vector x
!    x   -- input, vector containing the n values to be evaluated
!    res -- output, vector containing erf(x)
! Reference: Book "Computation of Special functions"
! Vectorized algorithm adapted from the Fortran scalar subroutines provided in
! http://jin.ece.illinois.edu/routines/routines.html
! Example:
! rs <- .Fortran('slow_erf',n=as.integer(1),x=as.double(1), res=as.double(0)); rs$res
! should provide 0.8427008, i.e., Erf[1]
subroutine slow_erf(n, x, res)
    implicit none
    integer,intent(in)                  :: n
    real(8),dimension(n),intent(in)     :: x
    real(8),dimension(n),intent(out)    :: res
    real(8),parameter                   :: eps = 1.0d-15
    real(8),parameter                   :: dpi = 1.0d0 / sqrt(3.141592653589793d0)
    real(8),dimension(n)                :: a0, sq, esq, er
    logical,dimension(n)                :: test1, test2
    integer                             :: i, k
    real(8)                             :: r, rs, mysq

    a0    = abs(x)
    test1 = (a0 < 3.5d0)
    test2 = (x  < 0.0d0)
    sq    = x * x
    esq   = exp(-sq)
    do i=1,n
        mysq = sq(i)
        if(test1(i)) then
            r  = x(i)
            rs = r
            do k = 1, 120
                r = r * mysq / ( k + 0.5d0 )
                rs = rs + r
                if( abs(r/rs) <= eps ) exit
            enddo
        else
            r  = 1.0d0 / x(i)
            rs = r
            do k = 1, 13
                r = -r * (k - 0.5d0) / mysq
                rs = rs + r
                if( abs(r/rs) <= eps ) exit
            enddo
        end if
        er(i) = rs
    enddo
    where(test1)
        res = 2.0d0 * dpi * esq * er
    elsewhere
        res = 1.0d0 - dpi * esq * er
        where(test2)
            res = -res
        end where
    end where
end subroutine

!Subroutine slow_erfi: computes the imaginary error function of a complex number,
!using a slow but accurate algorithm
! Arguments:
!    n   -- input, dimension of vector x
!    x   -- input, vector containing the n values to be evaluated
!    res -- output, vector containing erfi(x)
! Reference: Book "Computation of Special functions"
! Vectorized algorithm adapted from the Fortran scalar subroutines provided in
! http://jin.ece.illinois.edu/routines/routines.html
! Example:
! ris <- .Fortran('slow_erfi',n=as.integer(1),z=complex(imaginary=1),
!                 res=as.complex(0)); Im(ris$res)
! should provide 1.650426, i.e., Erfi[1]
subroutine slow_erfi(n, z, res)
    implicit none
    integer,intent(in)                  :: n
    complex(8),dimension(n),intent(in)  :: z
    complex(8),dimension(n),intent(out) :: res
    real(8),parameter                   :: eps = 1.0d-15
    real(8),parameter                   :: dpi = 1.0d0 / sqrt(3.141592653589793d0)
    real(8),dimension(n)                :: a0, esq
    complex(8),dimension(n)             :: er
    logical,dimension(n)                :: test1, test2
    complex(8),dimension(n)             :: z1, sq
    integer                             :: i, k
    complex(8)                          :: r, rs, mysq

    a0    = cdabs(z)
    test1 = (a0 <= 5.8d0)
    test2 = (real(z) < 0.0)
    where(test2)
        z1 = -z
    elsewhere
        z1 = z
    end where
    sq    = z1 * z1
    esq   = real(cdexp(-z*z))
    do i=1,n
        mysq = sq(i)
        if(test1(i)) then
            r  = z1(i)
            rs = r
            do k = 1, 120
                r  = r * mysq / (k + 0.5d0)
                rs = rs + r
                if(cdabs(r/rs) < eps) exit
            enddo
        else
            r  = 1.0d0 / z1(i)
            rs = r
            do k = 1, 13
                r  = -r * (k - 0.5d0) / mysq
                rs = rs + r
                if(cdabs(r/rs) < eps) exit
            enddo
        endif
        er(i) = rs
    enddo
    where(test1)
        res = 2.0d0 * dpi * esq * er
    elsewhere
        res = 1.0d0 - dpi * esq * er
    end where
    where(test2)
        res = -res
    end where
end subroutine

!* Wrapping subroutine, used to compute the Exponential Integral **************
subroutine expint(x, res)
    implicit none
    real(8), intent(in)    :: x
    real(8), intent(out)   :: res
    real(8)                :: r8_ei

    res = r8_ei(x)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TESTING, remove when done
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convolution_normal_continuous(n, xl, dxl, b0l, b1l, b2l, normct, m, v, x, y)
    implicit none
    integer,intent(in)                  :: n
    real(8),dimension(n),  intent(in)   :: xl
    real(8),dimension(n-1),intent(in)   :: dxl
    real(8),dimension(n-1),intent(in)   :: b0l, b1l, b2l
    real(8),intent(in)                  :: normct, m, v
    real(8),dimension(n),intent(out)    :: x, y
    real(8),dimension(n-1)              :: a0part, a1part, alpha0, alpha1, alpha2
    real(8),dimension(n-1)              :: xlow, xhigh, integrals
    integer                             :: i
    real(8)                             :: evaluation_point, c

    c      = - 0.5d0 * log( 2d0 * 3.141592653589793d0 * v) - log(normct)
    xlow   = xl(1:(n-1))
    xhigh  = xl(2:n)
    a0part = b0l - b1l * xlow/dxl + b2l * (xlow/dxl)**2 + c
    a1part = b1l/dxl - 2.0d0 * b2l * xlow / dxl**2
    alpha2 = -0.5d0 / v + b2l / dxl**2
    do i = 1, n
      evaluation_point = xl(i) - 3.5d0 + 7.0d0*dble(i-1)/dble(n-1)
      alpha0 = a0part - 0.5d0*(evaluation_point - m)**2/v
      alpha1 = a1part + (evaluation_point - m)/v
      call continuous_integrals(n-1, alpha0, alpha1, alpha2, xlow, xhigh, integrals)
      x(i) = evaluation_point
      y(i) = log(sum(integrals * dxl))
    enddo
end subroutine convolution_normal_continuous
