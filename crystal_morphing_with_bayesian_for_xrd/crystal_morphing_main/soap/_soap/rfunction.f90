
!####################################################################################################
! small functions
!####################################################################################################
subroutine get_ref_gvecs(num_j, ref_mvecs,bb, ref_gvecs)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  integer, dimension(1:num_j,1:3), intent(in)::  ref_mvecs
  real(8), dimension(1:3,1:3), intent(in)::  bb
  real(8), dimension(1:num_j,1:3), intent(out):: ref_gvecs
  integer, dimension(3):: mvec

  do j=1, num_j
    mvec = ref_mvecs(j,:)
    do ix=1, 3
      ref_gvecs(j,ix) = dot_product(mvec,bb(:,ix))
    end do
  end do
  return
end subroutine get_ref_gvecs

subroutine index_in_rcut(num_j, ref_gvecs,rcut, iis)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  real(8), dimension(1:num_j,1:3), intent(in)::  ref_gvecs
  real(8), intent(in)::  rcut
  real(8), dimension(3):: gvec
  integer, dimension(1:num_j), intent(out):: iis
  do j=1, num_j
    gvec = ref_gvecs(j,:)
    gr = dot_product(gvec,gvec)**0.5
    if (gr < rcut) then
      iis(j) = j-1
    else
      iis(j) = -1  !!! in case for r>rcut : insert -1. Delete from the array by postprocessing on the python side.Priority was given to ease of implementation.
    endif
  end do
  return
end subroutine index_in_rcut

subroutine get_thetas_phis(num_i, xyzis, ris, thetas, phis)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  dimension xyzi(1:3)
  real(8), dimension(1:num_i,1:3), intent(in):: xyzis
  real(8), dimension(1:num_i), intent(in):: ris
  real(8), dimension(1:num_i), intent(out):: thetas, phis
  do i=1, num_i
    xyzi = xyzis(i,1:3)
    ri = ris(i)
    if (ri == 0.0) then
      thetas(i) = 0.0
      phis(i) = 0.0
    else
      thetas(i) = dacos(xyzi(3)/ri)
      phis(i)   = datan2(xyzi(2), xyzi(1))
    end if
  enddo
  return
end subroutine get_thetas_phis

subroutine get_ylms(nlmax, num_i, xyzis, ris, ylms)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  complex(kind(0d0)):: ylm
  real(8), dimension(1:num_i):: thetas, phis

  real(8), dimension(1:num_i,1:3), intent(in):: xyzis
  real(8), dimension(1:num_i), intent(in):: ris
  complex(kind(0d0)), dimension(nlmax+1,2*nlmax+1, 1:num_i):: ylms

  call get_thetas_phis(num_i, xyzis, ris, thetas, phis)
  do l=0, nlmax
    do m=0, 2*nlmax
      if (abs(m-nlmax) > l) then
        ylms(l+1, m+1, :) = (0d0, 0d0)
      end if
      do i=1, num_i
        theta = thetas(i)
        phi = phis(i)
        call sph_harm (l, m-nlmax, theta, phi, ylm)
        ylms(l+1, m+1, i) = ylm
      end do !i
    end do !m
  end do !l
  return
end subroutine get_ylms
!####################################################################################################
! chig
!####################################################################################################
subroutine chig (num_j, rjs, gvec, ajs, sigma, vu, res)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  complex(kind(0d0)) img

  intent(in)  sigma, vu
  real(8), dimension(1:3), intent(in):: gvec
  real(8), dimension(1:num_j,1:3), intent(in):: rjs
  real(8), dimension(1:num_j), intent(in):: ajs
  double complex, intent(out):: res

  img=(0d0, 1d0)
  sigma2=sigma*sigma
  vu_inv=1d0/vu

  dot_gg = dot_product(gvec,gvec)
  res   = (0.0,0.0)
  do j=1, num_j
    dot_gr = dot_product(gvec,rjs(j,:))
    res    = res + ajs(j)*zexp(-1d0*img*dot_gr)
  end do
  res = res*vu_inv*exp(-0.5d0*dot_gg*sigma2)
  return
end subroutine chig

subroutine get_ref_chigs(num_j, num_i, rjs, ref_gvecs, ajs, sigma, vu, ress)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)

  intent(in)  sigma, vu
  real(8), dimension(1:num_i,1:3), intent(in)::  ref_gvecs
  real(8), dimension(1:num_j,1:3), intent(in):: rjs
  real(8), dimension(1:num_j), intent(in):: ajs
  double complex, dimension(1:num_i), intent(out):: ress

  real(8), dimension(1:3) :: gvec
  double complex:: res

  do i=1, num_i
    gvec=ref_gvecs(i,:)
    call chig(num_j, rjs, gvec, ajs, sigma, vu, res)
    ress(i) = res
  enddo  
end subroutine get_ref_chigs

!####################################################################################################
! ds_dr
!####################################################################################################
subroutine ds_dr (rjs, gvecs, ajs, num_i, num_j, sigma, vu, tol_F, ds_dr_val)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  dimension dchig_real(1:3)
  complex(kind(0d0)) img, f_val, exp_igr
  complex(kind(0d0)), dimension(1:3):: dchig

  intent(in)  sigma, vu, tol_F
  real(8), dimension(1:num_i,1:3), intent(in):: gvecs
  real(8), dimension(1:num_j,1:3), intent(in):: rjs
  real(8), dimension(1:num_j), intent(in):: ajs
  real(8), dimension(1:num_j,1:num_i,1:3), intent(out):: ds_dr_val

  img=(0d0, 1d0)
  sigma2=sigma*sigma

  do i=1, num_i
    f_val=(0.0,0.0)
    do k=1, num_j
      dot_gr = dot_product(gvecs(i,:),rjs(k,:))
      f_val = f_val + ajs(k)*zexp(-1d0*img*dot_gr)
    end do
    dot_gg = dot_product(gvecs(i,:),gvecs(i,:))
    exp_gg = exp(-0.5d0*dot_gg*sigma2)
    abs_f = dsqrt(real(f_val)**2d0 + aimag(f_val)**2d0)
    if (abs_f < tol_F) then
      do j=1, num_j
        dot_gr = dot_product(gvecs(i,:),rjs(j,:))
        factor = ajs(j) * exp_gg / vu
        dchig_real = factor*(-1d0)*gvecs(i,:)*sin(dot_gr)
        ds_dr_val(j,i,:) = dchig_real 
      enddo
    else
      do j=1, num_j
        dot_gr = dot_product(gvecs(i,:),rjs(j,:))
        factor = ajs(j) * exp_gg / vu
        exp_igr = zexp(img*dot_gr)*f_val
        dchig = factor*img*gvecs(i,:)*(exp_igr - conjg(exp_igr))/(2.0*abs_f)
        ds_dr_val(j,i,:) = real(dchig)
      enddo
    end if
  end do
  return
end subroutine ds_dr


!####################################################################################################
! dpp_ds_res
!####################################################################################################
subroutine dpp_ds_res (num_i, num_m, nlmchi, ndlmchi, dnlmchi_ds, dndlmchi_ds, res)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  complex(kind(0d0)) res1, res2

  double complex, dimension(1:num_m), intent(in):: nlmchi, ndlmchi
  double complex, dimension(1:num_i,1:num_m), intent(in):: dnlmchi_ds, dndlmchi_ds
  double complex, dimension(1:num_i), intent(out):: res

  do i=1, num_i
    res1   = dot_product(conjg(ndlmchi),dnlmchi_ds(i,:))
    res2   = dot_product(conjg(dndlmchi_ds(i,:)),nlmchi)
    res(i) = res1 + res2
  end do
  return
end subroutine dpp_ds_res


!####################################################################################################
! nlmchi_arr
!####################################################################################################
subroutine nlmchi_arr (nbasis, num_i, xyzis, ris, sis, r2values, nlmax, watom, nlmchi_val)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  complex(kind(0d0)):: comp_coef, ylm
  complex(kind(0d0)), dimension(nlmax+1,2*nlmax+1, 1:num_i):: ylms

  intent(in)  nbasis, nlmax, watom
  real(8), dimension(1:num_i,1:3), intent(in):: xyzis
  real(8), dimension(1:num_i), intent(in):: ris
  real(8), dimension(1:num_i), intent(in):: sis
  real(8), dimension(nbasis, nlmax+1,1:num_i), intent(in):: r2values
  double complex, dimension(nbasis,nlmax+1,2*nlmax+1), intent(out):: nlmchi_val

  pi=dacos(-1d0)
  coef = (2.0/pi)**0.5 / (watom)**3.0
  
  call get_ylms(nlmax, num_i, xyzis, ris, ylms)
  
  nlmchi_val=(0d0, 0d0)
  do n=1, nbasis
    do l=0, nlmax
      do m=0, 2*nlmax
        if (abs(m-nlmax) > l) then
          cycle
        end if
        do i=1, num_i
          ylm = ylms(l+1, m+1, i)
          comp_coef = sis(i) * coef * conjg(ylm) 
          nlmchi_val(n,l+1,m+1) = nlmchi_val(n,l+1,m+1) + comp_coef * r2values(n,l+1,i)
        end do !i
      end do !m
    end do !l
  end do !n
  return
end subroutine nlmchi_arr

!####################################################################################################
! dnlmchi_dt_arr
!####################################################################################################
subroutine dnlmchi_dt_arr (nbasis, num_i, bb, xyzis, ris, sis, r2values, r3values, mfrac, gaunt_1d, &
                       nlmax, watom, res)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  integer xi, eta, nlld
  integer, dimension(2):: llds
  real(8), dimension(1:3):: xyzi, bxi
  real(8), allocatable, dimension(:,:,:,:):: gaunt
  complex(kind(0d0)):: ylm=(0d0, 0d0)
  complex(kind(0d0)), dimension(1:3):: jlm=(0d0, 0d0)
  complex(kind(0d0)), dimension(nlmax+1, 2*nlmax+1, num_i, 3):: jlms
  complex(kind(0d0)), dimension(1:3):: term1, term2, termvec
  complex(kind(0d0)), dimension(nbasis,nlmax+1,2*nlmax+1,1:3,1:3):: termvecs
  complex(kind(0d0)), dimension(nlmax+1,2*nlmax+1, 1:num_i):: ylms
  real(8), dimension(1:num_i):: thetas, phis

  integer, intent(in):: nbasis, nlmax
  real(8), intent(in):: watom
  real(8), dimension(:), intent(in):: gaunt_1d
  real(8), dimension(1:num_i,1:3), intent(in):: mfrac
  real(8), dimension(1:3,1:3), intent(in):: bb
  real(8), dimension(1:num_i,1:3), intent(in):: xyzis
  real(8), dimension(1:num_i), intent(in):: ris, sis
  real(8), dimension(nbasis, nlmax+1,1:num_i), intent(in):: r2values, r3values
  double complex, dimension(nbasis,nlmax+1,1:3,1:3,2*nlmax+1), intent(out):: res
  
  allocate(gaunt(nlmax+1,2*nlmax+1,2,3))
  gaunt = reshape(gaunt_1d, shape(gaunt))

  pi=dacos(-1d0)
  coef = 1d0 / (watom)**5d0
  termvecs=(0d0, 0d0)
  res=(0d0, 0d0)

  !! for speedup
  call get_thetas_phis(num_i, xyzis, ris, thetas, phis)
  call get_ylms(nlmax, num_i, xyzis, ris, ylms)
  do l=0, nlmax
    if (l == 0) then
      nlld = 1
      llds = (/ 1, 0 /)
    else
      nlld = 2
      llds = (/ l-1, l+1 /)
    end if
    do m=0, 2*nlmax
      if (abs(m-nlmax) > l) then
        cycle
      end if
      do i=1, num_i
        theta = thetas(i)
        phi = phis(i)
        call jgaunt_2 (gaunt(l+1,m+1,:,:), nlld, llds, m-nlmax, theta, phi, jlm)
        jlms(l+1, m+1, i, :) = jlm(:)
      enddo
    enddo
  enddo

  do n=1, nbasis
    do l=0, nlmax
      if (l == 0) then
        nlld = 1
        llds = (/ 1, 0 /)
      else
        nlld = 2
        llds = (/ l-1, l+1 /)
      end if
      do m=0, 2*nlmax
        if (abs(m-nlmax) > l) then
          cycle
        end if
        do i=1, num_i
          xyzi = xyzis(i,1:3)
          ri = ris(i)
          theta = thetas(i)
          phi = phis(i)
          ylm = ylms(l+1, m+1, i)

          !!! speedup          
          jlm(:) = jlms(l+1, m+1, i, :)
          !!! original
          ! call jgaunt_2 (gaunt(l+1,m+1,:,:), nlld, llds, m-nlmax, theta, phi, jlm)

          term1 = -1d0*(2d0/pi)**5d-1 * conjg(ylm) * r2values(n,l+1,i) * xyzi
          term2 = (-1d0)**(m-nlmax)*2d0/(3d0)**5d-1 * r3values(n,l+1,i) * jlm

          do j=1, 3
            termvecs(n,l+1,m+1,j,:) = termvecs(n,l+1,m+1,j,:) + sis(i)*mfrac(i,j)*(term1+term2)
          end do !j
        end do !i
        do xi=1, 3
          do eta=1, 3
            bxi = bb(xi,:)
            termvec = termvecs(n,l+1,m+1,eta,:)
            res(n,l+1,xi,eta,m+1) = coef * dot_product(bxi, termvec)
            if (n==1 .and. l==0 .and. xi==1 .and. eta==1) then
            end if
          end do !eta
        end do !xi
      end do !m
    end do !l
  end do !n
  deallocate(gaunt)
  return
end subroutine dnlmchi_dt_arr

!####################################################################################################
! dnlmchi_ds_arr
!####################################################################################################
subroutine dnlmchi_ds_arr (nbasis, num_i, xyzis, ris, r2values, nlmax, watom, res)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  dimension xyzi(1:3)
  complex(kind(0d0)):: comp_coef, ylm
  complex(kind(0d0)), dimension(nlmax+1,2*nlmax+1, 1:num_i):: ylms
  real(8), dimension(1:num_i):: thetas, phis

  integer, intent(in):: nbasis, nlmax
  real(8), intent(in):: watom
  real(8), dimension(1:num_i,1:3), intent(in):: xyzis
  real(8), dimension(1:num_i), intent(in):: ris
  real(8), dimension(nbasis, nlmax+1,1:num_i), intent(in):: r2values
  double complex, dimension(nbasis,nlmax+1,1:num_i,2*nlmax+1), intent(out):: res
  ! complex(kind(0d0)):: nlmchi_val
  ! real(8), intent(out):: nlmchi_val_real

  pi=dacos(-1d0)
  coef=(2d0/pi)**5d-1 / (watom)**3d0
  res=(0d0, 0d0)

  call get_thetas_phis(num_i, xyzis, ris, thetas, phis)
  call get_ylms(nlmax, num_i, xyzis, ris, ylms)

  do n=1, nbasis
    do l=0, nlmax
      do m=0, 2*nlmax
        if (abs(m-nlmax) > l) then
          cycle
        end if
        do i=1, num_i
          xyzi = xyzis(i,1:3)
          ri = ris(i)
          theta = thetas(i)
          phi = phis(i)
          ylm = ylms(l+1, m+1, i)
          comp_coef = coef * conjg(ylm)
          res(n,l+1,i,m+1) = comp_coef * r2values(n,l+1,i)
        end do !i
      end do !m
    end do !l
  end do !n
  return
end subroutine dnlmchi_ds_arr

!####################################################################################################
! Spherical harmonics
!####################################################################################################
subroutine sph_harm (l,m_pm,theta,phi,ylm)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  complex(kind(0d0)):: factor,img

  intent(in) l, m_pm, theta, phi
  complex(kind(0d0)), intent(out):: ylm
  !Computes the associated Legendre polynomialPml(x).Heremandlare integers satisfying0≤m≤l, whilexlies in the range−1≤x≤1.
  !call assert(m >= 0, m <= l, abs(x) <= 1.0, "plgndr_s args")
  m=abs(m_pm)
  pi4=4d0*dacos(-1d0)
  img=(0d0, 1d0)
  x=dcos(theta)
  ! if(x < 1e-10) then
  !   x = 0d0
  ! end if
  pmm=1d0
  pll=1d0
  !ComputePmm.
  if (m > 0) then
     somx2=sqrt((1d0-x)*(1d0+x))
     pmm=fact(2*m)/2d0**m/fact(m)*somx2**m
     if (mod(m,2) == 1) pmm=-pmm
  end if
  if(l == m) then
     plgndr=pmm
  else
     pmmp1=x*dble(2*m+1)*pmm
     !ComputePmm+1.
     if (l == m+1) then
        plgndr=pmmp1
     else
        !ComputePml,l>m+1.
        do ll=m+2,l
           pll=(x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm)/dble(ll-m)
           pmm=pmmp1
           pmmp1=pll
        end do
        plgndr=pll
     end if
  end if
  factor = zexp(img*dble(m)*phi)*dsqrt(dble(2.0*l+ 1)*dble(fact(l-m))/pi4/dble(fact(l+m)))
  ylm = factor*plgndr

  if (m_pm < 0) then
    ylm = (-1d0)**m * conjg(ylm)
  end if

contains
  ! factorial
  ! :param  integer n
  ! :return integer fact
  integer(kind=8) recursive function fact(n)
    implicit none
    integer, intent(IN) :: n
    integer :: i

    fact = 1
    do i = 1, n
       fact = fact * i
    enddo
    return
  end function fact
end subroutine sph_harm

!####################################################################################################
! Jlm for fixed array size
!####################################################################################################
subroutine jgaunt_2 (gaunt,nlld,llds,mm,theta,phi,jlm)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  complex(kind(0d0)):: img, term1, term2, term3, ylm

  intent(in) nlld, mm, theta, phi
  integer, dimension(1:2), intent(in):: llds
  real(8), dimension(1:2,1:3), intent(in):: gaunt !! <= different from jgaunt
  complex(kind(0d0)), dimension(1:3), intent(out):: jlm

  img=(0d0, 1d0)
  jlm=(0d0, 0d0)

  do i=1, nlld
    call sph_harm(llds(i), mm+1, theta, phi, ylm)
    term1 = conjg(ylm) * gaunt(i,1)
    call sph_harm(llds(i), mm-1, theta, phi, ylm)
    term2 = conjg(ylm) * gaunt(i,2)
    call sph_harm(llds(i), mm, theta, phi, ylm)
    term3 = conjg(ylm) * gaunt(i,3)

    jlm(1) = jlm(1) + (term1 - term2)
    jlm(2) = jlm(2) + img*(term1 + term2)
    jlm(3) = jlm(3) + term3*2d0**5d-1
  end do

  ! write(*,*) llds(1), mm, theta, phi, ylm, gaunt(1,3)

  return
end subroutine jgaunt_2

!####################################################################################################
! Jlm
!####################################################################################################
subroutine jgaunt (gaunt,nlld,llds,mm,theta,phi,jlm)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  complex(kind(0d0)):: img, term1, term2, term3, ylm

  intent(in) nlld, mm, theta, phi
  integer, dimension(1:nlld), intent(in):: llds
  real(8), dimension(1:nlld,1:3), intent(in):: gaunt
  complex(kind(0d0)), dimension(1:3), intent(out):: jlm

  img=(0d0, 1d0)
  jlm=(0d0, 0d0)

  do i=1, nlld
    call sph_harm(llds(i), mm+1, theta, phi, ylm)
    term1 = conjg(ylm) * gaunt(i,1)
    call sph_harm(llds(i), mm-1, theta, phi, ylm)
    term2 = conjg(ylm) * gaunt(i,2)
    call sph_harm(llds(i), mm, theta, phi, ylm)
    term3 = conjg(ylm) * gaunt(i,3)

    jlm(1) = jlm(1) + (term1 - term2)
    jlm(2) = jlm(2) + img*(term1 + term2)
    jlm(3) = jlm(3) + term3*2d0**5d-1
  end do

  ! write(*,*) llds(1), mm, theta, phi, ylm, gaunt(1,3)

  return
end subroutine jgaunt

!####################################################################################################
! Radial basis set
!####################################################################################################
subroutine radial_basis_set (nbasis,nrad,rcut, watom, drad, rbasis, rgrid)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  dimension rgrid(1:nrad)
  dimension rbasis(1:nrad,1:nbasis)
  intent(in)  nbasis, nrad, rcut, watom
  intent(out) drad, rbasis, rgrid
  watom_tmp=5d-1/watom**2
  call construct_radial_basis_set (drad,nbasis,nrad,rbasis,rcut,rgrid,watom_tmp)
end subroutine radial_basis_set

!####################################################################################################
! Radial integration
!####################################################################################################
subroutine r2_arr(num_i, nbasis, ris, nlmax, rcut, watom, nrad, drad, rbasis, rgrid, r2values)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  intent(in)  nbasis, nlmax, rcut, watom, nrad, drad
  real(8), dimension(1:num_i), intent(in):: ris
  real(8), dimension(1:nrad,1:nbasis), intent(in):: rbasis
  real(8), dimension(1:nrad), intent(in):: rgrid
  real(8), dimension(1:nbasis, nlmax+1, 1:num_i), intent(out):: r2values
  do n=1, nbasis
    do l=0, nlmax
      do i=1, num_i
        ri = ris(i)
        call r2 (n,l, nbasis, nlmax, rcut, watom, nrad, ri, drad, rbasis, rgrid, r2value)
        r2values(n, l+1, i) = r2value
      end do !i
    enddo
  enddo
end subroutine r2_arr

subroutine r2 (ibasis,ll, nbasis, nlmax, rcut, watom, nrad, rjir, drad, rbasis, rgrid, clm)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  !parameter(nrad=30)
  dimension expr(1:nrad)
  dimension ril(1:nrad,0:(nlmax+1))
  dimension rr(1:nrad)
  dimension cr(1:nrad)
  ! double complex clm
  ! double complex, dimension(1:nrad) :: c

  intent(in)  nbasis, nlmax, rcut, watom, nrad, rjir, ibasis, ll, drad
  real(8), dimension(1:nrad,1:nbasis), intent(in):: rbasis
  real(8), dimension(1:nrad), intent(in):: rgrid
  intent(out) clm

  ! write(*,*) ibasis,ll, nbasis, nlmax
  ! write(*,*) rcut, watom, nrad, drad
  ! write(*,*) rbasis, rgrid
  watom_tmp=5d-1/watom**2

  ! call construct_radial_basis_set (drad,nbasis,nrad,rbasis,rcut,rgrid,watom_tmp)
  !! => out: drad, rbasis, rgrid

  call set_gaussian_r (expr,eps,nrad,pi,pi4,prefactor,rgrid,watom_tmp)
  !! => out: expr, eps, pi, pi4, prefactor

  call gaussian (expjir,rjir,watom_tmp)
  !! => out: expji

  call cutoff_bp (fji,pi,rcut,rjir)
  !! => fji

  do krad=1, nrad
    rr(krad)=2d0*watom_tmp*rjir*rgrid(krad)
  enddo


  call modified_spherical_Bessel_1st (nlmax,nrad,nrad,rr(1),ril(1,0))
  ! call modified_spherical_Bessel_1st (nlmax,nrad,nrad,rr,ril)
  ! check
  ! do krad=1, nrad
  !   write(*,*) rr(krad), ril(krad,0)
  ! enddo
  ! do krad=1, nrad
  !   write(*,*) rr(krad), ril(krad,1)
  ! enddo
  ! stop
  !! => ril

  kamn = ll
  do krad=1, nrad
    if (rjir == 0.0) then
      if (ll==0) then
        ril(krad,kamn) = 1.0
      else
        ril(krad,kamn) = 0.0
      endif
    endif
    cr(krad)=prefactor*fji*expjir*expr(krad)*ril(krad,kamn)
    ! write(*,*) krad*drad, cr(krad)
  enddo
  clm=0d0
  do irad=1, nrad
    clm=clm+rgrid(irad)**2*rbasis(irad,ibasis)*cr(irad)*drad
  enddo

  return
end subroutine r2

subroutine r3_arr(num_i, nbasis, ris, nlmax, rcut, watom, nrad, drad, rbasis, rgrid, r3values)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  intent(in)  nbasis, nlmax, rcut, watom, nrad, drad
  real(8), dimension(1:num_i), intent(in):: ris
  real(8), dimension(1:nrad,1:nbasis), intent(in):: rbasis
  real(8), dimension(1:nrad), intent(in):: rgrid
  real(8), dimension(1:nbasis, nlmax+1, 1:num_i), intent(out):: r3values
  do n=1, nbasis
    do l=0, nlmax
      do i=1, num_i
        ri = ris(i)
        call r3 (n,l, nbasis, nlmax, rcut, watom, nrad, ri, drad, rbasis, rgrid, r3value)
        r3values(n, l+1, i) = r3value
      end do !i
    enddo
  enddo
end subroutine r3_arr

subroutine r3 (ibasis,ll, nbasis,nlmax, rcut, watom, nrad, rjir, drad, rbasis, rgrid, clm)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  !parameter(nrad=30)
  dimension expr(1:nrad)
  dimension rgrid(1:nrad)
  dimension rbasis(1:nrad,1:nbasis)
  dimension ril(1:nrad,0:(nlmax+1))
  dimension rr(1:nrad)
  dimension cr(1:nrad)
  ! double complex clm
  ! double complex, dimension(1:nrad) :: c

  intent(in)  nbasis, nlmax, rcut, watom, nrad, rjir, ibasis, ll
  intent(out) clm

  watom_tmp=5d-1/watom**2

  ! call construct_radial_basis_set (drad,nbasis,nrad,rbasis,rcut,rgrid,watom_tmp)
  !! => out: drad, rbasis, rgrid

  call set_gaussian_r (expr,eps,nrad,pi,pi4,prefactor,rgrid,watom_tmp)
  !! => out: expr, eps, pi, pi4, prefactor

  call gaussian (expjir,rjir,watom_tmp)
  !! => out: expji

  call cutoff_bp (fji,pi,rcut,rjir)
  !! => fji

  do krad=1, nrad
    rr(krad)=2d0*watom_tmp*rjir*rgrid(krad)
  enddo

  ! call modified_spherical_Bessel_1st (nlmax,nrad,nrad,rr(1),ril(1,0))
  call modified_spherical_Bessel_1st (nlmax,nrad,nrad,rr(1),ril(1,0))
  !! => ril

  kamn = ll
  do krad=1, nrad
    if (rjir == 0.0) then
      if (ll==0) then
        ril(krad,kamn) = 1.0
      else
        ril(krad,kamn) = 0.0
      endif
    endif
    cr(krad)=prefactor*fji*expjir*expr(krad)*ril(krad,kamn)
  enddo
  clm=0d0
  do irad=1, nrad
    clm=clm+rgrid(irad)**3*rbasis(irad,ibasis)*cr(irad)*drad
  enddo

  return
end subroutine r3

!####################################################################################################
! Construct radial basis set
!####################################################################################################
subroutine construct_radial_basis_set (drad,nbasis,nrad,rbasis,rcut,rgrid,watom_tmp)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  dimension ik(1:nbasis)
  dimension rbasis(1:nrad,1:nbasis)
  dimension rcrad(1:nbasis)
  dimension rgaussian(1:nrad,1:nbasis)
  dimension rgrid(1:nrad)
  dimension s(1:nbasis,1:nbasis),u(1:nbasis,1:nbasis)
  dimension w(1:nbasis),wk(1:3*nbasis-1)

  intent(in) nbasis, nrad, rcut, watom_tmp
  intent(out) drad, rbasis, rgrid

  eps=1d6*epsilon(1d0)
  drad=rcut/dble(nrad)
  ! wbasis_tmp=5d-1/wbasis**2
  wbasis_tmp = watom_tmp
  do i=1, nbasis
    ! rcrad(i)=dble(i)/dble(nbasis)*rcut
    rcrad(i)=dble(i)/dble(nbasis)*rcut
  enddo
  do i=1, nrad
    rgrid(i)=dble(i)*drad
  enddo
  do i=1, nbasis
    do j=1, nrad
      rgaussian(j,i)=dexp(-wbasis_tmp*(rgrid(j)-rcrad(i))**2)
    enddo
  enddo
  !!! Kajita, trapezoidl rule
  do i=1, nbasis
    do j=1, nbasis
      !!! Kajita, trapezoidl rule
      s(j,i)=0d0
      do k=1, nrad
        if (k==1 .or. k==nrad) then
          coef = 0.5
        else
          coef = 1.0
        endif
        s(j,i)=s(j,i)+coef*rgaussian(k,j)*rgaussian(k,i)*rgrid(k)**2*drad
      enddo
      ! write(*,*) "trapezoidal, 1", s
      ! !!
      ! s(j,i)=0d0
      ! do k=1, nrad
      !   s(j,i)=s(j,i)+rgaussian(k,j)*rgaussian(k,i)*rgrid(k)**2*drad
      ! enddo
      ! write(*,*) "original, 1", s
    enddo
  enddo
  do i=1, nbasis
    do j=1, nbasis
      u(j,i)=s(j,i)
    enddo
  enddo
  call dsyev ('V','U',nbasis,u,nbasis,w,wk,3*nbasis-1,info)
  lflag=0
  do i=1, nbasis
    if(w(i) .lt. eps) then
      lflag=1
    endif
  enddo
  if(lflag .eq. 1) then
    call abend("sub.construct_radial_basis_set","1. Overcompleteness in radial basis set. Set smaller nbasis or wbasis")
  endif
  do i=1, nbasis
    do j=1, nbasis
      u(j,i)=s(j,i)
    enddo
  enddo
  call dpotrf ('U',nbasis,u,nbasis,info)
  do i=1, nbasis
    do j=i+1, nbasis, 1
      u(j,i)=0d0
    enddo
  enddo
  call dgetrf (nbasis,nbasis,u,nbasis,ik,info)
  call dgetri (nbasis,u,nbasis,ik,wk,nbasis,info)
  do i=1, nbasis
    do k=1, nrad
      rbasis(k,i)=0d0
    enddo
    do j=1, nbasis
      do k=1, nrad
        rbasis(k,i)=rbasis(k,i)+u(j,i)*rgaussian(k,j)
      enddo
    enddo
  enddo
  lflag=0
  do i=1, nbasis
    do j=1, nbasis
      !!! Kajita, trapezoidl rule
      s(j,i)=0d0
      do k=1, nrad
        if (k==1 .or. k==nrad) then
          coef = 0.5
        else
          coef = 1.0
        endif
        s(j,i)=s(j,i)+coef*rbasis(k,j)*rbasis(k,i)*rgrid(k)**2*drad
      enddo
      ! write(*,*) "trapezoidal, 2", s
      ! s(j,i)=0d0
      ! do k=1, nrad
      !   s(j,i)=s(j,i)+rbasis(k,j)*rbasis(k,i)*rgrid(k)**2*drad
      ! enddo
      ! write(*,*) "original, 2", s
      if(j .eq. i) then
        if(dabs(s(j,i)-1d0) .gt. eps) then
          lflag=1
        endif
      else
        if(dabs(s(j,i)) .gt. eps) then
          lflag=1
        endif
      endif
      if(lflag .eq. 1) then
        call abend("sub.construct_radial_basis_set","2. Overcompleteness in radial basis set. Set smaller nbasis or wbasis")
      endif
    enddo
  enddo

  ! write(*,*) u

end subroutine construct_radial_basis_set

!****************************************************************************************************
! Constants for radial basis set expansion
!****************************************************************************************************

subroutine set_gaussian_r (expr,eps,nrad,pi,pi4,prefactor,rgrid,watom_tmp)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  dimension expr(1:nrad)
  dimension rgrid(1:nrad)

  intent(in) nrad, rgrid, watom_tmp
  intent(out) expr, eps, pi, pi4, prefactor

  eps=1d-6
  pi=dacos(-1d0)
  pi4=4d0*dacos(-1d0)
  ! watom_tmp=5d-1/watom**2
  ! watom_tmp=1d0/watom**2
  ! prefactor=pi4*(watom_tmp/pi)**1.5
  prefactor = 1.0
  do irad=1, nrad
    expr(irad)=dexp(-watom_tmp*rgrid(irad)**2)
  enddo

end subroutine set_gaussian_r

!####################################################################################################
! Gaussian
!####################################################################################################

subroutine gaussian (expji,rji,w)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)

  intent(in) rji, w
  intent(out) expji

  expji=dexp(-w*rji**2)

  return
end subroutine gaussian

!####################################################################################################
! Behler-Parrinello cutoff function
!####################################################################################################

subroutine cutoff_bp (fji,pi,rcut,rji)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)

  intent(in) pi, rcut, rji
  intent(out) fji

  if(rji .le. rcut) then
    fji=5d-1*(dcos(pi*rji/rcut)+1d0)
    dfjidrji=5d-1*(-pi/rcut*dsin(pi*rji/rcut))
  else
    fji=0d0
    dfjidrji=0d0
  endif

  ! if(rji .le. rcut) then
  !   fji= 1d0
  !   dfjidrji=0.0
  ! else
  !   fji=1d0
  !   dfjidrji=0d0
  ! endif
  return
end subroutine  cutoff_bp

!####################################################################################################
! Modified spherical Bessel
!####################################################################################################

subroutine modified_spherical_Bessel_1st (lmax,nrad,mrad,r,ril)
  ! implicit real(kind(1d0)) (a-h,o-z)
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  dimension r(1:mrad)
  dimension ril(1:mrad,0:(lmax+1))

  intent(in) lmax, nrad, mrad, r
  intent(out) ril

  if(lmax .eq. 0) then
    do irad=1, nrad
      ril(irad,0)=dsinh(r(irad))/r(irad)
      ril(irad,1)=(r(irad)*dcosh(r(irad))-dsinh(r(irad)))/r(irad)**2
    enddo
  else if(lmax .eq. 1) then
    do irad=1, nrad
      ril(irad,0)=dsinh(r(irad))/r(irad)
      ril(irad,1)=(r(irad)*dcosh(r(irad))-dsinh(r(irad)))/r(irad)**2
      ril(irad,2)=ril(irad,0)-3d0*ril(irad,1)/r(irad)
    enddo
  else
    do irad=1, nrad
      ril(irad,0)=dsinh(r(irad))/r(irad)
      ril(irad,1)=(r(irad)*dcosh(r(irad))-dsinh(r(irad)))/r(irad)**2
      ril(irad,2)=ril(irad,0)-3d0*ril(irad,1)/r(irad)
    enddo
    do l=3, lmax+1
      do irad=1, nrad
        ril(irad,l)=ril(irad,l-2)-(2*(l-1)+1)/r(irad)*ril(irad,l-1)
      enddo
    enddo
  endif
  return
end subroutine modified_spherical_Bessel_1st

!####################################################################################################
! Output error message
!####################################################################################################

subroutine abend(subname,comment)
  character(len=*) subname,comment

  intent(in) subname, comment

  write(*,'(/5x,"!!! error in ",a," !!!"/5x,a)') subname(1:len_trim(subname)),comment(1:len_trim(comment))
  stop
end subroutine abend

