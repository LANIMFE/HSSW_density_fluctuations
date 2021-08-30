Module hsay_structure
  Use nrtype
  Use hs_structure
  Use quadratures
  Implicit None
  Contains

    !Function that calculates the direct correlation function of a Hard Sphere +
    !Attractive Yukawa System using Sharma-Sharma approximation Scheme and Verlet-Weiss
    !for the hard core part
    ! c(k,ϕ,T,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function c_hsay_shvw(phi,T,z,k) result(c)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8 :: c
      Real * 8 :: k_func,ckdum,chs
      Real * 8 :: k2,z2,k4
      chs = c_hs_vw(phi,k)
      !chs = c_hs_py(phi,k)
      k2 = k * k
      z2 = z ** 2
      if (k > 0.075d0) then
        c = ( k * cos(k) ) + ( z * sin(k) )
        c = c / k
      else
        k4 = k2 * k2
        c = (1.d0 + z) - ((3.d0 + z) / 6.d0) * k2 + ((5.d0 + z) / 120.d0) * k4
      end if
      c = 4.d0 * pi * c / (T * (k2 + z2))
      c = chs + c
      return
    End Function

    !Function that calculates the inverse structure factor of a Hard Sphere +
    !Attractive Yukawa System using Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! 1/S(k,ϕ,T,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function is_hsay_shvw(phi,T,z,k) result(is)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent (in) :: k !Value of the wave vector
      Real * 8 :: is
      Real * 8 :: c
      Real * 8 :: chs
      Real * 8 :: k2,z2,k4, phi_vw
      phi_vw = phi*(1.d0 - (phi / 16.d0))
      chs = c_hs_vw(phi,k)
      !chs = c_hs_py(phi,k)
      k2 = k * k
      z2 = z ** 2
      if (k > 0.075d0) then
        c = ( k * cos(k) ) + ( z * sin(k) )
        c = c / k
      else
        k4 = k2 * k2
        c = (1.d0 + z) - ((3.d0 + z) / 6.d0) * k2 + ((5.d0 + z) / 120.d0) * k4
      end if
      c = 4.d0 * pi * c / (T * (k2 + z2))
      is =  (phi_vw * chs) + (phi * c)
      is = 1.d0 - 6.d0 * is / pi
      return
    End Function

    !Function that calculates the structure factor of a Hard Sphere +
    !Attractive Yukawa System using Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function s_hsay_shvw(phi,T,z,k) result(s)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent (in) :: k !Value of the wave vector
      Real * 8 :: s
      s = 1.d0 / is_hsay_shvw(phi,T,z,k)
      return
    End Function



    !!!!!#Note: do not use any of the gsvw functions, there are bad implementations, instead use gsvw2
    !Function that calculates the structure factor of a Hard Sphere +
    !Attractive Yukawa System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function c_hsay_gsvw(phi,T,z,k,qn,qw) result(c)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: c
      Real * 8 :: cay,u,chs,r,w, dumsin,kr
      Real * 8 :: k2,z2,k4
      Integer :: i1
      cay = 0.d0
      do i1=1, size(qn)
        r = qn(i1)
        kr = k * r
        if (kr > 0.075d0 ) then
          w = qw(i1) * r * sin(k * r) / k
        else
          w = qw(i1) * r**2 * (1.d0 + ( (k**2) * r/6.d0) )
        end if
        u = - exp(z * (r-1.d0) )/r
        cay = cay + ( ( exp(-u/T) - 1.d0 ) * w )
      end do
      chs = c_hs_vw(phi,k)
      c = chs + cay
      return
    End Function

    !Function that calculates the inverse structure factor of a Hard Sphere +
    !Attractive Yukawa System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! 1/S(k,ϕ,T,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function is_hsay_gsvw(phi,T,z,k,qn,qw) result(is)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: is
      is = c_hsay_gsvw(phi,T,z,k,qn,qw)
      is = 1.d0 - (6.d0 * phi * is / pi)
      return
    End Function

    !Function that calculates the structure factor of a Hard Sphere +
    !Attractive Yukawa System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function s_hsay_gsvw(phi,T,z,k,qn,qw) result(s)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: s
      s = 1.d0 / is_hsay_gsvw(phi,T,z,k,qn,qw)
      return
    End Function

    Function Ts_hsay_shvw(phi,z) result(Ts)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8 :: Ts !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8 :: S_HS
      S_HS = s_hs_vw(phi,0.d0)
      Ts = 24.d0 * phi * (z + 1.d0)  * S_HS / (z*z)
      return
    End Function

    Function u_yuk(r,z) result(u)
      Implicit None
      Real * 8, Intent(in) :: r !interparticle distance
      Real * 8, Intent(in) :: z !Screening value defined as positive
      Real * 8 :: u !Energy of the interacion scaled with the yukawa energy constant
      Real * 8 :: exponent !dummy variable for the exponent
      exponent=-z*(r-1.d0)
      if (exponent > -36.d0) Then
        u = - exp(exponent) / r
      else
        u = 0.d0
      endif
      return
    End Function

    Subroutine c_pk_hsay_gsvw_aux(z,T,np,rp,rw)
      Use quadratures
      Implicit None
      Real * 8, Intent(in) :: z, T !Yukawa parameter and system Temperature
      Integer, Intent(in) :: np !number of points
      Real * 8, Dimension(np), Intent(out) :: rp, rw
      Real * 8, Dimension(size(rw)) :: u , x
      Real * 8, Dimension(3) :: q_params
      Real * 8 :: bu !exponential dummy
      Real * 8 :: rx !Transformation of r
      Real * 8 :: drxj !Jacobian of Transformation of r
      Real * 8 :: omrp !dummy var
      integer :: i1
      q_params(1) = 0.d0
      q_params(2) = 1.d0
      q_params(3) = DBLE(np)
      call quad_select("CLCU",q_params,rp,rw)
      do i1=1, size(rw)
        omrp= 1.d0-rp(i1)
        if (omrp>0.d0) then
          if (rp(i1) < 0.999999999999999d0) then
            rx = 1.d0 + (rp(i1)/omrp)
          else
            rx = 1.d16
          endif
          bu = u_yuk(rx,z) / T
          if (-bu > 1.d-16) then
            u(i1) = bu
          else
            u(i1) = 0.d0
          endif
          drxj = omrp**(-2)
          rw(i1) = (exp(-u(i1)) - 1.d0) * rw(i1) * rx * drxj
        else
          rw(i1)=0.d0
        end if
      end do
    End Subroutine

    Function c_pk_hsay_gsvw(k,rp,rw) result(cpk)
      Implicit None
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: rp  !Quadrature points
      Real * 8, Intent(in), dimension(:) :: rw  !Integrand weights (quadrature weights included)
      Real * 8, dimension(size(rp)) :: Integrand  !Integrand
      Real * 8 :: cpk !Perturbation Direct correlation function in k-space
      Real * 8 :: rx !radial distance
      Real * 8 :: omrp,rxk,k2,k4,rx3,rx5,sinrxk!dummy variables
      Integer :: i1
      k2 = k*k
      k4 = k2 * k2
      do i1=1, size(rp)
        omrp= 1.d0-rp(i1)
        if (rp(i1) < 0.999999999999999d0) then
          rx = 1.d0 + (rp(i1)/omrp)
        else
          rx = 1.d16
        endif
        rxk = rx * k
        if (rxk > 0.075d0) then
          sinrxk = sin(rxk) / k
        else
          rx3 = k2 * rx  * rx * rx / 6.d0
          rx5 = k4 * rx3 * rx * rx / 120.d0
          sinrxk = rx - rx3 + rx5
        endif
        Integrand(i1) = rw(i1) * sinrxk
      end do
      cpk = 4.d0 * pi * sum(Integrand)
      return
    End Function

    Function c_hsay_gsvw2(phi,k,rp,rw) result(c)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: rp  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: rw  !Quadrature weights
      Real * 8 :: c
      Integer :: i1
      c = c_hs_vw(phi,k) + c_pk_hsay_gsvw(k,rp,rw)
      return
    End Function

    Function is_hsay_gsvw2(phi,k,rp,rw) result(is)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: rp  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: rw  !Quadrature weights
      Real * 8 :: is, phi_vw
      phi_vw = phi*(1.d0 - (phi / 16.d0))
      is = phi_vw * c_hs_vw(phi,k) + phi * c_pk_hsay_gsvw(k,rp,rw)
      is = 1.d0 - (6.d0 * is / pi)
      return
    End Function

    Function s_hsay_gsvw2(phi,k,rp,rw) result(s)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: rp  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: rw  !Quadrature weights
      Real * 8 :: s
      s = 1.d0 / is_hsay_gsvw2(phi,k,rp,rw)
      return
    End Function

End Module
