Module hdsw_structure
  Use nrtype
  Use hd_structure
  Implicit None
  Contains




    !Function that calculates the structure factor of a Hard Disk +
    !Square Well System using Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Quadrature points vector "qn", the range should be qn E [1:lambda] where lambda is the well length
    !Quadrature weights "qw"
    Function c_hdsw_shvw(phi,T,z,k,qn,qw) result(c)
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
      Real * 8 :: dum1,ck_dum,y, j0
      cay = 0.d0
      !open(unit=221, file="c_test.dat",status="replace")
      Do i1=1, size(qn)
        r    = qn(i1)
        kr = k * r
        if (kr > 0.075d0) then
          j0   = BESSEl_JN(0, kr)
        else
          j0   = 1.d0 - (0.5d0*kr*kr)
        endif
        u = - 1.d
        w    = qw(i1) * r * j0
        cay = cay - (u*w)
        !write(221,*) r, u*w
      End Do
      !close(221)
      cay = cay / T
      c = ck_HD_ros(phi,k) + 2.d0*pi*cay
      return
    End Function

    Function is_hdsw_shvw(phi,T,z,k,qn,qw) result(is)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent (in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: is
      is =  phi * c_hday_shvw(phi,T,z,k,qn,qw)
      is = 1.d0 - 4.d0 * is / pi
      return
    End Function

    !Function that calculates the structure factor of a Hard Sphere +
    !Attractive Yukawa System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function s_hdsw_shvw(phi,T,z,k,qn,qw) result(s)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: s
      s = 1.d0 / is_hday_shvw(phi,T,z,k,qn,qw)
      return
    End Function



    !Function that calculates the structure factor of a Hard Disk +
    !Square Well System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Quadrature points vector "qn", the range should be qn E [1:lambda] where lambda is the well length
    !Quadrature weights "qw"
    Function c_hdsw_gsvw(phi,T,z,k,qn,qw) result(c)
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
      Real * 8 :: dum1,ck_dum,y, j0,dum
      cay = 0.d0
      do i1=1, size(qn)
        r   = qn(i1)
        kr  = k * r
        j0  = BESSEl_JN(0, kr)
        u   = - 1.d0
        w   = qw(i1) * r * j0
        dum = -u/T
        if (dum > 1.d-16 ) cay = cay + ( ( exp(dum) - 1.d0 ) * w )
      end do
      c = ck_HD_ros(phi,k) + 2.d0*pi*cay
      return
    End Function



    Function is_hdsw_gsvw(phi,T,z,k,qn,qw) result(is)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent (in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: is
      is =  phi * c_hday_gsvw(phi,T,z,k,qn,qw)
      is = 1.d0 - 4.d0 * is / pi
      return
    End Function



    !Function that calculates the structure factor of a Hard Sphere +
    !Attractive Yukawa System using Generalized Sharma-Sharma approximation Scheme and Verlet-Weiss
    ! S(k,ϕ,z)
    !wave-vector "k" for "k>=0"
    !density "ϕ" for "ϕ>0"
    !temperature "T" for "T>0"
    !Yukawa potential length parameter "z>0"
    Function s_hdsw_gsvw(phi,T,z,k,qn,qw) result(s)
      Implicit None
      Real * 8, Intent(in) :: phi !Volume fraction value
      Real * 8, Intent(in) :: T !Temperature
      Real * 8, Intent(in) :: z !Yukawa potential length parameter
      Real * 8, Intent(in) :: k !Value of the wave vector
      Real * 8, Intent(in), dimension(:) :: qn  !Quadrature nodes
      Real * 8, Intent(in), dimension(:) :: qw  !Quadrature weights
      Real * 8 :: s
      s = 1.d0 / is_hday_gsvw(phi,T,z,k,qn,qw)
      return
    End Function



End Module
