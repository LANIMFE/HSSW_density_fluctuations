Module heterogeneities_module
  Use nrtype
  Use numerical_recipes
  Use structure_function_selector
  !Use FFT_module
  Implicit None

  Complex * 16, Dimension(:), Allocatable :: C_d_rho_k,C_d_rho_r_u,C_d_rho_k_i
	!Complex * 16, Dimension(:,:), Allocatable :: C_d_rho_k
  Real * 8 :: mean_rho,idk,dkh
  !Integer :: grid_size
  Contains

Subroutine reading_intial_fluctuation(SDimen,box_size,mean_rho,grid,rho_box,np)
  Implicit None
  Integer, intent(in) :: SDimen,grid, np
	Real * 8, intent(in) :: box_size
  Real * 8, Dimension(:), intent(inout):: rho_box
  Real * 8, intent(in) :: mean_rho
  Real * 8 :: x,y,z
  Integer :: i1,i2,i3,i4,ix,iy,iz,ibox

  rho_box=0.d0
  Open(unit=192, File="config_sim_phi_0_1000.dat", Status="Old")
  if (SDimen == 2) then
    Do i1=1, np
      Read(192,*) x,y
      ix = int( x / box_size ) + 1
      iy = int( y / box_size )
      if (ix .le. grid .and. iy < grid ) then
        ibox=ix + (grid*iy)
        rho_box(ibox) = rho_box(ibox) + 1.d0
      endif
    End Do
  else if (SDimen == 3) then
    Read(192,*) x,y,z
    ix = int( x / box_size ) + 1
    iy = int( y / box_size )
    iz = int( z / box_size )
    if (ix .le. grid .and. iy < grid .and. iz < grid ) then
      ibox=ix + (grid*iy) + (grid*grid*iz)
      rho_box(ibox) = rho_box(ibox) + 1.d0
    endif
  endif
  Close(192)
	rho_box = rho_box / (box_size **SDimen)
  rho_box = rho_box - mean_rho
  print*, "Done reading initial fluctuation"
  Print*, "dimen=", sdimen
  print*, "ix=",ix, "iy=",iy
End Subroutine

Subroutine reading_intial_condition_old(SDimen,box_size,phi,grid,box_posx,box_posy,box_posz,rho_box)
  Implicit None
  Integer, intent(in) :: SDimen,grid
	Real * 8, intent(in) :: box_size
  Real * 8, Dimension(:), intent(inout), allocatable :: box_posx,box_posy,box_posz,rho_box
  Real * 8, intent(in) :: phi
  Integer :: i1,i2,i3,i4
  !grid=2**5
	Allocate (rho_box(grid**SDimen))
	Allocate (box_posx(grid**SDimen))
	Allocate (box_posy(grid**SDimen))
	Allocate (box_posz(grid**SDimen))
  mean_rho=2.d0*SDimen*phi / pi
  !box_size=50.d0/DBLE(grid)
  idk=box_size*grid/(2.d0*pi)
  dkh=1.d0/idk
  Open(unit=192, File="delta_n0_32x32x32_L50_kim.dat", Status="Old")
  Do i1=1, grid**SDimen
    Read(192,*) box_posx(i1),box_posy(i1),box_posz(i1),rho_box(i1)
  End Do
  Close(192)
	rho_box=rho_box+mean_rho
  print*, mean_rho
End Subroutine


Subroutine spatial_grid_creation(SDimen,box_size,grid,box_posx,box_posy,box_posz)
  Implicit None
  Integer, intent(in) :: SDimen,grid
	Real * 8, intent(in) :: box_size
	Real * 8, Dimension(:), intent(inout), allocatable :: box_posx,box_posy,box_posz
  Real * 8 :: pos_dum,hgrid
  Integer :: i1,ix,iy,iz,ibox,ixyz
  Logical :: conditional1

  !box_size=5.d0!*pi*sigma(1)!/DBLE(grid_size)
  idk=box_size*grid/(2.d0*pi)
  hgrid = box_size * grid/2
  print*, "dk=",1.d0/idk
  !box_size=pi*sigma(1)*2.d0*idk/DBLE(grid_size)
  Allocate (box_posx(grid**SDimen))
  Allocate (box_posy(grid**SDimen))
  If (SDimen==3) Allocate (box_posz(grid**SDimen))
  ix=0
  iy=0
  iz=0
  !pos_dum=DBLE(grid_size)/2.d0
  Open (unit=120,FILE="grid_test.dat",STATUS="Replace")
  If (SDimen==3) Then
    Do i1=1,grid**SDimen

      If (ix==grid) Then
        Write (120,*) " "
        ix=0
        iy=iy+1
        If (iy==grid) Then
          If(iz<grid)  iz=iz+1
          iy=0
        End If
      End If

      box_posx(i1)=(DBLE(ix))*box_size!-hgrid
      box_posy(i1)=(DBLE(iy))*box_size!-hgrid
  		box_posz(i1)=(DBLE(iz))*box_size!-hgrid

      Write (120,*) SNGL(box_posx(i1)),SNGL(box_posy(i1)),SNGL(box_posz(i1))
      ix=ix+1
    End do
  elseif(SDimen==2) Then
    Do i1=1,grid**SDimen

      If (ix==grid) Then
        Write (120,*) " "
        ix=0
        iy=iy+1
      End If

      box_posx(i1)=(DBLE(ix))*box_size!-hgrid
      box_posy(i1)=(DBLE(iy))*box_size!-hgrid

      Write (120,*) SNGL(box_posx(i1)),SNGL(box_posy(i1))
      ix=ix+1
    End do
  end if
  Close(120)
  print*, "writting test"
  Open (unit=121,FILE="grid_test_v2.dat",STATUS="Replace")
  If (SDimen==3) Then
    Do iz=0, grid-1
      Do iy=0,grid-1
        ibox=0
        ix=0
        Do while (ibox<2**(SDimen-1))
          ix=ix+1
          ixyz=ix + iy*grid + ( iz*grid **2 )
          if (ibox==0) then
            Write (121,*) SNGL(box_posx(ixyz)),SNGL(box_posy(ixyz)),SNGL(box_posz(ixyz))
            Write (121,*) SNGL(box_posx(ixyz)+box_size) , SNGL(box_posy(ixyz)) , SNGL(box_posz(ixyz))
          else if (ibox==1) then
            Write (121,*) SNGL(box_posx(ixyz)) , SNGL(box_posy(ixyz)+box_size) , SNGL(box_posz(ixyz))
            Write (121,*) SNGL(box_posx(ixyz)+box_size) , SNGL(box_posy(ixyz)+box_size) , SNGL(box_posz(ixyz))
          else if (ibox==2) then
            Write (121,*) SNGL(box_posx(ixyz)) , SNGL(box_posy(ixyz)) , SNGL(box_posz(ixyz)+box_size)
            Write (121,*) SNGL(box_posx(ixyz)+box_size) , SNGL(box_posy(ixyz)) , SNGL(box_posz(ixyz)+box_size)
          else if (ibox==3) then
            Write (121,*) SNGL(box_posx(ixyz)) , SNGL(box_posy(ixyz)+box_size) , SNGL(box_posz(ixyz)+box_size)
            Write (121,*) SNGL(box_posx(ixyz)+box_size) , SNGL(box_posy(ixyz)+box_size) , SNGL(box_posz(ixyz)+box_size)
          endif
          If (ix==grid) Then
            Write (121,*) " "
            ix = 0
            ibox=ibox+1
          End If
        End do
      End do
    Enddo
  elseif(SDimen==2) Then
    Do iy=0,grid-1
      ibox=0
      ix=0
      Do while (ibox<2**(SDimen-1))
        ix=ix+1
        ixyz=ix + iy*grid


        if (ibox==0) then
          Write (121,*) box_posx(ixyz),box_posy(ixyz)
          Write (121,*) box_posx(ixyz)+box_size,box_posy(ixyz)
        else if (ibox==1) then
          Write (121,*) box_posx(ixyz),box_posy(ixyz)+box_size
          Write (121,*) box_posx(ixyz)+box_size,box_posy(ixyz)+box_size
        endif

        If (ix==grid) Then
          Write (121,*) " "
          ix = 0
          ibox=ibox+1
        End If

      End do
    end do
  endif
  Close(121)
End Subroutine

Subroutine spatial_grid_creation_slab(SDimen,box_size,grid,box_posx,box_posy,box_posz)
  Implicit None
  Integer, intent(in) :: SDimen,grid
	Real * 8, intent(in) :: box_size
	Real * 8, Dimension(:), intent(inout), allocatable :: box_posx,box_posy,box_posz
  Real * 8 :: pos_dum,hgrid
  Integer :: i1,ix,iy,iz,ibox,ixyz
  Logical :: conditional1

  !box_size=5.d0!*pi*sigma(1)!/DBLE(grid_size)
  idk=box_size*grid/(2.d0*pi)
  hgrid = box_size * grid/2
  print*, "dk=",1.d0/idk
  !box_size=pi*sigma(1)*2.d0*idk/DBLE(grid_size)
  Allocate (box_posx(grid**SDimen))
  Allocate (box_posy(grid**SDimen))
  If (SDimen==3) Allocate (box_posz(grid**SDimen))
  ix=0
  iy=0
  iz=0
  !pos_dum=DBLE(grid_size)/2.d0
  Open (unit=120,FILE="grid_test.dat",STATUS="Replace")
  If (SDimen==3) Then
    Do i1=1,grid**SDimen

      If (ix==grid) Then
        Write (120,*) " "
        ix=0
        iy=iy+1
        If (iy==grid) Then
          If(iz<grid)  iz=iz+1
          iy=0
        End If
      End If

      box_posx(i1)=(DBLE(ix))*box_size!-hgrid
      box_posy(i1)=(DBLE(iy))*box_size!-hgrid
  		box_posz(i1)=(DBLE(iz))*box_size!-hgrid

      Write (120,*) SNGL(box_posx(i1)),SNGL(box_posy(i1)),SNGL(box_posz(i1))
      ix=ix+1
    End do
  elseif(SDimen==2) Then
    Do i1=1,grid**SDimen

      If (ix==grid) Then
        Write (120,*) " "
        ix=0
        iy=iy+1
      End If

      box_posx(i1)=(DBLE(ix))*box_size!-hgrid
      box_posy(i1)=(DBLE(iy))*box_size!-hgrid

      Write (120,*) SNGL(box_posx(i1)),SNGL(box_posy(i1))
      ix=ix+1
    End do
  end if
  Close(120)
  print*, "writting test"
  Open (unit=121,FILE="grid_test_v2.dat",STATUS="Replace")
    Do iy=0,grid-1
      ibox=0
      ix=0
      Do while (ibox<2)
        ix=ix+1
        ixyz=ix + iy*grid


        if (ibox==0) then
          Write (121,*) box_posx(ixyz),box_posy(ixyz)
          Write (121,*) box_posx(ixyz)+box_size,box_posy(ixyz)
        else if (ibox==1) then
          Write (121,*) box_posx(ixyz),box_posy(ixyz)+box_size
          Write (121,*) box_posx(ixyz)+box_size,box_posy(ixyz)+box_size
        endif

        If (ix==grid) Then
          Write (121,*) " "
          ix = 0
          ibox=ibox+1
        End If

      End do
    end do

  Close(121)
End Subroutine

Subroutine calc_equilibrium_density_deviation(srho,mean,box_size,SDimen,phi,grid,rho_box)
  Implicit None
  Integer, intent(in) :: SDimen,grid
  Real * 8, intent(in) :: phi,srho,mean
	Real * 8, intent(in) :: box_size
  Real * 8, Dimension(:), intent(out) :: rho_box
  Integer :: i1,ix,ibox,iy,iz,ixyz
  Integer * 4 :: iseed1
  Real * 8 :: sigma_rho,Dran
  Real * 8 :: X_0,rho_max
  Real * 8 :: r2,mrho,rcp,x,pdf,truncated_mean,new_mean,truncated_variance,new_variance
	Real * 8 :: nt_variance,yx,b,mean_realization
  iseed1=123
  print *, "sigma rho=", srho
  x_0=(grid/2)*box_size
	if (SDimen == 3) rcp = 0.64d0
	if (SDimen == 2) rcp = 0.80d0
  rho_max = rcp * 2.d0 * SDimen / pi
	mrho = phi *  2.d0 * SDimen / pi
  Print*, mrho,phi,srho,box_size
  Open (unit=130, File="density_test.dat",Status="Replace")
	Open (unit=230, File="density_pdf.dat",Status="Replace")
	write(230,*) "### x // pdf"
	x=0.d0
	Do while (x<rho_max)
		call truncated_normal_ab_pdf ( x, mrho, srho, 0.d0, rho_max, pdf )
		write(230,*) x,pdf
		x=x+1.d-5
	End Do
	close(230)
  ix=0

	call truncated_normal_ab_mean ( mean, srho, -mrho, rho_max-mrho, truncated_mean )
	call truncated_normal_ab_variance ( mean, srho, -mrho, rho_max-mrho, truncated_variance )
	print*, truncated_mean,truncated_variance


  Do i1=1,grid**SDimen
    !Dran=!gasdev_ran2_DP(iseed1)
    !Do while(DAbs(Dran*sigma_rho)>mean_rho)
      !Dran=gasdev_ran2_DP(iseed1)
    !End Do
    If (box_size <2.0d0 ) then
    !  x=r8_uniform_01 ( iseed1 )
    !  yx = b * (1.d0 - sqrt(1.d0-x))
    !  Dran = yx - mean
    rho_box(i1)=-mean
    if (r8_uniform_01 ( iseed1 ) .le. mean) rho_box(i1)=1.d0-mean
    else
      call truncated_normal_ab_sample( mean, srho, -mrho, rho_max-mrho, iseed1, Dran )
      rho_box(i1)= Dran
    endif

    !r2=(box_posx(i1))**2+(box_posy(i1))**2
    !rho_box(i1)=exp(-r2) !TEST FUNCTION GAUSSIAN
    Write(130,*) rho_box(i1)
    flush (130)
    ix=ix+1
    if (ix==grid) Then
      ix=0
      Write(130,*) " "
    endif
  End Do
  mean_realization=sum(rho_box)/dble(grid**SDimen)
  print*, "mean density of the realization=", mean_realization,mean_realization*pi/4.d0
  Close(130)
  Open (unit=131, File="density_test_v2.dat",Status="Replace")
  ix=0

  If (SDimen==2) Then
    Do iy=0, grid-1
      ibox=0
      Do while (ibox < 2**(SDimen-1))
        Do i1=1,grid
          ixyz=i1+iy*grid
          Write(131,*) SNGL(rho_box(ixyz))
          Write(131,*) SNGL(rho_box(ixyz))
          flush (131)
          ix=ix+1
          if (ix==grid) Then
            ix=0
            ibox=ibox+1
            Write(131,*) " "
          endif
        End Do
      End Do
    End Do
  else if (SDimen==3) Then
    Do iz=0, grid-1
      Do iy=0, grid-1
        ibox=0
        Do while (ibox < 2**(SDimen-1))
          Do i1=1,grid
            ixyz=i1+iy*grid+(iz*grid**2)
            Write(131,*) SNGL(rho_box(ixyz))
            Write(131,*) SNGL(rho_box(ixyz))
            flush (131)
            ix=ix+1
            if (ix==grid) Then
              ix=0
              ibox=ibox+1
              Write(131,*) " "
            endif
          End Do
        End Do
      End Do
    End Do
  End if
  Close(131)
  return
End Subroutine

Subroutine calc_equilibrium_density_deviation_slab(srho,mean,box_size,SDimen,phi,grid,rho_box)
  Implicit None
  Integer, intent(in) :: SDimen,grid
  Real * 8, intent(in) :: phi,srho,mean
	Real * 8, intent(in) :: box_size
  Real * 8, Dimension(:), intent(out) :: rho_box
  Integer :: i1,ix,ibox,iy,iz,ixyz
  Integer * 4 :: iseed1
  Real * 8 :: sigma_rho,Dran
  Real * 8 :: X_0,rho_max
  Real * 8 :: r2,mrho,rcp,x,pdf,truncated_mean,new_mean,truncated_variance,new_variance
	Real * 8 :: nt_variance,yx,b,mean_realization
  iseed1=123
  print *, "sigma rho=", srho
  x_0=(grid/2)*box_size
	if (SDimen == 3) rcp = 0.64d0
	if (SDimen == 2) rcp = 0.80d0
  rho_max = rcp * 2.d0 * SDimen / pi
	mrho = phi *  2.d0 * SDimen / pi
  Print*, mrho,phi,srho,box_size
  Open (unit=130, File="density_test.dat",Status="Replace")
	Open (unit=230, File="density_pdf.dat",Status="Replace")
	write(230,*) "### x // pdf"
	x=0.d0
	Do while (x<rho_max)
		call truncated_normal_ab_pdf ( x, mrho, srho, 0.d0, rho_max, pdf )
		write(230,*) x,pdf
		x=x+1.d-5
	End Do
	close(230)
  ix=0

	call truncated_normal_ab_mean ( mean, srho, -mrho, rho_max-mrho, truncated_mean )
	call truncated_normal_ab_variance ( mean, srho, -mrho, rho_max-mrho, truncated_variance )
	print*, truncated_mean,truncated_variance


  Do i1=1,grid**SDimen
    !Dran=!gasdev_ran2_DP(iseed1)
    !Do while(DAbs(Dran*sigma_rho)>mean_rho)
      !Dran=gasdev_ran2_DP(iseed1)
    !End Do
    If (box_size <2.0d0 ) then
    !  x=r8_uniform_01 ( iseed1 )
    !  yx = b * (1.d0 - sqrt(1.d0-x))
    !  Dran = yx - mean
    rho_box(i1)=-mean
    if (r8_uniform_01 ( iseed1 ) .le. mean) rho_box(i1)=1.d0-mean
    else
      call truncated_normal_ab_sample( mean, srho, -mrho, rho_max-mrho, iseed1, Dran )
      rho_box(i1)= Dran
    endif

    !r2=(box_posx(i1))**2+(box_posy(i1))**2
    !rho_box(i1)=exp(-r2) !TEST FUNCTION GAUSSIAN
    Write(130,*) rho_box(i1)
    flush (130)
    ix=ix+1
    if (ix==grid) Then
      ix=0
      Write(130,*) " "
    endif
  End Do
  mean_realization=sum(rho_box)/dble(grid**SDimen)
  print*, "mean density of the realization=", mean_realization,mean_realization*pi/4.d0
  Close(130)
  Open (unit=131, File="density_test_v2.dat",Status="Replace")
  ix=0
    Do iy=0, grid-1
      ibox=0
      Do while (ibox < 2)
        Do i1=1,grid
          ixyz=i1+iy*grid
          Write(131,*) SNGL(rho_box(ixyz))
          Write(131,*) SNGL(rho_box(ixyz))
          flush (131)
          ix=ix+1
          if (ix==grid) Then
            ix=0
            ibox=ibox+1
            Write(131,*) " "
          endif
        End Do
      End Do
    End Do
  Close(131)
  return
End Subroutine

Subroutine calc_k_grid_FT(SDimen,box_size,grid,box_poskx,box_posky,box_poskz,box_k2)
  Implicit None
  Integer, intent(in) :: SDimen,grid
	Real * 8, intent(in) :: box_size
  Integer :: i1,ix,iy,iz,i2,i3,idum,idumy,ibox,ixyz
  Integer, Dimension(SDimen) :: nn
  Real * 8 :: kx,ky,kz,x_0,dum1,dum2
  Real * 8, Dimension(:), Allocatable, intent (inout) :: box_poskx,box_posky,box_poskz,box_k2
  dkh=1.d0/idk
  Do i1=1,SDimen
    nn(i1)=grid
  End Do
  Allocate (box_poskx(grid**SDimen))
  Allocate (box_posky(grid**SDimen))
  If (SDimen==3) Allocate (box_poskz(grid**SDimen))
	Allocate (box_k2(grid**SDimen))
  Open(unit=140,File="k_grid_test.dat",Status="Replace")
  idum=0
  ix=idum
  iy=0
  iz=0
  x_0=DBLE((grid)/2)*box_size!+3*box_size
  Do i1=idum+1, grid**SDimen
    If (ix<=grid/2) kx=DBLE(ix)*dkh
    If (ix>grid/2) kx=-DBLE(grid-ix)*dkh
    If (iy<=grid/2) ky=DBLE(iy)*dkh
    If (iy>grid/2) ky=-DBLE(grid-iy)*dkh
    If (iz<=grid/2) kz=DBLE(iz)*dkh
    If (iz>grid/2) kz=-DBLE(grid-iz)*dkh
    box_poskx(i1) = kx
    box_posky(i1) = ky
		box_k2(i1)= kx*kx + ky*ky
    If (SDimen==3) then
      box_poskz(i1) = kz
			box_k2(i1)= box_k2(i1) + kz * kz
      Write(140,*) SNGL(kx),SNGL(ky),SNGL(kz)
    else if (SDimen==2) then
      Write(140,*) SNGL(kx),SNGL(ky)
    endif
    ix=ix+1
    If(ix==grid) Then
      Write(140,*) ""
      ix=0
      iy=iy+1
      If(iy==grid) Then
        iy=0
        If(SDimen>2.and.iz<=grid)iz=iz+1
      End If
    End If
    !print*, i1,"holis",iy,ix,dum1
  End Do
  print*, grid**SDimen
  Close(140)

  Open (unit=121,FILE="k_grid_test_v2.dat",STATUS="Replace")
  If (SDimen==3) Then
    Do iz=0, grid-1
      Do iy=0,grid-1
        ibox=0
        ix=0
        Do while (ibox<2**(SDimen-1))
          ix=ix+1
          ixyz=ix + iy*grid + (iz*grid **2 )
          if (ibox==0) then
            Write (121,*) SNGL(box_poskx(ixyz)),SNGL(box_posky(ixyz)),SNGL(box_poskz(ixyz))
            Write (121,*) SNGL(box_poskx(ixyz)+dkh) , SNGL(box_posky(ixyz)) , SNGL(box_poskz(ixyz))
          else if (ibox==1) then
            Write (121,*) SNGL(box_poskx(ixyz)) , SNGL(box_posky(ixyz)+dkh) , SNGL(box_poskz(ixyz))
            Write (121,*) SNGL(box_poskx(ixyz)+dkh) , SNGL(box_posky(ixyz)+dkh) , SNGL(box_poskz(ixyz))
          else if (ibox==2) then
            Write (121,*) SNGL(box_poskx(ixyz)) , SNGL(box_posky(ixyz)) , SNGL(box_poskz(ixyz)+dkh)
            Write (121,*) SNGL(box_poskx(ixyz)+dkh) , SNGL(box_posky(ixyz)) , SNGL(box_poskz(ixyz)+dkh)
          else if (ibox==3) then
            Write (121,*) SNGL(box_poskx(ixyz)) , SNGL(box_posky(ixyz)+dkh) , SNGL(box_poskz(ixyz)+dkh)
            Write (121,*) SNGL(box_poskx(ixyz)+dkh) , SNGL(box_posky(ixyz)+dkh) , SNGL(box_poskz(ixyz)+dkh)
          endif
          If (ix==grid) Then
            Write (121,*) " "
            ix = 0
            ibox=ibox+1
          End If
        End do
      End do
    Enddo
  elseif(SDimen==2) Then
    Do iy=0,grid-1
      ibox=0
      ix=0
      Do while (ibox<2**(SDimen-1))
        ix=ix+1
        ixyz=ix + iy*grid

        if (ibox==0) then
          Write (121,*) SNGL(box_poskx(ixyz)),SNGL(box_posky(ixyz))
          Write (121,*) SNGL(box_poskx(ixyz)+dkh),SNGL(box_posky(ixyz))
        else if (ibox==1) then
          Write (121,*) SNGL(box_poskx(ixyz)),SNGL(box_posky(ixyz)+dkh)
          Write (121,*) SNGL(box_poskx(ixyz)+dkh),SNGL(box_posky(ixyz)+dkh)
        endif

        If (ix==grid) Then
          Write (121,*) " "
          ix = 0
          ibox=ibox+1
        End If

      End do
    end do
  endif
  Close(121)
End Subroutine

Subroutine calc_k_grid_FT_slab(SDimen,box_size,grid,box_poskx,box_posky,box_poskz,box_k2)
  Implicit None
  Integer, intent(in) :: SDimen,grid
	Real * 8, intent(in) :: box_size
  Integer :: i1,ix,iy,iz,i2,i3,idum,idumy,ibox,ixyz
  Integer, Dimension(SDimen) :: nn
  Real * 8 :: kx,ky,kz,x_0,dum1,dum2
  Real * 8, Dimension(:), Allocatable, intent (inout) :: box_poskx,box_posky,box_poskz,box_k2
  dkh=1.d0/idk
  Do i1=1,SDimen
    nn(i1)=grid
  End Do
  Allocate (box_poskx(grid**SDimen))
  Allocate (box_posky(grid**SDimen))
  If (SDimen==3) Allocate (box_poskz(grid**SDimen))
	Allocate (box_k2(grid**SDimen))
  Open(unit=140,File="k_grid_test.dat",Status="Replace")
  idum=0
  ix=idum
  iy=0
  iz=0
  x_0=DBLE((grid)/2)*box_size!+3*box_size
  Do i1=idum+1, grid**SDimen
    If (ix<=grid/2) kx=DBLE(ix)*dkh
    If (ix>grid/2) kx=-DBLE(grid-ix)*dkh
    If (iy<=grid/2) ky=DBLE(iy)*dkh
    If (iy>grid/2) ky=-DBLE(grid-iy)*dkh
    If (iz<=grid/2) kz=DBLE(iz)*dkh
    If (iz>grid/2) kz=-DBLE(grid-iz)*dkh
    box_poskx(i1) = kx
    box_posky(i1) = ky
		box_k2(i1)= kx*kx + ky*ky
    If (SDimen==3) then
      box_poskz(i1) = kz
			box_k2(i1)= box_k2(i1) + kz * kz
      Write(140,*) SNGL(kx),SNGL(ky),SNGL(kz)
    else if (SDimen==2) then
      Write(140,*) SNGL(kx),SNGL(ky)
    endif
    ix=ix+1
    If(ix==grid) Then
      Write(140,*) ""
      ix=0
      iy=iy+1
      If(iy==grid) Then
        iy=0
        If(SDimen>2.and.iz<=grid)iz=iz+1
      End If
    End If
    !print*, i1,"holis",iy,ix,dum1
  End Do
  print*, grid**SDimen
  Close(140)

  Open (unit=121,FILE="k_grid_test_v2.dat",STATUS="Replace")

    Do iy=0,grid-1
      ibox=0
      ix=0
      Do while (ibox<2)
        ix=ix+1
        ixyz=ix + iy*grid

        if (ibox==0) then
          Write (121,*) SNGL(box_poskx(ixyz)),SNGL(box_posky(ixyz))
          Write (121,*) SNGL(box_poskx(ixyz)+dkh),SNGL(box_posky(ixyz))
        else if (ibox==1) then
          Write (121,*) SNGL(box_poskx(ixyz)),SNGL(box_posky(ixyz)+dkh)
          Write (121,*) SNGL(box_poskx(ixyz)+dkh),SNGL(box_posky(ixyz)+dkh)
        endif

        If (ix==grid) Then
          Write (121,*) " "
          ix = 0
          ibox=ibox+1
        End If

      End do
    end do
  Close(121)
End Subroutine

Subroutine calc_density_spatial_fluctuation_FT(SDimen,box_size,rho_box,FT_rho_box,grid)
  Implicit None
  Integer, intent(in) :: SDimen,grid
	Real * 8, intent(in) :: box_size
  Real * 8, Dimension(:), intent(in) :: rho_box
  Complex * 16, Dimension(:), intent(out) :: FT_rho_box
  Complex * 16, Dimension(:), allocatable :: IFT_rho_box
  Integer :: i1,ix,iy,iz,i2,i3,idum,idumy,ibox,ixyz
  Integer, Dimension(SDimen) :: nn
  Real * 8 :: kx,ky,kz,x_0,dum1,dum2
  Do i1=1,SDimen
    nn(i1)=grid
  End Do
  allocate(IFT_rho_box(grid**SDimen))
  Do i1=1, grid**SDimen
    FT_rho_box(i1)=CMPLX(rho_box(i1),0.d0)
    !C_d_rho_k(i1)=CMPLX(rho_box(i1)-mean_rho,0.d0)
  End Do
    print*, "calling FT"
	Call fourn_gather(FT_rho_box,nn,-1)
  FT_rho_box  = FT_rho_box  * (box_size**SDimen) !Renormalization of the FT

  IFT_rho_box=FT_rho_box
  Call fourn_gather(IFT_rho_box,nn,1)
  IFT_rho_box = IFT_rho_box * ((dkh/(2.d0*pi))**SDimen) !Renormalization of the IFT
  !Call fourn_gather(C_d_rho_k,nn,1)
  Open(unit=140,File="k_density_test.dat",Status="Replace")
  idum=0
  ix=idum
  iy=0
  iz=0

  x_0=DBLE((grid)/2)*box_size!+3*box_size
  Do i1=idum+1, grid**SDimen
    !*Cexp(CMPLX(0.d0,(kx+ky+kz)*x_0))
    dum1=Real(FT_rho_box(i1))
    dum2=Imag(FT_rho_box(i1))
  !  Write(140,*) kx,ky,kz,dum1,dum2
    Write(140,*) dum1,dum2,Real(IFT_rho_box(i1))
    ix=ix+1
    If(ix==grid) Then
      Write(140,*) ""
      ix=0
      iy=iy+1
      If(iy==grid) Then
        iy=0
        If(SDimen>2.and.iz<=grid)iz=iz+1
      End If
    End If
    !print*, i1,"holis",iy,ix,dum1
  End Do
  Close(140)
  Open (unit=131, File="k_density_test_v2.dat",Status="Replace")
  ix=0
  If (SDimen==2) Then
    Do iy=0, grid-1
      ibox=0
      Do while (ibox < 2**(SDimen-1))
        Do i1=1,grid
          ixyz=i1+iy*grid
          Write(131,*) SNGL(Real(FT_rho_box(ixyz)))
          Write(131,*) SNGL(Real(FT_rho_box(ixyz)))
          flush (131)
          ix=ix+1
          if (ix==grid) Then
            ix=0
            ibox=ibox+1
            Write(131,*) " "
          endif
        End Do
      End Do
    End Do
  else if (SDimen==3) Then
    Do iz=0, grid-1
      Do iy=0, grid-1
        ibox=0
        Do while (ibox < 2**(SDimen-1))
          Do i1=1,grid
            ixyz=i1+iy*grid+(iz*grid**2)
            Write(131,*) SNGL(Real(FT_rho_box(ixyz)))
            Write(131,*) SNGL(Real(FT_rho_box(ixyz)))
            flush (131)
            ix=ix+1
            if (ix==grid) Then
              ix=0
              ibox=ibox+1
              Write(131,*) " "
            endif
          End Do
        End Do
      End Do
    End Do
  End if
  Close(131)
End Subroutine

Subroutine calc_density_spatial_fluctuation_FT_slab(SDimen,box_size,rho_box,FT_rho_box,grid)
  Implicit None
  Integer, intent(in) :: SDimen,grid
	Real * 8, intent(in) :: box_size
  Real * 8, Dimension(:), intent(in) :: rho_box
  Complex * 16, Dimension(:), intent(out) :: FT_rho_box
  Complex * 16, Dimension(:), allocatable :: IFT_rho_box
  Integer :: i1,ix,iy,iz,i2,i3,idum,idumy,ibox,ixyz
  Integer, Dimension(SDimen) :: nn
  Real * 8 :: kx,ky,kz,x_0,dum1,dum2
  Do i1=1,SDimen
    nn(i1)=grid
  End Do
  allocate(IFT_rho_box(grid**SDimen))
  Do i1=1, grid**SDimen
    FT_rho_box(i1)=CMPLX(rho_box(i1),0.d0)
    !C_d_rho_k(i1)=CMPLX(rho_box(i1)-mean_rho,0.d0)
  End Do
    print*, "calling FT"
	Call fourn_gather(FT_rho_box,nn,-1)
  FT_rho_box  = FT_rho_box  * (box_size**SDimen) !Renormalization of the FT

  IFT_rho_box=FT_rho_box
  Call fourn_gather(IFT_rho_box,nn,1)
  IFT_rho_box = IFT_rho_box * ((dkh/(2.d0*pi))**SDimen) !Renormalization of the IFT
  !Call fourn_gather(C_d_rho_k,nn,1)
  Open(unit=140,File="k_density_test.dat",Status="Replace")
  idum=0
  ix=idum
  iy=0
  iz=0

  x_0=DBLE((grid)/2)*box_size!+3*box_size
  Do i1=idum+1, grid**SDimen
    !*Cexp(CMPLX(0.d0,(kx+ky+kz)*x_0))
    dum1=Real(FT_rho_box(i1))
    dum2=Imag(FT_rho_box(i1))
  !  Write(140,*) kx,ky,kz,dum1,dum2
    Write(140,*) dum1,dum2,Real(IFT_rho_box(i1))
    ix=ix+1
    If(ix==grid) Then
      Write(140,*) ""
      ix=0
      iy=iy+1
      If(iy==grid) Then
        iy=0
        If(SDimen>2.and.iz<=grid)iz=iz+1
      End If
    End If
    !print*, i1,"holis",iy,ix,dum1
  End Do
  Close(140)
  Open (unit=131, File="k_density_test_v2.dat",Status="Replace")
  ix=0
    Do iy=0, grid-1
      ibox=0
      Do while (ibox < 2)
        Do i1=1,grid
          ixyz=i1+iy*grid
          Write(131,*) SNGL(Real(FT_rho_box(ixyz)))
          Write(131,*) SNGL(Real(FT_rho_box(ixyz)))
          flush (131)
          ix=ix+1
          if (ix==grid) Then
            ix=0
            ibox=ibox+1
            Write(131,*) " "
          endif
        End Do
      End Do
    End Do
  Close(131)
End Subroutine

Subroutine calc_density_spatial_fluctuation_of_u(u_time,pre_prop_n,kd_rho_box0,SDimen,grid,drho_box_u)
  Implicit None
  Integer, intent(in) :: SDimen,grid
  Real * 8, Intent (in) :: u_time
	Real * 8, Dimension(:), Intent(in) :: pre_prop_n
	Complex * 16, Dimension(:), Intent(in) :: kd_rho_box0
	Real * 8, Dimension(:), Intent(out) :: drho_box_u
	Complex * 16, Dimension(:), Allocatable :: integrand
	Integer, Dimension(:), Allocatable :: nn
	Integer :: i1
  Allocate(nn(SDimen))
  Do i1=1, SDimen
    nn(i1)=grid
  End Do
  Allocate(integrand(grid**SDimen))
  Do i1=1, grid**SDimen
    !if( DEXP(u_time*pre_prop_n(i1)) < 1.d0 ) then
    !  integrand(i1)=0.d0
    !else
    !  integrand(i1)=kd_rho_box0(i1)
    !endif
		integrand(i1) = kd_rho_box0(i1)* DEXP(u_time*pre_prop_n(i1))
	End Do
	Call fourn_gather(integrand,nn,1)
	drho_box_u = Real(integrand) * ( (dkh/(2.d0*pi) )**SDimen )
	Deallocate(integrand)
End Subroutine

Subroutine writting_drho_box_v1(i_unit,drho_box_u,grid)
	Implicit None
	Integer, Intent(in) :: i_unit,grid
	Real * 8, Dimension(:), Intent(in) :: drho_box_u
	Integer :: i1,ix,sizerho
	ix=0
	sizerho=size(drho_box_u)
	Do i1=1, sizerho
    Write(i_unit,*) drho_box_u(i1)
    ix=ix+1
    If(ix==grid) Then
      Write(i_unit,*) ""
			ix=0
    End If
  End Do
End Subroutine

Subroutine writting_drho_box_v2(i_unit,drho_box_u,grid,SDimen)
	Implicit None
	Integer, Intent(in) :: i_unit,grid,SDimen
	Real * 8, Dimension(:), Intent(in) :: drho_box_u
	Integer :: i1,ix,sizerho,iy,iz,ibox,ixyz
	If (SDimen==3) Then
    Do iz=0, grid-1
      Do iy=0,grid-1
        ibox=0
        ix=0
        Do while (ibox<2**(SDimen-1))
          ix=ix+1
          ixyz=ix + iy*grid + (iz*grid **2)
          Write (i_unit,*) SNGL(drho_box_u(ixyz))
          Write (i_unit,*) SNGL(drho_box_u(ixyz))
          If (ix==grid) Then
            Write (i_unit,*) " "
            ix = 0
            ibox=ibox+1
          End If
        End do
      End do
    Enddo
  elseif(SDimen==2) Then
    Do iy=0,grid-1
      ibox=0
      ix=0
      Do while (ibox<2**(SDimen-1))
        ix=ix+1
        ixyz=ix + iy*grid
				Write (i_unit,*) SNGL(drho_box_u(ixyz))
				Write (i_unit,*) SNGL(drho_box_u(ixyz))
				If (ix==grid) Then
          Write (i_unit,*) ""
          ix = 0
          ibox=ibox+1
        End If
      End do
    end do
  endif
  Close(i_unit)
End Subroutine

Subroutine writting_drho_box_v2_slab(i_unit,drho_box_u,grid,SDimen)
	Implicit None
	Integer, Intent(in) :: i_unit,grid,SDimen
	Real * 8, Dimension(:), Intent(in) :: drho_box_u
	Integer :: i1,ix,sizerho,iy,iz,ibox,ixyz
  Do iy=0,grid-1
    ibox=0
    ix=0
    Do while (ibox<2)
      ix=ix+1
      ixyz=ix + iy*grid
			Write (i_unit,*) SNGL(drho_box_u(ixyz))
			Write (i_unit,*) SNGL(drho_box_u(ixyz))
			If (ix==grid) Then
        Write (i_unit,*) ""
        ix = 0
        ibox=ibox+1
      End If
    End do
  end do
  Close(i_unit)
End Subroutine



End Module
