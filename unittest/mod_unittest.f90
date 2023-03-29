module unittest
  use global_var, only: wp, stdout
  use fileIO, only: eamFile, Readfile, Atom, ReadData, region
  use computeUE, only: interpolate
  implicit none

  private
  public :: TestRho, TestF, TestZ

contains

  subroutine TestRho(filename)
    !! 测试势函数中对数组ρ的插值
    character(len=*), intent(in) :: filename    !< 势函数文件名称
    type(eamFile) :: eampot                     !< 包含了势函数所有信息和处理之后的数据的集合
    integer :: i                                !< 循环变量
    real(wp) :: rho                             !< 插值计算出来的ρ

    call Readfile(filename, eampot)

    do i = 0, eampot%Nr - 1
      rho = interpolate(eampot%dri, eampot%rx, eampot%rho, i * eampot%dr)
      !! 如果插值和数组对应位置的误差过大, 说明插值不正确
      if (abs(rho - eampot%Rho(i+1)) > 1e-10) stop 'Interpolation error when interpolate array of rho.'
    end do

    write(stdout, '(A)') 'Interpolation of Rho successful!'

  end subroutine TestRho

  subroutine TestF(filename)
    !! 测试势函数中对数组F的插值
    character(len=*), intent(in) :: filename  !< 势函数文件名称
    type(eamFile) :: eampot                   !< 包含了势函数所有信息和处理之后的数据的集合
    integer :: i                              !< 循环变量
    real(wp) :: F                             !< 插值计算出来的ρ

    call Readfile(filename, eampot)

    do i = 0, eampot%Nrho - 1
      F = interpolate(eampot%drhoi, eampot%rhox, eampot%F, i * eampot%drho)
      !! 如果插值和数组对应位置的误差过大, 说明插值不正确
      if (abs(F - eampot%F(i+1)) > 1e-10) stop 'Interpolation error when interpolate array of F.'
    end do

    write(stdout, '(A)') 'Interpolation of F successful!'
  end subroutine TestF

  subroutine TestZ(filename)
    !! 测试势函数中对数组Z的插值
    character(len=*), intent(in) :: filename    !< 势函数文件名称
    type(eamFile) :: eampot                     !< 包含了势函数所有信息和处理之后的数据的集合
    integer :: i                                !< 循环变量
    real(wp) :: Z                               !< 插值计算出来的ρ

    call Readfile(filename, eampot)

    do i = 0, eampot%Nr - 1
      Z = interpolate(eampot%dri, eampot%rx, eampot%Z, i * eampot%dr)
      !! 如果插值和数组对应位置的误差过大, 说明插值不正确
      if (abs(z - eampot%z(i+1)) > 1e-10) stop 'Interpolation error when interpolate array of Z.'
    end do

    write(stdout, '(A)') 'Interpolation of Z successful!'

  end subroutine TestZ

end module unittest
