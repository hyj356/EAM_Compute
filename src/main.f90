program main
  !! 该程序目前还不能正常计算模型势能, 因为其中某些原子
  !! 的电荷密度rho可能会超过势函数文件可以提供的上限
  use global_var, only: wp, stdout
  use fileIO, only: eamFile, Readfile, Atom, ReadData, region
  use computeUE, only: Pair, interpolate, Energy
  implicit none
  type(eamFile) :: cu_u3                !< 包含势函数的数据和处理之后的数据的类
  type(Atom), allocatable :: cuSC(:)    !< 单晶Cu的模型文件
  type(region)  :: box                  !< 仿真盒子
  real(wp) :: power

  call Readfile('Cu_u3.eam', cu_u3)
  call ReadData('cu.lmp', box, cuSC)
  call Energy(cu_u3, cuSC, power, box)    !! 根据势函数和模型文件计算势能

  write(stdout, *) power, ' ev'           !! 输出势能
  ! lammps的输出: -14160.0000090988 ev

end program main
 !! 输出插值结果查看, 以供debug
 ! write(stdout, '(ES9.2)') interpolate(cu_u3%dri, cu_u3%rx, cu_u3%Z, cu_u3%dr)                !! 插值计算有效电荷Z值, 正确输出: 1.08E+01
 ! write(stdout, '(ES9.2)') interpolate(cu_u3%drhoi, cu_u3%rhox, cu_u3%F, 2.d0 * cu_u3%drho)   !! 插值计算嵌入能F值, 正确输出: -5.23E-01
 ! write(stdout, '(ES9.2)') interpolate(cu_u3%dri, cu_u3%rx, cu_u3%rho, 3.615d0)               !! 插值计算电荷密度, 正确输出: 5.44E-05
 ! write(stdout, '(ES9.2)') Pair(cu_u3, cu_u3%dr)                                              !! 计算2体势能, 正确输出: 1.68E+05
 ! write(stdout, '(ES9.2)') cu_u3%cutoff             !! 截断半径
 ! write(stdout, *) maxval(cu_u3%rhox)
 !! 参考教程: 1. https://zhuanlan.zhihu.com/p/512798840?utm_source=wechat_session&utm_medium=social&s_r=0
 !!          2. https://docs.lammps.org/pair_eam.html
