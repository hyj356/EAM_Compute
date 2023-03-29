module computeUE
  use global_var, only: wp, stdout
  use fileIO, only: eamFile, Atom, region, AllocMem
  implicit none
  private
  public :: interpolate, Pair, Energy

contains

  pure function interpolate( dri, x, y, xi ) result(res)
    !! 这里因为是只需要对势函数进行插值, 而势函数的数组是确定的, 所以传入的数据有所不同
    real(wp), intent(in) ::  dri                  !< dri等于1/dr, 为了避免频繁的计算除法, 且dr可以早早从势函数文件中读取, 因此这里提前计算好
    real(wp), dimension(:), intent(in) :: x       !< 用于插值的数组, x取值
    real(wp), dimension(:), intent(in) :: y       !< 用于插值的数组, y取值
    real(wp), intent(in)               :: xi      !< 被插的点
    real(wp) :: res
    integer                            :: prev    !< 数组中位于xp之前的下标
    integer                            :: next    !< 数组中位于xp之后的下标

    !! 在实际运算的时候, 必须确保原子之间距离小于截断半径才会进行运算, 所以无需担心待插的xi超出范围
    !! 寻找传入距离对应的2个下标
    prev = int(xi * dri) + 1      !! 如果太靠近起点, 这里是0, 但是Fortran数组从1开始
    next = prev + 1               !! 所以这里要注意

    res = y(prev) + ((xi - x(prev)) * (y(next) - y(prev))) * dri    !! 进行内插

  end function interpolate

  pure function Pair(eamPot, r) result(res)
    !! 根据传入的2个原子之间的距离线性内插计算2体势
    type(eamFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: r             !< 2个原子之间的距离
    real(wp) :: EffectiveZ               !< 用于计算2体势φ的有效电荷值
    real(wp) :: res                      !< 计算对势值, 单位为ev

    !! 根据从势函数读取到的Z数组, 内插获得在r距离下的有效电荷值, 注意这里应该使用Z2R数组
    EffectiveZ = interpolate(eamPot%dri, eamPot%rx, eamPot%Z2R, r)
    res = EffectiveZ * EffectiveZ / r

  end function Pair

  pure function Prof(eampot, r) result(res)
    !! 此程序用于计算电子密度, 根据2个原子之间的距离线性内计算电子密度
    type(eamFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: r             !< 2个原子之间的距离
    real(wp) :: res

    res = interpolate(eamPot%dri, eamPot%rx, eamPot%Rho, r)

  end function Prof

  pure function Embed(eampot, rho) result(res)
    !! 根据传入的电子密度ρ(rho), 线性内插计算嵌入能F
    type(eamFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: rho           !< 中心原子与近邻原子形成的电子密度
    real(wp) :: res

    res = interpolate(eamPot%drhoi, eamPot%rhox, eamPot%F, rho)

  end function Embed

  subroutine BuildList(Atoms, cutoff2, box)
    !! 构建每个原子的邻居列表, 使用verlet方法
    type(Atom), intent(inout), allocatable :: Atoms(:) !< 存储原子坐标和邻居列表的ATOM类数组
    real(wp), intent(in) :: cutoff2                    !< 建立邻居列表的截断半径的平方
    type(region), intent(in)  :: box                   !< 仿真盒子, 包含了坐标上下限和三个维度长度的信息
    !real(wp), intent(in), dimension(3, 2) :: box      !< 存储了盒子的几个坐标极限的数组
    real(wp) :: dR                                     !< 两个原子之间的距离
    integer :: nAtoms                                  !< 原子个数
    integer :: i, j                                    !< 循环变量

    !! 初始化变量
    nAtoms = size(Atoms)

    !! 开始verlet循环构建邻居列表, 时间复杂度为N(N - 1)/2, N为原子数量
	!! 所以当原子数超过1W时不建议使用这种方法构建邻居列表, 而应该使用时间复
	!! 杂度为O(N)的Cell-list法.
    do i = 1, nAtoms - 1
      do j = i + 1, nAtoms
        dR = Distance(Atoms, i, j, box)
        if (dR < cutoff2) then
          Atoms(i)%slave = [Atoms(i)%slave, j]    !! Fortran2003特性, 可以动态修改可分配数组的元素, 但因涉及反复申请内存, 效率极低
          Atoms(j)%slave = [Atoms(j)%slave, i]    !! 当然也可以给每个原子提前分配一个足够大的比如400个邻居ID数组, 但是内存使用又会增加许多
												  !! 比较好的折中办法是使用指针和派生数据类型搭建链表
        end if
      end do
    end do

    !! 邻居列表构建完毕
  end subroutine BuildList

  pure function Distance(Atoms, i, j, box) result(res)
    !! 考虑周期性边界条件计算2个原子之间的距离
    type(Atom), intent(in), dimension(:) :: Atoms     !< 存储原子坐标和邻居列表的ATOM类数组
    integer, intent(in) :: i, j                       !< 2个原子之间的下标
    type(region), intent(in)  :: box                  !< 仿真盒子, 包含了坐标上下限和三个维度长度的信息
    real(wp) :: dx, dy, dz                            !< xyz方向上的差距, 考虑周期性边界条件
    real(wp) :: res                                   !< 考虑周期性边界条件下第i个与第j个原子之间的距离


    dx = Atoms(i)%x - Atoms(j)%x
    dx = dx - nint(dx / box%lx) * box%lx
    dy = Atoms(i)%y - Atoms(j)%y
    dy = dy - nint(dy / box%ly) * box%ly
    dz = Atoms(i)%z - Atoms(j)%z
    dz = dz - nint(dz / box%lz) * box%lz

    res = dx * dx + dy * dy + dz * dz

  end function Distance

  subroutine Energy(eampot, Atoms, power, box)
    type(eamFile), intent(in) :: eampot                 !< 一个包含了势函数所有数据和处理好数据之后的类
    type(Atom), intent(inout), allocatable :: Atoms(:)  !< 存储原子坐标和邻居列表的ATOM类数组
    real(wp), intent(out) :: power                      !< 输出的模型的能量
    type(region), intent(in)  :: box                    !< 仿真盒子, 包含了坐标上下限和三个维度长度的信息
    real(wp) :: phi     !< 两体势
    real(wp) :: rho     !< 电子密度
    !real(wp) :: emb     !< 嵌入势, 其实计算用不太上
    real(wp) :: dR      !< 两个原子之间的距离
    integer :: i, j

    !! 首先构建邻居列表
    call BuildList(Atoms, eampot%cutoff2, box)

    !! 初始化每一个变量
    phi = 0.d0; rho = 0.d0;
    power = 0.d0 ! emb = 0.d0;
    !! 第一层循环遍历每一个原子, 第二层循环遍历每一个原子的近邻
    do i = 1, size(Atoms)
      do j = 1, size(Atoms(i)%slave)
        dR = sqrt(Distance(Atoms, i, Atoms(i)%slave(j), box))
        phi = phi + pair(eampot, dR)
        rho = rho + prof(eampot, dR)
      end do
      power = power + Embed(eampot, rho)
      rho = 0.d0
    end do
    power = power + 0.5d0 * phi

  end subroutine Energy
end module
