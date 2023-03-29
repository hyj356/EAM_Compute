module fileIO
  use iso_fortran_env, only: si => int32, li => int64
  use global_var, only: wp, stdout, Bohr, Hartree
  implicit none
  private
  public :: eamFile, Readfile, ReadData, Atom, region, AllocMem

  type eamFile
    !! 一个用于存储eam类势函数文件中所有有效信息的类, 使用type进行管理,方便调用
    character(len=:), allocatable :: filename   !< 势函数文件名称
    real(wp), allocatable :: F(:)               !< 多体势, 也称为嵌入势的数组
    real(wp), allocatable :: Z(:)               !< 有效电荷数组
    real(wp), allocatable ::  Rho(:)            !< 电子密度, 与距离有关, 用于计算多体势
    real(wp), allocatable :: Z2R(:)             !< 处理之后的数组Z, 避免计算势能时多次计算乘法
    real(wp), allocatable :: rx(:)              !< 原子距离在截断半径上的等差数组
    real(wp), allocatable :: rhox(:)            !< 电荷密度的等差数组
    real(wp) :: dr                 !< r的取点间隔
    real(wp) :: dri                !< r的取点间隔的导数, 避免反复的除法运算
    real(wp) :: drho               !< 电子密度rho的取点间隔
    real(wp) :: drhoi              !< 电子密度rho的取点间隔
    real(wp) :: cutoff             !< eam势函数的截断半径
    real(wp) :: cutoff2            !< eam势函数的截断半径的平方
    real(wp) :: Amass              !< 原子质量
    real(wp) :: lattice            !< 晶格常数
    integer  :: Nr                 !< r的取点个数
    integer  :: Nrho               !< 电子密度rho的取点间隔
    integer  :: AtomicNumber       !< 原子序数
    character(len=3) :: latType    !< 晶格类型, 对于金属来说, 一般只有FCC, BCC, HCP三种晶格类型
  end type

  type Atom
    integer :: typeId   !< 原子的类型ID
    real(wp) :: x       !< 原子的x坐标
    real(wp) :: y       !< 原子的y坐标
    real(wp) :: z       !< 原子的z坐标
    integer, allocatable :: slave(:)    !< 原子的所有在截断半径以内的邻居
  contains
    procedure, public, pass(self) :: InitialList
  end type

  type region
    real(wp) :: xlo, xhi
    real(wp) :: ylo, yhi
    real(wp) :: zlo, zhi
    real(wp) :: lx, ly, lz
  end type

  interface AllocMem
    module procedure :: AllocMem_1d_dp, AllocMem_2d_dp, AllocMem_3d_dp
    module procedure :: AllocMem_1d_si, AllocMem_2d_si, AllocMem_3d_si
    module procedure :: AllocMem_1d_li, AllocMem_2d_li, AllocMem_3d_li
    module procedure :: AllocMem_1d_Atom
  end interface

contains
  subroutine Readfile(filename, eamPotential)
    !! 读取eam势函数文件到数组中
    character(len=*), intent(in) :: filename            !< 势函数文件名称
    type(eamFile), intent(out) :: eamPotential          !< 势函数文件的数据集合
    integer :: fileId             !< 势函数文件对应的通道ID
    integer :: io_Stat            !< 检测文件读取是否出错的整数
    logical :: isExist            !< 检测文件是否存在的逻辑变量
    integer :: i                  !< 循环变量

    !! 查询文件是否存在, 如果不存在就报错并终止程序
    inquire(file=filename, exist=isExist)
    if (.not. isExist) then
      write(stdout, '(A, A)') filename, " can't be found in this path."
      stop
    end if
    eamPotential%filename = filename
    open(newunit=fileId, file=filename, action='read')

    !! 空读取跳过前1行
    read(fileId, *, iostat=io_Stat)
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading 1st line of File.'
    end if

    !! 第二行其实有有用的信息, 分别是原子序数, 原子质量, 晶格常数, 晶格类型, 但是对于MC迭代来说没有用处
    !! 但是这里还是把相关程序写上以便未来调用
    read(fileId, *) eamPotential%AtomicNumber, eamPotential%Amass, eamPotential%lattice, eamPotential%latType
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading 2nd line of File.'
    end if

    !! 第三行包含了eam势函数文件中非常重要的5个参数, 具体含义在变量声明的时候已经说明
    read(fileId, *) eamPotential%Nrho, eamPotential%drho, eamPotential%Nr, eamPotential%dr, eamPotential%cutoff
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading 3rd line of File.'
    end if

    !! 根据读取到的变量分配内存
    call AllocMem(eamPotential%F, eamPotential%Nrho)
    call AllocMem(eamPotential%Z, eamPotential%Nr)
    call AllocMem(eamPotential%Rho, eamPotential%Nr)
    call AllocMem(eamPotential%Z2R, eamPotential%Nr)

    !! 之后从文件中读取电子密度, 电荷密度, 和多体势能数组, 并对其中一些数据进行预处理减少计算量
    !! 如2体势φ的计算: r · φ = 27.2 · 0.529 · Zi · Zj, 中Hartree = 27.2d0, Bohr = 0.529d0
    !! 对于单元素的eam]势函数来说, Zi = Zj.
    read(fileId, *) (eamPotential%F(i), i = 1, eamPotential%Nrho)
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading array of F.'
    end if
    read(fileId, *) (eamPotential%Z(i), i = 1, eamPotential%Nr)
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading array of Z.'
    end if
    read(fileId, *) (eamPotential%Rho(i), i = 1, eamPotential%Nr)
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading array of Rho.'
    end if

    !! 关闭文件
    close(fileId)

    !! 对读取到的数据进行初步处理, 避免后续的反复乘法或者除法运算
    eamPotential%Z2R =  eamPotential%Z * Bohr * Hartree
    call InitialArray(eamPotential%Rhox, eamPotential%rx, &
      eamPotential%dRho, eamPotential%dR, &
      eamPotential%Nrho, eamPotential%Nr)
    eamPotential%dri = 1.d0 / eamPotential%dr
    eamPotential%drhoi = 1.d0 / eamPotential%drho
    eamPotential%cutoff2 = eamPotential%cutoff * eamPotential%cutoff

  end subroutine Readfile

  subroutine InitialArray(Rhox, rx, dRho, dR, Nrho, Nr)
    !! 初始化对于Rho和根据cutoff和dR计算的数学公式
    real(wp), allocatable, dimension(:), intent(out) :: Rhox, rx  !< 分别对应电子密度的均分和原子距离的均分数组
    real(wp), intent(in) ::dRho, dR         !< 上述数组相邻元素之间的步长
    integer, intent(in) :: Nrho, Nr         !< 数组的大小
    integer :: i                            !< 循环变量

    call AllocMem(Rhox, Nrho)
    call AllocMem(rx, Nr)

    do concurrent(i = 0 : Nrho - 1)
      Rhox(i + 1) = i * dRho
    end do

    do concurrent(i = 0 : Nrho - 1)
      Rx(i + 1) = i * dR
    end do

  end subroutine InitialArray

  subroutine ReadData(filename, box, Atoms)
    character(len=*), intent(in) :: filename            !< 模型文件的名称
    type (region),intent(out) :: box        !< 记录仿真盒子尺寸的矩阵
    type(Atom), allocatable, intent(out) :: Atoms(:)    !< 存储原子坐标和邻居列表的ATOM类数组
    integer :: fileId             !< 模型文件对应的通道ID
    !integer :: io_Stat            !< 检测文件读取是否出错的整数
    integer :: i                  !< 循环变量
    logical :: isExist            !< 检测文件是否存在的逻辑变量
    integer :: nTypes             !< 模型文件中有多少种元素
    integer :: nAtoms             !< 模型中有多少个原子
    integer :: tempi              !< 用来跳过第一列原子ID的数据

    inquire(file=filename, exist=isExist)
    if (.not. isExist) then
      write(stdout, '(A, A)') filename, " can't be found in this path."
      stop
    end if
    open(newunit=fileId, file=filename, action='read')

    !! 跳过前2行注释行
    read(fileId, *)
    read(fileId, *)

    !! 读取有多少个原子和多少种元素, 并分配内存
    read(fileId, *) nAtoms
    read(fileId, *) nTypes
    read(fileId, *)
	if (nTypes /= 1) STOP "EAM potential File only have one element, but there are two or more element in atom data file!"
    allocate(Atoms(nAtoms))
    call Atoms%InitialList()    !! 初始化近邻列表的内存

    !! 读取盒子的几个极限
    read(fileId, *) box%xlo, box%xhi
    read(fileId, *) box%ylo, box%yhi
    read(fileId, *) box%zlo, box%zhi
    box%lx = box%xhi - box%xlo
    box%ly = box%yhi - box%ylo
    box%lz = box%zhi - box%zlo

    !! 连续跳过7行, 原子质量不会被读取因为在计算势能的时候不会用到
    do i = 1, 7
      read(fileId, *)
    end do

    !! 读取原子坐标
    do i = 1, nAtoms
      read(fileId, *) tempi, Atoms(i)%typeId, Atoms(i)%x, Atoms(i)%y, Atoms(i)%z
    end do

    close(fileId) !! 关闭文件

  end subroutine ReadData

  pure elemental subroutine InitialList(self)
    !! 初始化所有的slaveList为0数组
    class(Atom), intent(inout) :: self

    allocate(self%slave(0))
  end subroutine InitialList

  subroutine AllocMem_1d_dp(Array, dim)
    !! 对传入的双精度数组进行分配内存, 如果报错可以给出报错信息
    real(wp), intent(inout), allocatable :: Array(:)  !< 传入的双精度数组
    integer, intent(in) :: dim                        !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(dim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_1d_dp

  subroutine AllocMem_2d_dp(Array, xdim, ydim)
    !! 对传入的二维双精度数组进行分配内存, 如果报错可以给出报错信息
    real(wp), intent(inout), allocatable :: Array(:, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim                    !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_2d_dp

  subroutine AllocMem_3d_dp(Array, xdim, ydim, zdim)
    !! 对传入的三维双精度数组进行分配内存, 如果报错可以给出报错信息
    real(wp), intent(inout), allocatable :: Array(:, :, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim, zdim                 !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim, zdim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_3d_dp

  subroutine AllocMem_1d_si(Array, dim)
    !! 对传入的双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(si), intent(inout), allocatable :: Array(:)  !< 传入的双精度数组
    integer, intent(in) :: dim                        !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(dim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_1d_si

  subroutine AllocMem_2d_si(Array, xdim, ydim)
    !! 对传入的二维双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(si), intent(inout), allocatable :: Array(:, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim                    !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_2d_si

  subroutine AllocMem_3d_si(Array, xdim, ydim, zdim)
    !! 对传入的三维双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(si), intent(inout), allocatable :: Array(:, :, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim, zdim                 !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim, zdim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_3d_si

  subroutine AllocMem_1d_li(Array, dim)
    !! 对传入的双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(li), intent(inout), allocatable :: Array(:)  !< 传入的双精度数组
    integer, intent(in) :: dim                        !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(dim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_1d_li

  subroutine AllocMem_2d_li(Array, xdim, ydim)
    !! 对传入���二维双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(li), intent(inout), allocatable :: Array(:, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim                    !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_2d_li

  subroutine AllocMem_3d_li(Array, xdim, ydim, zdim)
    !! 对传入的三维双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(li), intent(inout), allocatable :: Array(:, :, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim, zdim                 !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim, zdim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_3d_li

  subroutine AllocMem_1d_Atom(Array, dim)
    !! 对传入的Atom类数组进行分配内存, 如果报错可以给出报错信息
    type(Atom), intent(inout), allocatable :: Array(:)  !< 传入的双精度数组
    integer, intent(in) :: dim                        !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(dim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if
  end subroutine AllocMem_1d_Atom
end module fileIO
