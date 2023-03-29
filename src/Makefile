FC = gfortran						# 编译器名字
TARGET = calculate			# 可执行程序名称
FFLAGS = -c -O2 -Wall -Wno-unused-function 		# 编译选项
OBJ = global_var.o mod_fileIO.o mod_computePE.o main.o	# 二进制中间文件名称, 其实就是源代码文件名称后缀改成.o
# FILES = global_var.f90 mod_fileIO.f90 mod_computePE.f90 main.f90
# SRC = $(wildcard *.f90)									# 所有的以.f90为后缀的文件
# OBJ = ${patsubst, %.f90, %.o, ${SRC}}		# 将SRC中所有以.f90结尾的字符串改成以.o结尾的
# TESTS = global_var.f90 mod_fileIO.f90 mod_computePE.f90 mod_unittest.f90 test.f90

${TARGET}: ${OBJ}
	${FC} -o $@  $^

%.o: %.f90
	${FC} ${FFLAGS} $< -o $@

.PHONY: clean

clean:
	rm -rf *.mod *.o ${TARGET}