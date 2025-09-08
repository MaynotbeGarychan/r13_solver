# makefile to build my umat
SOLVER_DIR=C:\LSDYNA\program
# build user define function and model
build_udf:
	cd ./udf && $(MAKE)

# build the solver
build_solver:
	cd ./src && $(MAKE) TARGET=lsdyna_cp.exe
send_solver:
	scp ./src/lsdyna_cp.exe $(SOLVER_DIR)
clean_solver:
	cd ./src && $(MAKE) clean TARGET=lsdyna_cp.exe

# build the solver for test
build_solver_test:
	cd ./src && $(MAKE) TARGET=lsdyna_cp_test.exe
send_solver_test:
	scp ./src/lsdyna_cp_test.exe $(SOLVER_DIR)
clean_solver_test:
	cd ./src && $(MAKE) clean TARGET=lsdyna_cp_test.exe

# FCC model: International Journal of Plasticity, J.Rossiter, 2010
# A new crystal plasticity scheme for explicit time integration codes
# to simulate deformation in 3D microstructures: Effects of strain path,
# strain rate and thermal softening on localized deformation in the aluminum
# alloy 5754 during simple shear, International Journal of Plasticity
build_fcc:
	cd ./src && $(MAKE) UDFMATOBJ=dyn21umats_fcc.obj TARGET=lsdyna_fcc.exe
send_fcc:
	scp ./src/lsdyna_fcc.exe $(SOLVER_DIR)
clean_fcc:
	cd ./src && $(MAKE) clean TARGET=lsdyna_fcc.exe
# BCC model: following the same appraoch as FCC
build_bcc:
	cd ./src && $(MAKE) UDFMATOBJ=dyn21umats_bcc.obj TARGET=lsdyna_bcc.exe
send_bcc:
	scp ./src/lsdyna_bcc.exe $(SOLVER_DIR)
clean_bcc:
	cd ./src && $(MAKE) clean TARGET=lsdyna_bcc.exe

# Dislocation based FCC: International Journal of Plasticity, M.G. Lee, 2010
# A dislocation density-based single crystal constitutive equation
build_fcc_dsl:
	cd ./src && $(MAKE) UDFMATOBJ=dyn21umats_fcc_dsl.obj TARGET=lsdyna_fcc_dsl.exe
send_fcc_dsl:
	scp ./src/lsdyna_fcc_dsl.exe $(SOLVER_DIR)
clean_fcc_dsl:
	cd ./src && $(MAKE) clean TARGET=lsdyna_fcc_dsl.exe

# clean function
#clean:
#	cd ./src && $(MAKE) clean
