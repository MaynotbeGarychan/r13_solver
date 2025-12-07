#
# Makefile for building lsdyna.exe
# on 'Windows X64
# SMPD
FC=ifort
FLAGS = -c -traceback  -4Yportlib -fp:precise -assume:byterecl,buffered_io  -QxSSE2 -w -WB \
-fpp -DPCWIN -DINTEL -DIA64 -DAMD64 -DWINX64 -DNET_SECURITY -DXEON64 -DLSTCODE -DLONGFORMAT 

TARGET=lsdyna.exe  
SMPS = -DOPENMP -Qopenmp -DDOUBLEXR 
SMPD = -DOPENMP -Qopenmp -DAUTODOUBLE -4R8 -4I8 

FFLAGS = $(FLAGS) $(SMPD)
LFLAGS = -F:1610612736 -link -force -nodefaultlib:msvcrt.lib  \
-ignore:4217 -ignore:4049 -nodefaultlib:vcomp.lib 
OBJS= couple2other_user.obj dyn21cnt.obj dyn21ueos.obj dyn21umats.obj  dyn21ushl.obj \
      dyn21utan.obj dyn21.obj dyn21tumat.obj  dyn21umatc.obj  dyn21umatv.obj \
      dyn21usld.obj  dynrfn_user.obj dyn21umat.obj init_dyn21.obj \
	  dyn21umats_cp_ph.obj dyn21umats_fcc_implicit.obj
UDFOBJS=udf_model_dislocation_vector.obj\
	udf_model_elasticity.obj udf_model_eshelby_tensor.obj\
	udf_model_hardening_dsl.obj udf_model_hardening.obj\
	udf_model_kinematics_configuration.obj udf_model_kinematics.obj\
	udf_model_rotation.obj\
	udf_model_schmid_tensor.obj\
	udf_model_slip_system_operation.obj udf_model_slip_system.obj\
	udf_model_slip.obj\
	udf_model_stress_update.obj\
	udf_model_orientation_init.obj\
	udf_math.obj udf_math_tranformation.obj udf_math_vector.obj udf_math_matrix.obj udf_math_array.obj\
	udf_math_distribution.obj udf_math_lusolve.obj\
	udf_strain_calculation.obj udf_stress_calculation.obj

LIBS= libdyna.lib libansys.lib shell32.lib user32.lib comctl32.lib comdlg32.lib Iphlpapi.lib \
gdi32.lib mkl_intel_ilp64.lib mkl_sequential.lib mkl_core.lib

$(TARGET): $(OBJS)
	$(FC) -w -q -o $(TARGET) $(OBJS) $(UDFOBJS) $(LIBS) $(LFLAGS)
$(OBJS) :nhisparm.inc
.F.obj:  
 	$(FC) $(FFLAGS) $*.F

clean:
	-if exist *.obj erase *.obj
#	-if exist del $(OBJS)
	-if exist $(TARGET) erase $(TARGET)

