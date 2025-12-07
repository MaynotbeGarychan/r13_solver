# prerequiste
#FC=gfortran
#FFLAGS=-g
#CLEAN=rm

FC=ifort
FFLAGS= 
CLEAN=erase
# object lists
OBJS=udf_model_dislocation_vector.obj\
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
	udf_math_distribution.obj\
	udf_strain_calculation.obj udf_stress_calculation.obj
user:
# model
#	cd build && $(CLEAN) *.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_dislocation_vector.f -o ./build/udf_model_dislocation_vector.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_elasticity.f -o ./build/udf_model_elasticity.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_eshelby_tensor.f -o ./build/udf_model_eshelby_tensor.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_hardening_dsl.f -o ./build/udf_model_hardening_dsl.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_hardening.f -o ./build/udf_model_hardening.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_kinematics_configuration.f -o ./build/udf_model_kinematics_configuration.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_kinematics.f -o ./build/udf_model_kinematics.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_rotation.f -o ./build/udf_model_rotation.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_schmid_tensor.f -o ./build/udf_model_schmid_tensor.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_slip_system_operation.f -o ./build/udf_model_slip_system_operation.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_slip_system.f -o ./build/udf_model_slip_system.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_slip.f -o ./build/udf_model_slip.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_stress_update.f -o ./build/udf_model_stress_update.obj
	$(FC) $(FFLAGS) -c ./model/udf_model_orientation_init.f -o ./build/udf_model_orientation_init.obj
# tool
	$(FC) $(FFLAGS) -c ./tool/udf_math.f -o ./build/udf_math.obj
	$(FC) $(FFLAGS) -c ./tool/udf_math_tranformation.f -o ./build/udf_math_tranformation.obj
	$(FC) $(FFLAGS) -c ./tool/udf_math_vector.f -o ./build/udf_math_vector.obj
	$(FC) $(FFLAGS) -c ./tool/udf_math_matrix.f -o ./build/udf_math_matrix.obj
	$(FC) $(FFLAGS) -c ./tool/udf_math_array.f -o ./build/udf_math_array.obj
	$(FC) $(FFLAGS) -c ./tool/udf_math_distribution.f -o ./build/udf_math_distribution.obj
	$(FC) $(FFLAGS) -c ./tool/udf_strain_calculation.f -o ./build/udf_strain_calculation.obj
	$(FC) $(FFLAGS) -c ./tool/udf_stress_calculation.f -o ./build/udf_stress_calculation.obj
	$(FC) $(FFLAGS) -c ./tool/udf_math_lusolve.f -o ./build/udf_math_lusolve.obj
# copy
	scp ./build/*.obj ../src

debug:
	$(FC) $(FFLAGS) -c ./test/udf_main.f -o ./build/udf_main.obj
	cd build && $(FC) $(FFLAGS) -o main.exe udf_main.obj $(OBJS)

clean:
	cd build && $(CLEAN) *.obj
