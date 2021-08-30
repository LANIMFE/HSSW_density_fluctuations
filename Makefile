#FFTW library
fftwdir = usr/local/lib
#Fortran compiler
FC  = gfortran
#Fortran Flags
FFlagsO  = -c -J./source_code/objects -I./source_code/objects -ffast-math -o3
FFlagsMain  = -J./source_code/objects -I./source_code/objects -ffast-math -o3
#Program Output:
Output = -o program.out
#Current directory
current_dir = $(shell pwd)
#DIRECTORIES:
#Main Program Directory:
main_program_DIR = ./source_code
#Object directory:
obj_dir = ./source_code/objects
#Numerical recipes and math code directory :
math_dir = ./source_code/math
#Structures directory:
structures_dir = ./source_code/structures
#Dynamics directory
Dynamics_dir = ./source_code/dynamics
#Systems directories:
##Hard_Sphere##
HS_dir = $(structures_dir)/Hard_Sphere
##Hard_Sphere_Square_Well##
HSSW_dir = $(structures_dir)/Hard_Sphere_Square_Well
##Hard_Sphere_Attractive_Yukawa##
HSAY_dir = $(structures_dir)/Hard_Sphere_Attractive_Yukawa
##Hard Disk##
HD_dir = $(structures_dir)/Hard_Disk
##Hard_Disk_Square_Well##
HDSW_dir = $(structures_dir)/Hard_Disk_Square_Well
##Hard_Disk_Attractive_Yukawa##
HDAY_dir = $(structures_dir)/Hard_Disk_Attractive_Yukawa
## Structures Files ##
Structures_files = structure_function_selector general_structural_functions
## Systems Files ##
System_files = hs_structure hssw_structure hsay_structure hd_structure hday_structure hdsw_structure
## Math Files and numerical recipes  ##
Math_files = quadratures numerical_recipes nrtype nrutil
## Dynamic Files ##
Dynamics_files = scgle nescgle heterogeneities_module
###########################
#Main Program Compilation#
###########################
#program compilation
program_dep = $(foreach System_files,$(System_files),$(obj_dir)/$(System_files).o)\
$(foreach Structures_files,$(Structures_files),$(obj_dir)/$(Structures_files).o)\
$(foreach Dynamics_files,$(Dynamics_files),$(obj_dir)/$(Dynamics_files).o)\
$(foreach Math_files,$(Math_files),$(obj_dir)/$(Math_files).o)

program.out: $(program_dep) $(main_program_DIR)/program.f03
	$(FC) $(FFlagsMain) $(program_dep) $(main_program_DIR)/program.f03 -o program.out

###############################
#Structures Files Compilation#
###############################
#Structure Function Selector#
structure_function_selector_dep = $(obj_dir)/hs_structure.o $(obj_dir)/hssw_structure.o $(obj_dir)/quadratures.o
$(obj_dir)/structure_function_selector.o: $(structures_dir)/structure_function_selector.f03 $(structure_function_selector_dep)
	$(FC) $(FFlagsO) $(structures_dir)/structure_function_selector.f03 -o $(obj_dir)/structure_function_selector.o
#General structural functions#
general_structural_functions_dep = $(obj_dir)/nrtype.o $(obj_dir)/hssw_structure.o $(obj_dir)/quadratures.o
$(obj_dir)/general_structural_functions.o: $(structures_dir)/general_structural_functions.f03 $(general_structural_functions_dep)
	$(FC) $(FFlagsO) $(structures_dir)/general_structural_functions.f03 -o $(obj_dir)/general_structural_functions.o
#Hard Sphere:
hs_dep = $(obj_dir)/nrtype.o
$(obj_dir)/hs_structure.o: $(HS_dir)/hs_structure.f03 $(hs_dep)
	$(FC) $(FFlagsO) $(HS_dir)/hs_structure.f03 -o $(obj_dir)/hs_structure.o
#Hard Disk:
hd_dep = $(obj_dir)/nrtype.o $(obj_dir)/quadratures.o
$(obj_dir)/hd_structure.o: $(HD_dir)/hd_structure.f03 $(hd_dep)
	$(FC) $(FFlagsO) $(HD_dir)/hd_structure.f03 -o $(obj_dir)/hd_structure.o
#Hard Sphere Square Well:
hssw_dep = $(obj_dir)/nrtype.o $(obj_dir)/hs_structure.o $(obj_dir)/quadratures.o
$(obj_dir)/hssw_structure.o: $(HSSW_dir)/hssw_structure.f03 $(hssw_dep)
	$(FC) $(FFlagsO) $(HSSW_dir)/hssw_structure.f03 -o $(obj_dir)/hssw_structure.o
#Hard Sphere Attractive Yukawa:
hsay_dep = $(obj_dir)/nrtype.o $(obj_dir)/hs_structure.o
$(obj_dir)/hsay_structure.o: $(HSAY_dir)/hsay_structure.f03 $(hsay_dep)
	$(FC) $(FFlagsO) $(HSAY_dir)/hsay_structure.f03 -o $(obj_dir)/hsay_structure.o
#Hard Disk Square Well:
hdsw_dep = $(obj_dir)/nrtype.o $(obj_dir)/hd_structure.o $(obj_dir)/quadratures.o
$(obj_dir)/hdsw_structure.o: $(HDSW_dir)/hdsw_structure.f03 $(hdsw_dep)
	$(FC) $(FFlagsO) $(HDSW_dir)/hdsw_structure.f03 -o $(obj_dir)/hdsw_structure.o
#Hard Disk Attractive Yukawa:
hday_dep = $(obj_dir)/nrtype.o $(obj_dir)/hd_structure.o
$(obj_dir)/hday_structure.o: $(HDAY_dir)/hday_structure.f03 $(hday_dep)
	$(FC) $(FFlagsO) $(HDAY_dir)/hday_structure.f03 -o $(obj_dir)/hday_structure.o

###############################
#Dynamics Files Compilation#
###############################
#SCGLE:
scgle_dep = $(obj_dir)/nrtype.o
$(obj_dir)/scgle.o: $(Dynamics_dir)/scgle.f03 $(scgle_dep)
	$(FC) $(FFlagsO) $(Dynamics_dir)/scgle.f03 -o $(obj_dir)/scgle.o
#NESCGLE:
nescgle_dep = $(obj_dir)/nrtype.o $(obj_dir)/scgle.o
$(obj_dir)/nescgle.o: $(Dynamics_dir)/nescgle.f03 $(nescgle_dep)
	$(FC) $(FFlagsO) $(Dynamics_dir)/nescgle.f03 -o $(obj_dir)/nescgle.o
#Heterogeneities
heterogeneities_module_dep= $(obj_dir)/scgle.o $(obj_dir)/numerical_recipes.o $(obj_dir)/nrtype.o $(obj_dir)/structure_function_selector.o
$(obj_dir)/heterogeneities_module.o: $(Dynamics_dir)/heterogeneities_module.f03 $(heterogeneities_module_dep)
	$(FC) $(FFlagsO) $(Dynamics_dir)/heterogeneities_module.f03 -o $(obj_dir)/heterogeneities_module.o
###############################
#Math Files Compilation#
###############################
#Quadratures:
$(obj_dir)/quadratures.o: $(math_dir)/quadratures.f03 $(obj_dir)/nrtype.o
	$(FC) $(FFlagsO) $(math_dir)/quadratures.f03 -o $(obj_dir)/quadratures.o
###############################
#Numerical Recipes compilation#
###############################
#numerical recipes
numerical_recipes_dep= $(obj_dir)/nrtype.o $(obj_dir)/nrutil.o
$(obj_dir)/numerical_recipes.o: $(math_dir)/numerical_recipes.f03 $(numerical_recipes_dep)
	$(FC) $(FFlagsO) $(numerical_recipes_dep) $(math_dir)/numerical_recipes.f03 -o $(obj_dir)/numerical_recipes.o
##nrutil
nrutil_dep= $(obj_dir)/nrtype.o
$(obj_dir)/nrutil.o: $(math_dir)/nrutil.f03 $(nrutil_dep)
	$(FC) $(FFlagsO) $(math_dir)/nrutil.f03 -o $(obj_dir)/nrutil.o
##nrtype
$(obj_dir)/nrtype.o: $(math_dir)/nrtype.f03
	$(FC) $(FFlagsO) $(math_dir)/nrtype.f03 -o $(obj_dir)/nrtype.o

#Special Functions:
#$(obj_dir)/fn.o: $(math_dir)/fn.f03
#	$(FC) $(FFlagsO) $(math_dir)/fn.f03 -o $(obj_dir)/fn.o
clean:
	rm ./source_code/objects/*.o
	rm ./source_code/objects/*.mod
	rm *.out
