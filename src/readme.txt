

How to make executable for User Defined Materials

	
        1.Open a Dos command prompt Windows and set compilers
          environment variables properly.

Start->Programs->Intel Parallel Studio XE 2019->Command Prompt
->Parallel Studio XE with Intel Compiler XE 2019->
Intel 64 Visual Studio 2019 mode

        2.Make the executable by default package.

Go to the working directory where you unzip User Defined
Materials package. Type nmake.exe under command line prompt.
See if you can generate the executable by default package.

        3.Modify user defined material subroutine

Modify related Fortran file e.g: dyn21.F with notepad or other
edit utility. Copy or insert your own code into dyn21.F.
	
	4.Make the executable

Type "nmake" and press Enter to Compile and link.
Ignore the warning messages caused by multiple-defined 
subroutines in the LIB files.

Compiler and version:

	1).Intel Fortran:

Intel Parallel Studio XE 2019

	2).Microsoft Visual C++:

Microsoft Visual C++ 2019 x64 cross tools

        5. Run the executable

The executable can be ran without an issue under building
Envronment. If you move it out to other system. You may 
get libiomp5md.dll was not found error. You can find 
libiomp5md.dll from Intel Fortran installation folder, 
then copy it together and put in the same folder with 
the executable.



 

