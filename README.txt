PROJECT: MINIGL

I. COMPILING:

Use 'make' command to compile both of the executables: minigl, minigl-nogl
If opengl is not available, use 'make minigl-nogl' to compile only the version without opengl viewer.

II. RUNNING:

Two different version of the program is available:
  minigl-nogl: This is the program that the grading script uses, this generates a output.png and diff.png files if the solution file is available
  minigl: Same as minigl-nogl, but it displays the rendered image (and the solution next to it, if available) in an opengl context.

  To run a test:
    ./minigl-nogl -i <test_filename>  #OR
    ./minigl -i <test_filename>
  e.g.
    ./minigl -i ./tests/00.txt