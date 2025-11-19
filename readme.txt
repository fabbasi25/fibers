fibers.jl contains the code that computes the fibers and bases in a given 4d reflexive polytope. It requires the files h11_bases.dat and toric_bases.txt to run. It takes two inputs: the paths to the input and output files.

The input file must be in palp processed format (run poly.x -v -d -g <input_file> <output_file> on the input file from Kreuzer-Skarke). Please have the input file name without dashes. 
However, if your input is broken into multiple files, you may name your input files as "<file_name>-<file_number>.txt" where the file number ranges from 1-N. If you do this, all output polytopes will be numbered. 
For example, if you have files v05-1.txt, v05-2.txt containing 1000 polytopes and 561 polytopes respectively, then the first polytope in the output file corresponding to v05-2.txt will be numbered 1001. You will have to change polytopes_per_file to the desired number though. 

Each line in the output file corresponds to one polytope. It is of the following format:
[[v, n], [h11, h21], number_of_fibers, [[singular_flag, base_number, h11_base, list_of_missing_rays, list_of_fiber_types]]],
where
v = number of vertices in the dual polytope
n = polytope index
h11, h21 = Hodge pair of the corresponding Calabi-Yau
number_of_fibers = total number of fibration structures found (with multiplicity)
singular_flag = 1 if the base is singular
                         0 if the base is smooth
base_number = the base index in the toric bases list from 1204.0283
h11_base = h11 of the base
list_of_missing_rays = [] if the base is smooth
                                  or a list of 0s and 1s where 0 corresponds to the missing rays if base is singular
list_of_fiber_types = [[fiber_type, fiber_multiplicity]] contains two element lists corresponding to the number of each type of fiber that was found.

For example, we could have the following line in the output:
[[5, 11], [38,2], 11, [[1, 212, 7, [1, 1, 1, 1, 1, 1, 0, 0, 1], [[1, 1]]], [0, 6414, 16, [], [[1, 1]]], [0, 550, 9, [], [[5, 3]]], [0, 4154, 14, [], [[10, 6]]]]],
Then, this is the 11th polytope in the list. The number of vertices in the dual polytope are 5. The Hodge numbers are 38, 2. There were 11 fibers found in total over three bases:

We have a singular base 212 with h11 = 7 with two rays missing. It has fiber 1 on it with multiplicity = 1.
We have base 6414 with h11 = 16. It also has a fiber 1 with multiplicity = 1.
We have base 4154 with h11 = 14. It has fiber 10 with multiplicity = 6.
