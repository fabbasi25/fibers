# Classifying Fibers and Bases in Toric Hypersurface Calabi-Yau Threefolds

This project implements algorithms to compute **elliptic fibrations of toric hypersurface Calabi–Yau threefolds** in the Kreuzer-Skarke database. It classifies possible fibers and their bases, and outputs the complete fibration structure of each polytope.

The code is written in **Julia** and is designed to efficiently handle large datasets from the Kreuzer-Skarke database.

---

## Features

- Computes fibers and bases for a given 4-dimensional reflexive polytope.  
- Handles singular and smooth bases, with identification of missing rays in singular cases.  
- Outputs fibration data and fiber multiplicities.  
- Can output basis transformations/implement basis transformations to put polytopes in `[fiber; base]` form.  

---

## Technical Overview

- **fibers.jl** – Main computation of fibers and bases.  
- **basis_change_matrix.jl** – Returns transformation matrices to convert polytopes to `[fiber; base]` basis.  
- **transformed_polytope.jl** – Produces all points in the polytope in the new basis.

**Input:**  
- 4D reflexive polytope files in PALP-processed format 
- `h11_bases.dat` and `toric_bases.txt`  
- Supports multiple input files with automatic polytope numbering.

**Output:**  
Each line corresponds to one polytope and contains:  
- Number of vertices  
- Polytope index  
- Hodge numbers  
- Number of fibers  
- Base details (singular flag, base index, Hodge numbers, missing rays)  
- Fiber types and multiplicities  

*(More details of output are included in the source files.)*

---

## Installation / Usage

1. Install **Julia**.  
2. Place the input files (`h11_bases.dat`, `toric_bases.txt`) in the project directory.  
3. Run the script from the command line by passing the input and output file paths as arguments:

```bash
julia fibers.jl path/to/input_file.txt path/to/output_file.txt
