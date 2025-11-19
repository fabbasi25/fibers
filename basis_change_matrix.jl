"""
Accompanied by Abbasi, Nally, Taylor --> 2511.10601 (Classifying Fibers and Bases in Toric Hypersurface Calabi-Yau Threefolds)

This code computes the BASIS CHANGE MATRIX corresponding to elliptic fibrations of Calabi-Yau threefolds. It takes as input a 4d lattice polytope and returns the
transformation matrix that transforms the polytope to a 4d polytope where the first two coordinates represent the fiber and the last two coordinates represent the
base. 

This requires the files h11_bases.dat and toric_bases.txt to run. It takes two inputs: the paths to the input and output files.

The input file must be in palp processed format (run poly.x -v -d -g <input_file> <output_file> on the input file from Kreuzer-
Skarke). Please have the input file name without dashes. However, if your input is broken into multiple files, you may name yo-
ur input files as "<file_name>-<file_number>.txt" where the file number ranges from 1-N. If you do this, all output polytopes 
will be numbered. For example, if you have files v05-1.txt, v05-2.txt containing 1000 polytopes and 561 polytopes respectively, 
then the first polytope in the output file corresponding to v05-2.txt will be numbered 1001. You will have to change 
polytopes_per_file to the desired number though. 

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
list_of_fiber_types = [[fiber_type, fiber_multiplicity, list_of_transformations]] contains three element lists corresponding to the number of each type 
                        of fiber that was found. 
list_of_transformations = [basis_change_matrix] contains the basis change matrices corresponding to a particular fiber_type over a base_number 
"""

function get(f::IOStream, points::Matrix{Int64}, vertices::Matrix{Int64})::Tuple{Int64, Int64, String, Bool}
    """
    Reads the palp processed file and stores information about points of the polytope and vertices of the dual polytope 
    """
    a = false
    h = 0
    while a == false
        s = readline(f)	
        if s == ""
            return 0, 0, "0", false
        end
        a = occursin("M", s) && occursin("H", s)

        h_ind = collect(findfirst("H", s))[1]
        end_ind = collect(findfirst("[", s))[1]

        h = s[h_ind+2:end_ind-2]
    end
    s = readline(f)
    v = parse(Int, split(s, " ")[2])
    #deal with annoying palp bug
    if v == 4
        v = parse(Int, split(s, " ")[1])
        for i = 1:v
            s = readline(f)
            vertices[:, i] = map(x -> parse(Int, x),
                                  split(s, r" +")[1:4])
        end
    else
        for i = 1:4
            s = readline(f)
            vertices[i,1:v] = map(x -> parse(Int, x),
                                  split(s, r" +")[2:v+1])
        end
    end
    s = readline(f)
    p = parse(Int, split(s, " ")[2])
    #deal with annoying palp bug (this one first for 9)
    if p == 4
        p = parse(Int, split(s, " ")[1])
        for i = 1:p
            s = readline(f)
            points[:, i] = map(x -> parse(Int, x),
                                  split(s, r" +")[1:4])
        end
    else
        for i = 1:4
            s = readline(f)
            points[i,1:p] = map(x -> parse(Int, x),
                                  split(s, r" +")[2:p+1])
        end
    end
    return p, v, h, true
end
    
function read_vector_of_vectors(filename::String)::Vector{Vector{Int64}}
    """
    Reads each line (which is of the form {-1, -3, -2}) from file and stores it in a vector. 
    """
    vecs = Vector{Vector{Int}}()  # Initialize empty vector of vectors
    buffer = ""  # Temporary storage for handling multiline entries
    
    open(filename, "r") do file
        for line in eachline(file)
            buffer *= line  # Append current line to buffer
            
            if occursin("}", line)  # Check if the line completes an entry
                clean_line = replace(buffer, "{" => "", "}" => "")  # Remove braces
                numbers = parse.(Int, split(clean_line, ","))
                push!(vecs, numbers)  # Append to list
                buffer = ""  # Reset buffer for the next entry
            end
        end
    end
    
    return vecs
end

using LinearAlgebra

function member(v::Vector{Int64}, list::Matrix{Int64}, size::Int64, multiple = 1)::Bool
    """ 
    Checks if a vector lives in a list
    (should be some more efficient implementation)
    """
    for i = 1:size
        if v == multiple*list[:, i]
            return true
        end
    end
    return false
end

function findshort(points::Matrix{Int64}, vertices::Matrix{Int64}, p::Int, v::Int, short::Array{Int64, 3}, s::Vector{Int64})
    """
    Identify points with small inner product with all vertices
    Adds points which have maximum dot product i with vertices of the dual polytope to short[:,:,i]
    s[i] = # of points with maximum dot product i with vertices of the dual  
    """
    s[:] = [0, 0, 0]
    for i = 1:p
        m = maximum([dot(points[:, i], vertices[:, j]) for j in 1:v])
        if m<4 && m != 0
            s[m] += 1
            short[:, s[m], m] = points[:, i]
        end
    end
end

function checkfibertype(s::Vector{Int64}, points::Int64)::Int64
    """
    Uses the cardinality of short lists + # of points to determine fiber type 
    """
    ftypes = Dict([0, 0, 0] => 0, [6, 3, 0] => 16, [8, 0, 0] => 15, [6, 2, 0] => 14, [4, 2, 2] => 13, [6, 1, 0] => 12, [4, 2, 1] => 11, 
    [6, 0, 0] => 7, [4, 1, 1] => 8, [4, 2, 0] => 9, [2, 2, 1] => 10, [2, 2, 1] => 6, [4, 1, 0] => 5, [2, 0, 2] => 4, [2, 2, 0] => 3, 
    [4, 0, 0] => 2, [0, 3, 0] => 1) 

    ftype = ftypes[s]

    if ftype == 6 || ftype == 10
        if points == 7 
            ftype = 10
        end
        if points == 6 
            ftype =  6
        end
    end
    return ftype
end

function getfibersubpolytope(points::Matrix{Int64}, vertices::Matrix{Int64}, p::Int64, v::Int64, x::Vector{Int64}, y::Vector{Int64})::Tuple{Int64, Int64, Matrix{Int64}}
    """
    Looks for subpolytope spanned by x and y
    Returns the fibertype, the number of points in the subpolytope and the full subpolytope 
    """
    subpoly = zeros(Int, 4, p)
    # sp = size of subpoly
    sp = 0
    # get hyperplane spanned by x and y 

    for i = 1:p
        z = points[:, i]

        dot_xy = dot(y, x)

        b = [dot(z, x), dot(z, y)]
        A = [dot(x, x) dot_xy; dot_xy dot(y, y)]
        try
            coeffs = inv(A)*b
            error = sum(abs.(coeffs[1]*x + coeffs[2]*y - z))
            in_plane = false

            if abs(error) < 1.0e-7
                in_plane = true 
            end 
            
            if in_plane
                # add to subpolytope 
                sp += 1
                subpoly[:, sp] = z[:]
            end 
        catch
            # in case matrix is singular
            println("Something has gone terribly wrong. Is your polytope reflexive?")
        end 
    end

    # find short for each point in the subpolytope contained in the plane spanned by x and y
    short_subpoly = zeros(Int, 4, sp, 3)
    s = zeros(Int, 3)
    findshort(subpoly, vertices, sp, v, short_subpoly, s)

    # compare to get fiber type
    return checkfibertype(s, sp), sp, subpoly
end

function bezcoeffs(a::Int64, b::Int64)::Tuple{Int64, Int64, Int64}
    """     
    Runs extended Euclidean/Bezout's lemma to calculate gcd(a, b), x, y such that a*x + b*y = gcd(a, b)
    """
    oldr, r = a, b
    olds, s = 1, 0
    oldt, t = 0, 1

    while r != 0
        q = div(oldr,r)  
        oldr, r = r, oldr - q*r
        olds, s = s, olds - q*s
        oldt, t = t, oldt - q*t
    end 

    thisgcd = oldr
    if thisgcd < 0 
        thisgcd, olds, oldt = -thisgcd, -olds, -oldt 
    end
    return thisgcd, olds, oldt 
end 

function project3d(a::Int64, b::Int64, c::Int64, d::Int64)::Vector{Int64}
    """
    Finds x, y, z, w such that ax + by + cz + dw = 1 assuming that there exists a GL(4, Z) transformation 
    that takes point [a, b, c, d] -> [1, b, c, d] 
    """
    if abs(a) == 1
        return [a, 0, 0, 0]
    end
    bcd = gcd(b,c,d) 

    # for a GL(4, Z) transform, "a" can only be scaled by +-1 
    if (1-a) % bcd == 0 
        a_scale = 1 
    elseif (1+a) % bcd == 0 
        a_scale = -1
    else
        return ["Failed", "Failed", "Failed", "Failed"]
    end 

    scale = div((1-a_scale*a), bcd)
    b_ = div(b, bcd)
    c_ = div(c, bcd)
    d_ = div(d, bcd)

    gcd1, x1, x2 = bezcoeffs(b_, c_)
    gcd2, y1, y2 = bezcoeffs(c_, d_)

    gcd3, z1, z2 = bezcoeffs(gcd1, gcd2)

    
    if (a_scale*a) + (b*scale*z1*x1) + c*(scale*(z1*x2+z2*y1))+d*(scale*z2*y2) != 1
        return ["Danger", "Danger", "Danger", "Danger"]
    end 

    return [a_scale, scale*z1*x1, scale*(z1*x2+z2*y1), scale*z2*y2]
end

function project2d(a::Int64, b::Int64, c::Int64)::Vector{Int64}
    """
    Finds x, y, z such that ax + by + cz = 1 assuming that there exists a GL(3, Z) transformation 
    that takes point [a, b, c] -> [1, b, c]
    """

    if abs(a) == 1
        return [a, 0, 0]
    end
    bc = gcd(b,c)

    if (1-a) % bc == 0  
        a_scale = 1 
    elseif (1+a) % bc == 0 
        a_scale = -1 
    else  
        return ["Failed", "Failed", "Failed"]
    end 

    scale = div((1- a_scale*a), bc)
    b_ = div(b, bc)
    c_ = div(c, bc)

    gcd1, x1, x2 = bezcoeffs(b_, c_)

    if (a_scale*a) + (b*scale*x1) + (c*scale*x2) != 1
        return ["Danger", "Danger", "Danger"]
    end 

    return [a_scale, scale*x1, scale*x2]
end

function transform_vertices(p::Int64, dim::Int64, poly, s::Vector{Int64}, r::Vector{Int64})::Matrix{Int64}
    """
    Given the scale factors s1, s2, s3, s4, it transforms first component of each point (a,b,c,d) 
    from a -> s1*a+s2*b+s3*c+s4*d
    Given the scale factors r1, r2, r3, it transforms ith component (i>1) of each point (a,b,c,d) 
    from b -> b + a*r1//c -> c + a*r2// d -> d + a*r3
    """
    proj = zeros(Int, dim, p)
    for i=1:p
        for j=1:dim
            proj[1,i] += poly[j,i]*s[j]
        end
    end 
    for i=1:p
        for j=1:dim-1
            proj[j+1,i] = poly[j+1, i] + proj[1,i]*r[j]
        end
    end  

    return proj
end 

function delzerofrombase(p::Int64, poly::Matrix{Int64})::Matrix{Int64}
    """    
    Removes points [0,0] from base that were completely projected 
    """
    newpoly = zeros(Int, 2, p)
    b = 0
    # b = size of base
    for i=1:p
        if poly[1, i] != 0 || poly[2, i] != 0
            b += 1
            newpoly[:, b] = poly[:, i]
        end
    end
    return newpoly[:, 1:b]
end 

function getbasesubpolytopeandtransform(p::Int64, sp::Int64, poly::Matrix{Int64}, subpoly::Matrix{Int64}, transforms::Array{Int64, 3}, fibnum::Int64)::Matrix{Int64}
    """
    Projects out the fiber and returns the base 
    """
    # p = size of poly; sp = size of subpoly/fiber 
    point1 = subpoly[:, 1]
    point2 = zeros(Int, 3)
    second_point = zeros(Int, 4)

    poly3d = zeros(Int, 3, p)
    poly2d = zeros(Int, 2, p)

    scalars1 = project3d(point1[1], point1[2], point1[3], point1[4])
    firstprojected = transform_vertices(p, 4, poly, scalars1, [-point1[2], -point1[3], -point1[4]])

    # result of projecting out one point 
    poly3d[:, :] .= firstprojected[2:4, 1:p]
    
    for i=2:sp
        # find a point in the fiber that was not already projected out
        second_point = transform_vertices(1, 4, subpoly[:, i], scalars1, [-point1[2], -point1[3], -point1[4]])
        point2 = second_point[2:4]
        if gcd(point2[1], point2[2], point2[3]) == 1
            break
        end
    end 

    scalars2 = project2d(point2[1], point2[2], point2[3])
    secondprojected = transform_vertices(p, 3, poly3d, scalars2, [-point2[2], -point2[3]])

    # result of projecting out two points (we have now projected out the fiber)
    poly2d[:, :] .= secondprojected[2:3, 1:p]

    # get the four transformations 
    # mat1 takes [a, b, c, d] -> [1, b, c, d] and mat2 takes [1, b, c, d] -> [1, 0, 0, 0]
    # mat3 takes [0, x, y, z] -> [0, 1, y, z] and mat4 takes [0, 1, y, z] -> [0, 1, 0, 0]
    mat1 = [
        scalars1[1]  scalars1[2]  scalars1[3]  scalars1[4]
        0 1 0 0 
        0 0 1 0 
        0 0 0 1
    ]

    mat2 = [
        1 0 0 0 
        -point1[2] 1 0 0 
        -point1[3] 0 1 0
        -point1[4] 0 0 1
    ]

    mat3 = [
        1 0 0 0 
        0 scalars2[1]  scalars2[2]  scalars2[3]
        0 0 1 0 
        0 0 0 1
    ]

    mat4 = [
        1 -second_point[1] 0 0 
        0 1 0 0 
        0 -second_point[3] 1 0 
        0 -second_point[4] 0 1
    ]

    # the full GL(4, Z) transform that projects the fiber is the matrix product of these four transformations
    full_transform = mat4*mat3*mat2*mat1
    transforms[:, :, fibnum] = full_transform

    # gets rid of extra points that were fully projected out
    return delzerofrombase(p, poly2d)
end

function get_sorted_base(base::Matrix{Int64})::Vector{Tuple{Int64, Int64}}
    """
    Get base in the format we want for canonical and compute_intersection_sequence
    i.e., change to a tuple and sort the rays by angles
    """
    rays = [(base[1, i], base[2, i]) for i in 1:size(base)[2]]
    # get rid of repetitions 
    rays_ = collect(Set(rays))
    # get rid of non-primitive rays
    p_rays = [(x[1], x[2]) for x in rays_ if gcd(x[1], x[2]) == 1]
    # sort rays
    s_base = sort(p_rays, by=r->((atan(r[2],r[1]))%(2*pi)))
    return s_base
end 

function check_singular(base::Vector{Tuple{Int64, Int64}})::Bool
    """
    Checks whether the base is singular by computing determinant of matrix formed by pairs of adjacent rays 
    """
    n = size(base)[1]
    for i=1:n
        if i == n 
            next = 1 
        else 
            next = i+1
        end 

        ray_matrix = [base[i][1] base[next][1]
                     base[i][2] base[next][2]]
        if abs(det(ray_matrix)-1) > 0.0001
            return true 
        end 
    end 
    return false 
end 

function convexhull(rays::Vector{Tuple{Int64, Int64}})::Vector{Vector{Any}}
    """
    Adds missing rays to singular bases to find corresponding smooth toric base
    """
    # find the size of polytope
    n = size(rays)[1]
    new_rays = [[(0, 0), 0] for _ in 1:2*n]
    index = 0

    for i in 1:n
        # Get three consecutive rays (with wrapping)
        prev = rays[mod1(i-1, n)]
        curr = rays[i]
        det = prev[1]*curr[2] - curr[1]*prev[2]

        index += 1
        if index > length(new_rays)
            resize!(new_rays, 2*length(new_rays))
        end 
        new_rays[index] = [prev, 1]

        # check if ray missing  
        if abs(det - 1) > 0.0001
            # iterate over the box containing the two rays 
            for x = min(prev[1], curr[1], 0):max(prev[1], curr[1], 0)
                for y = min(prev[2], curr[2], 0):max(prev[2], curr[2], 0)
                    extra = (x, y)
                    if !(extra in rays) & !(extra == (0, 0))
                        # check primitive
                        if gcd(x, y) == 1
                            # check whether it is under the line joining v and w                          
                            vec1 = extra .- prev
                            vec2 = curr .- extra
                            cp = cross([vec1[1], vec1[2], 0], [vec2[1], vec2[2], 0])

                            if cp[3] <= 0
                                # check whether between rays (using angles)
                                subrays = [prev, curr, extra]
                                sorted_ = sort(subrays, by=r->atan(r[2],r[1]))
                                                                
                                for i=1:3
                                    if sorted_[i] == prev 
                                        if sorted_[mod1(i+1, 3)] == extra 
                                            index += 1
                                            if index > length(new_rays)
                                                resize!(new_rays, 2*length(new_rays))
                                            end 
                                            
                                            new_rays[index] = [extra, 0]
                                        end
                                    end
                                end 
                            end
                        end
                    end
                end
            end
        end
    end
    return new_rays[1:index]
end

function canonical(intersections::Vector{Int64})::Vector{Int64}
    """
    Return the lexicographically first equivalent sequence under rotations and reflections.
    """
    n = length(intersections)
    min_sequence = intersections
    
    # Check all rotations
    for i in 1:n
        # Forward rotation
        rotated = circshift(intersections, -(i-1))
        if rotated < min_sequence
            min_sequence = copy(rotated)
        end
        
        # Reversed rotation
        reversed = reverse(rotated)
        if reversed < min_sequence
            min_sequence = copy(reversed)
        end
    end
    return min_sequence
end

function canonical_sing(intersections::Matrix{Int64})::Matrix{Int64}
    """
    Finds the canonical ordering of sequence if base is singular
    The format is different because this keeps track of which rays were missing in the singular base 
    """
    n = size(intersections)[2]
    min_sequence = intersections
    
    # Check all rotations
    for i in 1:n
        # Forward rotation
        rotated = circshift(intersections, (0,-(i-1)))
        if rotated[1,:] < min_sequence[1,:]
            min_sequence = copy(rotated)
        end
        if (rotated[1,:] == min_sequence[1,:]) && (rotated[2,:] < min_sequence[2,:])
            min_sequence = copy(rotated)
        end
        
        # Reversed rotation
        reversed = reverse(rotated, dims=2)
        if reversed[1,:] < min_sequence[1,:]
            min_sequence = copy(reversed)
        end
        if (reversed[1,:] == min_sequence[1,:]) && (reversed[2,:] < min_sequence[2,:])
            min_sequence = copy(reversed)
        end
    end
    
    return min_sequence
end


function compute_intersection_sequence(rays::Vector{Tuple{Int,Int}})::Vector{Int}
    """
    Compute intersection numbers from a sequence of ordered rays.
    For each ray, uses the relation c * r_curr + r_prev + r_next = 0
    where c is the intersection number we're solving for.
    Returns the sequence of intersection numbers.
    """
    n = size(rays)[1]
    intersections = Vector{Int}(undef, n)
    
    for i in 1:n
        # Get three consecutive rays (with wrapping)
        prev = rays[mod1(i-1, n)]
        curr = rays[i]
        next = rays[mod1(i+1, n)]

        # If curr = (x,y), then c * (x,y) + (x_prev,y_prev) + (x_next,y_next) = (0,0)
        # Taking either component gives us c
        if curr[1] != 0
            c = -(prev[1] + next[1]) ÷ curr[1]
        else
            c = -(prev[2] + next[2]) ÷ curr[2]
        end
        
        intersections[i] = c
    end
    
    return intersections
end

function compute_intersection_sequence_sing(rays::Vector{Vector{Any}})::Matrix{Int64}
    """
    Computes intersection sequence of toric base that corresponds to a singular base 
    The format is different because this keeps track of which rays were missing in the singular base 
    """
    n = size(rays)[1]
    intersections = zeros(Int, 2, n)
    
    for i in 1:n
        # Get three consecutive rays (with wrapping)
        prev = rays[mod1(i-1, n)][1]
        curr = rays[i][1]
        next = rays[mod1(i+1, n)][1]
        # If curr = (x,y), then c * (x,y) + (x_prev,y_prev) + (x_next,y_next) = (0,0)
        # Taking either component gives us c
        if curr[1] != 0
            c = -(prev[1] + next[1]) ÷ curr[1]
        else
            c = -(prev[2] + next[2]) ÷ curr[2]
        end
        
        intersections[:, i] = collect([c, rays[i][2]])
    end

    return intersections
end

function get_base_number(base::Vector{Int64})
    """
    Base number according to toric_base_list in 1204.0283
    """
    index = 0
    for b in toric_bases_list
        index += 1 
        if base == b
            return index 
        end
    end 
    return base
end

function get_fiber_and_base(points::Matrix{Int64}, vertices::Matrix{Int64}, p::Int64, v::Int64, short1::Vector{Int64}, short2::Vector{Int64}, fibnum::Int64, fibers::Array{Int64, 3}, fibsize::Vector{Int64}, bases::Dict{Any, Any}, transforms::Array{Int64, 3})::Int64
    """
    short1 and short2 are points x and y that span the fiber 
    This function saves the fiber subpolytope in fibers array and updates the counter in the bases dictionary 
    for number of subpolytopes of each type found corresponding to each base type. 
    """
    
    # call getfibersubpolytope to get the fibertype [ft], number of points in subpolytope [sp], and the subpolytope [subpoly] 
    ft, sp, subpoly = getfibersubpolytope(points, vertices, p, v, short1, short2)
    fibnum += 1

    if fibnum > size(fibers)[3]
        # resize arrays 
        fibers = cat(fibers, zeros(Int64, 4, p, 2*fibnum), dims=3)
        fibsize = cat(fibsize, zeros(Int64, 2*fibnum), dims=1)
        transforms = cat(transforms, zeros(Int64, 4, 4, 2*fibnum), dims=3)
    end

    fibers[:, :, fibnum] .= subpoly
    fibsize[fibnum] = sp

    # calls getbasesubpolytopeandtransform to compute the base subpolytope and then sorts it by angles 
    base = getbasesubpolytopeandtransform(p, sp, points, subpoly, transforms, fibnum)
    s_base = get_sorted_base(base)

    # check if base singular and call appropriate functions to get the intersection sequences 
    if check_singular(s_base)
        new_base = collect(Set(convexhull(s_base)))
        sorted = sort(new_base, by=r->atan(r[1][2], r[1][1]))
        int_seq = canonical_sing(compute_intersection_sequence_sing(sorted))
        base_number = int_seq
    else
        int_seq = canonical(compute_intersection_sequence(s_base))
        base_number = get_base_number(int_seq)
    end 
    
    already_found = keys(bases)
    if !(base_number in already_found)
        bases[base_number] = Dict(1 => [0,[],[]], 2 => [0,[],[]], 3 => [0,[],[]], 4 => [0,[],[]], 5 => [0,[],[]], 6 => [0,[],[]], 7 => [0,[],[]], 8 => [0,[],[]], 9 => [0,[],[]], 10  => [0,[],[]], 
        11 => [0,[],[]], 12 => [0,[],[]], 13 => [0,[],[]], 14 => [0,[],[]], 15 => [0,[],[]], 16 => [0,[],[]])
    end
    bases[base_number][ft][1] += 1
    push!(bases[base_number][ft][2], transforms[:,:, fibnum])
    push!(bases[base_number][ft][3], fibers[:, 1:sp, fibnum])

    return fibnum
end

function delzerofromdict(dict::Dict{Int64, Vector{Any}})
    """
    Remove keys with value 0 from dictionary
    """ 
    for i = 1:16 
        if dict[i] == [0, [], []] 
            delete!(dict, i)
        end
    end
end

function check_for_fiber(points::Matrix{Int64}, vertices::Matrix{Int64}, p::Int64, v::Int64, short::Array{Int64, 3}, s::Vector{Int64})::Tuple{Int64, Dict{Any, Any}}
    """
    Checks three conditions for existence of 2d reflexive subpolytope, calls get_fiber_and_base to compute the fiber 
    and base subpolytopes, and stores fibers and bases in arrays 
    """
    # dictionary to contain all bases 
    bases = Dict()

    # arrays to contain all fibers
    # fibsize keep track of size of each fiber  
    fibers = zeros(Int, 4, p, 500)
    transforms = zeros(Int, 4, 4, 500)
    fibsize = zeros(Int, 500)
    fibnum = 0

    for i = 1:s[1]
        if member(-short[:, i, 1], short[:,:, 1], s[1])
            for j = i+1:s[1]
                checked = false

                if (member(-short[:, j, 1], short[:,:, 1], s[1])
                    && short[:, j, 1] !=  -short[:, i, 1])

                    for m=1:fibnum
                        if (member(short[:, i, 1], fibers[:,:,m], fibsize[m]))  && (member(short[:, j, 1], fibers[:,:,m], fibsize[m]))
                            # if both points are in a fiber already found, skip 
                            checked = true 
                            break
                        end 
                    end

                    if !checked 
                        fibnum = get_fiber_and_base(points, vertices, p, v, short[:, i, 1], short[:, j, 1], fibnum, fibers, fibsize, bases, transforms)
                    end
                end
            end
        end
    end

    for i = 1:s[2]
        for j = i+1:s[2]
            checked = false
            if (member(-short[:, i, 2]-short[:, j, 2], short[:,:, 1], s[1])
                || member(-short[:, i, 2]-short[:, j, 2], short[:,:, 2], s[2]))

                for m=1:fibnum
                    if (member(short[:, i, 2], fibers[:,:,m], fibsize[m]))  && (member(short[:, j, 2], fibers[:,:,m], fibsize[m]))
                        checked = true 
                        break
                    end 
                end

                if !checked
                    fibnum = get_fiber_and_base(points, vertices, p, v, short[:, i, 2], short[:, j, 2], fibnum, fibers, fibsize, bases, transforms)
                end
            end
        end
    end

    for i = 1:s[3]
        for j = i+1:s[3]
            checked = false
            if member(-short[:, i, 3]-short[:, j, 3], short[:,:, 1], s[1], 2)

                for m=1:fibnum
                    if (member(short[:, i, 3], fibers[:,:,m], fibsize[m]))  && (member(short[:, j, 3], fibers[:,:,m], fibsize[m]))
                        checked = true 
                        break
                    end 
                end

                if !checked
                    fibnum = get_fiber_and_base(points, vertices, p, v, short[:, i, 3], short[:, j, 3], fibnum, fibers, fibsize, bases, transforms)
                end
            end
        end
    end

    for bn in keys(bases)
        # delete keys with value = 0 
        delzerofromdict(bases[bn])
    end

    return fibnum, bases
end

function reformat(bases::Dict{Any, Any})::Vector{Any}
    """
    Changes bases format from a dictionary to a vector
    """ 
    data = []  
    
    for base in keys(bases)
        # change the base code to have 1/0 for singular
        fibers = []
        for fiber in keys(bases[base])
            push!(fibers, [fiber, bases[base][fiber][1], bases[base][fiber][2], bases[base][fiber][3]])
        end 

        if (typeof(base) == Int)
            push!(data, [0, base, parse(Int, h11_bases[base]), [], fibers])         
        else
            base_number = get_base_number(base[1, :])
            push!(data, [1, base_number, parse(Int, h11_bases[base_number]), base[2, :], fibers]) 
        end
    end
    return data 
end 

function printbase(basedata::Any, out::IOStream)
    """
    Writes the base array to the output file 
    """
    if !isa(basedata, Int) && !isa(basedata, Matrix{Int})
        item_num = 0 
        write(out, "[")
        for item in basedata
            item_num += 1 
            if item_num > 1
                write(out, ", ")
            end 
            printbase(item, out)
        end 
        write(out, "]")
    elseif isa(basedata, Matrix{Int})
        write(out, "[")
        for i=1:size(basedata)[2]
            if i > 1
                write(out, ", ")
            end 
            write(out, string(basedata[1:4, i]))
        end 
        write(out, "]")
    else 
        write(out, string(basedata))
    end 
end 

function check(file::String, output::String, cumulative = "fiber-data")
    """
    Checks everything in a file, outputs to file-f, saves data
    """
    f = open(string(file))
    out = open(string(output), "a")
    data = open(cumulative, "a")
    points = zeros(Int, 4, 1000)
    vertices = zeros(Int, 4, 36)
    short = zeros(Int, 4, 1000, 3)
    s = zeros(Int, 3)
    # change this to actual number of polytopes per file 
    polytopes_per_file = 1000000
    total = 0; fibers = 0

    if (occursin("-", file))
        # corrects for polytope number when using multiple files in a series, e.g. v08-n.txt
        index = last(findall("-", file))[1]
        try
            filenum = file[index+1:last(findall(".", file))[1]-1]
            filenum = parse(Int, filenum)
            already_seen = polytopes_per_file*(filenum-1)
            total += already_seen 
        catch
            println("Error processing file. Please name your file as 'name-number.txt' or 'name.txt' without dashes.")
        end
    end 

    p, v, h, active = get(f, points, vertices)

    while active
        total += 1

        findshort(points, vertices, p, v, short, s)
        if p > 1000 || v > 36
            write(out, "error")
        end
        fibernum, bases = check_for_fiber(points, vertices, p, v, short, s)
        bases_output = reformat(bases)
        write(out, string("[[$v, $total], [$h], $fibernum, ["))
        base_ind = 0 

        # deals with brackets etc, ensures that output is in a nice format 
        for base in bases_output
            base_ind += 1 
            if base_ind > 1
                write(out, ", ")
            end 
            printbase(base, out)
        end        
        write(out, "]]", "\n")

        fibers += fibernum
        p, v, h, active = get(f, points, vertices)
    end

    close(f)
    close(out)
    write(data, string(file, " ", fibers, " ", total, " ", total-fibers, " ",  1.0-fibers/total, "\n"))
    close(data)
end

h11_bases = readlines("h11_bases.dat") 
toric_bases_list = read_vector_of_vectors("toric_bases.txt")
check(ARGS[1], ARGS[2])
