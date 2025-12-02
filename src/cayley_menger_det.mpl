libname := (libname, "/users/stuettgen/custom_maple_libraries"):

read("/users/stuettgen/Documents/Masterarbeit/Maple/Initialization_local.mpl"):
with(LinearAlgebra):

# Procedure to compute the Cayley-Menger Determinant and Volume of a Simplex
# Parameters:
#   n: Number of vertices (corners)
#   DistMatrix: An n x n symmetric matrix where Entry [i,j] is the edge length between vertex i and vertex j
CayleyMengerSimplex := proc(n, DistMatrix)
    local M, i, j, detM, volSq, factor;

    # 1. Construct the Cayley-Menger Matrix M
    # M is an (n+1)x(n+1) matrix.
    # The top-left nxn block contains squared distances.
    # The last row and column contain 1s, with a 0 in the bottom-right.
    
    M := Matrix(n + 1, n + 1);

    for i from 1 to n do
        for j from 1 to n do
            M[i, j] := DistMatrix[i, j]^2;
        end do;
        # Border the matrix with 1s
        M[i, n + 1] := 1;
        M[n + 1, i] := 1;
    end do;
    
    # Set the bottom-right corner to 0
    M[n + 1, n + 1] := 0;

    # 2. Compute the Determinant
    detM := Determinant(M);

    # 3. Compute Volume Squared
    # For a simplex with n vertices (dimension k = n-1):
    # V^2 = ((-1)^n / (2^(n-1) * ((n-1)!)^2)) * det(M)
    
    factor := (-1)^n / (2^(n - 1) * (factorial(n - 1))^2);
    volSq := factor * detM;

    # Return a record containing the Matrix, the Determinant, and the Volume Squared
    return Record(Matrix = M, Determinant = detM, VolumeSquared = volSq);
end proc;

# --- EXAMPLE 1: Triangle (n=3) ---
# This reproduces Heron's Formula
# Let sides be a, b, c.
# Distance Matrix structure:
# [0, c, b]
# [c, 0, a]
# [b, a, 0]

DistTri := Matrix([
    [0, c, b], 
    [c, 0, a], 
    [b, a, 0]
]);

ResultTri := CayleyMengerSimplex(3, DistTri);

print("--- Triangle (Heron's Formula) ---");
print("Determinant:", ResultTri:-Determinant);
print("Area Squared:", simplify(ResultTri:-VolumeSquared));


# --- EXAMPLE 2: Tetrahedron (n=4) ---
# Vertices 1,2,3,4 with pairwise edge lengths d_ij
# Distance Matrix must be symmetric with 0 on diagonal

DistTetra := Matrix([
    [0, d12, d13, d14],
    [d12, 0, d23, d24],
    [d13, d23, 0, d34],
    [d14, d24, d34, 0]
]);

ResultTetra := CayleyMengerSimplex(4, DistTetra);

print("--- Tetrahedron Volume ---");
print("Volume Squared:", ResultTetra:-VolumeSquared);