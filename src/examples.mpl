read("Initialization_SSH.mpl"):

with(SimplicialSurfaceEmbeddings):
with(ListTools):
with(combinat):

CombSimplSurfByVertsOfFaces := proc(facelist::listlist)
    local surf;

    ASSERT(map(f -> evalb(nops(f) = 3), facelist) = [true$nops(facelist)], "Faces needs to be a list of lists of length 3.");
    ASSERT({op(map(f -> op(f), facelist))} = {$1..max(facelist)}, "Faces must be numbered from 1 to n.");

    surf := NewSurface();
    surf[vertex_names] := [$1..max(facelist)];
    surf[faces] := facelist;
    surf[edges] := map(e -> [op(e)], MakeUnique(map(f -> op(map(e -> {op(e)}, choose(f, 2))), facelist)));

    return surf
end proc:

Plesken := CombSimplSurfByVertsOfFaces([[1, 2, 3], [2, 1, 4], [3, 5, 1], [5, 4, 1], [2, 6, 3], [6, 2, 4], [6, 7, 3], [7, 5, 3], [4, 8, 6], [8, 4, 5], [7, 6, 8], [7, 8, 5]]);
Octahedron := CombSimplSurfByVertsOfFaces([[1, 2, 5], [2, 3, 5], [3, 4, 5], [4, 1, 5], [1, 2, 6], [2, 3, 6], [3, 4, 6], [4, 1, 6]]);
Icosahedron := CombSimplSurfByVertsOfFaces([[1,2,9],[1,2,10],[1,5,6],[1,5,9],[1,6,10],[2,7,8],[2,7,9],[2,8,10],[3,4,11],[3,4,12],[3,5,6],[3,5,11],[3,6,12],[4,7,8],[4,7,11],[4,8,12],[5,9,11],[6,10,12],[7,9,11],[8,10,12]]);
TriangularPrism := CombSimplSurfByVertsOfFaces([[1,2,3], [1,3,4], [2,3,5], [1,2,6], [3,4,9], [3,9,5], [2,5,8], [2,8,6], [1,6,7], [1,7,4], [4,7,9], [6,7,8], [5,8,9], [7,8,9]]);
OpenTriangularPrism := CombSimplSurfByVertsOfFaces([[1,2,3], [1,3,4], [2,3,5], [1,2,6], [3,4,9], [3,9,5], [2,5,8], [2,8,10], [1,6,7], [1,7,4], [4,7,9], [10,7,8], [5,8,9], [7,8,9]]);