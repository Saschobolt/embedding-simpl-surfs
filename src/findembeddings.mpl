# Load the LinearAlgebra package for vector and matrix operations
with(LinearAlgebra):

# Load Simplicial Surface package
libname := (libname, "/users/stuettgen/custom_maple_libraries"):
read("/users/stuettgen/Documents/Masterarbeit/Maple/Initialization_local.mpl"):

# change Jupyter Renderer Outpu
# with(Jupyter):
# SetOutputRendererByType(equation, "text/plain"):
# SetOutputRendererByType(exprseq, "text/plain"):
# SetOutputRendererByType(Vector, "text/plain"):

with(ArrayTools):

# =============================================================================
# HELPER PROCEDURE 1: Trilaterate
# Calculates intersection in the canonical XY-plane system.
# =============================================================================
Trilaterate := proc(r1, r2, r3, d, i, j, assumptions:={})
    local x, y, z_sq, z, p1, p2;
    description "Calculates intersection for spheres in canonical alignment.";

    if d = 0 or j = 0 then
        # This case is handled by AlignTriangle's collinearity check,
        # but is included for robustness.
        error "Sphere centers cannot be collinear.";
    end if;

    x := evala((r1^2 - r2^2 + d^2) / (2 * d));
    userinfo(2, Trilaterate, "x := ", x);
    y := evala((r1^2 - r3^2 + i^2 + j^2 - 2 * i * x) / (2 * j));
    userinfo(2, Trilaterate, "y := ", y);
    z_sq := r1^2 - x^2 - y^2;
    userinfo(2, Trilaterate, "z_sq := ", z_sq);

    # if z_sq < 0 then
    #     return [];
    # elif z_sq = 0 then
    #     return [[x, y, 0.0]]; # Return as a list of one point
    # else
    z := evala(sqrt(z_sq));
    userinfo(2, Trilaterate, "z := ", z);
    p1 := Vector([x, y, z]);   # Point with positive z
    p2 := Vector([x, y, -z]);  # Point with negative z
    return [p1, p2];
    # end if;
end proc:

# =============================================================================
# HELPER PROCEDURE 2: AlignTriangle
# Computes the transformation to move sphere centers to the XY-plane.
# =============================================================================
AlignTriangle := proc(p1::Vector, p2::Vector, p3::Vector, assumptions:={}, check:=true)
    local v1, v2, u, v, w, R, T, p1_new, p2_new, p3_new, cross_prod_norm;
    description "Computes rigid transformation to align 3 points with XY plane.";

    T := p1;
    userinfo(2, AlignTriangle, `coordinates T[1] := `, T[1]);
    userinfo(2, AlignTriangle, `coordinates T[2] := `, T[2]);
    userinfo(2, AlignTriangle, `coordinates T[3] := `, T[3]);
    v1 := p2 - T;
    userinfo(2, AlignTriangle, `coordinates v1[1] := `, v1[1]);
    userinfo(2, AlignTriangle, `coordinates v1[2] := `, v1[2]);
    userinfo(2, AlignTriangle, `coordinates v1[3] := `, v1[3]);
    v2 := p3 - T;
    userinfo(2, AlignTriangle, `coordinates v2[1] := `, v2[1]);
    userinfo(2, AlignTriangle, `coordinates v2[2] := `, v2[2]);
    userinfo(2, AlignTriangle, `coordinates v2[3] := `, v2[3]);


    if check then 
        if Norm(v1, 2) = 0 then
            error "Sphere centers 1 and 2 are coincident.";
        end if;
        if Norm(CrossProduct(v1, v2), 2) = 0 then
            error "Sphere centers are collinear; cannot define a unique plane.";
        end if;
    end if;

    u := evala(Normalize(v1, 2)) assuming op(assumptions);
    userinfo(2, AlignTriangle, `u[1] := `, u[1]);
    userinfo(2, AlignTriangle, `u[2] := `, u[2]);
    userinfo(2, AlignTriangle, `u[3] := `, u[3]);
    w := evala(Normalize(CrossProduct(v1, v2), 2)) assuming op(assumptions);
    userinfo(2, AlignTriangle, `w[1] := `, w[1]);
    userinfo(2, AlignTriangle, `w[2] := `, w[2]);
    userinfo(2, AlignTriangle, `w[3] := `, w[3]);
    v := evala(CrossProduct(w, u)) assuming op(assumptions);
    userinfo(2, AlignTriangle, `v[1] := `, v[1]);
    userinfo(2, AlignTriangle, `v[2] := `, v[2]);
    userinfo(2, AlignTriangle, `v[3] := `, v[3]);


    R := evala(Matrix([Vector(u), Vector(v), Vector(w)])^+) assuming op(assumptions);

    userinfo(1, AlignTriangle, `Computed transformation matrix`);
    userinfo(2, AlignTriangle, `R := `, R);
    
    p1_new := ZeroVector(nops(p1));
    userinfo(2, AlignTriangle, `p1_new := `, p1_new);
    p2_new := evala(R . (p2 - T)) assuming op(assumptions);
    userinfo(2, AlignTriangle, `p2_new := `, p2_new);
    p3_new := evala(R . (p3 - T)) assuming op(assumptions);
    userinfo(2, AlignTriangle, `p3_new := `, p3_new);

    return [p1_new, p2_new, p3_new], R, T;
end proc:

# =============================================================================
#         Generalized Intersection of Three Arbitrary Spheres in 3D
# =============================================================================
#
# This script provides a primary function, IntersectThreeSpheres, which
# computes one of the two possible intersection points of three spheres defined
# by their centers (C1, C2, C3) and radii (r1, r2, r3).
#
# The method works as follows:
# 1. A rigid transformation (rotation + translation) is computed using the
#    helper function `AlignTriangle`. This transformation moves the sphere
#    centers to a simplified, canonical orientation on the XY-plane.
# 2. The intersection is calculated in this simplified system using the helper
#    function `Trilaterate`.
# 3. A `sign_choice` argument (+1 or -1) determines which of the two
#    intersection points (on either side of the triangle's plane) is chosen.
# 4. The chosen intersection point is then mapped back to the original
#    coordinate system using the inverse of the initial transformation.
#
# =============================================================================
IntersectThreeSpheres := proc(C1::Vector, C2::Vector, C3::Vector, r1, r2, r3, sign_choice::{+1,-1}, assumptions:={})
    local alignment_data, transformed_centers, R, T, d, i, j,
          local_points, selected_point_local, P_local, P_final;

    description "Computes the intersection of three arbitrary spheres.";

    # --- Step 1: Align the sphere centers to the XY plane ---
    # The try-catch block will handle collinearity errors from AlignTriangle.
    userinfo(1, IntersectThreeSpheres, `Computing transformed sphere centers.`);
    try
        transformed_centers, R, T := AlignTriangle(C1, C2, C3, assumptions);
    catch:
        error "Error during alignment: " || lasterror[2];
    end try;

    userinfo(1, IntersectThreeSpheres, `Found transformed sphere centers.`);

    # --- Step 2: Extract parameters for the canonical solver ---
    # d is the distance between transformed C1 and C2 (it's the x-coord of C2_new)
    d := transformed_centers[2][1];
    # i and j are the coordinates of the transformed C3
    i := transformed_centers[3][1];
    j := transformed_centers[3][2];

    userinfo(2, IntersectThreeSpheres, `d := `, d);
    userinfo(2, IntersectThreeSpheres, `i := `, i);
    userinfo(2, IntersectThreeSpheres, `j := `, j);

    # --- Step 3: Solve for the intersection in the simplified system ---
    userinfo(1, IntersectThreeSpheres, `Trilaterating in transformed coordinate system.`);
    local_points := Trilaterate(r1, r2, r3, d, i, j, assumptions);

    # If Trilaterate returned a string, it means no solution, so propagate it.
    if nops(local_points) = 0 then
        return;
    end if;

    # --- Step 4: Select the desired intersection point ---
    # `local_points` is a list of lists, e.g., [[x,y,z], [x,y,-z]]
    if nops(local_points) = 1 then
        # Tangent case: only one point, sign_choice doesn't matter.
        P_local := local_points[1];
    else
        # Two points exist. The first has z >= 0, the second has z <= 0.
        if sign_choice = +1 then
            P_local := local_points[1]; # Positive z in local system
        else
            P_local := local_points[2]; # Negative z in local system
        end if;
    end if;

    # --- Step 5: Apply the inverse transformation ---
    # To restore the point to its original coordinate system, we first apply
    # the inverse rotation (Transpose(R)) and then the inverse translation (+ T).
    userinfo(1, IntersectThreeSpheres, `Applying inverse transformation.`);
    P_final := Transpose(R) . P_local + T;
    userinfo(2, IntersectThreeSpheres, `P_final := `, P_final);
    # P_final := evala(P_final) assuming op(assumptions);
    # userinfo(2, IntersectThreeSpheres, `P_final after evala := `, P_final);

    return P_final;
end proc:

sqdist := proc(v::Vector, w::Vector) # square distance between two 3d-vectors v,w
    local x;
    x := v-w;
    return evala(add(map(entry -> entry^2, x)));
end proc:

sqnorm := proc(v::Vector) # return the square of the 2-norm of a ed-vector v
    return sqdist(v,Vector([seq(0, 1..nops(v))]))
end proc:


# =============================================================================
# PROCEDURE: Embed the next vertex
# Computes the candidates for the position of the next vertex, if embedded_coords are already given. 
# neighbor_indices are the indices of the neighbors of the next vertex in the already embedded list.
# distances include the distances of the next vertex to the neighbors in the neighbor_indices list.
# next_param_index is the index of the next parameter that is to be used, if the next vertex is underconstrained. Parameters are expected to be _t[1], _t[2], _t[3], ...
# sign_choice indicates the half space the next vertex is supposed to be located in wrt the position of neighbor_indices[[1..3]]. (three finger rule with vectors neighbor[1] -> neighbor[2], neighbor[1] -> neighbor[3])
# =============================================================================
EmbedNextVertex := proc(embedded_coords, neighbor_indices, distances, next_param_index, {sign_choice::{+1,-1} := +1})
    local new_vertex_index, new_vertex_coords, n, neighbor_points, new_embedding;
    
    description "Computes coordinates for a new vertex, including canonical placement for the first three.";

    if max(op(neighbor_indices)) > nops(embedded_coords) then
        error "neighbor indices not valid: more neighbor indices than embedded coordinates."
    end if;

    if not nops(neighbor_indices) = nops(distances) then
        error "number of neighbor indices and distances don't match."
    end if;

    new_vertex_index := nops(embedded_coords) + 1;

    n := nops(neighbor_indices);

    # --- PART 1: Canonical placement for the first three vertices ---
    if new_vertex_index = 1 then
        new_vertex_coords := Vector([0, 0, 0]);

        return [[new_vertex_coords]], next_param_index;

    elif new_vertex_index = 2 then
        if n <> 1 or neighbor_indices[1] <> 1 then error "Vertex 2 must be a neighbor of vertex 1 only."; end if;
        new_vertex_coords := Vector([distances[1], 0, 0]);
        new_embedding := [[op(embedded_coords), new_vertex_coords]];

        return new_embedding, next_param_index;

    elif new_vertex_index = 3 then
        local v, w, dvw, x, y;
        v := embedded_coords[1]; 
        w := embedded_coords[2];
        dvw := Norm(v - w, 2); # distance between v and w
        x := (dvw^2 + distances[1]^2 - distances[2]^2) / (2*dvw);
        y := sqrt(distances[2]^2 - x^2);
        new_vertex_coords := Vector([x, y, 0]);

        new_embedding := [[op(embedded_coords), new_vertex_coords]];

        return new_embedding, next_param_index;
    end if;

    # --- PART 2: General embedding logic for all subsequent vertices ---
    neighbor_points := embedded_coords[neighbor_indices];

    if n <= 2 then
        local p_idx, new_next_param_index, dists;
        new_next_param_index := next_param_index;
        
        for p_idx from 1 to nops(embedded_coords) do
            if not p_idx in neighbor_indices then
                neighbor_points := [op(neighbor_points), embedded_coords[p_idx]];
                dists := [op(distances), _t[new_next_param_index]];
                new_next_param_index := new_next_param_index + 1;
            end if;

            if nops(neighbor_points) = 3 then
                new_vertex_coords := IntersectThreeSpheres(neighbor_points[1], neighbor_points[2], neighbor_points[3], dists[1], dists[2], dists[3], sign_choice);

                new_embedding := [[op(embedded_coords), new_vertex_coords]];

                return new_embedding, new_next_param_index;
            end if;
        end do;
        
        if nops(neighbor_points) < 3 then error "Could not find enough unique anchor points."; end if;

    elif n = 3 then
            new_vertex_coords := IntersectThreeSpheres(neighbor_points[1], neighbor_points[2], neighbor_points[3], distances[1], distances[2], distances[3], sign_choice);

            new_embedding := [[op(embedded_coords), new_vertex_coords]];

            return new_embedding, next_param_index;

    else # n >= 4
        local remaining_points, remaining_dists, coords, vars, eqns, numer_eqns, denom_eqns, ext, s, u, sub, new_coords, neigh;

        new_vertex_coords := IntersectThreeSpheres(neighbor_points[1], neighbor_points[2], neighbor_points[3], distances[1], distances[2], distances[3], sign_choice);

        remaining_points := neighbor_points[4..nops(neighbor_points)];
        remaining_dists := distances[4..nops(distances)];

        coords := [op(embedded_coords), new_vertex_coords];

        # need to solve system of equations to ensure that distance to the remaining neighbors is as desired.
        vars := map(i -> _t[i], [$1..(next_param_index - 1)]);
        eqns := map(i -> sqdist(new_vertex_coords, remaining_points[i]) - remaining_dists[i]^2, [$1..nops(remaining_points)]);
        numer_eqns := evala(convert(numer(eqns), RootOf));
        denom_eqns := evala(convert(denom(eqns), RootOf));
        ext := [op(indets(numer_eqn, algext))];
        ext := sort(ext, key = (expr -> nops(indets(expr, algext))));
        ext := [op(map(expr -> [expr], ext))];
        userinfo(1, EmbedNextVertex, `vars := `, vars);
        userinfo(1, EmbedNextVertex, `eqns := `, eqns);
        userinfo(1, EmbedNextVertex, `numer_eqns := `, numer_eqns);
        userinfo(1, EmbedNextVertex, `denom_eqns := `, denom_eqns);
        userinfo(1, EmbedNextVertex, `ext := `, ext);

        if eqn = [0$nops(eqn)] then
            s := [[]];
        else
            s, u := `SimplicialSurfaceEmbeddings/solve_polynomial_system`(numer_eqns, vars, denom_eqns, ext);
            userinfo(1, EmbedNextVertex, `s := `, s);
        end if;

        new_embedding := [];

        for sub in s do
            try
                # print(evala(Simplify~(subs(sub, coordsNew))));
                new_coords := evala(subs(sub, coords));
                # print(evala(dist(coordsNew1[tip], coordsNew1[v])^2 - elng[Search(e, edgesSets)]^2) assuming real);
                # print(coordsNew1);
                # sanity check that all edge lengths are as desired.
                for neigh in neighbor_indices do
                    userinfo(1, EmbedNextVertex, `sqdist between`, neigh, `and new vertex `, nops(coords), `:`, evala(sqdist(new_coords[neigh], new_coords[nops(new_coords)])), `. Should be `, distances[neigh]^2, `(Difference : `, evala(sqdist(new_coords[neigh], new_coords[nops(new_coords)]) - distances[neigh]^2), `)`);
                    if not evala(sqdist(new_coords[neigh], new_coords[nops(new_coords)]) - distances[neigh]^2) = 0 then
                        next sub;
                    end if;
                end do;
                new_embedding := [op(new_embedding), new_coords];
            catch: next sub;
            end try;
        end do;
        
        return new_embedding, next_param_index;
    end if;
end proc:

# =============================================================================
# PROCEDURE: Perform plan
# Embed a graph by performing the steps indicated by the input. The list is a list of lists. Each element indicates one step of the plan.
# Each step is a list containing three elements: The neighbors of the vertex that is embedded in this step, the distances to the neighbors, the sign of the halfspace wrt the first three neighbor positions. (see EmbedNextVertex)
# =============================================================================
PerformPlan := proc(plan)
    local step, emb, embeddings, next_param_index, neighs, dists, sign_choice, new_embeddings, new_embeddings_this_emb;

    embeddings := [[]]; # list of list of vectors -> valid embeddings
    next_param_index := 1;

    for step in plan do
        # get relevant information from the step
        neighs := step[1];
        dists := step[2];
        sign_choice := step[3];

        userinfo(1, PerformPlan, `step:`, step);
        userinfo(2, PerformPlan, `neighs:`, neighs);
        userinfo(2, PerformPlan, `dists:`, dists);
        userinfo(2, PerformPlan, `sign_choice:`, sign_choice);
        userinfo(1, PerformPlan, `number of embeddings currently:`, nops(embeddings));
        userinfo(2, PerformPlan, `current embeddings:`, embeddings);

        new_embeddings := [];

        # perform the step for each embedding in embeddings.
        for emb in embeddings do
            userinfo(2, PerformPlan, `currently handling emb`, emb);
            new_embeddings_this_emb, next_param_index := EmbedNextVertex(emb, neighs, dists, next_param_index, sign_choice);
            userinfo(2, PerformPlan, `New embeddings found:`, new_embeddings_this_emb);
            new_embeddings := [op(new_embeddings), op(new_embeddings_this_emb)];
        end do;
        embeddings := new_embeddings;
    end do;

    return embeddings;
end proc:

# infolevel[EmbedNextVertex] := 1;
# infolevel[PerformPlan] := 1;
# waterbomb_cell := PerformPlan([[[], [], 1], [[1], [1], 1], [[1,2], [a,b], 1], [[1,3], [a,c], -1], [[1,4], [1,b],-1], [[1,5], [a,b], -1], [[1,2,5], [a,b,c], 1]]);