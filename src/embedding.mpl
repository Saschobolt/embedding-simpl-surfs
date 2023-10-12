with(ListTools):
with(ArrayTools):
with(SimplicialSurfaceEmbeddings):
with(combinat):
with(CodeGeneration):

dist := proc(v,w)
    # aux function that returns the distance of the lists v and w as vectors
  return sqrt~(add((v-w)^~2));
end proc:

indmax := proc(l)
    # returns index of the first appearance of the maximal value in list l
    local sol;

    member(max(l), l, sol);
    return sol
end proc:


indmin := proc(l)
    # returns index of the first appearance of the minimal value in list l
    local sol;

    member(min(l), l, sol);
    return sol
end proc:

indwhere := proc(condition, l)
  # l is a list, condition is a boolean evaluated function
  # returns index of the first entry of l that satisfies the condition and 0 if no such entries exist
  local sol, entry;

  if nops(select(condition, l)) = 0 then
    return 0
  end if;

  entry := select(condition, l)[1];
  member(entry, l, sol);
  return sol
end proc:


outgoingEdges := proc(surf, verts)
    ASSERT(evalb(map(v->evalb(v in Vertices(surf)), verts) = [true$nops(verts)]), "the verts need to be vertices of surf.");

    return select(e -> evalb(nops({op(e)} intersect {op(verts)}) = 1), Edges(surf))
end proc:

Triangle := proc(l::list)
    # returns the coordinates of a triangle with edge length defined by l as a list of lists.
    # if the vertices of the triangle are called a,b,c, then
    # l[1] = ||a-b||
    # l[2] = ||b-c||
    # l[3] = ||c-a||
    local first, second, third;

    first := [0,0,0];
    second := [l[1],0,0];
    third := subs([r = second[1], r[1] = l[3], r[2] = l[2]], [(r^2+r[1]^2-r[2]^2)/(2*r), sqrt(-r^4+2*r^2*r[1]^2+2*r^2*r[2]^2-r[1]^4+2*r[1]^2*r[2]^2-r[2]^4)/(2*r), 0]);

    return evala(convert([first, second, third], RootOf));
end proc:


Tetrahedron := proc(l::list)
    # return the coodrinates of a tetrahedron with edge length defined by l as a list of lists.
    # if the vertices of the tetrahedron are called a,b,c,d then 
    # l[1] = ||a-b||
    # l[2] = ||b-c||
    # l[3] = ||c-a||
    # l[4] = ||a-d||
    # l[5] = ||b-d||
    # l[6] = ||c-d||
    local first, second, third, fourth;

    first := [0,0,0];
    second := [l[1],0,0];
    third := subs([r = second[1], r[1] = l[3], r[2] = l[2]], [(r^2+r[1]^2-r[2]^2)/(2*r), sqrt(-r^4+2*r^2*r[1]^2+2*r^2*r[2]^2-r[1]^4+2*r[1]^2*r[2]^2-r[2]^4)/(2*r), 0]);
    fourth := subs([r = second[1], s = third[1], t = third[2], r[1] = l[4], r[2] = l[5], r[3] = l[6]], 
                    [(r^2+r[1]^2-r[2]^2)/(2*r), -(r^2*s-r*s^2-r*t^2-r*r[1]^2+r*r[3]^2+s*r[1]^2-s*r[2]^2)/(2*r*t), sqrt(-r^4*s^2-r^4*t^2+2*r^3*s^3+2*r^3*s*t^2+2*r^3*s*r[1]^2-2*r^3*s*r[3]^2-r^2*s^4-2*r^2*s^2*t^2-4*r^2*s^2*r[1]^2+2*r^2*s^2*r[2]^2+2*r^2*s^2*r[3]^2-r^2*t^4+2*r^2*t^2*r[2]^2+2*r^2*t^2*r[3]^2-r^2*r[1]^4+2*r^2*r[1]^2*r[3]^2-r^2*r[3]^4+2*r*s^3*r[1]^2-2*r*s^3*r[2]^2+2*r*s*t^2*r[1]^2-2*r*s*t^2*r[2]^2+2*r*s*r[1]^4-2*r*s*r[1]^2*r[2]^2-2*r*s*r[1]^2*r[3]^2+2*r*s*r[2]^2*r[3]^2-s^2*r[1]^4+2*s^2*r[1]^2*r[2]^2-s^2*r[2]^4-t^2*r[1]^4+2*t^2*r[1]^2*r[2]^2-t^2*r[2]^4)/(2*r*t)]);

    return evala(convert([first, second, third, fourth], RootOf));
end proc:

# Tetrahedron := proc(l::list)
#     a := l[1]:
#     b := l[2]:
#     c := l[3]:
#     d := l[4]:
#     e := l[5]:
#     t := l[6]:

#     tetra := evala([[0, 0, 0], [a, 0, 0], [1/2*(a^2-b^2+c^2)/a, 1/2*(-a^4+2*a^2*b^2+2*a^2*c^2-b^4+2*b^2*c^2-c^4)^(1/2)/a, 0], 
#                     [1/2*(a^2+d^2-e^2)/a, 1/2*(-a^4+2*a^2*d^2+2*a^2*e^2-d^4+2*d^2*e^2-e^4)^(1/2)*(1/(a^4-a^2*b^2-a^2*c^2-a^2*d^2-a^2*e^2+2*a^2*t^2-b^2*d^2+b^2*e^2+c^2*d^2-c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2))*(-a^4+a^2*b^2+a^2*c^2+a^2*d^2+a^2*e^2-2*a^2*t^2+b^2*d^2-b^2*e^2-c^2*d^2+c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2))-1)/a/(1/(a^4-a^2*b^2-a^2*c^2-a^2*d^2-a^2*e^2+2*a^2*t^2-b^2*d^2+b^2*e^2+c^2*d^2-c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2))*(-a^4+a^2*b^2+a^2*c^2+a^2*d^2+a^2*e^2-2*a^2*t^2+b^2*d^2-b^2*e^2-c^2*d^2+c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2))+1), -(-a^4+2*a^2*d^2+2*a^2*e^2-d^4+2*d^2*e^2-e^4)^(1/2)/a/(a^4-a^2*b^2-a^2*c^2-a^2*d^2-a^2*e^2+2*a^2*t^2-b^2*d^2+b^2*e^2+c^2*d^2-c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2))*((a^4-a^2*b^2-a^2*c^2-a^2*d^2-a^2*e^2+2*a^2*t^2-b^2*d^2+b^2*e^2+c^2*d^2-c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2))*(-a^4+a^2*b^2+a^2*c^2+a^2*d^2+a^2*e^2-2*a^2*t^2+b^2*d^2-b^2*e^2-c^2*d^2+c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2)))^(1/2)/(1/(a^4-a^2*b^2-a^2*c^2-a^2*d^2-a^2*e^2+2*a^2*t^2-b^2*d^2+b^2*e^2+c^2*d^2-c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2))*(-a^4+a^2*b^2+a^2*c^2+a^2*d^2+a^2*e^2-2*a^2*t^2+b^2*d^2-b^2*e^2-c^2*d^2+c^2*e^2+(-(a+b+c)*(a-b+c)*(b-c+a)*(-c+a-b))^(1/2)*(-(a+d+e)*(a-d+e)*(d-e+a)*(-e+a-d))^(1/2))+1)]]);

#     return evala(convert(tetra, RootOf));
# end proc:


dofThroughVertex := proc(v, surf, S)
    local edgesToS;

    if v in S then
        return 0
    end if; 

    edgesToS := select(e -> evalb(v in e and (e[1] in S or e[2] in S)), Edges(surf));
    return max(0, 3 - nops(edgesToS))
end proc:





EmbeddingPlan := proc(surf, recursionDepths := 1, startFace := [])
    local S, Q, assemblyPlan, feasableVerts, criticalIndex, determinedVerts, undeterminedVerts, newAdjacentVertsInS, anchor, keyfunction;

    if startFace = [] then
        S := Faces(surf)[indmax(map(f -> nops(outgoingEdges(surf, f)), Faces(surf)))]; # embedded vertices in order
    else
        ASSERT({op(startFace)} in map(f -> {op(f)}, Faces(surf)), cat(startFace, " isn't a face of the provided surface."));
        S := startFace;
    end if;
    Q := []; # order in which to add Vertices to S

    assemblyPlan := [[S, [], []]]; # strategy to embed vertices. First entry of each entry is added vertices, second new dofs represented by one vertex, third is vertices which distances to solve for.

    # aux function to fill queue in step 1: returns the number of edges connecting v to the subset S of the vertices of surf.
    keyfunction := proc(v, surf, S, futureSteps := 1)
        local steps;

        steps := min(futureSteps, nops(Vertices(surf)) - nops(S) - 1);

        if futureSteps = 0 then
            return [dofThroughVertex(v, surf, S)]
        end if;    

        return [dofThroughVertex(v, surf, S), op(sort(map(w -> keyfunction(w, surf, [op(S), v], steps-1), [op({op(Vertices(surf))} minus {op(S)} minus {v})] ), key = (x -> add(x)))) ]
    end proc:

    while not(nops(S) = nops(Vertices(surf))) do
        # step 1: fill queue
        feasableVerts := select(v -> evalb(dofThroughVertex(v, surf, S) <= 1), [op({op(Vertices(surf))} minus {op(S)})]); # vertices that are connected to S via at least 2 edges
        feasableVerts := sort(feasableVerts, key=(v -> keyfunction(v, surf, S, recursionDepths))); # sort those vertices by cost
        
        criticalIndex := indwhere(v -> evalb(dofThroughVertex(v, surf, S) = 1), feasableVerts); # all vertices with lower index can be embedded without dofs, the ones with higher index can be embedded with 1 additional dof
        if criticalIndex = 1 then # all relevant vertices have 2 edges connecting them to S
            determinedVerts := []; # vertices that can be embedded without dofs
            undeterminedVerts := feasableVerts; # vertices that can be embedded with one dof
        elif criticalIndex = 0 then # all relevant vertices have 3 or more edges connecting them to S
            determinedVerts := feasableVerts;
            undeterminedVerts := [];
        else
            determinedVerts := feasableVerts[[$1..(criticalIndex - 1)]];
            undeterminedVerts := feasableVerts[[$(criticalIndex - 1)..nops(feasableVerts)]];
        end if;
        
        Q := [op(determinedVerts), op(Q), op(undeterminedVerts)];
        
        # step 2: add first element of Q to S
        ASSERT(nops(Q) >= 1, "Q is empty, but not all elements of surf have been embedded.");
        newAdjacentVertsInS := map(e -> op({op(e)} minus {Q[1]}), select(e -> evalb(Q[1] in e), outgoingEdges(surf, S))); # adjacent verts of Q[1] which are in S
        S := [op(S), Q[1]];

        # step 3.1: if added element is an undetermined vertex, a new edge needs to be added.
        if nops(newAdjacentVertsInS) = 2 then
            anchor := select(v -> evalb(not v in newAdjacentVertsInS), S[[$1..3]])[1];
            assemblyPlan := [op(assemblyPlan), [[Q[1]], [anchor], newAdjacentVertsInS]];

        # step 3.2: if added element has >=3 adjacent vertices in S, no new edges need to be added
        elif nops(newAdjacentVertsInS) >= 3 then
            assemblyPlan := [op(assemblyPlan), [[Q[1]], [], newAdjacentVertsInS]];
        end if;

        # step 4: remove first element of Q
        Q := subsop(1 = NULL, Q);

        # step 5: pop first element of Q until the first is not in S
        while Q[1] in S do
            Q := subsop(1 = NULL, Q);

            if Q = [] then 
                return assemblyPlan
            end if;
        end do;
    end do;


    return assemblyPlan;
end proc:


PasteTetrahedra := proc(coordinates1, coordinates2, pasteVertices1, pasteVertices2)
    local pair, pasteCoordinates1, pasteCoordinates2, B, C, EidB, EidC, CidE, isometry;

    # returns the coordinates of tetrahedron 2 that is pasted to tetrahedron 1. pasteVertices1 are pasted to pasteVertices2
    for pair in choose(3, 2) do
        ASSERT(dist(coordinates1[pasteVertices1[pair[1]]], coordinates1[pasteVertices1[pair[2]]]) = dist(coordinates2[pasteVertices2[pair[1]]], coordinates2[pasteVertices2[pair[2]]]), "Tetrahedra cannot be pasted.");
    end do;

    pasteCoordinates1 := [coordinates1[pasteVertices1[1]], coordinates1[pasteVertices1[2]], coordinates1[pasteVertices1[3]]]; # coordinates of face of tetrahedron 1 to be pasted
    pasteCoordinates2 := [coordinates2[pasteVertices2[1]], coordinates2[pasteVertices2[2]], coordinates2[pasteVertices2[3]]]; # coordinates of face of tetrahedron 2 to be pasted

    B := [pasteCoordinates1[2] - pasteCoordinates1[1], pasteCoordinates1[3]- pasteCoordinates1[1], CrossProduct(pasteCoordinates1[2] - pasteCoordinates1[1], pasteCoordinates1[3]- pasteCoordinates1[1])];
    C := [pasteCoordinates2[2] - pasteCoordinates2[1], pasteCoordinates2[3]- pasteCoordinates2[1], CrossProduct(pasteCoordinates2[2] - pasteCoordinates2[1], pasteCoordinates2[3]- pasteCoordinates2[1])];

    EidB := Matrix(B)^+;
    EidC := Matrix(C)^+;
    CidE := CidE^(-1);

    # isometry to apply to coordinates of tetrahedron 2 for the tetras to align on pasteVertices1 and pasteVertices2 resp.
    isometry := x -> convert((EidB.CidE.Vector(x)), list); 

    return map(v -> isometry(v), coordinates2);
end proc:


MoveTetraBase := proc(coordinatesTriangle, coordinatesTetrahedron, sign := 1)
    # proc that returns the coordinates of a the tetrahedron defined by coordinatesTetrahedron, which was moved so that the first 3 coordinates are the coordinates in corrdinatedTriangle
    local pair, B, C, EidB, EidC, CidE, isometry;

    for pair in choose(3, 2) do
        ASSERT(dist(coordinatesTriangle[pair[1]], coordinatesTriangle[pair[2]]) = dist(coordinatesTetrahedron[pair[1]], coordinatesTetrahedron[pair[2]]), "New base coordinates are not valid.");
    end do;

    B := [coordinatesTriangle[2] - coordinatesTriangle[1], coordinatesTriangle[3] - coordinatesTriangle[1], sign* convert(CrossProduct(coordinatesTriangle[2] - coordinatesTriangle[1], coordinatesTriangle[3] - coordinatesTriangle[1]), list)];
    C := [coordinatesTetrahedron[2] - coordinatesTetrahedron[1], coordinatesTetrahedron[3] - coordinatesTetrahedron[1], convert(CrossProduct(coordinatesTetrahedron[2] - coordinatesTetrahedron[1], coordinatesTetrahedron[3] - coordinatesTetrahedron[1]), list)];


    EidB := Matrix(3, 3, B)^+;
    EidC := Matrix(3, 3, C)^+;
    CidE := EidC^(-1);

    # isometry to apply to coordinates of tetrahedron 2 for the tetras to align on pasteVertices1 and pasteVertices2 resp.
    isometry := x -> convert((EidB.CidE.Vector(3, x - coordinatesTetrahedron[1]), list)) + coordinatesTriangle[1]; 

    return evala(map(v -> isometry(v), coordinatesTetrahedron));
end proc:


###### Backup
# FindEmbeddings := proc(surf, elng::list := [1$nops(Edges(surf))])
#     # find all vertex faithful embeddings of the simplicial surface surf with edge lengths defined by elng.

#     if nops(Vertices(surf)) = 3 then 
#         return Triangle(elng);
#     end if;

#     if nops(Vertices(surf)) = 4 then 
#         return Tetrahedron(elng);
#     end if;

#     edgesSets := map(e -> {op(e)}, Edges(surf)); # edges of surf as sets

#     embeddingPlan := EmbeddingPlan(surf);
#     nvars := 0; # aux variable to name variables

#     embeddings := [[[0,0,0] $ nops(Vertices(surf))]]; # list of embedded vertex coordinates

#     edgeLengthTetrahedron := proc(tetra::list, coords, step)
#         # calculates the edge lengths that need to be passed to the Tetrahedron proc.
#         # the tetrahedron is given through tetra: the first three entries are the indices of the base of the tetrahedron in Vertices(surf), the 4th is the tetrahedron tip. 
#         local l := [0$6];
#         local varsthisstep := nvars;
#         for local enum in Enumerate([[tetra[1], tetra[2]], [tetra[2], tetra[3]], [tetra[3], tetra[1]], [tetra[1], tetra[4]], [tetra[2], tetra[4]], [tetra[3], tetra[4]]]) do
#             local iter := enum[1];
#             local e := enum[2];
#             local edgeIndex := Search({op(e)}, edgesSets);
            
#             if edgeIndex = 0 then
#                 if e[1] in [op(step[3]), op(step[2])] and e[2] in [op(step[3]), op(step[2])] then
#                     l[iter] := dist(coords[e[1]], coords[e[2]]);
#                 else
#                     varsthisstep := varsthisstep + 1;
#                     l[iter] := _t[varsthisstep];
#                 end if;
#             else
#                 l[iter] := elng[edgeIndex];
#             end if;
#         end do;

#         return l;
#     end proc;

#     # embed first 4 vertices
#     step := embeddingPlan[2];
#     base := [op(step[2]), op(step[3])];
#     tip := step[1][1];
#     l := edgeLengthTetrahedron([op(base), tip], embeddings[1], step);
#     tetra := Tetrahedron(l);

#     embeddings[1][base[1]] := tetra[1];
#     embeddings[1][base[2]] := tetra[2];
#     embeddings[1][base[3]] := tetra[3];
#     embeddings[1][tip] := tetra[4];

#     nvars := nvars + nops(step[2]);

#     performStep := proc(step, coords, m)
#         # proc to perform one step of the embeddingplan with starting embedding defined by coords

#         # no equations need to be solved
#         if nops(step[3]) <= 3 then
#             base := [op(step[2]), op(step[3])];
#             tip := step[1][1];
#             l := evala(edgeLengthTetrahedron([op(base), tip], coords, step));
#             try
#                 tetra := Tetrahedron(l);
#             catch: error "FAIL"
#             end try;

#             for sign in [1, -1] do
#                 try
#                     coordsNew := evala(subsop(tip = MoveTetraBase(map(v -> coords[v], base), tetra, sign)[4], coords));
#                     Threads:-Mutex:-Lock(m);
#                     Append(newEmbeddings, coordsNew);
#                     Threads:-Mutex:-Unlock(m);
#                 end try;
#             end do;
#         # need to solve equations
#         elif nops(step[3]) > 3 then
#             base := step[3][1..3];
#             tip := step[1][1];
#             l := evala(edgeLengthTetrahedron([op(base), tip], coords, step));
#             try
#                 tetra := Tetrahedron(l);
#             catch: error "FAIL"
#             end try;

#             for sign in [1, -1] do
#                 try
#                     coordsNew := evala(subsop(tip = MoveTetraBase(map(v -> coords[v], base), tetra, sign)[4], coords));
                    
#                     # solve equations
#                     eqn := [];
#                     for v in step[3][4..nops(step[3])] do
#                         e := {tip, v};
#                         eqn := [op(eqn), evala(dist(coordsNew[tip], coordsNew[v])^2 - elng[Search(e, edgesSets)]^2)];
#                     end do;
#                     numer_eqn := evala(numer(eqn));
#                     denom_eqn := evala(denom(eqn));
#                     # print(nops(eqn));
#                     print(cat("numer_eqn = ", numer_eqn));
#                     print(cat("denom_eqn = ", denom_eqn));
#                     # print(cat("extension = ", indets(numer_eqn, algext)));
#                     # rel, inv_su, v := `SimplicialSurfaceEmbeddings/defining_ideal_products`([op(map(expr -> [expr], indets(numer_eqn, algext)))]);
#                     ext := [op(indets(numer_eqn, algext))];
#                     ext := sort(ext, key = (expr -> nops(indets(expr, algext))));
#                     print(cat("ext = ", [op(map(expr -> [expr], ext))]));

#                     s, u := `SimplicialSurfaceEmbeddings/solve_polynomial_system`(numer_eqn, map(i -> _t[i], [$1..nvars]), denom_eqn, [op(map(expr -> [expr], ext))]);
#                     print(cat("s = ", s));

#                     if eqn = [0$nops(eqn)] then
#                         s := [[]];
#                     end if;

#                     for sub in s do
#                         try
#                             # print(evala(Simplify~(subs(sub, coordsNew))));
#                             coordsNew1 := evala(subs(sub, coordsNew));
#                             # print(evala(dist(coordsNew1[tip], coordsNew1[v])^2 - elng[Search(e, edgesSets)]^2) assuming real);
#                             # print(coordsNew1);
#                             if not evala(dist(coordsNew1[tip], coordsNew1[v])^2 - elng[Search(e, edgesSets)]^2) = 0 then
#                                 next sub;
#                             end if;
#                             Threads:-Mutex:-Lock(m);
#                             Append(newEmbeddings, coordsNew1);
#                             Threads:-Mutex:-Unlock(m);
#                         # catch: next sub;
#                         end try;
#                     end do;
#                 # catch: next sign;
#                 end try;
#             end do;

#         end if;
#     end proc;

#     # other vertices
#     for step in embeddingPlan[3..nops(embeddingPlan)] do
#         print(cat("number embeddings = ", nops(embeddings)));
#         print(step);
#         newEmbeddings := Array([]);
#         nvars := nvars + nops(step[2]);

#         m := Threads:-Mutex:-Create();

#         for coords in embeddings do
#             try
#                 performStep(step, coords, m);
#             catch "FAIL": next coords;
#             end try;
#         end do;
#         # Threads:-Task:-Start(null, Task=[performStep, step, embeddings[1], m]);
#         # Threads:-Task:-Start(null, map(coords -> Task=[performStep, step, coords, m], embeddings));
        
#         embeddings := MakeUnique(convert(newEmbeddings, list));
#     end do;

#     return embeddings
# end proc:


FindEmbeddings := proc(surf, recursionDepth := 1, elng::list := [1$nops(Edges(surf))])
    # find all vertex faithful embeddings of the simplicial surface surf with edge lengths defined by elng.

    if nops(Vertices(surf)) = 3 then 
        return Triangle(elng);
    end if;

    if nops(Vertices(surf)) = 4 then 
        return Tetrahedron(elng);
    end if;

    edgesSets := map(e -> {op(e)}, Edges(surf)); # edges of surf as sets

    embeddingPlan := EmbeddingPlan(surf, recursionDepth);
    print(embeddingPlan);
    nvars := 0; # aux variable to name variables

    embeddings := [[[0,0,0] $ nops(Vertices(surf))]]; # list of embedded vertex coordinates

    edgeLengthTetrahedron := proc(tetra::list, coords, step, nvars)
        # calculates the edge lengths that need to be passed to the Tetrahedron proc.
        # the tetrahedron is given through tetra: the first three entries are the indices of the base of the tetrahedron in Vertices(surf), the 4th is the tetrahedron tip. 
        local l := [0$6];
        local varsthisstep := nvars;
        for local enum in Enumerate([[tetra[1], tetra[2]], [tetra[2], tetra[3]], [tetra[3], tetra[1]], [tetra[1], tetra[4]], [tetra[2], tetra[4]], [tetra[3], tetra[4]]]) do
            local iter := enum[1];
            local e := enum[2];
            local edgeIndex := Search({op(e)}, edgesSets);
            
            if edgeIndex = 0 then
                if e[1] in [op(step[3]), op(step[2])] and e[2] in [op(step[3]), op(step[2])] then
                    l[iter] := dist(coords[e[1]], coords[e[2]]);
                else
                    varsthisstep := varsthisstep + 1;
                    l[iter] := _t[varsthisstep];
                end if;
            else
                l[iter] := elng[edgeIndex];
            end if;
        end do;

        return l;
    end proc;

    # embed first 4 vertices
    step := embeddingPlan[2];
    base := [op(step[2]), op(step[3])];
    tip := step[1][1];
    l := edgeLengthTetrahedron([op(base), tip], embeddings[1], step, nvars);
    tetra := Tetrahedron(l);

    embeddings[1][base[1]] := tetra[1];
    embeddings[1][base[2]] := tetra[2];
    embeddings[1][base[3]] := tetra[3];
    embeddings[1][tip] := tetra[4];

    nvars := nvars + nops(step[2]);

    performStep := proc(step, coords, sign, nvars)
        # proc to perform one step of the embeddingplan with starting embedding defined by coords; sign determines, in which half space the tip of the attached tetrahedron lies.

        newEmbeddings := [];
        
        # no equations need to be solved
        if nops(step[3]) <= 3 then
            base := [op(step[2]), op(step[3])];
            tip := step[1][1];
            l := evala(edgeLengthTetrahedron([op(base), tip], coords, step, nvars));
            try
                tetra := Tetrahedron(l);
            catch: error "FAIL"
            end try;

            try
                coordsNew := evala(subsop(tip = MoveTetraBase(map(v -> coords[v], base), tetra, sign)[4], coords));
                newEmbeddings := [op(newEmbeddings, coordsNew)];
            end try;

        # need to solve equations
        elif nops(step[3]) > 3 then
            base := step[3][1..3];
            tip := step[1][1];
            l := evala(edgeLengthTetrahedron([op(base), tip], coords, step, navrs));
            try
                tetra := Tetrahedron(l);
            catch: error "FAIL"
            end try;

            try
                coordsNew := evala(subsop(tip = MoveTetraBase(map(v -> coords[v], base), tetra, sign)[4], coords));
                
                # solve equations
                eqn := [];
                for v in step[3][4..nops(step[3])] do
                    e := {tip, v};
                    eqn := [op(eqn), evala(dist(coordsNew[tip], coordsNew[v])^2 - elng[Search(e, edgesSets)]^2)];
                end do;
                numer_eqn := evala(numer(eqn));
                denom_eqn := evala(denom(eqn));
                # print(nops(eqn));
                # print(cat("numer_eqn = ", numer_eqn));
                # print(cat("denom_eqn = ", denom_eqn));
                # print(cat("extension = ", indets(numer_eqn, algext)));
                # rel, inv_su, v := `SimplicialSurfaceEmbeddings/defining_ideal_products`([op(map(expr -> [expr], indets(numer_eqn, algext)))]);
                ext := [op(indets(numer_eqn, algext))];
                ext := sort(ext, key = (expr -> nops(indets(expr, algext))));
                # print(cat("ext = ", [op(map(expr -> [expr], ext))]));

                s, u := `SimplicialSurfaceEmbeddings/solve_polynomial_system`(numer_eqn, map(i -> _t[i], [$1..nvars]), denom_eqn, [op(map(expr -> [expr], ext))]);
                # print(cat("s = ", s));

                if eqn = [0$nops(eqn)] then
                    s := [[]];
                end if;

                for sub in s do
                    try
                        # print(evala(Simplify~(subs(sub, coordsNew))));
                        coordsNew1 := evala(subs(sub, coordsNew));
                        # print(evala(dist(coordsNew1[tip], coordsNew1[v])^2 - elng[Search(e, edgesSets)]^2) assuming real);
                        # print(coordsNew1);
                        if not evala(dist(coordsNew1[tip], coordsNew1[v])^2 - elng[Search(e, edgesSets)]^2) = 0 then
                            next sub;
                        end if;
                        newEmbeddings := [op(newEmbeddings), coordsNew1];
                    # catch: next sub;
                    end try;
                end do;
            # catch: next sign;
            end try;

        end if;

        return newEmbeddings;
    end proc;

    local Server := proc(step, embeddings, nvars)
        # task server for grid parallel programming. https://de.maplesoft.com/support/help/maple/view.aspx?path=ProgrammingGuide/Chapter15
        local msg;

        global updatedEmbeddings := [];

        for local coordinates in embeddings do
            for local sign in [1,-1] do
                # get a work request
                msg := Grid:-Receive();
                i := msg[1];
                # send out work
                Grid:-Send(i, [step, coordinates, sign, nvars]);

                # if we receive new embeddings, store them
                if msg[2] <> [] then
                    updatedEmbeddings := [op(updatedEmbeddings), op(msg[2])];
                end if;
            end do;
        end do;

        # print("all data has been sent");
        # all data has been sent. Receive last messages
        for j from 1 to Grid:-NumNodes()-1 do
            msg := Grid:-Receive();
            if msg[2] <> [] then
                    updatedEmbeddings := [op(updatedEmbeddings), op(msg[2])];
            end if;
            # send terminate message to the nodes
            i := msg[1];
            Grid:-Send(i, -1);
        end do;

        return updatedEmbeddings;
    end proc;

    local Client := proc(i)
        local newEmbeddings;
        # send initial data request
        Grid:-Send(0, [i, []]);

        do
            # wait for reply
            # printf("Node %d of %d is waiting for reply.\n", i, Grid:-NumNodes()-1);
            local msg := Grid:-Receive(0);

            # if terminate message, exit the loop
            if msg = -1 then 
                # printf("Shutting down Node %d of %d.\n", i, Grid:-NumNodes()-1);
                break;
            end if;

            # perform the step given the parameters in msg
            # printf("Node %d of %d is performing calculations.\n", i, Grid:-NumNodes()-1);
            local step := msg[1];
            local coordinates := msg[2];
            local sign := msg[3];
            local nvars := msg[4];

            try
                newEmbeddings := performStep(step, coordinates, sign, nvars);
            catch: 
                newEmbeddings := [];
            end try;

            # send solution 
            Grid:-Send(0, [i, newEmbeddings]);
            # printf("Node %d of %d sent results.\n", i, Grid:-NumNodes()-1);
        end do;

        return NULL
    end proc;

    local GridFunction := proc()
        with(SimplicialSurfaceEmbeddings):
        local i := Grid:-MyNode();

        if i = 0 then
            Server(step, embeddings, nvars);
        else
            Client(i);
        end if;
    end proc;

    # other vertices
    for step in embeddingPlan[3..nops(embeddingPlan)] do
        print(cat("number embeddings = ", nops(embeddings)));
        print(step);
        nvars := nvars + nops(step[2]);

        updatedEmbeddings := Grid:-Launch(GridFunction, 
                    numnodes=min(13, nops(embeddings)*2+1), 
                    imports=['performStep', 'Server', 'Client', 'Tetrahedron', 'MoveTetraBase', 'edgeLengthTetrahedron', 'dist', 
                            'edgesSets', 'step', 'embeddings', 'nvars']);
                            
        embeddings := MakeUnique(updatedEmbeddings);
        if nops(embeddings) = 0 then
            return []
        end if;
    end do;

    return embeddings
end proc: