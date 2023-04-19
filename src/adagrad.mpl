read("embedding.mpl");

# adapted gradient descent
adagrad := proc(objective, gradient, startingPoint, {stepSize:= 0.5, maxSteps:=100000, eps := 1e-8, tolerance := 1e-30})
    local nvars, squareGradSums, solution, step, grad, alpha;

    Digits := Digits + 20;

    nvars := nops(startingPoint);
    squareGradSums := [0$nvars];

    solution := startingPoint;

    for step in $1..maxSteps do
        if step mod 100 = 0 then
            print(abs(objective(solution)) );
        end if;
        
        if abs(objective(solution)) < tolerance then
            print(step);
            return solution
        end if;

        grad := gradient(solution);
        squareGradSums := squareGradSums + grad^~2;
        alpha := stepSize /~ (sqrt~(squareGradSums) +~ eps);
        solution := solution - (alpha *~ grad);
    end do;

    return solution;
end proc:

randadagrad := proc(objective, gradient, lowerBounds, upperBounds, {stepSize:= 0.5, maxSteps := 10000, eps := 1e-8, tolerance := 1e-12})
    local startingPoint;

    ASSERT(nops(lowerBounds) = nops(upperBounds), "lowerBounds and upperBounds need to have the same number of entries");

    startingPoint := map(i -> rand(evalf(lowerBounds[i])..evalf(upperBounds[i]))(), [$1..nops(lowerBounds)]);
    return adagrad(objective, gradient, startingPoint, stepSize, maxSteps, eps, tolerance)
end proc:

# objective := proc(x,y)
#     return x^2 + y^2;
# end proc:

# gradient := proc(x,y)
#     return [2*x, 2*y];
# end proc: 

SimplSurfEmbeddingByAdagrad := proc(surf, edgelengths);
    local objective, gradDistConstraint, gradient;

    ASSERT(nops(edgeLengths) = nops(Edges(surf)), "every edge has to have a length attached to it.");

    objective := proc(coordinates)
        local sol, enum, i, e, v, w;

        sol := 0;
        for enum in Enumerate(surf[edges]) do
            # add distance conditions
            i := enum[1]; # index
            e := enum[2]; # edge

            v := coordinates[((e[1]-1)*3 + 1)..(e[1]*3)]; # first vertex of edge e
            w := coordinates[((e[2]-1)*3 + 1)..(e[2]*3)]; # second vertex of edge e
            sol := sol + evalf((dist(v,w) - edgelengths[i]))^2 # ^2 because then the only possibility for the sum to be zero over the reals is if every summand is zero
        end do;

        # factor out rigid motions
        sol := sol + (dist(coordinates[1..3], [0,0,0]))^2;
        sol := sol + (dist(coordinates[4..5], [0,0]))^2;
        sol := sol + (dist([coordinates[7]], [0]))^2;

        return sol;
    end proc:

    gradDistConstraint := map(var -> diff((dist([x,y,z], [a,b,c]) - r)^2, var), [x,y,z,a,b,c]); # gradient of (dist([x,y,z], [a,b,c]) - r)^2

    gradient := proc(coordinates)
        local grad, enum, i, e, v, w, gradThisEdge;

        grad := [0$nops(coordinates)];

        # gradients of distance constraints
        for enum in Enumerate(surf[edges]) do
            i := enum[1]; # index
            e := enum[2]; # edge

            v := coordinates[((e[1]-1)*3 + 1)..(e[1]*3)]; # first vertex of edge e
            w := coordinates[((e[2]-1)*3 + 1)..(e[2]*3)]; # second vertex of edge e

            gradThisEdge := evalf(subs([x = v[1], y = v[2], z =v[3], a = w[1], b = w[2], c = w[3], r = edgelengths[i]], gradDistConstraint));
            gradThisEdge := [0$((e[1]-1)*3), op(gradThisEdge[1..3]), 0$(e[2]-e[1]-1)*3, op(gradThisEdge[4..6]), 0$((nops(Vertices(surf))-e[2])*3)];
            grad := grad + gradThisEdge;
        end do;

        # constraints that arise from eliminating rigid motions
        grad[1] := grad[1] + 2*coordinates[1];
        grad[2] := grad[2] + 2*coordinates[2];
        grad[3] := grad[3] + 2*coordinates[3];
        grad[4] := grad[4] + 2*coordinates[4];
        grad[5] := grad[5] + 2*coordinates[5];
        grad[7] := grad[7] + 2*coordinates[7];

        return grad;
    end proc:

    return randadagrad(objective, gradient, [0$(3*nops(Vertices(surf)))], [3$(3*nops(Vertices(surf)))])
end proc: