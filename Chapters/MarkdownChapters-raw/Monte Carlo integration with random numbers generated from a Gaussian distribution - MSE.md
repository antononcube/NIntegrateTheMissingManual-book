Because of the comment by @ciao I think it is a good idea to give a solution/answer using `MonteCarloRule` and `PointGenerator`.

`PointGenerator` is an object that expects a sub-value function definition for "MakePoint" with the signature:

    "MakePoint"[dim_, totalNumberOfPoints_, i_, wprec_]

where the arguments are: 

1. dim - integration dimension, 
2. totalNumberOfPoints - desired total number of points per rule application,
3. i - point index,
4. wprec - working precision.

Note that "MakePoint" provides very fine granularity, it makes only one point. For optimization purposes custom definitions (like the one given below) have to provide some sort of optimization for generating a bulk of points. The consequence of such an approach is that if `AdaptiveMonteCarlo` is used then the sampling points of each application of `MonteCarloRule` are going to be rescaled versions of the same set of points.

Here are definitions that make `PointGenerator` work with a general distribution argument.

    ClearAll[MCAbscGenerator, MCDistributionPoints];
    MCDistributionPoints[distr_, dim_, totalNumberOfPoints_] :=
      MCDistributionPoints[distr, dim, totalNumberOfPoints] =
       Block[{points = RandomVariate[distr, {totalNumberOfPoints, dim}]},
        Transpose[
         Map[Rescale[#, {Min[#], Max[#]}, {0, 1}] &, 
          Transpose[points]]]];
    MCAbscGenerator["MakePoint"[dim_, totalNumberOfPoints_, n_, wprec_]] :=
        MCDistributionPoints[NormalDistribution[], dim, 
        totalNumberOfPoints][[n]];
    MCAbscGenerator[distr_][
      "MakePoint"[dim_, totalNumberOfPoints_, i_, wprec_]] :=
     MCDistributionPoints[distr, dim, totalNumberOfPoints][[i]];

(The definitions are probably not very didactic, but I think they provide good general properties: signatures and funcitonality.)

Here is how we apply the custom point generator with `MonteCarloRule`:

    NIntegrate[x^2 + y^3, {x, 0, 1}, {y, 0, 1}, PrecisionGoal -> 2,
       Method -> {"MonteCarlo", 
       Method -> {"MonteCarloRule", 
         "PointGenerator" -> 
          MCAbscGenerator[TruncatedDistribution[{-1, 1}, NormalDistribution[]]]}}]

    (* 0.561876 *)

The real value is `7/12`:

    Integrate[x^2 + y^3, {x, 0, 1}, {y, 0, 1}] // N
    
    (* 0.583333 *)

Let us visualize the sampling points for different truncated versions of the Normal Distribution and different Monte-Carlo methods. Notice that the integral gets more underestimated for the wider truncated distribution ranges.

    Needs["Integration`NIntegrateUtilities`"]
    
    ds = {TruncatedDistribution[{-1, 1}, NormalDistribution[]], 
       TruncatedDistribution[{-2, 2}, NormalDistribution[]], 
       TruncatedDistribution[{-3, 3}, NormalDistribution[]], NormalDistribution[]};
    Grid[Riffle[Table[{Style[distr, Blue], SpanFromLeft}, {distr, ds}], 
      Table[(est =
         NIntegrate[x^2 + y^3, {x, 0, 1}, {y, 0, 1}, PrecisionGoal -> 2, 
          Method -> {meth, 
            Method -> {"MonteCarloRule", 
              "PointGenerator" -> MCAbscGenerator[distr]}}];
        gr = NIntegrateSamplingPoints@
          NIntegrate[x^2 + y^3, {x, 0, 1}, {y, 0, 1}, PrecisionGoal -> 2, 
           Method -> {meth, 
             Method -> {"MonteCarloRule", 
               "PointGenerator" -> MCAbscGenerator[distr]}}];
        Labeled[Append[gr, Frame -> True], 
         Grid[{{"method:", meth}, {"estimate:", est}}]]
        ), {distr, ds}, {meth, {"MonteCarlo", "AdaptiveMonteCarlo"}}]], 
     Dividers -> {None, {{True, False}}}]

[![enter image description here][1]][1]


  [1]: http://i.stack.imgur.com/IyMX6.png
