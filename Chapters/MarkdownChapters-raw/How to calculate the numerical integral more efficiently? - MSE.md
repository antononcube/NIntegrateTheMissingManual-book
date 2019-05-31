## General 

When given an array of integrands `NIntegrate` is run separately over each array element. That is not necessary though, the core `NIntegrate` integration strategies can work with any integrands as long the error estimates are real numbers. With the array-of-functions integration rule used below around 10 to 100 times speed-up is achieved.

The motivation for implementing [`ArrayOfFunctionsRule`](https://github.com/antononcube/MathematicaForPrediction/blob/master/Misc/ArrayOfFunctionsRule.m) is to provide a significant speed-up for integrands that are arrays of functions. That is achived by evaluating all functions with the same integration rule abscissas and weights.

As mentioned in the comments these posts/answers: ["NIntegrate over a list of functions"](http://mathematica.stackexchange.com/a/81436/34008), ["How to avoid repetitive calculation when doing numerical integral?"](http://mathematica.stackexchange.com/a/120737/34008). 

For more details how rules like `ArrayOfFunctionsRule` are implemented see ["How to implement custom integration rules for use by NIntegrate?"](http://mathematica.stackexchange.com/q/118324/34008).

## Performance comparison 

The code below uses the definitions in the question.

Load the code for integration rule `ArrayOfFunctionsRule`:

    Import["https://raw.githubusercontent.com/antononcube/\
    MathematicaForPrediction/master/Misc/ArrayOfFunctionsRule.m"]

Convert the integrand matrix into matrices of functions:

    I1exprFS = 
      Map[Function[{fx}, Function[Evaluate[fx /. x -> #]]], I1expr, {2}];
   
Call `NIntegrate` with the new rule:

    res1 = Map[
        NIntegrate[1, {x, 0, 1}, 
          Method -> {"GlobalAdaptive", "SingularityHandler" -> None, 
            Method -> {ArrayOfFunctionsRule, "Functions" -> #1}}] &, 
        I1exprFS]; // AbsoluteTiming

    (* {0.058653, Null} *)

Compare with the standard `NIntegrate` call:

    I1 = NIntegrate[I1expr, {x, 0, 1}]; // AbsoluteTiming    

    (* {0.503172, Null} *)

We see that `ArrayOfFunctionsRule` provides 10 times speed-up (for the functions defined in the question.)   

Verify agreement of the results:

    Norm[res1 - I1, 2]

    (* 4.39568*10^-7 *)

See the options of `ArrayOfFunctionsRule`. With the option "ErrorsNormFunction" different norms can be used to compute the integration errors.

Note that the rule has to be used with a strategy specification that has the option "SingularityHandler" -> None, and that the rule does not work correctly with ranges that have infinity.

### Larger matrices

This makes a larger matrix (based on the integrands in the question):

    funcsExpr = I1expr;
    funcsExpr = Table[funcsExpr, {100}, {10}];
    funcsExpr = Partition[Flatten[funcsExpr], 100];
    Dimensions[funcsExpr]    
    (* {490, 100} *)

Converting to functions:

    funcs = 
      Map[Function[{fx}, Function[Evaluate[fx /. x -> #]]], 
       funcsExpr, {2}];

This integration is ~ 40 times faster than the default options one:

    res1 = 
       Map[NIntegrate[1, {x, 0, 1}, 
          Method -> {"GlobalAdaptive", "SingularityHandler" -> None, 
            Method -> {ArrayOfFunctionsRule, "Functions" -> #1}}] &, 
        funcs]; // AbsoluteTiming
    
    (* {13.8379, Null} *)

Integration with `NIntegrate`'s default method:

    I1 = NIntegrate[funcsExpr, {x, 0, 1}]; // AbsoluteTiming
    
    (* {538.901, Null} *)

Adherence verification:

    Norm[res1 - I1, 2]    
    (* 6.31739*10^-6 *)

### Parallel computation

Using `ParallelMap` on 4 processors we get ~100 times speed-up.

    pres1 = 
       ParallelMap[
        NIntegrate[1, {x, 0, 1}, 
          Method -> {"GlobalAdaptive", "SingularityHandler" -> None, 
            Method -> {ArrayOfFunctionsRule, "Functions" -> #1}}] &, 
        funcs, DistributedContexts -> All]; // AbsoluteTiming
    (* {5.24539, Null} *)
    
    Norm[pres1 - I1, 2]  
    (* 6.31739*10^-6 *)
