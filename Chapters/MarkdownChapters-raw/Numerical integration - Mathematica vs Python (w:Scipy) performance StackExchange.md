**General comments**

First, if you plan to use multi-dimensional integrals it is better to test with multi-dimensional integrals not with one dimensional ones. One might think that the test in the question is an appropriate one if multi-dimensional integration is done by the integrator in a recursive manner. This seems to be case for scipy.integrate.nquad (see [scipy.integrate.nquad.html][1]), but it is not for NIntegrate. NIntegrate constructs and utilizes proper multi-dimensional integration rules and/or strategies.

Second, I do not think this is a test from which we can make general conclusions for the speed of a numerical integrator. The integral is too specific: an odd function over (-Infinity, Infinity). (Evaluates to zero.) I assume it is chosen with the specific research to be undertaken in mind.

Third, for very high dimensions the more useful integration strategies are (quite) different than the useful integration strategies in low dimensions. The precision and accuracy goals sought after are much smaller. These observations make the selected test less relevant.

Fourth, NIntegrate plays very well with the optimization functions in *Mathematica*. I would assume you would be better off using *Mathematica* than Python, but I do not have much experience with NumPy and SciPy.

**More technically**

it is better to call the integration routine multiple times in order to get a better timing estimate. I wanted to modify and run both the *Mathematica* and Python tests like this but I found the installation of NymPy and SciPy to be too much work. For example:

    def intfun3(ntimes): 
    mu = 0 
    sig = 1 

    def npdf(x): 
        return 1.0 / np.sqrt( 2 * np.pi * sig**2 ) * np.exp( - (x - mu)**2 / (2 * sig**2) ) 

    tic = time.time()
    for i in range(1,ntimes) :
        out = scipy.integrate.quad( lambda x : x * npdf(x), -np.inf, np.inf)
    toc = time.time() 

    print out 
    print toc - tic 

    intfun3(1000)

 
We can get NIntegrate to do the test around 5-6 times faster (on my laptop with *Mathematica* 10.2) by providing options settings that correspond to the default integration parameters arguments of scipy.integrate.quad. (I have read the descriptions of the parameters in [scipy.integrate.quad.html][2] ). 

Here are the original and the modified tests:
   

    testpdf[\[Mu]_, \[Sigma]_, x_] := 1/Sqrt[2*Pi*`[Sigma]^2]*Exp[-((x - \[Mu])^2/(2*\[Sigma]^2))];`

    n = 1000;
    res = Do[NIntegrate[
         x*testpdf[0, 1, x], {x, -Infinity, Infinity}], {n}] // 
       AbsoluteTiming;
    res[[1]]/n

    (* Out[521]= 0.00485096 *)

    n = 1000;
    res = Do[NIntegrate[x*testpdf[0, 1, x], {x, -Infinity, Infinity}, 
         PrecisionGoal -> 8, AccuracyGoal -> 8, 
         Method -> {"DoubleExponential", 
           "SymbolicProcessing" -> 0}], {n}] // AbsoluteTiming;
    res[[1]]/n

    (* Out[533]= 0.00090782 *)

Using the option "SymbolicProcessing"->0 prevents NIntegrate to do symbolic preprocessing. (See ["SymbolicProcessing"][3].) For the integral we are discussing, with the default option settings NIntegrate detects it is an odd function over (-Infinity,Infinity) and integrates only over (0,Infinity) as a numerical check. See ["EvenOddSubdivision"][4]

The settings PrecisionGoal->8, AccuracyGoal->8 correspond to "epsabs=1.49e-08, epsrel=1.49e-08" in scipy.integrate.quad.html . Using the method ["DoubleExponential"][5] corresponds to the description "If one of the integration limits is infinite, then a Fourier integral is computed[...]" in [scipy.integrate.quad.html][2] .

Note that when using the option "SymbolicProcessing"->0 NIntegrates gives warnings that the integral does not converge quickly enough with the message:

*NIntegrate::slwcon :  "Numerical integration converging too slowly; suspect one of the following: singularity, value of the integration is 0, highly oscillatory integrand, or WorkingPrecision too small. "*


  [1]: http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.nquad.html
  [2]: http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
  [3]: https://reference.wolfram.com/language/tutorial/NIntegrateIntegrationStrategies.html#188031681
  [4]: https://reference.wolfram.com/language/tutorial/NIntegrateIntegrationStrategies.html#172715735
  [5]: https://reference.wolfram.com/language/tutorial/NIntegrateIntegrationStrategies.html#526196975
