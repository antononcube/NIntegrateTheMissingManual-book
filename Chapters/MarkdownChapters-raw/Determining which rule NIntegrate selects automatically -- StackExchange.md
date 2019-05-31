This question comes up often enough. See this discussion at [community.wolfram.com][1] : [Integration method used in NIntegrate][2] , and the notebook [Finding the applied NIntegrate methods][3] attached to my second response in the discussion.

That notebook contains examples of usage of the undocumented function `NIntegrateSamplingPoints` and `NIntegrate`'s option `IntegrationMonitor`.

**The integral in the question**

For the integral in the question with `NIntegrateSamplingPoints` we get kind of a boring picture because of the infinite region. Taking logs of the sampling points might be more informative:

    Needs["Integration`NIntegrateUtilities`"]

    gr = NIntegrateSamplingPoints@
       NIntegrate[Sin[x]^2 Sin[1000 x]^2/x^(5/2), {x, 0, Infinity}, 
        Method -> {"LevinRule", 
          Method -> {"GaussKronrodRule", "Points" -> 11}}];

    Graphics[gr[[1]] /. 
      Point[{x_?NumericQ, y_?NumericQ}] :> 
       Point[{Log[10, x + 10^-12], y}], Frame -> True, 
     FrameLabel -> {"lg(sampling points)", "evaluation order"}, AspectRatio -> 1/1.5]

[![enter image description here][4]][4]

The plot shows the order of evaluation of the sampling points.

Using `IntegrationMonitor` we can see the application of the integrand over the regions derived with the `LevinRule` method:

    t = Reap[NIntegrate[Sin[x]^2 Sin[1000 x]^2/x^(5/2), {x, 0, Infinity}, 
        Method -> {"LevinRule", 
          Method -> {"GaussKronrodRule", "Points" -> 11}}, 
        PrecisionGoal -> 1.5, 
        IntegrationMonitor -> (Sow[
            Map[{#1@"Integrand", #1@"Boundaries", #1@"Integral", #1@
                "Error"} &, #1]] &)]];
    res = t[[1]];
    t = t[[2, 1]];
    t

[![enter image description here][5]][5]


**UPDATE**

*(I have had the code below for some time, but I have been hesitant to share it because of several concerns. It is somewhat hard to interpret its results and coming up with it does require some internal knowledge of NIntegrate's development. After dedicated discussions with/about the NIntegrate method tracing code at [WTC 2015][6] it seems that it is better to show and describe it.)*

We can trace `NIntegrate`'s method initialization by manipulating the top level initialization function implementations. The basic idea is to take the down values and up values of the initialization functions of NIntegrate's methods that have the form 

    Block[{v___}, b_CompoundExpression]

and replace them with 

    Block[{res = Block[{v}, b]}, Print[res]; res]

When an `NIntegrate` command is executed we will see a printed trace of the methods that are initialized.

Here is the tracing code:

    symbNames = Names["NIntegrate`*"];
    symbNames = 
      Append[Pick[symbNames, 
        StringMatchQ[
         symbNames, (__ ~~ "Rule") | (__ ~~ 
            "Global" | "Local" | "MonteCarlo" | "Principal" | "Levin" | 
             "Osc" ~~ ___)]], "NIntegrate`AutomaticStrategy"];
    symbs = ToExpression[#] & /@ symbNames;
    dvs = DownValues /@ symbs;
    uvs = UpValues /@ symbs;
    Unprotect /@ symbs;
    dvsNew = MapThread[
       With[{s = #2},
         DownValues[s] = 
          ReplaceAll[#1, 
           HoldPattern[
             a_ :> b___] :> (a :> (Print["DownValue call for: ", 
                Style[s, Red]]; b))]] &, {dvs, symbs, symbNames}];
    uvsNew = MapThread[
       With[{s = #2},
         UpValues[s] =
          ReplaceAll[#1,
           HoldPattern[Block[vars_, CompoundExpression[b___]]] :>
            Block[{res = Block[vars, CompoundExpression[b]]}, 
             Print["UpValue call for: ", Style[s, Blue], 
              Style[" ::\n", Blue], res]; res]
           ]
         ] &, {uvs, symbs, symbNames}];

In that code using `Pick` I have reduced the number of `NIntegrate` context symbols being traced. Of course, if it is desired, the down/up values of the full list of the `NIntegrate` context symbols can be manipulated for tracing.

Let us look at tracing examples.

Here is a numerical integration computation with an automatically selected method:

[![enter image description here][7]][7]

The basic objects of `NIntegrate` are integration regions. Each region has its own integration function and integration rule. In order to interpret the printed trace it is helpful to know that `NIntegate` uses the [software design patterns][8] [Strategy][9], [Composite][10], [Decorator][11], and others. 

`NIntegrate`'s symbolic pre-processors are using Decorator, and we can see in the trace that the outcome of AutomaticStrategy are pre-processor symbols wrapped around the method "GlobalAdaptive". The method "GlobalAdaptive" uses a Gauss-Kronrod rule, which after its initialization is treated as a general rule. (I.e. a list of abscissas, integral estimate weights, and approximation error weights.) 

The method GlobalAdaptive is going to be directly used if symbolic processing is prevented:

[![enter image description here][12]][12]

Here is a numerical integration computation with a specifically selected pre-processing method:

[![enter image description here][13]][13]


**UPDATE 2 (`IntegrationMonitor` methods)**

*(Thanks to [Michael E2](http://mathematica.stackexchange.com/users/4999/michael-e2) for prompting these explanations.)*

Each integration strategy of `NIntegrate` creates and manipulates a collection of integration regions. Each integration region can have its own integrand and/or integration rule. `NIntegrate`'s main integration strategy "GlobalAdaptive" keeps the integration regions in a heap according to their error. The sum of the integral estimates of all regions make the global integral estimate. The sum of the integral errors make the global error. If the global error is larger than the desired tollerance "GlobalAdaptive" splits the region with the largest error estimate into two regions and applies the integration rule. If too many splittings have been done then a singularity handler is applied over the last region split.

At each step of an integration strategy the option `IntegrationMonitor` obtains as an argument the list of integration regions used in that step. Below is a table that shows methods that can be applied to each integration region in that list.

[![enter image description here][14]][14]

And here is (another) example of the application of those methods:

    iRegionMethods = {"Axis", "Boundaries", "Dimension", "Error", 
      "GetRule", "Integral", "Integrand", "WorkingPrecision"}; res = 
     Reap@NIntegrate[x^2 y^2, {x, 0, 4}, {y, 0, 2}, PrecisionGoal -> 1.1, 
       Method -> "AdaptiveMonteCarlo", 
       IntegrationMonitor :> 
        Function[{iregs}, 
         Sow[Association /@ 
           Transpose[
            Map[Thread[# -> Through[iregs[#]]] &, iRegionMethods]]]]];
    Dataset[Flatten[res[[2]]]]

[![enter image description here][15]][15]


  [1]: http://community.wolfram.com
  [2]: http://community.wolfram.com/groups/-/m/t/314124
  [3]: http://community.wolfram.com/groups?p_auth=nbRdQI9w&p_p_id=19&p_p_lifecycle=1&p_p_state=exclusive&p_p_mode=view&p_p_col_id=column-1&p_p_col_pos=1&p_p_col_count=3&_19_struts_action=%2Fmessage_boards%2Fget_message_attachment&_19_messageId=370695&_19_attachment=Finding%20the%20applied%20NIntegrate%20methods.nb
  [4]: http://i.stack.imgur.com/S2IUV.png
  [5]: http://i.stack.imgur.com/LWrw7.png
  [6]: http://www.wolfram.com/events/technology-conference/2015/
  [7]: http://i.stack.imgur.com/kVg91.png
  [8]: https://en.wikipedia.org/wiki/Software_design_pattern
  [9]: https://en.wikipedia.org/wiki/Strategy_pattern
  [10]: https://en.wikipedia.org/wiki/Composite_pattern
  [11]: https://en.wikipedia.org/wiki/Decorator_pattern
  [12]: http://i.stack.imgur.com/UC6yi.png
  [13]: http://i.stack.imgur.com/I2lhu.png
  [14]: http://i.stack.imgur.com/Fe1cY.png
  [15]: http://i.stack.imgur.com/h6GBd.png
