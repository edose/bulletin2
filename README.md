#  bulletin2
Demonstration code for an online replacement and extension of AAVSO's venerable LPV Bulletins.

 **Motivation:** 
 
 For all the 379 LPV (long-period variable) stars in Bulletin 2018, we want to predict the following:
 * minimum and maximum dates and their approximate V-mags during next calendar year,
 * approximate V-mags ON DEMAND, for any star and date, presumably from a future web tool, and
 * ideally, bonus predictions (dates, V-mags) for the entire next year of any star's:
    - extremely fast mag changes, and
    - lightcurve fine structure.
    
 All the above is intended to help observers better define the light curve over the next year.
 
The lightcurve model adopted here is second-order Fourier with first-order linear drift, which results 
in 6 fitted parameters. More recent observations are weighted more heavily as they are certainly more 
relevant (per obs) to next-year predictions.

Because the LPV period is a highly non-linear parameter, we fit it indirectly and approximately.
Even so, initial experiments indicate that when a LPV's period deviates from the historical
by more than perhaps 5%, this can be detected. If this result holds, each such star should be
marked as a special target for its next maximum and minimum.

---

 **An advantage of online tools: frequent model updates:** 

A severe requirement imposed on this demonstration is: for predictions over the year 2019,
only data through November 30, 2018 can be used in fitting. This is similar to the real-world
requirement when previous Bulletins were produced once annually (for the following calendar year).

However, this restriction is far more severe than would be encountered by an online tool,
because any online tool would have the advantage of using more recent data.

For example, previous Bulletins were produced annually, and thus needed to project minima and maxima
up to 13 months in advance. However, an online tool would have direct access to recent
observations, and if supporting least-squares fits were run only monthly, they would then need
to project magnitudes up to 2 months in advance; if fits were run more often,
the projections would be even less distant.

So the present experiment is required to project magnitudes into far more distant dates than
any online tool would need to do. Thus, if the current projection approach can work, a less
restrictive approach should certainly work well to support an online Bulletin tool.

---

 **Current fit-to-predict process, such as it is *in media res*: **
 
* **Collect required data.** At present, that comprises: 
    - the entire 2018 LPV bulletin in table form, and
    - observation data from the past 5 LPV cycles (periods).

  I don't think we need more than that for now.
  
* **Perform the fit to past observation data.** The fit model is described above, having 6 linear 
parameters; they are: 
    - V mag at JD0 (reference date),
    - linear V mag drift per year (e.g., negative for slowly brighening star),
    - First Fourier term sin(phase),
    - First Fourier term cos(phase),
    - Second Fourier term sin(2*phase),
    - Second Fourier term cos(2*phase),

  where phase is (JD - JD0) / period.
  
* **For each star, change the period if needed.** Unfortunately, the 2018 Bulletin cycle periods 
  sometimes aren't adequate to predict 2015-2018 (training data set) *or* 2019 (evaluation set)
  magnitudes or lightcurve shapes. 
  For some stars the Bulletin periods may be wrong, whereas for others the LPVs' periods may 
  actually be changing over a few years--not unknown.
  For our purposes, it doesn't matter: even a small period error can accumulate to 0.5-1 magnitude
  error when projected forward a year, and that's outside my acceptance.
  
  And also unfortunately, the period is a very non-linear parameter. Rather than take on non-linear
  regression and all its stability and accuracy headaches, it's much better just to re-run the 
  regression with modified periods and see what happens. It's far easier to program, and in fact even
  repeated linear regressions run much faster than one non-linear fit. So modifying the period in 
  search of better fits is what we do. If the pattern of fits indicates that a star's Bulletin
  period is off, that simply renders that star a more valuable target for the next year or two, or at
  least for the next few maxima.



* **For a comprehensive web document (like the 2018 Bulletin):** For each star, use its best model
fit to predict its V magnitude for every day over the next year. Really. It's blindingly fast on any
modern PC or server, even for hundreds of stars. Yay, linear models.

  From the projected light curves, pick out each star's minima and maxima dates over the next year, 
  together with their approximate V magnitudes. Add date ranges of special fine structure in the
  projected light curves, too.
  Brightly mark stars in probable need of observations, whether from lack of previous coverage or
  from need to refine their period durations. 
   
  Report all this in a customized document that observers can download, most
  profitably limiting it to their preferred: latitude and viewing altitude, magnitude range,
  and event type or observation need (minima, maxima). The idea is an on-demand, customized
  Bulletin with 

* **For a real-time web tool:**  Imagine: a talented observer--whether visual, CCD, or other--wants
to know "what's up": what LPVs are available for him to observe profitably. 

  On the AAVSO website, 
the new tool takes his location, telescope, size, sky brightness, star-type preferences
and other data (possibly part of his profile already), and delivers a listing of targets to go 
observe tonight, indeed **right now**, together with days since last observed, times that each
star is available *tonight*, estimated (approximate) magnitude range recommended observing 
frequency, etc etc. The observer can change his mind and resubmit. 

---

I'm continuing work to refine and test this approach and its breadth of application, 
but I see no technical obstacles at all. 
Certainly not to the math--that's already working.

When launched, the system could run its own model and projection updates, generate customized
documents, and give real-time observer guidance with *zero* real-time staff effort. 
Their manual intervention would be required only when a star acts out of model, 
and while the investigation into that star could be time-consuming, 
(1) in the aggregate it's bound to be a lot less burden than previous efforts to generate
annual Bulletins by hand, and 
(2) when something goes wrong, it's probably because a star is not behaving as expected, and... 

**_those_ are exactly the stars most worth investigating anyway.** 

That is, in addition to serving the observers directly, the proposed system would simply
be acting as an automated screen to flag those abberant stars, 
rather than waiting for someone to notice their behavior, perhaps only by accident.
The partially complete system contained in this code repository is an effort to do better.


                      in progress November 2019
                           Eric Dose, Albuquerque, NM USA
                           in support of the American Association of Variable Star Observers
                                          
