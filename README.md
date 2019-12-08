#  bulletin2
Demonstration code for an online resumption and extension of the venerable LPV Bulletins
previously published by the [American Association of Variable Star Observers (AAVSO)](https://aavso.org).

AAVSO has a strong interest in informing LPV observers, so that they can make the most
scientifically useful observations in line with their observing preferences.
There are two obvious approaches: 
1. Directly replace the 
previous annual LPV Bulletins. Observers download this once per 
year and used it for that entire calendar year. 
2. Provide a new online tool for LPV observers to access
through a regular browser, whenever they want. 
Ideally, they would receive information only on target stars suited to their 
geographic location and observing abilities or preferences. 

Both approaches require that the master list of LPVs have each target star's 
V / Visible magnitude predicted in advance, even if only approximately. 
But the two approaches' mag-prediction demands differ 
enough that each needed its own demonstration. 
Each "static" Bulletin document needs long-term (~15 months') predictions, whereas
an online tool could have its predictions refreshed monthly, weekly, or even conceivably
in real time, using all the data received to that moment for the freshest magnitude 
predictions. 

The two automated prediction strategies are quite different.
Through simulation with existing database data, I have demonstrated that both approaches 
are feasible. 

**Both demonstrations are included in this repo.** 

###Motivation
 
 For all the 379 LPV (long-period variable) stars in Bulletin 2018, we want to predict the following:
 * minimum and maximum dates and their approximate V-mags during next calendar year,
 * approximate V-mags ON DEMAND, for any star and date, presumably from a future web tool, and
 * ideally, bonus predictions (dates, V-mags) for the entire next year of any star's:
    - extremely fast mag changes, and
    - lightcurve fine structure.
    
All the above is intended to help observers better define the light curve over the next year.
 
###The lightcurve model

Both demonstrations use exactly the same foundation numerical model for
the persistent shape of LPV lightcurves, 
namely a **second-order Fourier** with first-order linear drift,
thus comprising 6 fitted parameters. 
We weight more recent observations more heavily 
as they are certainly more relevant (per obs) to next-year predictions.

A star's **LPV period** (on which that star's Fourier series is based) 
amounts to an additional parameter--unless one fixes it to constant value,
which proved impossible for two reasons: 
1. even the best sources of 
LPV period information are unreliable for perhaps 20% of the listed stars,
and 
2. LPV periods do change over time scale of years, so that listed
periods need to be adjusted for stable Fourier fits. 
Worse, period is an extremely non-linear parameter, 
so we must fit it indirectly and approximately.

Experiments indicate that when a LPV's period 
deviates from the historical by more than perhaps 5%, 
this can be detected.
As a side-benefit :each such deviant LPV star can and should be
marked as a special target for observation at (at least) its next maximum and minimum.
---
---

####Demonstration 1 of 2: Make an ANNUAL, STATIC LPV Bulletin table:

_Code file util.py. Primary function: make_new_bulletin()._ 

This successfully demonstrated the automated production of 
direct replacements of previous LPV Bulletins.

Previously done manually, the supporting predictions were too laborious
to continue.
We've demonstrated here that they can be automated with satisfactory
precision and accuracy for the vast majority of LPV target stars
listed in the latest 2018 LPV Bulletin. 

Requirements imposed on the predictions supporting this demo were
extremely severe:
For the demo predictions over the year 2019, 
only data through November 30, 2018 was used in fitting.
This is similar to the previous real-world requirement when 
previous Bulletins were produced once annually during December 
for the following calendar year, actually through February of the
next year still.

The prediction process for annual, static Bulletins:

* **Collect required data.** At present, that comprises: 
    - the entire 2018 LPV bulletin in table form, and
    - observation data from the past 5 LPV cycles (periods).
 
* **Model (fit to) the past observation data.** The 
  fit model is described above, having 6 linear 
   parameters; they are: 
    - V mag at JD0 (reference date),
    - linear V mag drift per year (e.g., negative for slowly brighening star),
    - First Fourier term sin(phase),
    - First Fourier term cos(phase),
    - Second Fourier term sin(2*phase),
    - Second Fourier term cos(2*phase),

  where phase = (JD - JD0) / period.
  
* **Adjust each star's model period if needed.** 
  Sometimes even the best documented periods 
  aren't adequate to predict 2015-2018 (training data set)
  *or* 2019 (evaluation set) magnitudes or lightcurve shapes. 
  For some stars the Bulletin periods may be wrong, whereas for others the LPVs' periods may 
  actually be changing over a few years.
  For our purposes, the cause doesn't matter: 
  even a 2-5% period error can accumulate to 0.5-1 magnitude
  error when projected forward a year, and that's 
  outside usefulness to the observers relying on the predictions.
  
  And also unfortunately, the period is a very non-linear parameter. Rather than take on non-linear
  regression and all its stability and accuracy headaches, it's much better just to re-run the 
  regression with modified periods and see what happens. It's far easier to program, and in fact even
  repeated linear regressions run much faster than one non-linear fit. So modifying the period in 
  search of better fits is what we do. If the pattern of fits indicates that a star's Bulletin
  period is off, that simply renders that star a more valuable target for the next year or two, or at
  least for the next few maxima.

* **Make the predictions:** For each star, use its best model
  fit to predict its V magnitude for every day over the next year. 
  Really. It's blindingly fast on any modern PC or server, 
  even for hundreds of stars. Yay, linear models.

  From the projected light curves, pick out each star's minima and maxima dates over the next year, 
  together with their approximate V magnitudes. Add date ranges of special fine structure in the
  projected light curves, too.
  Brightly mark stars in probable need of observations, whether from lack of previous coverage or
  from need to refine their period durations. 
   
  Report all this in a customized document that observers can download, most
  profitably limiting it to their preferred: latitude and viewing altitude, magnitude range,
  and event type or observation need (minima, maxima). The idea is an on-demand, customized
  Bulletin with 

* **Construct the static Bulletin document:**
    For the CSV (comma-separated variables) file format, this 
    reduces to one carefully constructed line of python code
    calling a pandas DataFrame method. The CSV file is for all practical
    purposes identical in function and format to the
    2018 LPV Bulletin CSV still available on the AAVSO website.

    Extending this to a more human-readable HTML format would be
    a routine exercise in reformatting the same data, and so 
    was not necessary to this demo. 

_This first demonstration was posted here and signaled to 
AAVSO mid-November 2019._

---

####Demonstration 2 of 2: Make SHORT-TERM predictions of LPV
 V/Vis magnitudes ON DEMAND:

_Code file util.py. Primary code: make_short_pred_demo()._ 

This successfully demonstrated that an automated short-term
prediction engine can deliver satisfactory V/Vis magnitude
guidance on demand.

The next step would be to adapt this modeling and prediction code
to drive the design and launch of a new web page 
to guide AAVSO LPV observers' planning of their 
next month's, next week's, or even same night's LPV observations,
as suitable to their own geographic location, 
sky altitude restrictions, preferred magnitude range, etc. 
 
The prediction process for on-demand predictions (for web tool): 

Requirements for these shorter-term predictions (perhaps from 2 weeks
to 2 months) are much easier statistically--but to access these
statistical gains requires additional computation and code.

* **Model by adapting the long-term model (above):** 

   The model for these short-term predictions starts with the same model 
   as for the long-term Bulletin demo (see above).
   But the short-term model required two modifications:
 
   - We averaged the best-fit period and the previously published
   period to yield the period used in the fit and prediction. 
   This is a form of statistical regularization, and
   was found to help stabilize the fits, especially where few
   recent observations existed or where the period was actually
   changing.
    
   - We adjusted the prediction for account for 
   recent mag shifts too recent for the long-term model to include.
   Specifically, we add to the long-term prediction an offset 
   calculated in this way: for each of the 20 most recent observations,
   compute the difference between the observed magnitude and the 
   (long-term model) predicted magnitude, then use the median
   (not mean, in case of outliers)
  
   The basis data set for each short-term model could not extend
   more recently than the date of the prediction--one cannot use 
   future data. This can never happen when real predictions are used
   in an online tool, but a demo that models past predictions
   must be careful to exclude any data past the targeted prediction 
   date.

* **Make the predictions:** 

   Thus, for a given short-term prediction date, we restricted
   the basis data set to observations at least 30 days old. 
   This is much more restrictive than would be required in a real
   online tool, which would have access to within the past week.
   But in case the predictions can only be run monthly (by a 
   server cron job, e.g.), we enforced this basis-data restriction
   without exception. 

   For each of the 379 stars in the 2018 LPV Bulletin, a V/Vis
   magnitude prediction was made for 15 dates in 2019 and early 2020, 
   spaced by 30 days. This comprised about 5,700 predicted magnitudes. 
   The basis-data restrictions describe above were applied without
   exception to ensure that this demo reflects actual performance
   of any future online tools based on this engine.

* **Make plots from the predictions:** 

   Then for each star, one plot was made which overlay the 
   predicted magnitudes (red dots) over all observation data
   in V, Vis, and TG. With these 379 plots, one can compare at
   a glance the realistically predicted magnitudes with the since-
   observed magnitudes through most of 2019.
   
   All plots were grouped 12 to a page and included in a PDF
   document delivered to AAVSO December 7, 2019.    

   These short-term predictions tracked actual 2019 observations very
   well for the vast majority of stars. There were perhaps a dozen
   stars for which the predictions were unsatisfactory; this seems
   to happen when stars change behavior (which is unavoidable, but
   represents valuable opportunities for future observations),
   or when there have been very few recent observations (which is
   unfortunate, but then an online tool is under consideration
   to remedy exactly that).  

---
---

##In conclusion:

Imagine: a talented observer--whether visual, CCD, or other--wants
to know "what's up" tonight: what LPVs are available for him to
observe profitably, with the most scientific value.

On the AAVSO website, 
the new tool takes his location, telescope, size, sky brightness, star-type preferences
and other data (possibly stored in his profile already), 
and delivers a listing of targets to go observe tonight, 
indeed **right now**, together with days since last observed, times that each
star is available *tonight*, estimated (approximate) magnitude 
range recommended observing frequency, etc etc. 
The observer can also change his mind, resubmit, get a refreshed observation roster. 

I recommend this be done. The fit and prediction engines have been demonstrated
to work at least as well as needed, and the process is entirely automated--with
the exception of stars that have very few previous observations or which have
changed their physical behavior. Those deviant situations will always require
human intervention, **_but also represent the greatest opportunities_**.

These demos are complete, their results are submitted, and 
I have no more work planned on them.

Should it be decided to restart the LPV Bulletins either as static, periodic
documents or as an online tool, I stand ready to adapt the computational
prediction engines described in this repo, 
working with others who would handle server aspects and web-page front end design.



                      work completed December 7, 2019
                           Eric Dose, Albuquerque, NM USA
                           in support of the American Association of Variable Star Observers
                                          
