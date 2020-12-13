#  bulletin2
Python code for an online resumption 2020- of the venerable LPV (Long-period variable) Bulletins
published through 2018 by the
[American Association of Variable Star Observers (AAVSO)](https://aavso.org).

AAVSO has a strong interest in informing LPV observers, so that they can make the most
scientifically useful observations that satisfy their observing preferences.
There are two obvious approaches: 
1. Directly replace the 
previous annual LPV Bulletins. Observers download this once per 
year and used it for that entire calendar year. 
2. Provide a new online tool for LPV observers to access
through a regular browser, whenever they want. 
Ideally, they would receive information only on target stars suited to their 
geographic location and observing abilities or preferences. 

**_The code in this repo satisfy 1. above, providing CSV files_** to observers of long-period variable
stars. These CSV files include the same 379 LPV stars that were included in the last (2018) Bulletin
produced by AAVSO staff. These new CSV files are semi-colon delimited and may be imported directly
into Excel, OpenOffice, etc, where all observers may easily screen and sort the stars to
emphasize those in line with their own geographic locations, telescope properties, etc.

_The previously proposed online tool for same-night LPV predictions and user-screening and -sorting
has proven too involved for the present. It will probably not be pursued in 2020._

I (Eric Dose, AAVSO member and active LPV observer) plan to generate these CSV-format Bulletins 
twice per year and to post them in the AAVSO forums, LPV section. 
The first Bulletin was posted there March 21 2020; it covers March-December 2020. 
Based on feedback received from observers, I plan to update the software and produce the next
Bulletin in June 2020 to cover June 2020-March 2021. 
If there is sufficient demand, and especially if coverage improves for the historically 
least-observed stars, then I might increase the frequency of publication to quarterly.

###  Motivation
 
 We want to provide to active LPV observers the following, 
 for all the 379 LPV (long-period variable) stars in Bulletin 2018, 
 in a new Bulletin LPV editions:
 * base information about each star, including RA and Declination, best magnitude-range estimates, and
 best Period estimates;
 * minimum and maximum dates during next few months;
 * number of observations of each star during the previous 12 months; and
 * magnitudes in V filters for each minimum and maximum date.
    
All the above is intended to help observers better define the best estimated of each LPV
light curve over the next few months.
 
### The lightcurve model

For these new Bulletins, the fundamental numerical model for the persistent shape of LPV
lightcurves is a **second-order Fourier** with first-order linear drift.
This model therefore comprises 6 fitted parameters:
* constant base (average) V magnitude,
* constant first derivative (dV/dt) of the base magnitude,
* sine and cosine terms at the fundamental frequency (1/Period), and
* sine and cosine terms at twice the fundamental.

We do not fit the LPV fundamental Period as a separate parameter.
Attempts to do so gave unreliable results for two reasons:
1. Period is an extremely non-linear parameter, which meshes very poorly with the other 6 parameters,
all of which are perfectly linear in fitting magnitudes; and
2. LPV periods do change over time scale of years, so that listed
periods need to be adjusted for stable Fourier fits. 
Worse, period is an extremely non-linear parameter, 
so we must fit it indirectly and approximately.

Sometimes even the best documented periods 
aren't adequate to predict magnitudes or lightcurve shapes. 
For some stars the Bulletin periods may be wrong, whereas for others the LPVs' periods may 
actually be changing over a few years.
For our purposes, the cause doesn't matter: 
even a 2-5% period error can accumulate to 0.5-1 magnitude
error when projected forward a year, and that's 
outside usefulness to the observers relying on the predictions.

Instead, we iterate the linear fits above with several potential values of the Period, clustered
around the historical Period as read from the latest AAVSO-produced 2018 Bulletin. 
This approach has proven numerically stable, and it provides very reasonable estimates of
the _current_ Period.

### Data used to determine model parameters

For each LPV star in the Bulletin, 
we first collect from AAVSO's AID all recent data in V-like passbands: V, Visual, and TG, where
recent is defined here as ending one week before the modeling run (to minimize selection bias in
very recently reported data), and beginning 5 nominal (2018 Bulletin) Periods before the end date.
Within the model, we weight more recent observations somewhat more heavily, 
as they are certainly more relevant (per obs) to next-year predictions.

The 379 LPVs vary extremely widely in their observation coverage: some have a few dozen reported
magnitudes, others have thousands. 
This means that the certainty in modeling (and thus in Bulletin mag predictions) varies widely
as well. 
The uncertain prediction of low-coverage LPV lightcurves results from that low coverage and is
not a feature of the current numerical model--indeed, projecting magnitudes from very limited
historical data would be even more uncertain by any method, and certainly would be worse than
trying to project via visual inspection, e.g.

### Workflow for production of new Bulletins

This "bulletin2" repo contains all the code needed to produce new CSV-format LPV bulletins.

The workflow comprises running several functions in a specific order. I run them from within
the Python Console facility in the IDE PyCharm, where they are imported via the statement

### `import bulletin2.util as u`

so that each function call begins with "u." as seen below.

---

**SHORTCUT**: the function call:

### `u.do_it_all(new_bulletin_start, n_bulletin_months)`

for example,

### `u.do_it_all(20201201, 13)`

will do everything. (You'll still need to sanity-check the results.)

---

Otherwise, you'll proceed (as do_it_all() does) with these sequential steps.

### `dates = u.calc_dates('20200301', 9)`

where the first parameter is the desired beginning date (yyyymmdd) of the new Bulletin, 
and the second
parameter is the desired length of the new Bulletin in whole calendar months.

Computes the Bulletin months, the end date for historical data used to
fit the numerical model (the start data depends on each LPV's nominal period), and the
months used to compute N(obs), the number of historical observations within the latest year.

These controlling data are written to a dictionary (compound variable) "dates" 
used by the succeeding steps.

### `u.capture_vsx_data(dates)`

Collects and organizes all the data needed for modeling and prediction.
Specifically this function:
* calculates the begin and end dates for each LPV star
* downloads each star's historical data from AAVSO's AID 
(AAVSO International Database) via its VSX (Variable Star Index) API,
* counts the previous year's V-like observations and writes all stars' N(obs) data
to a CSV file "df_nobs.csv", and
* stores all the downloaded observations (about 300,000) for all the stars into 
one enormous CSV file "df_obs_csv".

### `df_fit_results = u.process_all_stars(dates)`

This applies the numerical model and makes future magnitude predictions for all the 379
LPV stars. The output is written to a large pandas dataframe "df_fit_results" for use by the next step.

### `u.make_new_bulletin(dates, df_fit_results)`
### `u.make_new_bulletin(dates, df_fit_results, include_magnitudes=False)`


This function actually makes the TWO new Bulletin CSV files: the first CSV file
includes estimated min & max magnitudes, and the second CSV file does not include
estimated min & max magnitudes but only their dates (for observers who prefer that).
The CSV files are named, for example, LPVbulletin_2020-03.csv and 
LPVbulletin_2020-03_nomags.csv.

This function uses the fit results from u.process_all_stars(), 
finds the minima and maxima (and optionally the estimated magnitudes at those dates, 
organizes the results into a large table, and 
writes the table into the final CSV file. 

Again, the CSV file will include a magnitude estimate for each min and max 
only if include_magnitudes is
set to True. If set to False, the min and max entries will have no magnitude estimates,
as in the 2018 and earlier Bulletins.

You'll want to do a sanity check on a few star entries before posting. 
Especially, you'll want to manually remove from the CSV files any stars that have fewer than
maybe 10-12 observations over the past 5 years. 

The versions: 
* CCD observers will need the predicted magnitudes to optimal decide exposure times. 
* Some visual observers may prefer to plan observations without advance knowledge of the
star's expected brightness.

## Project status:

I (Eric Dose) am done with this project.
 
I have posted my last bulletins (Dec 2020-Dec 2021).
Anyone at AAVSO is welcome to clone this repo, edit it as they see fit,
and start generating Bulletins.
If you do so, all future credit as well as responsibility for future Bulletin postings is 
yours alone.

Best of luck in all your observing--may your skies be clear.


                      updated December 13, 2020
                           Eric Dose, Albuquerque, NM USA
                           in support of the American Association of Variable Star Observers
                                          
