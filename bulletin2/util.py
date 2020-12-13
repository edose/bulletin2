# Python core packages (need no separate installation):
import functools
from datetime import datetime, timezone, timedelta
import os
from math import sin, cos, pi, sqrt, floor, isnan
from collections import OrderedDict
from statistics import median

# External packages requiring separate installation:
import pandas as pd  # as standard as it gets.
import statsmodels.api as sm  # should be v.0.10 or higher.

# External packages probably not needed for online tools:
import urllib.request  # only for downloading of AID data only via web; prob. not needed by server tools.
import matplotlib.pyplot as plt  # only needed for my demo plots, v. prob. not for server tools.

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

DATA_DIRECTORY = 'C:/Dev/bulletin2/data'
BULLETIN2018_FILENAME = 'Bulletin2018.csv'
DF_OBS_FILENAME = 'df_obs.csv'
DF_NOBS_FILENAME = 'df_nobs.csv'
DATES_FILENAME_TEMPLATE = 'Dates_[date].txt'               # [date] to be replaced by '2020-03', e.g.
NEW_BULLETIN_FILENAME_TEMPLATE = 'LPVbulletin_[date].csv'  # "

NEW_BULLETIN_HEADER_STUB = ['Name', 'RA.hour', 'RA.min', 'RA.sec', 'DECL.deg', 'DECL.min', 'DECL.sec',
                            'Period', 'Range', 'N(obs)']
DF_BULLETIN_COLUMNS_STUB = ['star_id', 'ra_hour', 'ra_min', 'ra_sec', 'decl_deg', 'decl_min', 'decl_sec',
                            'period', 'range', 'n_obs']
MONTH_ABBRVS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
MIN_MAX_WINDOW_DAYS = 11  # window width to detect min or max mag, in days.
MIN_MAX_MARGIN_DAYS = 10  # must be the only min or max within this no. of days to be recognized.
MAX_BULLETIN_MONTHS = 14

DAYS_PER_YEAR = 365.25
HTTP_OK_CODE = 200  # "OK. The request has succeeded."
FIT_PERIODS_TO_CAPTURE = 7  # LPV periods to capture from VSX, counting back from FIT_END_DATE.
FIT_PERIODS = 5  # LPV periods to actually use during fit.
VSX_OBSERVATIONS_HEADER = 'https://www.aavso.org/vsx/index.php?view=api.delim'
VSX_DELIMITER = '@@@'  # NB: ',' fails as obsName values already have a comma.
MIN_MAG = 0.1  # magnitudes must be greater than this numeric value (mostly to remove zeroes).
MAX_MAG = 20.0  # maximum reasonable magnitude numeric value (subject to later revision for meter scopes).
RISING_ENTRY = 'rising'  # for bulletin cells with rising brightness (decreasing mag).
FADING_ENTRY = 'fading'  # for bulletin cells with fading brightness (increasing mag).
CONSTANT_ENTRY = 'const.' # for bulletin cells with constant brightness (probably never happens).



PERIOD_BEST_SHIFT_TO_FIT_ONLY_SHIFT = 0.6  # fraction of fit period shift to actually use (regularization).

# =====================================================================================================
#
#     bulletin2
#
#     Demonstration code for an online replacement and extension of AAVSO's venerable LPV Bulletins.
#
#     Motivation: For all the 379 LPV (long-period variable) stars in Bulletin 2018, we want
#     to predict the following:
#         * minimum and maximum dates and their approximate V-mags during next calendar year,
#         * approximate V-mags ON DEMAND, for any star and date, presumably from a future web tool, and
#         * ideally, bonus predictions (dates, V-mags) for the entire next year of any star's:
#                   * extremely fast mag changes, and
#                   * lightcurve fine structure
#                in order to help observers better define the light curve over the next year.
#
#    The lightcurve model adopted here is second-order Fourier with first-order linear drift,
#        which results in 6 fitted parameters. More recent observations are weighted more heavily as
#        they are certainly more relevant (per obs) to next-year predictions.

#    Because the LPV period is a highly non-linear parameter, we fit it indirectly and approximately.
#        Even so, initial experiments indicate that when a LPV's period deviates from the historical
#        by more than perhaps 5%, this can be detected. If this result holds, each such star should be
#        marked as a special target for its next maximum and minimum.
#
#                                         ---------------
#
#    Note: A severe requirement imposed on this demonstration is: for predictions over the year 2019,
#        only data through November 30, 2018 can be used in fitting. This is similar to the real-world
#        requirement when previous Bulletins were produced once annually (for the following calendar year).
#
#        However, this restriction is far more severe than would be encountered by an online tool,
#        because any online tool would have the advantage of using more recent data.
#
#    For example, previous Bulletins were produced annually, and thus needed to project minima and maxima
#        up to 13 months in advance. However, an online tool would have direct access to recent
#        observations, and if supporting least-squares fits were run only monthly, they would then need
#        to project magnitudes up to 2 months in advance; if fits were run more often,
#        the projections would be even less distant.
#
#    So the present experiment is required to project magnitudes into far more distant dates than
#        any online tool would need to do. Thus, if the current projection approach can work, a less
#        restrictive approach should certainly work well to support an online Bulletin tool.
#
#
#                                           offered by
#                                           Eric Dose, Albuquerque, NM USA
#                                           in support of the American Assoc. of Variable Star Observers
#
# =====================================================================================================
#
#   Workflow [version 1.01, 2020-03-21]:
#
#   * First, ensure old Bulletin2018.csv is present (mostly to have list of target stars).
#   1. dates = u.calc_dates(new_bulletin_start='20200301', n_bulletin_months=MAX_BULLETIN_MONTHS)  # dict
#   2. u.capture_vsx_data(dates)  # writes df_obs.csv and df_nobs.csv.
#   3. df_fit_results = u.process_all_stars(dates)
#   4a. u.make_new_bulletin(dates, df_fit_results)  # to include magnitudes
#   4b. u.make_new_bulletin(dates, df_fit_results, include_magnitudes=False)  # to omit magnitudes
#
#   or just run: u.do_it_all('202003', 10)
#
# =====================================================================================================
# =====================================================================================================


def do_it_all(new_bulletin_start, n_bulletin_months):
    """
    :param new_bulletin_start: date to start next bulletin, e.g. '20200301' or '202003'. [string or int]
    :param n_bulletin_months: number of months to include in bulletin, <= MAX_BULLETIN_MONTHS. [int]
    """
    dates = calc_dates(new_bulletin_start, n_bulletin_months)
    capture_vsx_data(dates)
    df_fit_results = process_all_stars(dates)
    make_new_bulletin(dates, df_fit_results, include_magnitudes=True)
    make_new_bulletin(dates, df_fit_results, include_magnitudes=False)


def calc_dates(new_bulletin_start, n_bulletin_months=MAX_BULLETIN_MONTHS):
    """ From 2 key values, compute all needed date ranges for next bulletin,
        write the results to dates.csv.
    :param new_bulletin_start: date to start next bulletin, e.g. '20200301' or '202003'. [string or int]
    :param n_bulletin_months: number of months to include in bulletin, <= MAX_BULLETIN_MONTHS. [int]
    :return: dictionary of dates etc needed for downstream functions. [py dict]
    """
    # Parse new_bulletin_start (always 1st of month):
    start = str(new_bulletin_start)  # ensure str.
    if len(start) not in [6, 8]:
        print(' >>>>> ERROR: new_bulletin_start of', new_bulletin_start,
              'is', start, ', s/b like "20200301".')
        exit(0)
    if not(1 <= n_bulletin_months <= MAX_BULLETIN_MONTHS):
        print(' >>>>> ERROR: n_bulletin_months of', n_bulletin_months,
              'is', str(n_bulletin_months), ', but must be within 1 to',
              MAX_BULLETIN_MONTHS, 'inclusive.')
        exit(0)
    year_start = int(start[0:4])
    month_start = int(start[4:6])
    date_string = '-'.join([start[0:4], start[4:6]])  # '2020-03', e.g. for use in filenames.

    # Make utc_bulletin_start (always 1st of month) & validate it:
    utc_bulletin_start = datetime(year_start, month_start, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
    utc_now = datetime.now(timezone.utc)
    if abs(utc_bulletin_start - utc_now) > timedelta(days=(2 * DAYS_PER_YEAR) + 1):
        print(' >>>>> ERROR: new_bulletin_start of', start,
              'is not within 2 years of the present but must be.')
        exit(0)

    # Compute utc_bulletin_end (always 1st of following month):
    month_end_raw = month_start + n_bulletin_months
    years_increment, month_end = divmod(month_end_raw, 12)
    year_end = year_start + years_increment
    utc_bulletin_end = datetime(year_end, month_end, 1, 0, 0, 0).replace(tzinfo=timezone.utc)

    # Compute all other dates:
    fit_end_date = utc_now - timedelta(days=7)  # Last date of historical date used in lightcurve fits.
    utc0 = fit_end_date  # Phase-zero during fit (~ arbitrary).
    nobs_end = datetime(fit_end_date.year, fit_end_date.month, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
    nobs_start = datetime(nobs_end.year - 1, nobs_end.month, 1, 0, 0, 0).replace(tzinfo=timezone.utc)

    # Make and return dictionary of all dates (datetimes in UTC) and related data:
    bulletin_filename = NEW_BULLETIN_FILENAME_TEMPLATE.replace('[date]', date_string)
    dates = dict([('bulletin_start', utc_bulletin_start),
                  ('bulletin_end', utc_bulletin_end),
                  ('bulletin_filename', bulletin_filename),  # e.g., 'bulletin_202003.csv'.
                  ('n_bulletin_months', n_bulletin_months),
                  ('fit_end', fit_end_date),
                  ('utc0', utc0),
                  ('jd0', jd_from_datetime_utc(utc0)),
                  ('nobs_start', nobs_start),
                  ('nobs_end', nobs_end)])
    lines = ['NEW BULLETIN will cover ' + str(n_bulletin_months) + ' months from ' +
             str(utc_bulletin_start.date()) + ' through ' +
             str((utc_bulletin_end - timedelta(days=1)).date()),
             '   Fit on data through ' + str(fit_end_date.date()),
             '   Obs count for period ' + str(nobs_start.date()) + ' through ' +
             str((nobs_end - timedelta(days=1)).date()),
             '   Output to file:' + os.path.join(DATA_DIRECTORY, bulletin_filename).replace('\\', '/')]
    lines = '\n'.join(lines)

    dates_filename = DATES_FILENAME_TEMPLATE.replace('[date]', date_string)
    fullpath = os.path.join(DATA_DIRECTORY, dates_filename)
    with open(fullpath, 'w') as f:
        f.writelines(lines)
    print('\n', lines, '\nNext: u.capture_vsx_data(dates)')
    return dates


@functools.lru_cache(maxsize=128, typed=False)
def get_vsx_obs(star_id, jd_start, jd_end=None):
    """
    Downloads observations from AAVSO's webobs/VSX for ONE star, returns dataframe of results.
       If star not in AAVSO's AID, return a dataframe with no rows.
       Return columns: target_name, date_string, filter, observer, jd, mag, error.
    :param star_id: the STAR id (not the fov's name) [string].
    :param jd_start: optional Julian date [float].
    :param jd_end: optional JD [float].
    :return: DataFrame containing data for 1 star; 1 row per observation downloaded,
        (or None if there was some problem).
    """
    # Construct required URL:
    parm_ident = '&ident=' + star_id.replace("+", "%2B").replace(" ", "+")  # make syntax safe within URL.
    if jd_end is None:
        jd_end = jd_now()
    parm_tojd = '&tojd=' + '{:20.5f}'.format(jd_end).strip()
    parm_fromjd = '&fromjd=' + '{:20.5f}'.format(jd_start).strip()
    parm_delimiter = '&delimiter=' + VSX_DELIMITER
    url = VSX_OBSERVATIONS_HEADER + parm_ident + parm_tojd + parm_fromjd + parm_delimiter

    # Get & parse data from URL:
    byte_text = urllib.request.urlopen(url)
    text = [line.decode('utf-8') for line in byte_text]
    data = [line.split(VSX_DELIMITER) for line in text]
    if len(data) == 0:
        return None
    this_dict = dict()
    column_names = []
    for column_name in data[0]:
        column_names.append(column_name.strip())
        this_dict[column_name.strip()] = []
    for row in data[1:]:
        for i_col, column_name in enumerate(column_names):
            this_dict[column_name].append(row[i_col].strip())
    df_star_obs = pd.DataFrame(this_dict)
    df_star_obs = df_star_obs.set_index('obsID', drop=False)
    return df_star_obs


def capture_vsx_data(dates):
    """  From VSX API, capture test data from VSX API, make a dataframe and write it to .csv file.
            Dataframe is df_obs; with one row per (screened & weighted) observation.
        Write these dataframes locally as .csv files, for use in testing magnitude-projection model(s).
        The idea is to get this data just one time and then store it locally for repeated use,
        saving us testing time, and saving AAVSO server time.
        Further, the number of all observations on each star over the period NOBS_START_DATE to
            NOBS_END_DATE is written to one line of a new dataframe, which is then written to .csv file.
    :param dates: dict of dates etc
    :return: [None] ... rather, df_obs dataframe is written to [DF_OBS_FILENAME] in [DATA_DIRECTORY].
                       Also, a df_nobs dataframe is written to [DF_NOBS_FILENAME] in [DATA_DIRECTORY].
    """
    df_obs = pd.DataFrame()  # master dataframe of all obs and metadata, to be written to local .csv file.
    nobs_dict_list = []

    star_ids = get_bulletin_star_ids()
    df_bulletin = get_df_bulletin()
    jd_df_nobs_start = jd_from_datetime_utc(dates['nobs_start'])  # these are the same for all stars.
    jd_df_nobs_end = jd_from_datetime_utc(dates['nobs_end'])      # "

    for star_id in star_ids:
        # Ensure enough dates to cover both df_obs and df_nobs:
        # ... where df_obs has all the observations, and df_nobs exists only to count nobs in past year.
        period = float(df_bulletin.loc[star_id, 'PERIOD'])  # in days.
        jd_df_obs_start = jd_from_datetime_utc(dates['fit_end']) - FIT_PERIODS_TO_CAPTURE * period
        jd_vsx_start = min(jd_df_obs_start, jd_df_nobs_start)  # jd_df_obs will almost always be the min.
        jd_vsx_end = jd_from_datetime_utc(dates['fit_end'])
        df_vsx = get_vsx_obs(star_id, jd_start=jd_vsx_start, jd_end=jd_vsx_end)

        # Calculate nobs for this star_id, save it as a dict, store dict in list.
        jd_all_obs = pd.Series([float(jd) for jd in df_vsx['JD']])
        keep_for_df_nobs = (jd_df_nobs_start <= jd_all_obs) & (jd_all_obs <= jd_df_nobs_end)
        nobs = sum(keep_for_df_nobs)
        nobs_dict = {'star_id': star_id, 'nobs': nobs}
        nobs_dict_list.append(nobs_dict)

        # Screen received observation dataframe, add weights column, add df to end of df_obs.
        keep_for_df_obs = list(jd_all_obs >= jd_df_obs_start)
        df_vsx_obs = df_vsx[keep_for_df_obs]  # this will usually keep all obs
        df_screened = screen_star_obs(df_vsx_obs)
        df_this_star = add_obs_weights(df_screened)
        df_obs = pd.concat([df_obs, df_this_star])
        print(star_id.ljust(10),
              str(len(df_this_star)), 'added of',
              str(len(df_vsx)), 'downloaded in JD range',
              str(int(jd_vsx_start)), 'to', str(int(jd_vsx_end)),
              'for running df_obs running total of ', str(len(df_obs)) + '.')

    df_obs_csv_fullpath = os.path.join(DATA_DIRECTORY, DF_OBS_FILENAME).replace('\\', '/')
    df_obs.to_csv(df_obs_csv_fullpath, sep=';', quotechar='"', encoding='UTF-8',
                  quoting=2, index=False)  # quoting=2-->quotes around non-numerics.
    print('   df_obs dataframe written to', df_obs_csv_fullpath)

    df_nobs = pd.DataFrame(nobs_dict_list)  # number of observations in NOBS date range, one row per star.
    df_nobs = reorder_df_columns(df_nobs, left_column_list=['star_id'])
    df_nobs = df_nobs.set_index('star_id', drop=False)
    df_nobs_csv_fullpath = os.path.join(DATA_DIRECTORY, DF_NOBS_FILENAME).replace('\\', '/')
    df_nobs.to_csv(df_nobs_csv_fullpath, sep=';', quotechar='"', encoding='UTF-8',
                   quoting=2, index=False)  # quoting=2-->quotes around non-numerics.
    print('   df_nobs dataframe written to', df_nobs_csv_fullpath)
    print('\nNext: df_fit_results = u.process_all_stars(dates).')


def screen_star_obs(df):
    """  Screen observations (in df_obs) for quality and relevance to V-band magnitude projection.
         Remove observations (df_obs rows) that do not qualify for use in projecting future magnitudes.
         Requirements to keep an observation (row) include:
            Band (filter) being one of: V (CCD V filter), 'Vis.' (visual obs), or 'TG' (V-like);
            Having a reasonable magnitude that is neither null nor a less-than report,
            Validation being either full or provisional.
    :param df: this star's unscreened observations (probably df_obs) [pandas DataFrame].
    :return: this star's screened observations [smaller pandas DataFrame].
    """
    mag_not_null = ~ df['mag'].isnull()
    df = df[mag_not_null]
    mag_value_ok = pd.Series([MIN_MAG <= float(mag) <= MAX_MAG for mag in df['mag']], index=df.index)
    band_is_v_like = df['band'].isin(['V', 'Vis.', 'TG'])
    not_fainter_than = df['fainterThan'] == '0'
    obstype_ok = df['obsType'].isin(['Visual', 'CCD', 'DSLR', 'VISDIG'])
    validation_ok = df['val'].isin(['V', 'Z'])

    obs_to_keep = mag_value_ok & band_is_v_like & not_fainter_than & obstype_ok & validation_ok
    df = df[obs_to_keep]
    return df


def add_obs_weights(df):
    """  Double the fitting weight for each transformed CCD V observation, relative to other observations.
    :param df: this star's screened observations [pandas DataFrame].
    :return: Dataframe with added column 'weight' but otherwise unchanged [pandas DataFrame].
    """
    weights = [2.0 if tr == '1' else 1.0 for tr in df['transformed']]
    df.loc[:, 'weight'] = weights
    return df


def make_df_x(dates, jd_list, fit_period):
    """  Given a list of JDs and other input data, deliver a dataframe ready to use as X-value
         (independent variable) input to a regression function; one column per variable, one row per obs.
    :param dates: dict of dates and related data. [py dict]
    :param jd_list: list of Julian date values to use in fitting observed magnitudes or in predicting
               magnitudes for new JDs [list of floats].
    :param fit_period: the LPV period in days to *assume* in the fit (may differ from Bulletin P) [float].
    :return: small dataframe with required X-value data [pandas Dataframe].
    """
    # jd0 is the Reference Julian Date to use in constructing dataframe [float].
    #    Two effects:
    #       (1) the fitted v_const parameter value will represent the projected mag at this date, and
    #       (2) all phases will use this date as to establish zero phase.
    jd0 = dates['jd0']
    df_x = pd.DataFrame()
    df_x['v_const'] = len(jd_list) * [1.0]  # ca. the mean V mag
    df_x['dv_year'] = [(jd - jd0) / DAYS_PER_YEAR for jd in jd_list]  # linear drift in V mag per year.
    phases = [(jd - jd0) / fit_period for jd in jd_list]  # numerically negative, since jd < JD0, usually.
    df_x['sin1'] = [sin((2 * pi) * 1.0 * ph) for ph in phases]  # first Fourier term
    df_x['cos1'] = [cos((2 * pi) * 1.0 * ph) for ph in phases]  # "
    df_x['sin2'] = [sin((2 * pi) * 2.0 * ph) for ph in phases]  # second Fourier term
    df_x['cos2'] = [cos((2 * pi) * 2.0 * ph) for ph in phases]  # "
    return df_x, jd0


def fit_one_star(dates, star_id, df_obs=None, df_bulletin=None, period_factor=1.0):
    """ The kernel & fitting engine.
        Arranges input data and performs weighted least-squares (via statsmodels WLS function),
        returns the entire WLS results object, as well as the JD0 reference JD used
        (as required for future interpretation of the fit results).
    :param dates: dict of dates and related data. [py dict]
    :param star_id: the STAR id, as NAME in Bulletin 2018 [string].
    :param df_obs:this star's fully screened observations [pandas DataFrame].
    :param df_bulletin: data from LPV Bulletin (probably via get_df_bulletin()) [pandas Dataframe].
    :param period_factor: the factor by which to multiply Bulletin's period when performing fit [float].
    :return: (result, jd0), where result is a RegressionResult object from python's statsmodel package,
             and jd0 is the *actual* jd0 reference JD needed to interpret or use those results [2-tuple].
    """
    # First, ensure that we have the required observation and star data:
    jd0 = dates['jd0']
    if df_obs is None:
        df_obs = get_local_df_obs()
    if df_bulletin is None:
        df_bulletin = get_df_bulletin()

    # Prepare dataframe for this star:
    df = df_obs[df_obs['starName'] == star_id]
    df = df[['JD', 'band', 'mag', 'obsID', 'obsType', 'weight']]
    df['JD'] = [float(jd) for jd in df['JD']]  # because strings are inherited from VSX.
    df['mag'] = [float(mag) for mag in df['mag']]  # because strings are inherited from VSX.
    df = df.set_index('obsID', drop=False)
    df = df.sort_values(by='JD')

    # Get variable-star period from 2018 Bulletin, then set the period actually used in this fit:
    bulletin_period = float(df_bulletin.loc[star_id, 'PERIOD'])
    this_fit_period = period_factor * bulletin_period  # because best period may differ from bulletin's.

    # Prepare dataframes for weighted multivariate regression fit:
    df_jd = pd.DataFrame()
    df_jd['JD'] = df.copy()['JD']
    df_jd.index = df.index
    df_x, jd0 = make_df_x(dates, df['JD'], this_fit_period)
    df_x.index = df.index
    df_y = pd.DataFrame()
    df_y['y'] = df.copy()['mag']
    df_y.index = df.index
    df_wt = pd.DataFrame()
    df_wt['weight'] = df.copy()['weight']
    df_wt.index = df.index

    # Select only obs from fit end date or before :
    jd_end = jd_from_datetime_utc(dates['fit_end'])
    jd_start = jd_end - FIT_PERIODS * this_fit_period
    to_keep = (df['JD'] >= jd_start) & (df['JD'] < jd_end)
    df_jd = df_jd[to_keep]
    df_x = df_x[to_keep]
    df_y = df_y[to_keep]
    df_wt = df_wt[to_keep]
    # In case we later need to fit consistent shifts between observed magnitudes in V, Vis., etc:
    # df_fit['dvis'] = [1.0 if obsType == 'Vis.' else 0.0
    #     for obsType in df_fit['obsType']]  # pseudo-categorical.
    # df_fit['dtg'] = [1.0 if obsType == 'TG' else 0.0 for obsType in df_fit['obsType']]  # "
    # df_fit['dvisdig'] = [1.0 if obsType == 'VISDIG' else 0.0 for obsType in df_fit['obsType']]  # "

    # For each observation falling within most recent LPV period, double that obs' pre-existing weight:
    days_before_jd_end = jd_end - df_jd['JD']
    weights = [2.0 * wt if db <= this_fit_period else 1.0 * wt
               for (wt, db) in zip(df_wt['weight'], days_before_jd_end)]

    result = sm.WLS(df_y, df_x, weights).fit()  # do the fit.
    # print(result.summary())
    return result, jd0


def process_one_star(dates, star_id, df_obs=None, df_bulletin=None, quiet=False):
    """  Fully process (comprehensive fit) one star.
         Calls fit_one_star() multiple times with various trial periods, to:
             (1) better fit most recent obs, if Bulletin period does not apply to them exactly, and
             (2) signal AAVSO and future online Bulletin facility that period appears to be changing.
         Returns a python dictionary with all the data that will be needed to  construct one row
             (representing this star's fitted behaviour) in a summary dataframe.
    :param dates: dict of dates and related data. [py dict]
    :param star_id: the STAR id, as NAME in Bulletin 2018 [string].
    :param df_obs:this star's fully screened observations [pandas DataFrame].
    :param df_bulletin: data from LPV Bulletin (probably via get_df_bulletin()) [pandas Dataframe].
    :param quiet: True=suppress printed updates; False=print them (e.g. when running manually) [boolean].
    :return: a dictionary of all relevant results of a weighted LS fit on one star [python dict object],
             e.g., result_dict['params'] contains fitted parameter values, and results_dict['fit_result']
             contains the entire RegressionResult object needed (with JD0) to perform predictions from
             the fitted WLS model of magnitudes from any JD, or list of JDs.
    """
    jd0 = dates['jd0']
    if star_id not in df_bulletin.index:
        print(' >>>>> Star', star_id, 'is absent from the 2018 Bulletin.')
        return None
    bulletin_period = float(df_bulletin.loc[star_id, 'PERIOD'])
    fit_factors = [0.9, 1.0, 1.1]  # these must be evenly spaced.
    all_results = []

    # Perform fit for one star at a time:
    for fit_factor in fit_factors:
        result, jd0 = fit_one_star(dates, star_id, df_obs, df_bulletin, period_factor=fit_factor)
        all_results.append(result)
    if not quiet:
        for i, result in enumerate(all_results):
            print(i, fit_factors[i], all_results[i].rsquared)

    # Find the best LPV period (with highest R-squared) if possible:
    r2_1, r2_2, r2_3 = tuple([all_results[i].rsquared for i in range(3)])
    if r2_2 > (r2_1 + r2_3) / 2:  # if quadratic through 3 points is concave downward:
        # Lagrange-interpolate for best period:
        x1, x2, x3 = tuple(fit_factors)
        y1, y2, y3 = tuple(all_results[i].rsquared for i in range(3))
        denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
        a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
        b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom
        # c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom
        fit_factor = -b / (2 * a)  # x-value at interpolated maximum.
        if x1 <= fit_factor <= x3:
            # Period seems to be within range, so we will use the interpolated period:
            fit_period = fit_factor * bulletin_period
            period_sign = '~'
            fit_result, jd0 = fit_one_star(dates, star_id, df_obs, df_bulletin, period_factor=fit_factor)
            fit_r_squared = fit_result.rsquared
        elif fit_factor < x1:
            # Period seems to be below range:
            fit_period = fit_factors[0] * bulletin_period
            period_sign = '<'
            fit_r_squared = all_results[0].rsquared
            fit_result = all_results[0]
        else:
            # Period seems to be above range:
            fit_period = fit_factors[2] * bulletin_period
            period_sign = '>'
            fit_r_squared = all_results[2].rsquared
            fit_result = all_results[2]
    else:
        # Here, R^2 vs period appears concave up, so just take the highest R^2, mark with uncertainty sign:
        r2_list = [all_results[i].rsquared for i in range(3)]
        i_max = r2_list.index(max(r2_list))
        fit_period = fit_factors[i_max] * bulletin_period
        period_sign = '?'
        fit_r_squared = all_results[i_max].rsquared
        fit_result = all_results[i_max]

    best_period = bulletin_period + PERIOD_BEST_SHIFT_TO_FIT_ONLY_SHIFT * (fit_period - bulletin_period)

    if not quiet:
        print(star_id.ljust(16), '   Best: ',
              '   factor' + period_sign + '{0:.3f}'.format(best_period / bulletin_period),
              '   P(best)=' + '{0:.1f}'.format(best_period),
              '   R^2=' + '{0:.3f}'.format(fit_r_squared))

    # Make dictionary of results to return to calling function:
    amplitude_1 = 2.0 * sqrt(fit_result.params['sin1']**2 + fit_result.params['cos1']**2)
    result_dict = {'star_id': star_id,
                   'nobs': int(fit_result.nobs),
                   'period_sign': period_sign,
                   'bulletin_period': bulletin_period,
                   'period_factor': best_period / bulletin_period,
                   'best_period': best_period,
                   'jd0': jd0,
                   'r_squared': fit_r_squared,
                   'amplitude_1': amplitude_1,
                   'params': fit_result.params,
                   'params_se': fit_result.bse,
                   'condition_number': fit_result.condition_number,
                   'se': fit_result.mse_resid,
                   'fit_result': fit_result}
    return result_dict


def process_all_stars(dates=None):
    """  Calls process_one_star() for each Bulletin star, then makes a dataframe holding all fit data.
    :param dates: dictionary of dates and related data. [py dict]
    :return: df_fit_results, the dataframe of all fit data; this is enough to perform model predictions
                 for any Bulletin star at any desired Julian Date, future or past [pandas DataFrame].
    """
    df_obs = get_local_df_obs()      # dataframe of all stars' fully screened observations.
    df_bulletin = get_df_bulletin()  # dataframe of previous LPV Bulletin (via get_df_bulletin():
    star_ids = df_bulletin['NAME']
    if dates is None:
        print(' >>>>> ERROR: dates dictionary must be given.')
        exit(0)

    result_dict_list = []
    for i, star_id in enumerate(star_ids):
        result_dict = process_one_star(dates, star_id, df_obs, df_bulletin, quiet=True)
        result_dict_list.append(result_dict)
        print('{0:5d}'.format(i + 1), '  ',
              star_id.ljust(8),
              '{0:4d}'.format(int(result_dict['nobs'])),  # which is the fit nobs, not df_obs nobs.
              '  r2=' + '{0:.3f}'.format(result_dict['r_squared']),
              '  se=' + '{0:.2f}'.format(result_dict['se']), 'mag')
    df_fit_results = pd.DataFrame(result_dict_list).set_index('star_id', drop=False)
    print('\nNext: u.make_new_bulletin(dates, df_fit_results).')
    return df_fit_results


def make_new_bulletin(dates, df_fit_results, include_magnitudes=True):
    """  Project 2019 daily magnitudes for each 2018 Bulletin star,
         construct a 2019 .csv file in same form as 2018 Bulletin .csv.
    :param dates: dict of dates etc
    :param df_fit_results: comprehensive data from best fit for all LPV stars,
               from process_all_stars() [very large pandas DataFrame].
    :param include_magnitudes: True iff predicted magnitude to be written with each min or max. [boolean]
    :return: [None] rather, writes a .csv file containing new Bulletin data
                 in form very like that of 2018 bulletin's .csv file.
    """
    # Define all dates required to make this bulletin:
    bulletin_start_date = dates['bulletin_start']
    bulletin_end_date = dates['bulletin_end']
    predict_start_date = bulletin_start_date - timedelta(days=MIN_MAX_MARGIN_DAYS + 1)
    predict_end_date = bulletin_end_date + timedelta(days=MIN_MAX_MARGIN_DAYS + 1)
    n_dates_to_predict = (predict_end_date - predict_start_date).days + 1
    dates_to_predict = [predict_start_date + timedelta(days=i) for i in range(n_dates_to_predict)]
    jds_to_predict = [jd_from_datetime_utc(date) for date in dates_to_predict]  # don't be clever here.

    # Identify month boundaries, including at end of bulletin timespan:
    month_start_dates = []
    bulletin_start_year = bulletin_start_date.year
    for i_month in range(dates['n_bulletin_months'] + 1):  # + 1 to capture the bulletin end date, too.
        naive_month = dates['bulletin_start'].month + i_month
        month_start_month = ((naive_month - 1) % 12) + 1
        year_increment = (naive_month - 1) // 12
        month_start_date = datetime(bulletin_start_year + year_increment, month_start_month, 1,
                                    0, 0, 0).replace(tzinfo=timezone.utc)
        month_start_dates.append(month_start_date)

    # Make header_dict for header row in new bulletin (inelegant, but sturdy):
    month_columns = [MONTH_ABBRVS[d.month - 1] + '.' + str(d.year) for d in month_start_dates[:-1]]
    df_bulletin_columns = DF_BULLETIN_COLUMNS_STUB + month_columns
    bulletin_header_labels = NEW_BULLETIN_HEADER_STUB + month_columns
    if len(df_bulletin_columns) != len(bulletin_header_labels):
        print(' >>>>> ERROR: df_bulletin has', str(len(df_bulletin_columns)),
              'columns, but bulletin header has', str(len(bulletin_header_labels)),
              'labels (must be equal).')
        exit(0)

    df_bulletin = get_df_bulletin()
    df_nobs = get_local_df_nobs()
    header_dict = dict(zip(df_bulletin_columns, bulletin_header_labels))  # header mimicks a row of data.
    star_dict_list = [header_dict]  # initiate CSV file data with this first (header) row.
    star_ids = df_bulletin['NAME']

    # For each star, predict mag at all days, find each month's status (max, min, etc):
    dates_to_test = [d for d in dates_to_predict if bulletin_start_date <= d < bulletin_end_date]
    print('Testing min/max for dates:', dates_to_test[0].date(), 'through', dates_to_test[-1].date())
    for i_star, star_id in enumerate(star_ids):
        print('{0:5d}'.format(i_star + 1), '  ', star_id)
        period = df_fit_results.loc[star_id, 'best_period']
        # jd0 = df_fit_results.loc[star_id, 'jd0']
        df_x, _ = make_df_x(dates, jds_to_predict, period)
        fit_result = df_fit_results.loc[star_id, 'fit_result']  # a RegressionResults object (large).
        predicted_mags = fit_result.predict(df_x)
        predicted_mag_series = pd.Series(data=predicted_mags.values, index=dates_to_predict)
        rolling_max = predicted_mag_series.rolling(window=MIN_MAX_WINDOW_DAYS, center=True).max()
        rolling_min = predicted_mag_series.rolling(window=MIN_MAX_WINDOW_DAYS, center=True).min()

        # Search all dates for minima and maxima, build list of them:
        min_max_list = []  # will be list of dicts, in date order.
        days_since_last_min_max = None  # None if none yet found, # of days if one has been found.
        for this_date in dates_to_test:
            this_mag = predicted_mag_series[this_date]
            if days_since_last_min_max is None:
                need_to_test_this_date = True
            else:
                need_to_test_this_date = (days_since_last_min_max >= MIN_MAX_MARGIN_DAYS)
            if need_to_test_this_date:
                if isnan(rolling_min[this_date]):
                    print(' >>>>> WARNING: no rolling_min value for date: ', this_date.date())
                else:
                    if this_mag >= rolling_max[this_date]:
                        min_max_dict = {'event': 'min',  # max magnitude is MIN BRIGHTNESS.
                                        'date': this_date,
                                        'mag_pred': this_mag}
                        min_max_list.append(min_max_dict)
                        days_since_last_min_max = 0
                if isnan(rolling_max[this_date]):
                    print(' >>>>> WARNING: no rolling_max value for date: ', this_date.date())
                else:
                    if this_mag <= rolling_min[this_date]:
                        min_max_dict = {'event': 'max',  # min magnitude is MAX BRIGHTNESS.
                                        'date': this_date,
                                        'mag_pred': this_mag}
                        min_max_list.append(min_max_dict)
                        days_since_last_min_max = 0
            if days_since_last_min_max is not None:
                days_since_last_min_max += 1

        # Store this star's minima and maxima behavior in month bins:
        n_months = dates['n_bulletin_months']
        behavior = pd.Series(data=n_months * ['NA'], index=month_start_dates[:n_months])
        for d in min_max_list:
            this_month_start_date = datetime(d['date'].year,
                                             d['date'].month, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
            date_string = '(' + str(d['date'].day) + ')'
            mag_string = ' {0:.1f}'.format(d['mag_pred']) if include_magnitudes is True else ''
            if d['event'] == 'max':
                behavior[this_month_start_date] = 'MAX' + date_string + mag_string
            if d['event'] == 'min':
                behavior[this_month_start_date] = 'MIN' + date_string + mag_string

        # Fill star's behavior for other months' (with rising, fading, or const.):
        for i_month in range(dates['n_bulletin_months']):
            if behavior[i_month] == 'NA':
                month_start_date = month_start_dates[i_month]
                month_end_date = month_start_dates[i_month + 1]
                i_start = (month_start_date - predict_start_date).days
                i_end = (month_end_date - predict_start_date).days
                if predicted_mags[i_end] > predicted_mags[i_start]:
                    behavior[i_month] = FADING_ENTRY  # when brightness is fading, mag is increasing.
                elif predicted_mags[i_end] < predicted_mags[i_start]:
                    behavior[i_month] = RISING_ENTRY  # when brightness is rising, mag is decreasing.
                else:
                    behavior[i_month] = CONSTANT_ENTRY  # (probably never occurs).

        # Now, make star_dict for this star (to become a row in df_new_bulletin), and append it to list:
        star_dict = {'star_id': star_id,
                     'ra_hour': df_bulletin.loc[star_id, 'RA.HOUR'],
                     'ra_min': df_bulletin.loc[star_id, 'RA.MIN'],
                     'ra_sec': df_bulletin.loc[star_id, 'RA.SEC'],
                     'decl_deg': df_bulletin.loc[star_id, 'DECL.DEG'],
                     'decl_min': df_bulletin.loc[star_id, 'DECL.MIN'],
                     'decl_sec': df_bulletin.loc[star_id, 'DECL.SEC'],
                     'period': '{0:.1f}'.format(period),
                     'range': '<' + '{0:.1f}'.format(min(predicted_mags)) + '-' +
                              '{0:.1f}'.format(max(predicted_mags)) + '>',
                     'n_obs': df_nobs.loc[star_id, 'nobs']
                     }
        for i_month in range(dates['n_bulletin_months']):
            # this_key = 'month_' + str(i_month + 1)
            this_key = month_columns[i_month]
            this_behavior = behavior[i_month]
            star_dict[this_key] = this_behavior
        star_dict_list.append(star_dict)

    # Make the dataframe:
    df_new_bulletin = pd.DataFrame(star_dict_list)
    # reorder_df_columns(df_new_bulletin, df_bulletin_columns)
    df_new_bulletin.index = df_new_bulletin['star_id'].values
    # df_new_bulletin = df_new_bulletin.set_index('star_id', drop=False)

    # Write the dataframe to demo .csv file:
    if include_magnitudes is True:
        bulletin_filename = dates['bulletin_filename']  # e.g., 'LPVbulletin_2020-03.csv'
    else:
        parts = dates['bulletin_filename'].rsplit('.', maxsplit=1) + ['']
        bulletin_filename = parts[0] + '_nomags.' + parts[1]  # e.g., 'LPVbulletin_2020-03_nomags.csv'
    fullpath = os.path.join(DATA_DIRECTORY, bulletin_filename)
    df_new_bulletin.to_csv(fullpath, sep=';', quotechar='"', encoding='UTF-8',
                           header=False, quoting=2,
                           index=False)  # quoting=2-->quotes around non-numerics.


def make_short_prediction(dates, star_id, jd_preds=None, df_bulletin=None, df_obs=None):
    """  Make a short prediction for one star at one JD, using both long fit lightcurve and recent
         observations to estimate recent mag shift from that fit lightcurve.
    :param dates: dictionary of dates and related data. [py dict]
    :param star_id: the STAR id, as NAME in Bulletin 2018 [string].
    :param jd_preds: Julian Dates for which short prediction is desired [float, or list of floats].
    :param df_bulletin: data from LPV Bulletin (probably via get_df_bulletin()).
               If absent or None, this function will download data from text file. [pandas Dataframe].
    :param df_obs: dataframes with observations from VSX--must cover all needed ranges.
               If absent or None, this function will download data (slow) [pandas Dataframe].
    :return: best short predictions for star at jd_preds [pandas Series with jd_preds as index].
    """
    # Ensure prediction JD list and bulletin data are ready:
    if jd_preds is None:
        jd_preds = make_jd_preds(datetime(2019, 1, 1).replace(tzinfo=timezone.utc),
                                 datetime(2020, 3, 1).replace(tzinfo=timezone.utc),
                                 30.0)
    elif not isinstance(jd_preds, list):
        jd_preds = [jd_preds]

    star_id = star_id.upper()
    if df_bulletin is None:
        df_bulletin = get_df_bulletin()
    period = float(df_bulletin.loc[star_id, 'PERIOD'])

    # Ensure df_all_obs is ready (fitting shall not use obs more recent than most recent date to predict):
    if df_obs is None:
        earliest_jd_pred, latest_jd_pred = min(jd_preds), max(jd_preds)
        jd_download_start = earliest_jd_pred - FIT_PERIODS_TO_CAPTURE * period
        df_all_obs = get_vsx_obs(star_id, jd_download_start, latest_jd_pred)
    else:
        df_all_obs = df_obs.copy()
    df_all_obs['JD'] = [float(jd) for jd in df_all_obs['JD']]
    df_all_obs = df_all_obs.sort_values(by='JD')
    df_all_obs = add_obs_weights(screen_star_obs(df_all_obs))

    # Make prediction for each date in jd_preds:
    # We must use a loop, because each jd_pred will be fit on a different df_obs (subset of df_all_obs).
    best_mags = []
    jd0 = dates['jd0']
    for jd_pred in jd_preds:
        # Select obs on which to fit:
        latest_jd_fit = jd_pred - 30.0
        earliest_jd_fit = latest_jd_fit - FIT_PERIODS_TO_CAPTURE * period
        to_use_in_fit = [(jd >= earliest_jd_fit) and (jd <= latest_jd_fit) for jd in df_all_obs['JD']]
        df_obs_fit = df_all_obs[to_use_in_fit]

        # Make long prediction for this jd_pred:
        pred_result_dict = process_one_star(star_id, df_obs_fit, df_bulletin, jd0, quiet=True)
        df_x, jd0 = make_df_x([jd_pred], pred_result_dict['best_period'], jd0)
        fit_result = pred_result_dict['fit_result']  # a RegressionResults object (large).
        long_prediction = fit_result.predict(df_x)[0]

        # Calculate mean offset of recent obs from long prediction mags:
        df_most_recent_obs = df_obs_fit[-20:]  # 20 most recent obs used in long prediction fit (above).
        if len(df_most_recent_obs) <= 0:
            # print(' >>>>> make_short_prediction(): No recent obs for', star_id, '... mean offset -> 0.')
            median_offset = 0
        else:
            df_x, jd0 = make_df_x(df_most_recent_obs['JD'], pred_result_dict['best_period'], jd0)
            predicted_mags = fit_result.predict(df_x)  # same fit model as for long prediction above.
            offsets = [(float(o) - p) for (o, p) in zip(df_most_recent_obs['mag'], predicted_mags)]
            median_offset = median(offsets)
            # print('median offset =', median_offset)
        best_short_prediction = long_prediction + median_offset
        best_mags.append(best_short_prediction)

    best_mags = pd.Series(data=best_mags, index=jd_preds)
    return best_mags


def make_short_pred_demo(jd_preds=None):
    """
    :param jd_preds: For all Bulletin 2018 stars, make mag short predictions at series of JDs, render plots.
    """
    # Set up base data:
    df_bulletin = get_df_bulletin()
    star_ids = get_bulletin_star_ids()
    if jd_preds is None:
        jd_preds = make_jd_preds(datetime(2019, 1, 1).replace(tzinfo=timezone.utc),
                                 datetime(2020, 3, 1).replace(tzinfo=timezone.utc),
                                 30.0)
    earliest_jd_pred, latest_jd_pred = min(jd_preds), max(jd_preds)
    jd_floor = 1000.0 * floor(earliest_jd_pred / 1000.0)
    jd_preds_to_plot = [jd - jd_floor for jd in jd_preds]
    label_x = 'JD - ' + str(int(jd_floor))

    # For initial test & debugging:
    df_bulletin = df_bulletin
    star_ids = star_ids

    n_plot_columns, n_plot_rows = 4, 3
    i_plot_row, i_plot_column = 0, 0
    n_plots_completed = 0
    i_first_plot = 1

    for star_id in star_ids:
        # Get df_all_obs here, because we'll need it for plotting below:
        period = float(df_bulletin.loc[star_id, 'PERIOD'])
        jd_download_start = earliest_jd_pred - FIT_PERIODS_TO_CAPTURE * period
        df_vsx_obs = get_vsx_obs(star_id, jd_download_start, latest_jd_pred)
        df_screened = screen_star_obs(df_vsx_obs)
        df_all_obs = add_obs_weights(df_screened)
        df_all_obs['JD'] = [float(jd) for jd in df_all_obs['JD']]
        in_plot_range = (df_all_obs['JD'] >= earliest_jd_pred) & (df_all_obs['JD'] <= latest_jd_pred)
        df_plot_obs = df_all_obs[in_plot_range].copy()
        jd_obs_to_plot = [jd - jd_floor for jd in df_plot_obs['JD']]

        # Compute short-term predictions at each jd_preds value, each with its own base data for fit:
        short_preds = make_short_prediction(star_id.upper(), jd_preds, df_bulletin, df_all_obs)

        # Start new figure page if needed:
        if i_plot_column == 0 and i_plot_row == 0:
            fig, axes = plt.subplots(ncols=n_plot_columns, nrows=n_plot_rows,
                                     figsize=(16, 10))  # (width, height) in "inches"

        # Add plot to overall figure page:
        ax = axes[i_plot_row, i_plot_column]
        ax.set_title(star_id.upper(), y=0.89)
        ax.set_xlabel(label_x, labelpad=0)
        ax.set_ylabel('Mag V/Vis.', labelpad=0)

        # Plot experimental observations and (short) predicted mags:
        obs_mags = [float(mag) for mag in df_plot_obs['mag']]
        ax.scatter(x=jd_obs_to_plot, y=obs_mags, alpha=0.4, color='black', zorder=+500)
        ax.scatter(x=jd_preds_to_plot, y=short_preds, alpha=0.9, color='red', zorder=+1000)

        # Finish plot (y-scale inverted, per custom of plotting mags brighter=upward):
        all_mags = list(obs_mags) + list(short_preds)
        max_mag = max(all_mags)
        min_mag = min(all_mags)
        mag_margin = 0.06 * (max_mag - min_mag)
        ax.set_ylim(max_mag + mag_margin, min_mag - 1.4 * mag_margin)
        ax.grid(True, color='lightgray', zorder=-1000)
        n_plots_completed += 1
        print(star_id.upper(), 'plotted.')

        # Prepare column and row indices for the next plot:
        i_plot_column += 1
        if i_plot_column >= n_plot_columns:
            i_plot_column, i_plot_row = 0, i_plot_row + 1
        if i_plot_row >= n_plot_rows:
            i_plot_row = 0

        # Flush page of plots if page full or if all plots have now been made:
        page_finished = i_plot_column == 0 and i_plot_row == 0 and n_plots_completed != 0
        all_plots_finished = (n_plots_completed == len(star_ids))
        if page_finished or all_plots_finished:
            fig.tight_layout(rect=(0, 0, 1, 0.925))
            fig.subplots_adjust(left=0.06, bottom=0.075, right=0.94, top=0.875, wspace=0.3, hspace=0.3)
            title_text = 'Demo for online LPV tool: 2018 Bulletin stars ' + \
                         str(i_first_plot) + ' through ' + str(n_plots_completed)
            fig.suptitle(title_text, color='dimgray', fontsize=18)
            fig.canvas.set_window_title(title_text)
            plt.show()
            fullpath = 'C:/Astro/AAVSO/Bulletin web project 2019/Demo ' +\
                       '{0:03d}'.format(i_first_plot) + '.png'
            plt.savefig(fullpath, orientation='landscape')
            i_plot_row, i_plot_column = 0, 0  # for next plot page (if any).
            i_first_plot = n_plots_completed + 1  # "
            print('Page plotted:', str(n_plots_completed), 'plots done.')


SUPPORT_________________________________________________ = 0


def get_bulletin_star_ids(path=None):
    """  Read file Bulletin2018.csv and extract star names only.
         (faster than reading in entire df_bulletin Dataframe for 'NAME' column, if names alone are needed)
    :param path: exactly where to find file Bulletin2018.csv (default is usually correct) [string].
    :return: all star names in AAVSO Bulletin 2018, order left unchanged from the csv file [list of strings]
    """
    if path is None:
        path = os.path.join(DATA_DIRECTORY, BULLETIN2018_FILENAME)
    with open(path) as f:
        lines = f.readlines()
    possible_star_ids = [line.split(',')[0].strip().upper() for line in lines]
    star_ids = [n for n in possible_star_ids if n not in ['NAME', '']]
    return star_ids


def get_df_bulletin(path=None):
    """  Delivers dataframe df_bulletin containing all data in a .csv file (usually Bulletin2018.csv).
    :param path: exactly where to find file Bulletin2018.csv (default is usually correct) [string].
    :return: df_bulletin [pandas Dataframe]
    """
    if path is None:
        path = os.path.join(DATA_DIRECTORY, BULLETIN2018_FILENAME)
    df_bulletin = pd.read_csv(path, sep=',', encoding="UTF-8", na_filter=False)
    df_bulletin['NAME'] = [name.strip().upper() for name in df_bulletin['NAME']]
    keep_rows = [name not in ['NAME', ''] for name in df_bulletin['NAME']]
    df_bulletin = df_bulletin[keep_rows]
    df_bulletin = df_bulletin.set_index('NAME', drop=False)
    return df_bulletin


def make_jd_preds(first_datetime, last_datetime, increment_days):
    """ Makes Julian Date list for user in make_short_prediction().
        If number of increments between first and last dates is non-integral, count from last_datetime.
    :param first_datetime: datetime to begin series [python Datetime object].
    :param last_datetime: datetime to end series [python Datetime object].
    :param increment_days: number of days between jd_preds [float].
    :return: list of Julian Dates [list of floats].
    """
    jd_first = jd_from_datetime_utc(first_datetime)
    jd_last = jd_from_datetime_utc(last_datetime)

    if jd_last < jd_first:
        return None
    elif jd_last == jd_first:
        return [jd_first]

    n_increments = int(floor((jd_last - jd_first) / increment_days))
    jd_preds = [jd_last - i * increment_days
                for i in range(n_increments + 1)][::-1]  # sorted by incr jds.
    return jd_preds


def get_local_df_nobs():
    """ Make df_nobs from locally stored .csv file. Avoids need to query VSX for number-of-observation
        data when they are needed (probably only when constructing a new demo bulletin).
    :return: dataframe of number of VSX observations of all types for each star_id [pandas DataFrame].
    """
    fullpath = os.path.join(DATA_DIRECTORY, DF_NOBS_FILENAME)
    if os.path.exists(fullpath):
        df_nobs = pd.read_csv(fullpath, sep=';', encoding="UTF-8", na_filter=False)
        df_nobs = df_nobs.set_index('star_id', drop=False)
    else:
        df_nobs = None
    return df_nobs


def get_local_df_obs():
    """  Make df_obs from locally stored .csv file. Avoids need to query VSX for observation data every
         time they are needed.
    :return: dataframe of VSX observation data stored in local .csv file [pandas Dataframe].
    """
    fullpath = os.path.join(DATA_DIRECTORY, DF_OBS_FILENAME)
    if os.path.exists(fullpath):
        df_obs = pd.read_csv(fullpath, sep=';', encoding="UTF-8", na_filter=False)
        df_obs = df_obs.set_index('obsID', drop=False)
    else:
        df_obs = None
    return df_obs


def jd_from_datetime_utc(datetime_utc=None):
    """  Converts a UTC datetime to Julian date. Imported from photrix (E. Dose).
    :param datetime_utc: date and time (in UTC) to convert [python datetime object]
    :return: Julian date corresponding to date and time [float].
    """
    if datetime_utc is None:
        return None
    datetime_j2000 = datetime(2000, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
    jd_j2000 = 2451544.5
    seconds_since_j2000 = (datetime_utc - datetime_j2000).total_seconds()
    return jd_j2000 + seconds_since_j2000 / (24*3600)


def datetime_utc_from_jd(jd=None):
    """  Converts a Julian Date to UTC datetime. Imported from photrix (E. Dose).
    :param jd: Julian date to be converted [float].
    :return: UTC datetime from Julian Date [python datetime object].
    """
    if jd is None:
        return datetime.now(timezone.utc)
    datetime_j2000 = datetime(2000, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)
    jd_j2000 = 2451544.5
    seconds_since_j2000 = 24 * 3600 * (jd - jd_j2000)
    return datetime_j2000 + timedelta(seconds=seconds_since_j2000)


def jd_now():
    """  Returns Julian date of moment this function is called. Imported from photrix (E. Dose).
    :return: Julian date for immediate present per system clock [float].
    """
    return jd_from_datetime_utc(datetime.now(timezone.utc))


def reorder_df_columns(df, left_column_list=None, right_column_list=None):
    """ Push desired columns to left edge of dataframe, others to right of dataframe.
           Other columns are left in between, in their original order.
           Largely cosmetic, but can be useful for debugging.
    :param df: any dataframe [pandas Dataframe]
    :param left_column_list: list of column names to go left [list of strings].
    :param right_column_list: list of column names to go right [list of strings].
    :return: dataframe with reordered columns [pandas DataFrame].
    """
    if left_column_list is None:
        left_column_list = []
    if right_column_list is None:
        right_column_list = []
    new_column_order = left_column_list +\
                       [col_name for col_name in df.columns
                        if col_name not in (left_column_list + right_column_list)] +\
                       right_column_list
    df = df[new_column_order]
    return df