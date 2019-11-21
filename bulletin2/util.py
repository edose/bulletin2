import functools
from datetime import datetime, timezone, timedelta
import os
from math import sin, cos, pi, sqrt
# import numpy as np
import pandas as pd
# import requests
# from bs4 import BeautifulSoup
import urllib.request
import statsmodels.api as sm

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

DATA_DIRECTORY = 'C:/Dev/bulletin2/data'
BULLETIN2018_FILENAME = 'Bulletin2018.csv'
DF_OBS_FILENAME = 'df_obs.csv'

DAYS_PER_YEAR = 365.25
HTTP_OK_CODE = 200  # "OK. The request has succeeded."
JD0 = 2458484.5  # Reference Julian Date for lightcurve phases. Jan 1 2019 00:00 utc
FIT_END_DATE = datetime(2018, 12, 1, 0, 0, 0).replace(tzinfo=timezone.utc)  # Dec 1 2018.
PREDICT_START_DATE = datetime(2019, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)  # Jan 1 2019.
PREDICT_END_DATE = datetime(2020, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc)  # Jan 1 2020.
FIT_PERIODS_TO_CAPTURE = 7  # LPV periods to capture from VSX, counting back from FIT_END_DATE.
FIT_PERIODS = 5
VSX_OBSERVATIONS_HEADER = 'https://www.aavso.org/vsx/index.php?view=api.delim'
VSX_DELIMITER = '@@@'  # NB: ',' fails as obsName values already have a comma.
MIN_MAG = 0.1  # magnitudes must be greater than this numeric value (mostly to remove zeroes).
MAX_MAG = 20.0  # maximum reasonable magnitude numeric value (subject to later revision for meter scopes).


@functools.lru_cache(maxsize=128, typed=False)
def get_vsx_obs(star_id, jd_start, jd_end=None):
    """
    Downloads observations from AAVSO's webobs for ONE star, returns dataframe of results.
       If star not in AAVSO's webobs site, return a dataframe with no rows.
       Columns: target_name, date_string, filter, observer, jd, mag, error.
    :param star_id: the STAR id (not the fov's name).
    :param jd_start: optional Julian date.
    :param jd_end: optional JD.
    :return: DataFrame containing data for 1 star, 1 row per observation downloaded,
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


def capture_vsx_data():
    """  From VSX API, capture test data from VSX API, make 2 dataframes:
            df_obs, with one row per (screened & weighted) observation, and
            df_stars, with one row per star, holding summary data for that star.
        Write these dataframes locally as .csv files, for use in testing magnitude-projection model(s).
        The idea is to get this data one time, holding it locally in csv files,
        thus saving us testing time, and saving AAVSO server time.
    :return:
    """
    df_obs = pd.DataFrame()  # master dataframe of all obs and metadata, to be written to local .csv file.
    # df_stars = pd.DataFrame()  # summary data, one row per star (we may end up not needing this).

    star_names = get_bulletin_star_names()
    df_bulletin = get_df_bulletin()
    for name in star_names:
        period = float(df_bulletin.loc[name, 'PERIOD'])
        jd_end = jd_now()  # collect all obs through 2019 as well (2019 obs to evaluate predictions).
        jd_start = jd_from_datetime_utc(FIT_END_DATE) - FIT_PERIODS_TO_CAPTURE * period
        df_vsx = get_vsx_obs(name, jd_start=jd_start, jd_end=jd_end)
        print('   ', str(len(df_vsx)), 'obs downloaded from VSX for JD range:',
              str(int(jd_start)), 'to', str(int(jd_end)) + ':')
        df_screened = screen_star_obs(df_vsx)
        df_this_star = add_obs_weights(df_screened)
        df_obs = pd.concat([df_obs, df_this_star])
        print('   ', str(len(df_this_star)), 'obs added for df_obs running total of',
              str(len(df_obs)) + '.')

    csv_fullpath = os.path.join(DATA_DIRECTORY, DF_OBS_FILENAME)
    df_obs.to_csv(csv_fullpath, sep=';', quotechar='"', encoding='UTF-8',
                  quoting=2, index=False)  # quoting=2-->quotes around non-numerics.
    print('df_obs dataframe written to', csv_fullpath)


def screen_star_obs(df):
    """  Screen observations for quality and relevance to V-band magnitude projection.
    :param df: this star's unscreened observations [pandas DataFrame].
    :return: this star's screened observations [pandas DataFrame].
    """
    mag_not_null = ~ df['mag'].isnull()
    df = df[mag_not_null]
    mag_value_ok = pd.Series([MIN_MAG <= float(mag) <= MAX_MAG for mag in df['mag']], index=df.index)
    band_not_v_like = df['band'].isin(['V', 'Vis.', 'TG'])
    not_fainter_than = df['fainterThan'] == '0'
    obstype_ok = df['obsType'].isin(['Visual', 'CCD', 'DSLR', 'VISDIG'])
    validation_ok = df['val'].isin(['V', 'Z'])
    obs_to_keep = mag_value_ok & band_not_v_like & not_fainter_than & obstype_ok & validation_ok
    df = df[obs_to_keep]
    return df


def add_obs_weights(df):
    """  Double the weights for all transformed (CCD) observations,
         relative to all other observations.
    :param df: this star's screened observations [pandas DataFrame]
    :return: same dataframe with added column 'weight' [pandas DataFrame].
    """
    weights = [2.0 if tr == '1' else 1.0 for tr in df['transformed']]
    # Next line removed, as validation seems to happen only for Visual data (why??),
    #    and thus columns 'val' and 'transformed' are overly correlated.
    # weights = [0.5 * wt if val == 'Z' else wt for (wt, val) in zip(weights, df['val'])]
    df.loc[:, 'weight'] = weights
    return df


def get_bulletin_star_names(path=None):
    """  Read file Bulletin2018.csv and extract star names only.
    :param path: exactly where to find file Bulletin2018.csv [string].
    :return: all star names in AAVSO Bulletin 2018, order left unchanged from the csv file.
    """
    if path is None:
        path = os.path.join(DATA_DIRECTORY, BULLETIN2018_FILENAME)
    with open(path) as f:
        lines = f.readlines()
    raw_star_names = [line.split(',')[0].strip().upper() for line in lines]
    star_names = [n for n in raw_star_names if n not in ['NAME', '']]
    return star_names


def get_df_bulletin(path=None):
    if path is None:
        path = os.path.join(DATA_DIRECTORY, BULLETIN2018_FILENAME)
    df_bulletin = pd.read_csv(path, sep=',', encoding="UTF-8", na_filter=False)
    df_bulletin = df_bulletin[df_bulletin['NAME'] != 'NAME']
    df_bulletin = df_bulletin.set_index('NAME', drop=False)
    return df_bulletin


def make_df_x(jd_list, fit_period, jd0=JD0):
    df_x = pd.DataFrame()
    df_x['v_const'] = len(jd_list) * [1.0]  # ca. the mean V mag
    df_x['dv_year'] = [(jd - jd0) / DAYS_PER_YEAR for jd in jd_list]  # linear drift in V mag per year.
    phases = [(jd - jd0) / fit_period for jd in jd_list]  # numerically negative, since jd < JD0, usually.
    df_x['sin1'] = [sin((2 * pi) * 1.0 * ph) for ph in phases]  # first Fourier term
    df_x['cos1'] = [cos((2 * pi) * 1.0 * ph) for ph in phases]  # "
    df_x['sin2'] = [sin((2 * pi) * 2.0 * ph) for ph in phases]  # second Fourier term
    df_x['cos2'] = [cos((2 * pi) * 2.0 * ph) for ph in phases]  # "
    return df_x, jd0


def fit_one_star(star_name, df_obs=None, df_bulletin=None, period_factor=1.0, jd0=JD0):
    # Ensure that we have the required observation and star data:
    if df_obs is None:
        df_obs = get_local_df_obs()
    if df_bulletin is None:
        df_bulletin = get_df_bulletin()

    # Prepare dataframe for this star:
    df = df_obs[df_obs['starName'] == star_name]
    df = df[['JD', 'band', 'mag', 'obsID', 'obsType', 'weight']]
    df['JD'] = [float(jd) for jd in df['JD']]  # because strings are inherited from VSX.
    df['mag'] = [float(mag) for mag in df['mag']]  # because strings are inherited from VSX.
    df = df.set_index('obsID', drop=False)
    df = df.sort_values(by='JD')

    # Get variable-star period from 2018 Bulletin, then set the period actually used in this fit:
    bulletin_period = float(df_bulletin.loc[star_name, 'PERIOD'])
    this_fit_period = period_factor * bulletin_period  # because best period may differ from bulletin's.

    # Prepare dataframes for weighted multivariate regression fit:
    df_jd = pd.DataFrame()
    df_jd['JD'] = df.copy()['JD']
    df_jd.index = df.index
    df_x, jd0 = make_df_x(df['JD'], this_fit_period)
    df_x.index = df.index
    df_y = pd.DataFrame()
    df_y['y'] = df.copy()['mag']
    df_y.index = df.index
    df_wt = pd.DataFrame()
    df_wt['weight'] = df.copy()['weight']
    df_wt.index = df.index

    # Select only obs Dec 1 2018 or before, :
    jd_end = jd_from_datetime_utc(FIT_END_DATE)
    jd_start = jd_end - FIT_PERIODS * this_fit_period
    to_keep = (df['JD'] >= jd_start) & (df['JD'] < jd_end)
    df_jd = df_jd[to_keep]
    df_x = df_x[to_keep]
    df_y = df_y[to_keep]
    df_wt = df_wt[to_keep]
    # The following commented-out lines are reserved in case we need to fit shifts between V, Vis., etc:
    # df_fit['dvis'] = [1.0 if obsType == 'Vis.' else 0.0
    #     for obsType in df_fit['obsType']]  # pseudo-categorical.
    # df_fit['dtg'] = [1.0 if obsType == 'TG' else 0.0 for obsType in df_fit['obsType']]  # "
    # df_fit['dvisdig'] = [1.0 if obsType == 'VISDIG' else 0.0 for obsType in df_fit['obsType']]  # "

    # For each obs within most recent period, double its pre-existing weight:
    days_before_jd_end = jd_end - df_jd['JD']
    weights = [2.0 * wt if db <= this_fit_period else 1.0 * wt
               for (wt, db) in zip(df_wt['weight'], days_before_jd_end)]

    # Do fit:
    # for weighted regression we'll want sm.WLS rather than OLS (ordinary least squares).
    # result = sm.OLS(df_y, sm.add_constant(df_x)).fit()
    result = sm.WLS(df_y, df_x, weights).fit()

    # Calculate first-order amplitude:
    # amplitude_1 = sqrt(result.params[2]**2 + result.params[3]**2)

    # print(result.summary())

    # Predict V magnitude for periodic dates:
    # jd_start = jd_from_datetime_utc(
    #     datetime(2019, 1, 1, 0, 0, 0).replace(tzinfo=timezone.utc))  # Jan 1 2019  00:00 utc
    # jds_to_predict = [jd_start + 10 * i for i in range(37)]  # 10-day spacing through all of 2019.
    # df_x_predict = make_df_x(jds_to_predict, this_fit_period)
    # prediction = result.predict(df_x_predict)
    # for jd, pred in zip(jds_to_predict, prediction):
    #     print(jd, pred)
    # return list(zip(jds_to_predict, prediction))

    return result, jd0  # to examine for data to store in df_fit_results.


def process_one_star(star_name, df_obs=None, df_bulletin=None, jd0=JD0, quiet=False):
    if star_name not in df_bulletin.index:
        print(' >>>>> Star', star_name, 'is not in the 2018 Bulletin.')
        return None
    bulletin_period = float(df_bulletin.loc[star_name, 'PERIOD'])
    fit_factors = [0.9, 1.0, 1.1]  # these must be evenly spaced.
    all_results = []
    for fit_factor in fit_factors:
        result, jd0 = fit_one_star(star_name, df_obs, df_bulletin, period_factor=fit_factor, jd0=jd0)
        all_results.append(result)
    if not quiet:
        for i, result in enumerate(all_results):
            print(i, fit_factors[i], all_results[i].rsquared)
    r2_1, r2_2, r2_3 = tuple([all_results[i].rsquared for i in range(3)])
    if r2_2 > (r2_1 + r2_3) / 2:  # if quadratic through 3 points is concave downward:
        # Lagrange-interpolate for best period:
        x1, x2, x3 = tuple(fit_factors)
        y1, y2, y3 = tuple(all_results[i].rsquared for i in range(3))
        denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
        a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
        b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom
        # c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom
        best_fit_factor = -b / (2 * a)  # x-value at interpolated maximum.
        if x1 <= best_fit_factor <= x3:
            # Period seems to be within range, so we will use the interpolated period:
            best_period = best_fit_factor * bulletin_period
            period_sign = '~'
            best_result, jd0 = fit_one_star(star_name, df_obs, df_bulletin,
                                            period_factor=best_fit_factor, jd0=jd0)
            best_r_squared = best_result.rsquared
        elif best_fit_factor < x1:
            # Period seems to be below range:
            best_period = fit_factors[0] * bulletin_period
            period_sign = '<'
            best_r_squared = all_results[0].rsquared
            best_result = all_results[0]
        else:
            # Period seems to be above range:
            best_period = fit_factors[2] * bulletin_period
            period_sign = '>'
            best_r_squared = all_results[2].rsquared
            best_result = all_results[2]
    else:
        # R^2 curve is concave up, so simply take the highest R^2 and mark with sign of uncertainty:
        r2_list = [all_results[i].rsquared for i in range(3)]
        i_max = r2_list.index(max(r2_list))
        best_period = fit_factors[i_max] * bulletin_period
        period_sign = '?'
        best_r_squared = all_results[i_max].rsquared
        best_result = all_results[i_max]

    if not quiet:
        print(star_name.ljust(16), '   Best: ',
              '   factor' + period_sign + '{0:.3f}'.format(best_period / bulletin_period),
              '   P=' + '{0:.1f}'.format(best_period),
              '   R^2=' + '{0:.3f}'.format(best_r_squared))

    # Make results to return to calling function:
    amplitude_1 = 2.0 * sqrt(best_result.params['sin1']**2 + best_result.params['cos1']**2)
    result_dict = {'star_name': star_name,
                   'nobs': best_result.nobs,
                   'period_sign': period_sign,
                   'bulletin_period': bulletin_period,
                   'period_factor': best_period / bulletin_period,
                   'best_period': best_period,
                   'jd0': jd0,
                   'r_squared': best_r_squared,
                   'amplitude_1': amplitude_1,
                   'params': best_result.params,
                   'params_se': best_result.bse,
                   'condition_number': best_result.condition_number,
                   'se': best_result.mse_resid,
                   'fit_result': best_result}
    return result_dict


def process_all_stars(df_obs=None, df_bulletin=None, jd0=JD0):
    if df_obs is None:
        df_obs = get_local_df_obs()
    if df_bulletin is None:
        df_bulletin = get_df_bulletin()
    result_dict_list = []
    star_names = df_bulletin['NAME']

    for star_name in star_names[0:5]:
        result_dict = process_one_star(star_name, df_obs, df_bulletin, jd0=JD0, quiet=True)
        result_dict_list.append(result_dict)
        print(star_name.ljust(8),
              '{0:.3f}'.format(result_dict['r_squared']),
              '{0:.3f}'.format(result_dict['se']))

    df_fit_results = pd.DataFrame(result_dict_list).set_index('star_name', drop=False)
    return df_fit_results




SUPPORT_________________________________________________ = 0


def get_local_df_obs():
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


DETRITUS________________________________________ = 0

# def get_obs_one_star(star_id, num_obs=200, jd_start=EARLIEST_FIRST_JD, jd_end=None):
#     """  Return a dataframe with all WebObs observations for 1 star, 1 row per observation downloaded.
#          Adapted from photrix.get_aavso_raw_table().
#          If no observations, return dataframe with no rows.
#     :param star_id: the STAR id [string].
#     :param jd_start: optional Julian date, default = January 1 2015. [float].
#     :param jd_end: optional JD, default = now [float].
#     :return: WebObs observation data for 1 star, 1 row per obs [pandas Dataframe].
#     """
#     safe_star_name = star_id.replace("+", "%2B").replace(" ", "+")
#     url = "https://www.aavso.org/apps/webobs/results/?star=" + safe_star_name + "&obs_types=vis+ccd+dslr"
#     if jd_start is not None:
#         url += "&start=" + str(jd_start)
#     if jd_end is not None:
#         url += "&end=" + str(jd_end)
#     # TODO: Try to use requests Session objects, for performance.
#     # print('get_aavso_webobs_raw_table() >' + url + '<')
#     r = requests.get(url)
#     obs_list = []
#     if r.status_code == HTTP_OK_CODE:
#         soup = BeautifulSoup(r.text, 'html.parser')
#         obs_lines = soup.find_all('tr', class_='obs')  # NB: "class_" not "class" (reserved).
#         for line in obs_lines:
#             cells = line.find_all('td')
#             cell_strings = [cell.text for cell in cells]
#             obs_list.append(cell_strings)
#     df = pd.DataFrame(obs_list)
#     return df
