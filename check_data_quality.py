## script checking data quality of raw data. dd august 8, 2012

from datetime import datetime, timedelta
from time import mktime
from calendar import timegm

from tables import *
import matplotlib.pyplot as plt
from numpy import average, linspace, sqrt, exp, nan
from scipy import optimize, ndimage

from hisparc.publicdb import download_data

DATAFILE = '50x_2012.h5'
# The datafile required is an HDF5 Data file as obtained with the download_data routine in the hisparc.publicdb library without downloading the blobs.
OUTPUTFILE = 'Logfile.h5'
# The outputfile is an HDF5 Data file using the class defined below

TIMESTEP = timedelta(hours=4) # This timestep determines the lenght of each block to be analyzed.
NR_BLOCKS = 24 / (TIMESTEP.seconds / 3600) # Used later on because as reference the average of a whole day is used.
STEPSIZE = 200 # Number of bins used in the histograms of the fit function
MARG_PEAKS = 25 # margin accepted in number of peaks, % of reference value
MARG_PLS_HT = 50 # margin accepted in shape analysis of pulse height histogram
MARG_PLS_INT = 50 # margin accepted in shape analysis of pulse integral histogram
PH_MPV_LOW = 120 # 76 mV
PH_MPV_HIGH = 439 # 250 mV
PI_MPV_LOW = 1200 # 760 mVns
PI_MPV_HIGH = 4390 # 2500 mVns


class Logs(IsDescription):
    ##Class used to create table LOGS in OUTPUTFILE.
    ref_date = Time32Col() # date of reference used in analyzing
    start = Time32Col() # start of the period diagnosed
    useful = BoolCol() # 1 for useful, 0 for conspicuous
    nr_ev = Int16Col(shape=4) # Number of events indexed to 100
    ph_diff = Int16Col(shape=4) # Events per ADC value of histogram, perct
    pi_diff = Int16Col(shape=4) # of ref period, accumulated over histogr.
    ph_mpv = Float32Col(shape=4) # Peak location of histogram in ADC counts
    pi_mpv = Float32Col(shape=4)


def first_download(station): # Not used in the current script; could be integrated later on
    """ Downloads data from HDF5

    Including data for referencedate (if necessary)
    """
    group = '/s%s' % station

    with openFile(filepath, 'a') as datafile:
        if group not in datafile:
            download_data(datafile, group, station, START, STOP)

            nextday = REFDATE + timedelta(days=1)
            if REFDATE < START or nextday > STOP:
                download_data(datafile, '/s%s' % station, station, REFDATE, nextday)


def make_timestamp(date):
    """ Turns date between brackets into timestamp

    Needed for reading data
    """
    timetuple = date.utctimetuple()
    timestamp = timegm(timetuple)
    return timestamp


def count_peaks(peaks, index):
    """ Counts the number of peaks in the group 'events' of DATAFILE

    Returns this number as an integer
    """
    peaks = peaks[:, index]
    cleaned_peaks = peaks.compress(peaks > 0)
    total = cleaned_peaks.sum()
    return total


def ref_peaks():
    """ Returns number of peaks for the selected referencedate

    Average for all day of reference for each scintillator.
    Returns refs = tuple of 4 integers
    """
    date = REFDATE
    one_day = timedelta(days=1)

    t0 = make_timestamp(date)
    t1 = make_timestamp(date + one_day)
    peaks = events.readWhere('(t0 <= timestamp) & (timestamp < t1)')['n_peaks']

    refa = max(1, (count_peaks(peaks, 0) / NR_BLOCKS))
    refb = max(1, (count_peaks(peaks, 1) / NR_BLOCKS))
    refc = max(1, (count_peaks(peaks, 2) / NR_BLOCKS))
    refd = max(1, (count_peaks(peaks, 3) / NR_BLOCKS))

    refs = [refa, refb, refc, refd]
    return refs


def compare_peaks(date, ref_pks):
    """Determines if number of peaks is within range of reference number

        divides the total numbers for each scintillator (and each TIMESTEP)
        by those of REFDATE.
        Output:peakdiff, tuple of deviations per scintillator, % of amounts on REFDATE
    """

    t0 = make_timestamp(date)
    t1 = make_timestamp(date + TIMESTEP)

    peaks = events.readWhere('(t0 <= timestamp) & (timestamp < t1)')['n_peaks']

    compa = (100 * (count_peaks(peaks, 0)) / ref_pks[0])
    compb = (100 * (count_peaks(peaks, 1)) / ref_pks[1])
    compc = (100 * (count_peaks(peaks, 2)) / ref_pks[2])
    compd = (100 * (count_peaks(peaks, 3)) / ref_pks[3])

    peakdiff = [compa, compb, compc, compd]
    return peakdiff


def get_ph(date):
    """Get the pulseheights for the specified time interval TIMESTEP

     returns
    """

    t0 = make_timestamp(date)
    t1 = make_timestamp(date + TIMESTEP)

    # Get the pulseheights in the specified time interval

    try:
        ph = events.readWhere('(t0 <= timestamp) & (timestamp < t1)')['pulseheights']
    except NoSuchNodeError:
        ph = array([[0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.float) # If there is no data all ph will be set to zero
        print 'No Such Node, ph is set to 0'

    # Rearrange the data

    ph_fit = zip(*ph)
    return ph, ph_fit


def get_pi(date):
    """Get the pulseheightintegral for the specified time interval TIMESTEP

    """
    t0 = make_timestamp(date)
    t1 = make_timestamp(date + TIMESTEP)

    # Get the pulseheightintegrals in the specified time interval

    try:
        pi = events.readWhere('(t0 <= timestamp) & (timestamp < t1)')['integrals']
    except NoSuchNodeError:
        pi = array([[0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.float) # If there is no data all pi will be set to zero
        print 'No Such Node, pi is set to 0'

    # Rearrange the data

    pi_fit = zip(*pi)
    return pi, pi_fit


def find_MPV(fit, conv=1):
    """Determine the ADC value of the MPV for each detector

    For pulseheights use conv = 1
    For pulseintegral values use conv = 10

    Returns a tuple of MPV values

    """
    MPV1, MPV2, MPV3, MPV4 = [], [], [], []

    for i in range(4): ##4 DETECTOREN
        try:
            y, bins, patches = plt.hist(fit[i], bins=linspace(0, conv * 2000 / 0.57, STEPSIZE), label="Data")

            x = bins[:-1] + .5 * (bins[1] - bins[0])

            # A first approximation of the position of the MPV and the minimum to the left is made. Used for fit bounds and initial guess of the parameters.
            guess_max_x = ndimage.extrema(y.compress((conv * 150 <= x) & (x < conv * 410)))[3][0] + len(y.compress(x < conv * 150))
            guess_max_x = x[guess_max_x]
            guess_min_x = ndimage.extrema(y.compress((conv * 53 <= x) & (x < guess_max_x)))[2][0] + len(y.compress(x < conv * 53))
            guess_min_x = max([conv * 120, x[guess_min_x]])
            guess_max_y = ndimage.extrema(y.compress((conv * 150 <= x) & (x < conv * 410)))[1]
            guess_min_y = ndimage.extrema(y.compress((conv * 53 <= x) & (x < guess_max_x)))[0]

            g0 = guess_max_y
            g1 = guess_max_x

            # The fit range. Since the right side of the ph histogram is most similar to a Gauss function, we do not want to set bound 2 too small.
            bound1 = guess_min_x + (guess_max_x - guess_min_x) / 4

            if (guess_max_x - guess_min_x) <= conv * 50:
                bound2 = guess_max_x + (guess_max_x - guess_min_x) * 1.5
            else:
                bound2 = guess_max_x + (guess_max_x - guess_min_x)

            # Fit a Gaussian to the peak.
            f = lambda x, a, b, c: a * exp(-((x - b) ** 2) / c)
            x2 = x.compress((bound1 <= x) & (x < bound2))
            y2 = y.compress((bound1 <= x) & (x < bound2))
            popt, pcov = optimize.curve_fit(f, x2, y2, (g0, g1, (guess_max_x - guess_min_x) / 2), sigma=sqrt(y2))

            # Find the x value of the peak:
            peak = popt[1]
            print ('The MPV of the pulseheight of detector %d lies at %d ADC' %
                   (i + 1, peak))

        # Safety nets: prevent the code from terminating if an error occurs
        except IndexError:
            peak = nan
            print 'Apparently there is no data'
        except RuntimeError:
            peak = nan
            print 'Unable to make a good fit'
        except TypeError:
            peak = nan
            print 'Not enough data to make a proper fit' # When it uses less than three bins to make the fit you get this error
        except ValueError:
            peak = nan
            print 'Value Error Occured'

        if i == 0:
            MIP1 = peak
        elif i == 1:
            MIP2 = peak
        elif i == 2:
            MIP3 = peak
        elif i == 3:
            MIP4 = peak

    mpv = (MIP1, MIP2, MIP3, MIP4)
    return mpv


def phdiff(date, ph_ref):
    """Checks if the Pulse Height histogram differs in shape from the one on REFDATE

    """

    ph_check = get_pulse(date, field='pulseheights')

    errors = []

    try:
        if len(ph_check) > 0:
            diff1_av = pulse_check(ph_check[0], ph_ref[0])
            diff2_av = pulse_check(ph_check[1], ph_ref[1])
            diff3_av = pulse_check(ph_check[2], ph_ref[2])
            diff4_av = pulse_check(ph_check[3], ph_ref[3])
        else:
            diff1_av = -1
            diff2_av = -1
            diff3_av = -1
            diff4_av = -1
    except:
        diff1_av = -1
        diff2_av = -1
        diff3_av = -1
        diff4_av = -1

    diff_av = [diff1_av, diff2_av, diff3_av, diff4_av]

    return diff_av


def pidiff(date, pi_ref):
    """Checks if the Pulse Integral histogram differs in shape from the one on REFDATE

    """

    pi_check = get_pulse(date, field='integrals')

    errors = []

    try:
        if len(pi_check) > 0:
            diff1_av = pulse_check(pi_check[0], pi_ref[0])
            diff2_av = pulse_check(pi_check[1], pi_ref[1])
            diff3_av = pulse_check(pi_check[2], pi_ref[2])
            diff4_av = pulse_check(pi_check[3], pi_ref[3])
        else:
            diff1_av = -1
            diff2_av = -1
            diff3_av = -1
            diff4_av = -1
    except:
        diff1_av = -1
        diff2_av = -1
        diff3_av = -1
        diff4_av = -1

    diff_av = [diff1_av, diff2_av, diff3_av, diff4_av]

    return diff_av


def get_pulse(date, ref=0, field='pulseheights'):
    """ Get pulseintegral or pulseheights data and averages for one hour

    """
    try:
        t0 = timegm(date.utctimetuple())
        if ref == 0:
            t1 = date + TIMESTEP
            hours = TIMESTEP.seconds / 3600
        else:
            t1 = date + timedelta(hours=24)
            hours = 24
        # The reference data is always based on 24 hours; so the reference data is not influenced by a possible day / night difference
        t1 = timegm(t1.utctimetuple())
        p = events.readWhere('(timestamp >= t0) & (timestamp < t1)')[field]
        p1, p2, p3, p4 = [], [], [], []

        if p[0, 0] != -1:
            p1 = plt.hist(p[:,0], bins=45, range=[877, 8770], histtype='step', log='True', color='k')
            p1 = p1[0] / float(hours)
        if p[0, 1] != -1:
            p2 = plt.hist(p[:, 1], bins=45, range=[877, 8770], histtype='step', log='True', color='r')
            p2 = p2[0] / float(hours)
        if p[0, 2] != -1:
            p3 = plt.hist(p[:, 2], bins=45, range=[877, 8770], histtype='step', log='True', color='g')
            p3 = p3[0] / float(hours)
        if p[0, 3] != -1:
            p4 = plt.hist(p[:, 3], bins=45, range=[877, 8770], histtype='step', log='True', color='c')
            p4 = p4[0] / float(hours)

        pulses = [p1, p2, p3, p4]
    except:
        print 'This date contains no data; \n', (date)
        pulses = [-1, -1, -1, -1]

    return pulses


def pulse_check(check, ref):

    if min(check) >= 0:
        pi = check
        diff= []
        for i in range(45):
            diff.append(100 * (abs(pi[i] - ref[i]) / ref[i]))
        diff_av = average(diff)
    else:
        if ref == []:
            diff_av = 0
        else:
            diff_av = -1

    return diff_av


if __name__ == '__main__':
    stationlist = [501, 502, 503, 504, 505, 506, 509]

    # For the analysis of the stations in the current station list three different reference dates were used as seen below.
    # The list can further be filled and if necessary new refdates have to be included.
    for i in stationlist:
        STATION = i
        if i == 503:
            REFDATE = datetime(2012, 2, 2)   #Date used to get reference values
        elif i == 505:
            REFDATE = datetime(2012, 4, 28)   #Date used to get reference values
        else:
            REFDATE = datetime(2012, 1, 1)   #Date used to get reference values
        START = REFDATE #Date to start the analysis; this is now the reference date itself, but this could be changed.
        STOP = datetime(2012, 7, 15) #Date to stop, will not be analyzed

        data = openFile(DATAFILE, 'r')
        events = data.getNode('/s%d' % STATION, 'events')

        output = openFile(OUTPUTFILE, 'a')

        try:
            table = output.getNode('/s%d' % STATION, 'Logs')
        except NoSuchNodeError:
            group = output.createGroup('/', 's%d' % STATION)
            table = output.createTable(group, 'Logs', Logs)
        logs = table.row

        refpks = ref_peaks()
        ph_ref = getph(REFDATE, ref=1)
        pi_ref = getpi(REFDATE, ref=1)

        date1 = START

        while date1 < STOP:
            print(date1)
            ph, ph_fit  = get_ph(date1)
            pi, pi_fit  = get_pi(date1)
            logs['start'] = mktime(date1.timetuple())
            logs['nr_ev'] = compare_peaks(date1, refpks)
            logs['ph_diff'] = phdiff(date1, ph_ref)
            logs['pi_diff'] = pidiff(date1, pi_ref)
            logs['ph_mpv'] = find_MPV(ph_fit)
            logs['pi_mpv'] = find_MPV(pi_fit, conv=10)
            logs['ref_date'] = mktime(REFDATE.timetuple())
            if (min(logs['nr_ev']) < (100 - MARG_PEAKS)) or (max(logs['nr_ev']) > (100 + MARG_PEAKS)):
                logs['useful'] = False
            elif (max(logs['ph_diff']) > MARG_PLS_HT) or (min(logs['ph_diff']) < 0):
                logs['useful'] = False
            elif (max(logs['pi_diff']) > MARG_PLS_INT) or (min(logs['pi_diff']) < 0):
                logs['useful'] = False
            elif (max(logs['ph_mpv']) > PH_MPV_HIGH) or (min(logs['ph_mpv']) < PH_MPV_LOW):
                logs['useful'] = False
            elif (max(logs['pi_mpv']) > PI_MPV_HIGH) or (min(logs['pi_mpv']) < PI_MPV_LOW):
                logs['useful'] = False
            else:
                logs['useful'] = True
            logs.append()
            date1 = date1 + TIMESTEP
            plt.close('all')
            # This closes the figures made bij the 'hist' functions. Otherwise the memory is slowly filled  until python crashes on a memory error
            print(STATION)

        data.close()
        output.close()
