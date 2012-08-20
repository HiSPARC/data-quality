## script checking data quality of raw data. dd august 8, 2012

## import

#import datetime
from datetime import datetime, timedelta
#import time
from time import mktime
from calendar import timegm
#import tables
from tables import *
#from hisparc import publicdb
from hisparc.publicdb import download_data
import matplotlib.pyplot as plt
#from numpy import *
from numpy import average, linspace, sqrt, exp, nan
from scipy import optimize, ndimage

# define globals

DATAFILE = '50x_2012.h5'
# The datafile required is an HDF5 Data file as obtained with the download_data routine in the hisparc.publicdb library without downloading the blobs.
OUTPUTFILE = 'Logfile.h5'
# The outputfile is an HDF5 Data file using the class defined below

stationlist = [501, 502, 503, 504, 505, 506, 509]

# For the analysis of the stations in the current station list three differen reference dates were used as seen below.
# The list can further be filled and if necessary new refdates have to be included.
for i in stationlist:
    STATION = i
    if i == 503:
        REFDATE = datetime(2012,2,2)   #Date used to get reference values
    elif i == 505:
        REFDATE = datetime(2012,4,28)   #Date used to get reference values
    else:
        REFDATE = datetime(2012,1,1)   #Date used to get reference values
    START = REFDATE #Date to start the analysis; this is now the reference date itself, but this could be changed.
    STOP = datetime(2012,7,15) #Date to stop, will not be analyzed
    TIMESTEP = timedelta(hours=4) # This timestep determines the lenght of each block to be analyzed.
    NR_BLOCKS = 24 / (TIMESTEP.seconds / 3600) # Used later on because as reference the average of a whole day is used.
    STEPSIZE = 200 #Number of bins used in the histograms of the fit function
    MARG_PEAKS = 25   #margin accepted in number of peaks, % of reference value
    MARG_PLS_HT = 50   #margin accepted in shape analysis of pulse height histogram
    MARG_PLS_INT = 50  #margin accepted in shape analysis of pulse integral histogram
    PH_MPV_LOW = 120 # 76 mV
    PH_MPV_HIGH = 439 # 250 mV
    PI_MPV_LOW = 1200 # 760 mVns
    PI_MPV_HIGH = 4390 # 2500mVns

    # define class

    class Logs(IsDescription):
        ##Class used to create table LOGS in OUTPUTFILE.
        ref_date = Time32Col()   #date of reference used in analyzing
        start = Time32Col()  #start of the period diagnosed
        useful = BoolCol()   #1 for useful, 0 for conspicuous
        nr_ev = Int16Col(shape=4)     #Number of events indexed to 100
        ph_diff = Int16Col(shape=4)   #Events per ADC value of histogram, perct
        pi_diff = Int16Col(shape=4)   # of ref period, accumulated over histogr.
        ph_mpv = Float32Col(shape=4)   #Peak location of histogram in ADC counts
        pi_mpv = Float32Col(shape=4)

    # define routines

    def first_download(): # Not used in the current script; could be integrated later on
        """downloads data from HDF5

            including data for referencedate (if necessary)
        """
        group = '/s%s' % STATION

        with openFile(DATAFILE, 'a') as datafile:
            if group not in datafile:
                download_data(datafile, group, STATION, START, STOP)

                nextday = REFDATE + timedelta(days = 1)
                if REFDATE < START or nextday > STOP:
                    download_data(datafile, '/s%s' %(STATION), STATION, REFDATE, nextday)


    def make_timestamp(date):
        """turns date between brackets into timestamp,

            needed for reading data
        """
        timetuple = date.utctimetuple()
        timestamp = timegm(timetuple)
        return timestamp


    def count_peaks(peaks, index):
        """Counts the number of peaks in the group 'events' of DATAFILE

        Returns this number as an integer
        """
        peaks = peaks[:,index]
        cleaned_peaks = peaks.compress(peaks > 0)
        total = cleaned_peaks.sum()
        return total


    def ref_peaks():
        """Returns number of peaks for the selected referencedate

            Average for all day of reference for each scintillator.
            Returns refs = tuple of 4 integers
        """
        date = REFDATE
        one_day = timedelta(days=1)

        t0 = make_timestamp(date)
        t1 = make_timestamp(date + one_day)
        peaks = events.readWhere('(t0 <= timestamp) & (timestamp < t1)')['n_peaks']

        refa = max(1,(count_peaks(peaks, 0)/NR_BLOCKS))
        refb = max(1,(count_peaks(peaks, 1)/NR_BLOCKS))
        refc = max(1,(count_peaks(peaks, 2)/NR_BLOCKS))
        refd = max(1,(count_peaks(peaks, 3)/NR_BLOCKS))

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

        compa = (100*(count_peaks(peaks, 0)) / ref_pks[0])
        compb = (100*(count_peaks(peaks, 1)) / ref_pks[1])
        compc = (100*(count_peaks(peaks, 2)) / ref_pks[2])
        compd = (100*(count_peaks(peaks, 3)) / ref_pks[3])

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
            ph = array([[0,0,0,0],[0,0,0,0]], dtype=np.float) # If there is no data all ph will be set to zero
            print 'No Such Node, ph is set to 0'

        # Rearrange the data

        ph_fit = zip(*ph)
        return ph, ph_fit


    def find_MPV(ph, ph_fit):
        """Determine the ADC value of the MPV for each detector

        return of tuple of MPV values
        """

        MPV1, MPV2, MPV3, MPV4 = [], [], [], []

        for i in range(4): ##4 DETECTOREN
            try:
                ph0 = ph_fit[i] #First fetches the pulseheight of the first detector and then that of the second.

                y, bins, patches = plt.hist(ph0, bins=linspace(0, 2000/0.57, STEPSIZE),  label="Data")

                x = bins[:-1] + .5 * (bins[1] - bins[0])

                # A first approximation of the position of the MPV and the minimum to the left is made. Used for fit bounds and initial guess of the parameters.
                guess_max_x = ndimage.extrema(y.compress((150<= x) & (x < 410)))[3][0] + len(y.compress(x < 150))
                guess_max_x = x[guess_max_x]
                guess_min_x = ndimage.extrema(y.compress((53 <= x) & (x < guess_max_x)))[2][0] + len(y.compress(x < 53))
                guess_min_x = max([120,x[guess_min_x]])
                guess_max_y = ndimage.extrema(y.compress((150 <= x) & (x < 410)))[1]
                guess_min_y = ndimage.extrema(y.compress((53 <= x) & (x < guess_max_x)))[0]

                g0 = guess_max_y
                g1 = guess_max_x

                # The fit range. Since the right side of the ph histogram is most similar to a Gauss function, we do not want to set bound 2 too small.
                bound1 = guess_min_x + (guess_max_x-guess_min_x)/4

                if (guess_max_x-guess_min_x) <= 50:
                    bound2 = guess_max_x + (guess_max_x-guess_min_x)*1.5
                else:
                    bound2 = guess_max_x + (guess_max_x-guess_min_x)

                # Fit a Gaussian to the peak.
                f = lambda x,a,b,c: a*exp(-((x-b)**2)/c)
                x2 = x.compress((bound1 <= x) & (x < bound2))
                y2 = y.compress((bound1 <= x) & (x < bound2))
                popt, pcov = optimize.curve_fit(f, x2, y2, (g0,g1,(guess_max_x-guess_min_x)/2), sigma=sqrt(y2))

                # Find the x value of the Peak:
                peak = popt[1]
                if i == 0:
                    MIP1 = peak
                elif i == 1:
                    MIP2 = peak
                elif i == 2:
                    MIP3 = peak
                elif i == 3:
                    MIP4 = peak
                print 'The MPV of detector',i + 1, 'lies at', peak, 'ADC'

            # Safety nets: prevent the code from terminating if an error occurs
            except IndexError:
                if i == 0:
                    MIP1 = nan
                if i == 1:
                    MIP2 = nan
                if i == 2:
                    MIP3 = nan
                if i == 3:
                    MIP4 = nan
                print 'Apparently there is no data'
            except RuntimeError:
                if i == 0:
                    MIP1 = nan
                if i == 1:
                    MIP2 = nan
                if i == 2:
                    MIP3 = nan
                if i == 3:
                    MIP4 = nan
                print 'Unable to make a good fit'
            except TypeError:
                if i == 0:
                    MIP1 = nan
                if i == 1:
                    MIP2 = nan
                if i == 2:
                    MIP3 = nan
                if i == 3:
                    MIP4 = nan
                print 'Not enough data to make a proper fit' # When it uses less than three bins to make the fit you get this error
            except ValueError:
                if i == 0:
                    MIP1 = nan
                if i == 1:
                    MIP2 = nan
                if i == 2:
                    MIP3 = nan
                if i == 3:
                    MIP4 = nan
                print 'Value Error Occured'

        mpv = (MIP1, MIP2, MIP3, MIP4)
        return mpv

    def get_pi(date):
        """Get the pulseheightintegral for the specified time interval TIMESTEP

        """
        t0 = make_timestamp(date)
        t1 = make_timestamp(date + TIMESTEP)

        # Get the pulseheightintegrals in the specified time interval

        try:
            pi = events.readWhere('(t0 <= timestamp) & (timestamp < t1)')['integrals']
        except NoSuchNodeError:
            pi = array([[0,0,0,0],[0,0,0,0]], dtype=np.float) # If there is no data all pi will be set to zero
            print 'No Such Node, pi is set to 0'

        # Rearrange the data

        pi_fit = zip(*pi)
        return pi, pi_fit


    def find_MPVINT(pi, pi_fit):
        """Determine the ADC value of the MPV for each detector

        return of tuple of MPV values
        """

        MPV1, MPV2, MPV3, MPV4 = [], [], [], []

        for i in range(4): ##4 DETECTOREN
            try:
                pi0 = pi_fit[i] #First fetches the pulseheight of the first detector and then that of the second.

                y, bins, patches = plt.hist(pi0, bins=linspace(0, 20000/0.57, STEPSIZE),  label="Data")

                x = bins[:-1] + .5 * (bins[1] - bins[0])

                # A first approximation of the position of the MPV and the minimum to the left is made. Used for fit bounds and initial guess of the parameters.
                guess_max_x = ndimage.extrema(y.compress((1500<= x) & (x < 4100)))[3][0] + len(y.compress(x < 1500))
                guess_max_x = x[guess_max_x]
                guess_min_x = ndimage.extrema(y.compress((530 <= x) & (x < guess_max_x)))[2][0] + len(y.compress(x < 530))
                guess_min_x = max([1200,x[guess_min_x]])
                guess_max_y = ndimage.extrema(y.compress((1500 <= x) & (x < 4100)))[1]
                guess_min_y = ndimage.extrema(y.compress((530 <= x) & (x < guess_max_x)))[0]

                g0 = guess_max_y
                g1 = guess_max_x

                # The fit range. Since the right side of the ph histogram is most similar to a Gauss function, we do not want to set bound 2 too small.
                bound1 = guess_min_x + (guess_max_x-guess_min_x)/4

                if (guess_max_x-guess_min_x) <= 500:
                    bound2 = guess_max_x + (guess_max_x-guess_min_x)*1.5
                else:
                    bound2 = guess_max_x + (guess_max_x-guess_min_x)

                # Fit a Gaussian to the peak.
                f = lambda x,a,b,c: a*exp(-((x-b)**2)/c)
                x2 = x.compress((bound1 <= x) & (x < bound2))
                x2 = x.compress((bound1 <= x) & (x < bound2))
                y2 = y.compress((bound1 <= x) & (x < bound2))
                popt, pcov = optimize.curve_fit(f, x2, y2, (g0,g1,(guess_max_x-guess_min_x)/2), sigma=sqrt(y2))

                # Find the x value of the Peak:
                peak = popt[1]
                if i == 0:
                    MIP1 = peak
                elif i == 1:
                    MIP2 = peak
                elif i == 2:
                    MIP3 = peak
                elif i == 3:
                    MIP4 = peak
                print 'The pulse integral MPV of detector',i + 1, 'lies at', peak, 'ADC'

            # Safety nets: prevent the code from terminating if an error occurs
            except IndexError:
                if i == 0:
                    MIP1 = nan
                if i == 1:
                    MIP2 = nan
                if i == 2:
                    MIP3 = nan
                if i == 3:
                    MIP4 = nan
                print 'Apparently there is no data'
            except RuntimeError:
                if i == 0:
                    MIP1 = nan
                if i == 1:
                    MIP2 = nan
                if i == 2:
                    MIP3 = nan
                if i == 3:
                    MIP4 = nan
                print 'Unable to make a good fit'
            except TypeError:
                if i == 0:
                    MIP1 = nan
                if i == 1:
                    MIP2 = nan
                if i == 2:
                    MIP3 = nan
                if i == 3:
                    MIP4 = nan
                print 'Not enough data to make a proper fit' # When it uses less than three bins to make the fit you get this error
            except ValueError:
                if i == 0:
                    MIP1 = nan
                if i == 1:
                    MIP2 = nan
                if i == 2:
                    MIP3 = nan
                if i == 3:
                    MIP4 = nan
                print 'Value Error Occured'

        mpvint = (MIP1, MIP2, MIP3, MIP4)
        return mpvint


    def getph(date, ref = 0):
        """ Takes the pulseheight data from the datafile and averages for one hour.

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
            ph = events.readWhere('(timestamp >= t0) & (timestamp < t1)')['pulseheights']

            if ph[0,0] == -1:
                ph1 = []
            else:
                ph1 = plt.hist(ph[:,0], bins = 45, range = [88,877], histtype = 'step', log = 'True', color = 'k')
                ph1 = ph1[0] / float(hours)

            if ph[0,1] == -1:
                ph2 = []
            else:
                ph2 = plt.hist(ph[:,1], bins = 45, range = [88,877], histtype = 'step', log = 'True', color = 'r')
                ph2 = ph2[0] / float(hours)

            if ph[0,2] == -1:
                ph3 = []
            else:
                ph3 = plt.hist(ph[:,2], bins = 45, range = [88,877], histtype = 'step', log = 'True', color = 'g')
                ph3 = ph3[0] / float(hours)

            if ph[0,3] == -1:
                ph4 = []
            else:
                ph4 = plt.hist(ph[:,3], bins = 45, range = [88,877], histtype = 'step', log = 'True', color = 'c')
                ph4 = ph4[0] / float(hours)

            pulseheights = [ph1, ph2, ph3, ph4]
        except:
            print 'This date contains no data; \n', (date)
            pulseheights = [-1, -1, -1, -1]

        return pulseheights

    def phdiff(date):
        """Checks if the Pulse Height histogram differs in shape from the one on REFDATE

        """

        ph_check = getph(date)

        errors = []

        try:
            if len(ph_check) > 0:
                if min(ph_check[0]) >= 0:
                    ph1 = ph_check[0]
                    diff= []
                    for i in range(45):
                        diff.append(100*(abs(ph1[i] - ph1_ref[i])/ ph1_ref[i]))
                    diff1_av = average(diff)
                else:
                    if ph1_ref == []:
                        diff1_av = 0
                    else:
                        diff1_av = -1

                if min(ph_check[1]) >= 0:
                    ph2 = ph_check[1]
                    diff= []
                    for i in range(45):
                        diff.append(100*(abs(ph2[i] - ph2_ref[i])/ ph2_ref[i]))
                    diff2_av = average(diff)
                else:
                    if ph2_ref == []:
                        diff2_av = 0
                    else:
                        diff2_av = -1

                if min(ph_check[2]) >= 0:
                    ph3 = ph_check[2]
                    diff= []
                    for i in range(45):
                        diff.append(100*(abs(ph3[i] - ph3_ref[i])/ ph3_ref[i]))
                    diff3_av = average(diff)
                else:
                    if ph3_ref == []:
                        diff3_av = 0
                    else:
                        diff3_av = -1

                if min(ph_check[3]) >= 0:
                    ph4 = ph_check[3]
                    diff= []
                    for i in range(45):
                        diff.append(100*(abs(ph4[i] - ph4_ref[i])/ ph4_ref[i]))
                    diff4_av = average(diff)
                else:
                    if ph4_ref == []:
                        diff4_av = 0
                    else:
                        diff4_av = -1

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

    def getpi(date, ref = 0):
         """ Takes the pulse integral data from the datafile and averages for one hour.

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
            pi = events.readWhere('(timestamp >= t0) & (timestamp < t1)')['integrals']

            if pi[0,0] == -1:
                pi1 = []
            else:
                pi1 = plt.hist(pi[:,0], bins = 45, range = [877,8770], histtype = 'step', log = 'True', color = 'k')
                pi1 = pi1[0] / float(hours)

            if pi[0,1] == -1:
                pi2 = []
            else:
                pi2 = plt.hist(pi[:,1], bins = 45, range = [877,8770], histtype = 'step', log = 'True', color = 'r')
                pi2 = pi2[0] / float(hours)

            if pi[0,2] == -1:
                pi3 = []
            else:
                pi3 = plt.hist(pi[:,2], bins = 45, range = [877,8770], histtype = 'step', log = 'True', color = 'g')
                pi3 = pi3[0] / float(hours)

            if pi[0,3] == -1:
                pi4 = []
            else:
                pi4 = plt.hist(pi[:,3], bins = 45, range = [877,8770], histtype = 'step', log = 'True', color = 'c')
                pi4 = pi4[0] / float(hours)

            pulseintegrals = [pi1, pi2, pi3, pi4]
        except:
            print 'This date contains no data; \n', (date)
            pulseintegrals = [-1, -1, -1, -1]

        return pulseintegrals

    def pidiff(date):
        """Checks if the Pulse Integral histogram differs in shape from the one on REFDATE

        """

        pi_check = getpi(date)

        errors = []

        try:
            if len(pi_check) > 0:
                if min(pi_check[0]) >= 0:
                    pi1 = pi_check[0]
                    diff= []
                    for i in range(45):
                        diff.append(100*(abs(pi1[i] - pi1_ref[i])/ pi1_ref[i]))
                    diff1_av = average(diff)
                else:
                    if pi_ref == []:
                        diff1_av = 0
                    else:
                        diff1_av = -1

                if min(pi_check[1]) >= 0:
                    pi2 = pi_check[1]
                    diff= []
                    for i in range(45):
                        diff.append(100*(abs(pi2[i] - pi2_ref[i])/ pi2_ref[i]))
                    diff2_av = average(diff)
                else:
                    if pi2_ref == []:
                        diff2_av = 0
                    else:
                        diff2_av = -1

                if min(pi_check[2]) >= 0:
                    pi3 = pi_check[2]
                    diff= []
                    for i in range(45):
                        diff.append(100*(abs(pi3[i] - pi3_ref[i])/ pi3_ref[i]))
                    diff3_av = average(diff)
                else:
                    if pi3_ref == []:
                        diff3_av = 0
                    else:
                        diff3_av = -1

                if min(pi_check[3]) >= 0:
                    pi4 = pi_check[3]
                    diff= []
                    for i in range(45):
                        diff.append(100*(abs(pi4[i] - pi4_ref[i])/ pi4_ref[i]))
                    diff4_av = average(diff)
                else:
                    if pi4_ref == []:
                        diff4_av = 0
                    else:
                        diff4_av = -1

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
    # main program

    if __name__ == '__main__':
        data = openFile(DATAFILE, 'r')
        events = data.getNode('/s%s' %STATION, 'events')

        output = openFile(OUTPUTFILE, 'a')

        try:
            table = output.getNode('/s%s' %STATION, 'Logs')
        except NoSuchNodeError:
            group = output.createGroup('/', 's%s' %STATION)
            table = output.createTable(group, 'Logs', Logs)
        logs = table.row

        refpks = ref_peaks()
        ph_ref = getph(REFDATE, ref = 1)
        ph1_ref = ph_ref[0]
        ph2_ref = ph_ref[1]
        ph3_ref = ph_ref[2]
        ph4_ref = ph_ref[3]
        pi_ref = getpi(REFDATE, ref = 1)
        pi1_ref = pi_ref[0]
        pi2_ref = pi_ref[1]
        pi3_ref = pi_ref[2]
        pi4_ref = pi_ref[3]

        date1 = START

        while date1 < STOP:
            print(date1)
            ph, ph_fit  = get_ph(date1)
            pi, pi_fit  = get_pi(date1)
            logs['start'] = mktime(date1.timetuple())
            logs['nr_ev'] = compare_peaks(date1, refpks)
            logs['ph_diff'] = phdiff(date1)
            logs['pi_diff'] = pidiff(date1)
            logs['ph_mpv'] = find_MPV(ph, ph_fit)
            logs['pi_mpv'] = find_MPVINT(pi, pi_fit)
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
