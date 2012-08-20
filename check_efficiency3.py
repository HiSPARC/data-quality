"""check_efficiency determines the ration of time a station is emitting useful data.

uses logfile.h5 files as input
in globals, set the station number and the period to be analyzed
returns efficiency as an integer percentage of time.
"""

# import modules
from tables import openFile, NoSuchNodeError
import csv

# define globals
DATAFILE = 'Logfile.h5'
OUTPUTFILE = 'Efficiency_Sciencepark_2012_2.dat'
STATION = [501, 502, 503, 504, 505, 506, 509]

# define routines

def efficiency(station):
    """determines the ratio of the number of useful periods to the total number of periods

        returns this ratio as an integer percentage
    """
#    station = STATION

    try:
        table = data.getNode('/s%s' %station, 'Logs')
        useful = table.col('useful')
        nr_ev = table.col('nr_ev')
        ph_diff = table.col('ph_diff')
        pi_diff = table.col('pi_diff')
        ph_mpv = table.col('ph_mpv')
        pi_mpv = table.col('pi_mpv')
        
        up = 0
        count_ev = 0
        count_phdiff = 0
        count_pidiff = 0
        count_phmpv = 0
        count_pimpv = 0
        
        for i in range(len(nr_ev)):
        
            if sum(nr_ev[i]) != 0:
                up = up + 1
                       
                if (min(nr_ev[i]) >= 75) and (max(nr_ev[i]) <= 125):
                    count_ev += 1
                
                if max(ph_diff[i]) < 50:
                    count_phdiff += 1
                
                if max(pi_diff[i]) < 50:
                    count_pidiff += 1
                
                if (min(ph_mpv[i]) > 120) and (max(ph_mpv[i]) < 439):
                        count_phmpv += 1
                
                if (min(pi_mpv[i]) > 1200) and (max(pi_mpv[i]) < 4390):
                        count_pimpv += 1
                                     
        uptime = 100*up / len(nr_ev)
        efficiency = 100*useful.sum()/up
        nr_ev_check = 100*count_ev / up
        ph_diff_check = 100*count_phdiff / up
        pi_diff_check = 100*count_pidiff / up
        ph_mpv_check = 100*count_phmpv / up
        pi_mpv_check = 100*count_pimpv / up
    except:
        NoSuchNodeError
        efficiency = 0
        print 'there are no data for station', station

    return uptime, efficiency, nr_ev_check, ph_diff_check, pi_diff_check, ph_mpv_check, pi_mpv_check

# main program

if __name__ == '__main__':

    efflist = []
    with openFile(DATAFILE, 'r') as data:

        for i in range(len(STATION)):
            st = STATION[i]
            uptime, eff, nr_ev_check, ph_diff_check, pi_diff_check, ph_mpv_check, pi_mpv_check = efficiency(st)
            efflist.append([st,uptime, eff, nr_ev_check, ph_diff_check, pi_diff_check, ph_mpv_check, pi_mpv_check])
            
            print 'the efficiency of station ', st, ' is ', eff, ('%')
    with open(OUTPUTFILE, 'w') as outputfile:
            file = csv.writer(outputfile, delimiter = '\t')
            outputfile.write ('st\tup\teff\tev\tphd\tpid\tphp\tpip\n')
            file.writerows(efflist)