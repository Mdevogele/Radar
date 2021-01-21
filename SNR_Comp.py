#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 14:33:15 2020

@author: mdevogel
"""

import numpy as np


Latitude = {'AO' : 18.34350}

System = {'AO leg': {'Tsys'     : 23,
                     'Power'    : 800,
                     'Gain'     : 9.011,
                     'T2R'      : 8,
                     'ZA'       : 19,
                     'Tau_min'  : 2,
                     'Freq'     : 2.380e9,
                     'Diameter' : 304.8,
                     'Type'     :'Single'},
          
          'NGAO': {  'Tsys'       : 25,
                     'Power'      : 5,
                     'N Dishes'   : 1000,
                     'D Dishes'   : 9,
                     'Eff Dishes' : 0.75,
                     'Gain'       : False,
                     'T2R'        : 1,
                     'ZA'         : 48,
                     'Tau_min'    : 0,
                     'Freq'       : 5e9,
                     'Type'       : 'Array'},
          
# Tau min is the Minimum integration time for Arecibo, even if RTT is
# less than $switch_A, unless RTT is less than 4 s
# At Arecibo, we usually allowed for at least two seconds of recording data
#   Recording for less than 2 s didn't seem to work          
          
          
          
          'Goldstone': {'Tsys': False,
                     'Power' : False,
                     'Gain' : False,
                     'T2R' : 2,
                     'ZA' : 80,
                     'Lat' : False}}






def SNR_Comp(RTT,diam,PERIOD = False, Site = 'AO', Sys = 'AO leg', ZA = False, a=False,Dec=False,l=False,Power=False,t=False,D=False,m=False,b=False,s=False,B=False):

    # Based on pearl program from Mike Nolan and Sean Marshall
    
    # Get the latitude from user input or list of known sites    
    if type(Site) == float:
        Lat = Site
    else:
        Lat = Latitude[Site]
    
    
    # Get instrument configuration from user dictionnary or list of known instruments
    if type(Sys) == dict:
        Instrument = Sys
    else:
        Instrument = System[Sys]
        
    # Allows the user to override the default ZA
    if ZA:
        Instrument['ZA'] = ZA
        
    # If no declination is given for an object, we are assuming that the object 
    # possesses the same declination as the latitude of the considered site
    if not Dec:
        Dec = Lat
        
    # If no period is provided, assume a default rotation period based on the size of the object
        
    if not PERIOD:
        Diam2Per = {
            5000:4,
            1000:3,
            500:2.5,
            140:2.1,
            100:1.5,
            80:1,
            50:0.5,
            40:0.25,
            30:0.15,
            10:0.1}
        
        period=Diam2Per[diam]*3600
    else:
        period = PERIOD*3600
     
    # allows the user to override the default transmitter power for the selectred instrument
             
    if Power:
        Instrument['Power'] = Power
        
    
    radius = diam * 0.5;    


        
    # Even if the asteroid's properties are known perfectly, these SNR calculations
    # are only accurate to about 20%, due to uncertainties in the telescope calibrations
    
    PI  = np.pi
    deg  = PI/180.0
    au   = 149597870700.0
    Clight  = 299792458.0
    kB     = 1.380649e-23
    
    # From 2019 SI redefinition (draft of February 6, 2019)
    #   https://www.bipm.org/en/measurement-units/rev-si/
    
    hour   = 3600.0
    #use constant SPD    => 86400.0;
    sdd    = 2.0*PI/7.292115e-5 # Sidereal day, in seconds
    Jy     = 1.0e-26; # Jy = 10^-26 W/(m^2*Hz)
    KJy    = 2.0e26*kB;
    # 1 K/Jy is equivalent to an effective area of 2*k_B*(1 K)/(1 Jy) ~ 2761 m^2
    #   http://www.cv.nrao.edu/course/astr534/Equations.html
    
    # SNR = G_tx*P_tx*A_rx*sigma*sqrt(tau/B)/[k_B*T_rx*(4*pi*Delta^2)^2]
    # 4*pi*A_eff = G*lambda^2; A_eff = G*lambda^2/(4*pi); G = 4*pi*A_eff/lambda^2
    # If monostatic:
    # SNR = G*P_tx*A*sigma*sqrt(tau/B)/[k_B*T_sys*(4*pi*Delta^2)^2]
    #     = G^2*P_tx*lambda^2*sigma*sqrt(tau/B)/[4*pi*k_B*T_sys*(4*pi*Delta^2)^2]
    #     = 4*pi*P_tx*A^2*sigma*sqrt(tau/B)/[k_B*T_sys*lambda^2*(4*pi*Delta^2)^2]
    # G^2 = SNR*[4*pi*k_B*T_sys*(4*pi*Delta^2)^2]/[P_tx*lambda^2*sigma*sqrt(tau/B)]
    # A^2 = SNR*[k_B*T_sys*lambda^2*(4*pi*Delta^2)^2]/[4*pi*P_tx*sigma*sqrt(tau/B)]    
    
    
    
    #   -a alb: Use specified radar albedo
    #   -d dec: Use specified declination (degrees) instead of full track
    #   -l lat: Specifies sub-radar latitude (in degrees; defaults to zero otherwise)
    #   -p power: Use specified power (in kW) for all transmitters
    #   -t time: Use time (s) as the max coherent integration time
    #   -D: Distance (in au) given instead of round-trip time
    #   -m: Size is given as H instead of D
    #   -b: Just print bandwidth
    #   -s: Just print S/N per run at AO
    # ";
    # # -B: Don't subtract Tx-Rx switch time
    # # -P period to specify period in hours (defaults to 2.1 hours otherwise)
        
        
    ZA_max = Instrument['ZA']
        
        
    if B:  # B for bistatic
    	#$RTT_min_ = 0.;
        switch  = 0.
        tau_min = 0.
    else:
        switch  = Instrument['T2R']
        tau_min = Instrument['Tau_min']
    
    if D:
    	dist = RTT * au; # Distance to target (in meters)
    	RTT = 2.0 * dist / Clight
    else:
        RTT = RTT # Round-trip time (seconds)
        dist = RTT * 0.5 * Clight
    
    
    if m:  # If given absolute magnitude, instead of diameter
        x = diam*(0.984) - 17.4
    	#print "  H magnitude $diam corresponds to diameter of ";
        print("H magnitude {:.3f} corresponds to diameter of {%.1f} m\n".format(diam,10 ** ((diam*(0.984) - 17.4)) / -5.) * 1000)
        diam = 10 ** (x / -5.) * 1000
    
    if a: 
    	alb = a
    else:
        alb = 0.10 # Radar albedo of asteroid; 10% by default
        
    # Assuming same radar albedo at all frequencies, since we usually don't know it anyway
    
    # Names of variables should be such that the dish and/or frequency band is always unambiguous
    # Notation for transmitters:

    # Define frequencies (all in Hz)

    freq  = Instrument['Freq']
    Lambda  = Clight / freq

    
    if Instrument['Type'] == 'Single':
        Ae = Instrument['Gain']*KJy 
        # Diameters of the dish (in meters)
        DA   = Instrument['Diameter'] 
    if Instrument['Type'] == 'Array':
        N_Dishes = Instrument['N Dishes']
        D_Dishes = Instrument['D Dishes']
        ap_eff = Instrument['Eff Dishes']
        DA = np.sqrt((D_Dishes/2)**2*N_Dishes**2)
        Ae = 0.25 * np.pi * DA**2 * ap_eff

    Lat = Latitude[Site]
        
    if (abs(d - Lat) > AO_ZA_max):
        quit("ERROR: Declination is beyond Arecibo's limits\n")


    Lat = Latitude[Site]
    # From http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
    #   sin(altitude) = cos(ZA) = sin(delta)*sin(latitude) + cos(delta)*cos(latitude)*cos(H)
    #   cos(delta)*cos(lat)*cos(H) = cos(ZA) - sin(delta)*sin(lat)
    #   cos(H) = [cos(ZA) - sin(delta)*sin(lat)]/[cos(delta)*cos(lat)]
    #   Range of hour angles: 2*arccos[[cos(ZA) - sin(delta)*sin(lat)]/[cos(delta)*cos(lat)]]
    HArange = np.arccos((np.cos(AO_ZA_max*deg) - np.sin(d*deg)*np.sin(Lat*deg)) / (np.cos(d*deg)*np.cos(Lat*deg)))
    #   Half-range of hour angles, in radians
    ttime = (HArange/PI)*sdd
   	# Note that this does NOT account for the minimum zenith angle of 1.06 deg (keyhole)
   	#   http://www.naic.edu/science/generalinfo_set.htm
    
    print("Tracking time for Arecibo: {:.1f} minutes\n".format(ttime/60.0))
    

    Tsys = Instrument['Tsys'] # S-narrow, http://www.naic.edu/~astro/RXstatus/rcvrtabz.shtml
    
    
    if p:
    	#print "Setting transmitter power to $opt_p kW\n";
    	Pt = p*1000.0 # Convert from kW to W
    else:
        Pt = Instrument['Power']
    
    xsec = alb * PI * radius**2
    
    # P_rx = [G_tx*P_tx/(4*pi*Delta^2)]*[alb*pi*radius^2]*A_rx/(4*pi*Delta^2)
    #      = [G_tx*P_tx*A_rx*alb*pi*radius^2]/(4*pi*Delta^2)^2
    #      = [4*pi*A_tx*P_tx*A_rx*alb*pi*(radius/lambda)^2]/(4*pi*Delta^2)^2
    #      = [4*pi^2*A_tx*P_tx*A_rx*alb*(radius/lambda)^2]/(4*pi*Delta^2)^2
    #      = [P_tx*A_tx*A_rx*alb*(radius/lambda)^2]/(4*Delta^4)
    #      = (P_tx*A_tx*A_rx*alb*radius^2)/(4*lambda^2*Delta^4)
    # 4*pi*A_e = G*lambda^2; G_tx = 4*pi*A_tx/lambda^2

    Pr = Pt* Ae**2 * alb*radius**2 / (4 * Lambda**2  * dist**4); # Arecibo S-band monostatic (pre-Maria)
    
    if l:
    	bwfac = np.cos(l * deg)
    else: 
        bwfac = 1.0
    
    #$bw = 2     * (2 * PI * $radius) / ($per * 3600) * 2   / Clight * $freqAS * $bwfac;
    #      2way     circum                  time       FW
    # c = lambda*f; lambda = c/f
    # Bandwidth: B = 4*pi*D*cos(sub-radar latitude)/(lambda*P) = 4*pi*D*f*cos(srl)/(c*P)

    bwcoef = (4 * np.pi * diam * bwfac) / (Clight * period)
    bw  = bwcoef * freq

    
    if b:
    	print("{%3.1f}".format(bw))
            
    # N_rms = k*T*B/sqrt(B*tau) = k*T*sqrt(B/tau)
    # Raw noise power (not divided by sqrt(B*tau)), in W
    Pn  = kB * Tsys  * bw;
    
    if RTT >= (switch + tau_min):
        # If RTT >= 10.0 seconds, integration time would be at least 2.0 seconds
        tau = RTT - switch
        run = 2.0*RTT # Total time needed for one transmit-receive cycle (run)
    else:
        if RTT >= 2.0*tau_min:
            # If RTT is between 4.0 (2*$tau_min_A) and 10.0 ($switch_A + $tau_min_A)
            #   seconds, set integration time to 2.0 seconds.
            tau = 1.0*tau_min
        else:
            # If RTT is less than 4.0 seconds, set integration time to half the RTT.
            # This is arbitary, but I want integration time to increase
            #   monotonically (or at least stay constant) as RTT increases.
            tau = 0.5*RTT
        print("WARNING: Resetting Arecibo integration time to {:4.2f} seconds\n".format(tau_A))
        run = RTT + tau + switch

    
    if tau > 0.0:
    	SNRrun = Pr / (Pn / np.sqrt(bw*tau))
    else:
        # This should not happen, based on how I am now (2019-10-11) handling $tau_A
    	SNRrun = 0
    
    
    print("RTT {:.1f} s\n".format(RTT))
    if RTT < (switch + tau_min):
        print("    Arecibo: integration time {:.1f} s, cycle time {:.1f} s\n".format(tau, run))
        
    print("  For radar albedo of {:.2f}, x-section = {:d} m^2 = {:.3f} km^2\n".format(alb, int(xsec), xsec * 1.e-6))
    
    if s:
    	print("%3.1f".format(SNRrun))    
    
    
    
    if RTT < ttime:
        runs = ttime / run
    else:
        runs = 0.
        # Can't do monostatic observations if RTT is too long
    
    
    SNR_Total = SNRrun * np.sqrt(runs)
        
    
    return SNR_Total
    


if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description='Extract spectrum')

    parser.add_argument('-a', help='',
                        default = 'alb: Use specified radar albedo')
    parser.add_argument('-d',
                        help='dec: Use specified declination (degrees) instead of full track',
                        default = "")
    parser.add_argument('-l',
                        help='lat: Specifies sub-radar latitude (in degrees; defaults to zero otherwise)',
                        action="store_true")    
    parser.add_argument('-p',
                        help='power: Use specified power (in kW) for all transmitters',
                        default = 'SpecOut.spec')    
    parser.add_argument('-t', help='time: Use time (s) as the max coherent integration time',
                        action="store_true")
    parser.add_argument('-D',
                        help='Distance (in au) given instead of round-trip time',
                        action="store_true")
    parser.add_argument('-m',
                        help='Size is given as H instead of D',
                        action="store_true")    
    parser.add_argument('-b',
                        help='Just print bandwidth',
                        action="store_true")    

    parser.add_argument('-s',
                        help='Just print S/N per run at AO',
                        action="store_true")    

    parser.add_argument('-B',
                        help="Don't subtract Tx-Rx switch time",
                        action="store_true")  



    parser.add_argument('RTT', help='Round trip time in (s)',
                        nargs='+')

    parser.add_argument('DIAMETER', help='Diamter in (m)',
                        nargs='+')
    
    parser.add_argument('PERIOD', help='Period in (h)',
                        nargs='+')    
    

    args = parser.parse_args()
    
    a = args.a
    d = args.d
    l = args.l
    p = args.p
    t = args.t
    D = args.D
    m = args.m
    b = args.b
    s = args.s
    B = args.B
    

    RTT = args.RTT
    DIAMETER = args.DIAMETER
    PERIOD = args.PERIOD

    
    SNR_Comp(RTT,DIAMETER,PERIOD,a=a,d=d,l=l,p=p,t=t,D=D,m=m,b=b,s=s,B=B)
    pass



