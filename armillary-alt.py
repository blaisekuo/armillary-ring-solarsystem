from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, AltAz
from astropy.coordinates import get_body_barycentric, get_body, get_moon

import math

# install astroplan

#conda install -c astropy astroplan

#to use Observer.sun_set_time

from astroplan import Observer
#from astroplan import download_IERS_A
#download_IERS_A()
from pytz import timezone
import datetime

#pip install timezonefinder
import timezonefinder

#imports for the plot

import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

# pick a date
dt = Time("2019-11-30")

#pick a location
#example
#location for Cebu City
#lat long for Bai hotel
locname = 'Cebu City'
lat = 10.3247
lng = 123.9345

def getalt(dt,locname,lat,lng):
    
    #get time zone
    tf = timezonefinder.TimezoneFinder()
    timezone_str = tf.certain_timezone_at(lat=lat,lng=lng)

    #set location for observer object
    loc = EarthLocation.from_geodetic(lng,lat)
    
    #setup observer object
    observer = Observer(name=locname,
               location=loc,
               timezone=timezone(timezone_str))
    
    #get the sunset time for the given day, given a date string. 
    # for our given location

    sunsettime=observer.sun_set_time(dt, which='next')
    
    #alternative, can also use nautical twilight as a criteria
    #get nautical twilight time for the given day, given a date string. 
    #nautical twilight is when the sun is 12 degrees below the horizon
    # for our given location

    #nautitime=observer.twilight_evening_nautical(dt, which='next')

    #convert time to datetime so that we can manipulate the timezone
    sunset_dt=sunsettime.to_datetime()
    
    #add utc timezone
    sunset_utc=sunset_dt.replace(tzinfo=timezone('utc'))
    
    #apply our local time zone
    localtime=sunset_utc.astimezone(timezone(timezone_str))
    
    #set time
    t=sunsettime
    
    #setup a dictionary for the visible planets, sun and moon
    bodies = {
        "sun" : '',
        "moon" : '',
        "mercury" : '',
        "venus" : '',
        "mars" : '',
        "jupiter" : '',
        "saturn" : ''
    }
    
    #set up a dictionary for the symbols for the bodies
    bodies_names = {
        "sun" : '\u2609',
        #"sun" : '\u1f31e',
        "moon" : '\u263d',
        #"moon" : '\u1f31b',
        "mercury" : '\u263f',
        "venus" : '\u2640',
        "mars" : '\u2642',
        "jupiter" : '\u2643',
        "saturn" : '\u2644'
    }
    
    #get the alt az position of the bodies for our given time and location
    #we will use the alt number

    with solar_system_ephemeris.set('builtin'):
        for key in bodies:
            bodies[key] = get_body(key,t,loc).transform_to(AltAz(obstime=t,location=loc))
            
    fig = plt.figure(figsize=(10,10))


    ax = fig.add_subplot(111, polar=True)


    ax.set_theta_direction(-1)
    ax.set_thetamin(180)
    ax.set_thetamax(360)
    ax.set_xticks(np.pi/180. * np.linspace(180,  360, 5))
    ax.set_xticklabels(['S', 'SW', 'W', 'NW',  'N'])
    ax.set_rorigin(-90)

    ax.set_title("Naked Eye planets, Sun and Moon at sunset on " + str(localtime) + " in " + locname)

    for key, values in bodies.items():
        az_rad=values.az.wrap_at(180 * u.deg).radian
        if key == 'sun':
            plotplanet=ax.scatter(az_rad, values.alt, marker='o'.format(bodies_names[key]), alpha=0.3, s=700, color='orange')
            plotplanet.set_label(key)
        else:
            plotplanet=ax.scatter(az_rad, values.alt, marker=r"$ {} $".format(bodies_names[key]), s=250)
            altstr=(values.alt).to_string(precision=0)
            plotplanet.set_label(key + " " + altstr)
        
    r = np.zeros(360)
    theta = np.arange(0,2 * np.pi,2*np.pi/360)

    ax.plot(theta,r, color='orange')

    ax.legend(loc='right',labelspacing=1.2,borderpad=1,bbox_to_anchor=(1.4, .55))

    plt.show()