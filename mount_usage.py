#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 18:32:46 2023

@author: edombek
"""

from mount import Mount
import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time

location = coord.EarthLocation(lat=20 * u.deg, lon=30 * u.deg, height=60.0 * u.m)

#m = Mount(device_name='Telescope Simulator') # indi_simulator_telescope
m = Mount(use_simulation=True) # indi_eqmod_telescope
m.connect()
m.set_Location_And_Time(location)

obs_time = Time.now()
target = coord.get_body('sun', obs_time, location, ephemeris='jpl')

m.goto(target, wait = True)
print(m.getAltAz())

target = coord.AltAz(alt = 30*u.deg, az = 20*u.deg)
m.goto(target, wait = True)
m.setTracking(False)