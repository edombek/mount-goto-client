#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 16:42:14 2023

@author: edombek
"""

import time
import sys
import datetime

import PyIndi
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u

class IndiClient(PyIndi.BaseClient):
    def __init__(self):
        super(IndiClient, self).__init__()
        self.isconnected = False

    def newMessage(self, d, m):
        print("Message for " + d.getDeviceName() + ":" + d.messageQueue(m))

    def serverConnected(self):
        self.isconnected = True
        print(
            "Server connected (" + self.getHost() + ":" + str(self.getPort()) + ")"
        )

    def serverDisconnected(self, code):
        self.isconnected = False
        print(
            "Server disconnected (exit code = "
            + str(code)
            + ","
            + str(self.getHost())
            + ":"
            + str(self.getPort())
            + ")"
        )

class Mount():
    def __init__(self, host="localhost", port=7624, device_name="EQMod Mount",
                 use_simulation=False, simprop='SIMULATION'):
        self.indiclient = IndiClient()
        self.indiclient.setServer(host, port)
        self.indiclient.watchDevice(device_name)
        self.device_name = device_name
        self.device = None
        self.use_simulation = use_simulation
        self.simprop = simprop
        
        self.numberFailures = set()
        self.switchFailures = set()
        self.textFailures = set()
        
    def getNumberWithRetry(self, prop, ntry=5, delay=0.2, as_prop=True):
        while ntry > 0:
            if as_prop:
                pr = self.device.getProperty(prop)
                p = pr.getNumber()
            else:
                p = self.device.getNumber(prop)
            if type(p) == PyIndi.PropertyViewNumber: #p.isValid():
                self.numberFailures.discard(prop)
                return p
            if prop in self.numberFailures:
                return None
            print("Unable to get number property " + prop + ", retrying")
            ntry -= 1
            time.sleep(delay)
        print("Unable to get number property " + prop + ", marking as failed")
        self.numberFailures.add(prop)
        return None


    def getSwitchWithRetry(self, prop, ntry=5, delay=0.2, as_prop=True):
        while ntry > 0:
            if as_prop:
                pr = self.device.getProperty(prop)
                p = pr.getSwitch()
            else:
                p = self.device.getSwitch(prop)
            if type(p) == PyIndi.PropertyViewSwitch: #p.isValid():
                self.switchFailures.discard(prop)
                return p
            if prop in self.switchFailures:
                return None
            print("Unable to get switch property " + prop + ", retrying")
            ntry -= 1
            time.sleep(delay)
        print("Unable to get switch property " + prop + ", marking as failed")
        self.switchFailures.add(prop)
        return None


    def getTextWithRetry(self, prop, ntry=5, delay=0.2, as_prop=True):
        while ntry > 0:
            if as_prop:
                pr = self.device.getProperty(prop)
                p = pr.getText()
            else:
                p = self.device.getText(prop)
            if type(p) == PyIndi.PropertyViewText: #p.isValid():
                self.textFailures.discard(prop)
                return p
            if prop in self.textFailures:
                return None
            print("Unable to get text property " + prop + ", retrying")
            ntry -= 1
            time.sleep(delay)
        print("Unable to get text property " + prop + ", marking as failed")
        self.textFailures.add(prop)
        return None
            
    def connect(self):
        print(f'Connecting server {self.indiclient.getHost()}:{self.indiclient.getPort()}')
        serverconnected = self.indiclient.connectServer()
        while not (serverconnected):
            print('No indiserver running')
            time.sleep(2)
            serverconnected = self.indiclient.connectServer()
        print('Connected to indiserver')
        
        device = self.indiclient.getDevice(self.device_name)
        while not (device):
            print(f'Trying to get device {self.device_name}')
            time.sleep(0.5)
            device = self.indiclient.getDevice(self.device_name)
        self.device = device
        print(f'Got device {self.device_name}')
        
        if not (device.isConnected()):
            if self.use_simulation:
                print("setting " + self.simprop + " On")
                device_sim = self.getSwitchWithRetry(self.simprop)
                if device_sim:
                    device_sim[0].setState(PyIndi.ISS_ON)  # the "ENABLE" switch
                    device_sim[1].setState(PyIndi.ISS_OFF)  # the "DISABLE" switch
                    self.indiclient.sendNewProperty(device_sim)
                else:
                    print('Simulation not work, exiting')
                    sys.exit(1)

        if not (device.isConnected()):
            device_connect = device.getSwitch("CONNECTION")
            while not (device_connect):
                print("Trying to connect device " + self.device_name)
                time.sleep(0.5)
                device_connect = device.getSwitch("CONNECTION")
        if not (device.isConnected()):
            device_connect[0].setState(PyIndi.ISS_ON)  # the "CONNECT" switch
            device_connect[1].setState(PyIndi.ISS_OFF)  # the "DISCONNECT" switch
            self.indiclient.sendNewProperty(device_connect)
        while not (device.isConnected()):
            time.sleep(0.2)
        print("Device " + self.device_name + " connected")
    
    def set_Location_And_Time(self, location: coord.earth.EarthLocation, utc_time=None):
        self.location = location
        if not utc_time: utc_time = datetime.datetime.now(datetime.timezone.utc)
        p = self.getTextWithRetry("TIME_UTC")
        p[0].setText(utc_time.isoformat())
        p[1].setText(str(int(utc_time.astimezone().utcoffset().total_seconds()//2600)))
        self.indiclient.sendNewProperty(p)
        
        p = self.getNumberWithRetry("GEOGRAPHIC_COORD")
        p[0].setValue(location.lat.to(u.deg).value)
        p[1].setValue(location.lon.to(u.deg).value)
        try: p[2].setValue(location.height.to(u.m).value)
        except: pass
        self.indiclient.sendNewProperty(p)
    
    def setTracking(self, guide):
        p = self.getSwitchWithRetry("TELESCOPE_TRACK_STATE")
        if guide:
            p[0].setState(PyIndi.ISS_ON)
            p[1].setState(PyIndi.ISS_OFF)
        else:
            p[0].setState(PyIndi.ISS_OFF)
            p[1].setState(PyIndi.ISS_ON)
        self.indiclient.sendNewProperty(p)
    
    def goto(self, target, Jnow = True, wait = False):
        obs_time = Time.now()
        if type(target) == coord.AltAz:
            target = coord.AltAz(alt=target.alt, az=target.az,
                                 obstime=Time.now() if not target.obstime else target.obstime,
                                 location=self.location if not target.location else target.location)
            target = target.transform_to(self.location.get_gcrs(obs_time))
        if Jnow:
            target = coord.FK5(ra = target.ra, dec = target.dec).transform_to(coord.FK5(equinox=obs_time))
            par = 'EQUATORIAL_EOD_COORD'
        else:
            par = 'EQUATORIAL_COORD'
        
        p = self.getNumberWithRetry(par)
        p[0].setValue(target.ra.to(u.hourangle).value)
        p[1].setValue(target.dec.to(u.deg).value)
        pcs = self.getSwitchWithRetry("ON_COORD_SET")
        pcs[0].setState(PyIndi.ISS_ON)
        pcs[1].setState(PyIndi.ISS_OFF)
        pcs[2].setState(PyIndi.ISS_OFF)
        self.indiclient.sendNewProperty(pcs)
        self.indiclient.sendNewProperty(p)
        if wait:
            print('Slewing')
            while p.getState() == PyIndi.IPS_BUSY:
                time.sleep(0.5)
                print('.', end='')
    
    def getAltAz(self):
        p = self.getNumberWithRetry('HORIZONTAL_COORD')
        if 'HORIZONTAL_COORD' in self.numberFailures:
            t = self.getJnow()
            return t.transform_to(coord.AltAz(obstime=Time.now(), location=self.location))
        return coord.AltAz(alt=p[0].getValue()*u.deg, az=p[1].getValue()*u.deg)
    
    def getJnow(self):
        p = self.getNumberWithRetry('EQUATORIAL_EOD_COORD')
        return coord.FK5(ra=p[0].getValue()*u.hourangle, dec=p[1].getValue()*u.deg, equinox=Time.now())