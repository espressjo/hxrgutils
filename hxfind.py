#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 15:29:18 2022

Utility module to find HxRG-ENG file on the computer.

@author: Jonathan St-Antoine
"""

from datetime import datetime,timedelta
from os import listdir
from natsort import natsorted
from os.path import join
from astropy.io.fits import getheader


class hxfile():
    """
    findlast
    ========
    Module to find HxRG-ENG file on a computer. 
    
    Main methods
    ------------
    *find   find all the files with uid #, return a list of files
    
    Search Path
    -----------
    Default search path is /opt/HxRG-ENG/data. You can
    change the default path when you initiallize the class. I.e.,
    
    .. code-block:: python
       :emphasize-lines: 1
       :linenos:
           
        from engfind import findlast 
        
        fl = findlast(p="/home/data")
        ls = fl.find(uid=20220504123602)
        print(ls)
    
    Usage
    -----
    .. code-block:: python
       :emphasize-lines: 1
       :linenos:
           
        from engfind import findlast 
        
        fl = findlast() #initialize the class. 
        
        ls = fl.find("20220504123602") # you can give the uid with an int or str
        
        print(ls) #ls will contain the file in ascending order i.e., ls[0]="../H4RG_R01_R01.fits"
       
    
    """
    def __init__(self,p = '/opt/spip/data'):
        self.path = p
    def find_dir(self,uid):
        """
        Unique identification number is given for each
        ramps taken with HxRG-ENG. Given an uid this 
        function returns a directory that contains some
        files with the specified UID in the header

        Parameters
        ----------
        uid : STR,INT
            Unique identification number.

        Returns
        -------
        directory : STR
            return an array of files with fullpath.

        """
        allfiles = [];
        uid = str(uid)
        start,stop = self._init_dt(uid)#find uid 1 hour before and 1 hour after given uid
        ls = [f for f in listdir(self.path) if f.isnumeric()]
        ls = [d for d in ls if all([start<int(d),stop>int(d),len(d)==len(str(start))])]
        for d in ls:
            for f in listdir(join(self.path,d)):
                if '.fits' not in f:
                    continue 
                h = getheader(join(join(self.path,d),f))
                if 'SEQID' not in h:
                    continue
                if uid in h['SEQID']:
                    return join(self.path,d)
        return ''
       
    def find(self,uid):
        """
        Unique identification number is given for each
        ramps taken with HxRG-ENG. Given an uid this 
        function returns a list of all the file with
        uid in the header.

        Parameters
        ----------
        uid : STR,INT
            Unique identification number.

        Returns
        -------
        allfiles : array[str]
            return an array of files with fullpath.

        """
        allfiles = [];
        if not isinstance(uid,str):
            uid = str(uid)
        start,stop = self._init_dt(uid)#find uid 1 hour before and 1 hour after given uid
        ls = [f for f in listdir(self.path) if f.isnumeric()]
        ls = [d for d in ls if all([start<int(d),stop>int(d),len(d)==len(str(start))])]
        for d in ls:
            for f in listdir(join(self.path,d)):
                if '.fits' not in f:
                    continue 
                h = getheader(join(join(self.path,d),f))
                if 'SEQID' not in h:
                    continue
                if uid in h['SEQID']:
                    allfiles.append(join(join(self.path,d),f))
        if len(allfiles)!=0:
            allfiles = natsorted(allfiles)
        return allfiles
    
    def _init_dt(self,uid):
        '''
        Return a folder name 1hour before the UID and 1 hour after UID

        Parameters
        ----------
        uid : STR
            Unique identidication number.

        Returns
        -------
        INT,INT
            Return yyymmddhhmmss(- 1hour),yyymmddhhmmss(+ 1hour) .

        '''
        print(uid)
        if not isinstance(uid,str):
            uid=str(uid)
        y = int(uid[:4])
        m = int(uid[4:6])
        d = int(uid[6:8])
        H = int(uid[8:10])
        M = int(uid[10:12])
        S = int(uid[12:14])
        dt = datetime(y,m,d,H,M,S)
        delta = timedelta(hours=1)
        dt1 = dt-delta
        dt2 = dt+delta
        start = dt1.strftime("%Y%m%d%H%M%S")
        stop = dt2.strftime("%Y%m%d%H%M%S")
        return int(start),int(stop)
if '__main__' in __name__:
    fl = hxfile()
    print(fl.find('20220421152140'))
