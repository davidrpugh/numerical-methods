from __future__ import division
import zipfile 
from urllib import urlopen 
from StringIO import StringIO 

import pandas as pd

def get_pwt_data_old(version=71, date='11302012', extract=False):
    """
    Load older (i.e., < 8.0) versions of the Penn World Tables data. If no local
    copy of the data is found, then a copy will be downloaded from the web.

    Arguments:
 
        version: (int) Version number for PWT data. Default is 71 (which is the 
                 most recent version).
        date:    (str) Date of the version update of the PWT data. Format is 
                 'mmddyyyy'. Default is '11302012' (which is the most recent 
                 update of version 71 of PWT).
        extract: (boolean) Whether or not you wish to save local copy of the 
                 PWT data in your working directory. Default is False.

    Returns:

        pwt:    A Pandas Panel object containing the Penn World Tables 
                Data along with the Solow residuals.

    TODO: Convert time index to datetime object.

    """        
    # first check for a local copy of PWT
    try: 
        path = 'pwt' + str(version) + '_wo_country_names_wo_g_vars.csv'
        pwt = pd.read_csv(path, index_col=['year', 'isocode'])

    # otherwise, download the appropriate zip file 
    except IOError:  
        url = ('http://pwt.econ.upenn.edu/Downloads/pwt' + str(version) + 
               '/pwt' + str(version) +'_' + date + 'version.zip') 
        archive = zipfile.ZipFile(StringIO(urlopen(url).read()), 'r') 

        # to extract or not to extract...
        tmp_file = 'pwt' + str(version) + '_wo_country_names_wo_g_vars.csv'
        
        if extract == True:
            archive.extractall()
            pwt = pd.read_csv(tmp_file, index_col=['year', 'isocode'])
        else:
            pwt = archive.read(tmp_file)
            pwt = pd.read_csv(StringIO(pwt), index_col=['year', 'isocode'])         
    
    # convert to Pandas Panel object
    pwt = pwt.to_panel()
    
    return pwt

def get_pwt_data(version=80):
    """
    Downloads the Penn World Tables database, including all of the necessary
    STATA programs and DO files for replicating it, into your working directory.

    Arguments:
 
        version:  (int) Version number for PWT data. Default is 80 (which is the 
                  most recent version).

    """        
    try:
        # first download the PWT data
        url = ('http://www.rug.nl/research/GGDC/data/pwt/V' + str(version) +
               '/pwt' + str(version) + '.zip') 
        data_archive = zipfile.ZipFile(StringIO(urlopen(url).read()), 'r') 

        # next download the PWT program files
        url = ('http://www.rug.nl/research/GGDC/data/pwt/V' + str(version) +
               '/Programs.zip') 
        programs_archive = zipfile.ZipFile(StringIO(urlopen(url).read()), 'r') 

        # extract the zip archives
        data_archive.extractall()
        programs_archive.extractall()
        
        print 'PWT data and program files successfully downloaded!'
        
    except:
        
        print 'Download failed! Perhaps you are not connected to the internet?'

def load_pwt_data(version=80, deltas=False):
    """
    Load the Penn World Tables data as a Pandas Panel object. Function expects
    a local copy of the PWT data file in the working directory. If no local copy
    exists, then one will be downloaded automatically (assuming that you have a
    working internet connection!).

    Arguments:
 
        version:  (int) Version number for PWT data. Default is 80 (which is the 
                  most recent version).
                  
        deltas:   (boolean) Whether or not you wish to load the data on 
                  depreciation rates (which is included in a separate .dta 
                  file). Default is False.
    Returns:

        pwt:    A Pandas Panel object containing the Penn World Tables data.
        
    TODO: Work out a way to merge Pandas Panel objects.
    
    """        
    # first check for a local copy of PWT
    try: 
        data = pd.read_stata('pwt' + str(version) + '.dta')
                
    # otherwise, 
    except IOError:  
        get_pwt_data(version)
        data = pd.read_stata('pwt' + str(version) + '.dta')
     
    # convert to a panel object   
    data.index =[data['year'], data['countrycode']]
    data = data.to_panel()
        
    if deltas == True:
        
        # load the separate file containing the depreciation rates
        dep_rates = pd.read_stata('depreciation_rates.dta')
        
        # convert to panel 
        dep_rates.index = [dep_rates['year'], dep_rates['countrycode']]
        dep_rates = dep_rates.to_panel()

        return [data, dep_rates]
        
    else:
        
        return data
        
    
