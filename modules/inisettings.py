import configparser as cfgParser
import os


#------------------------------------------------------------------------------

#https://python101.pythonlibrary.org/chapter14_config_parser.html
configFileName='pixek.ini'

def setDefaultIni():
    
    if os.path.exists(configFileName):
        return
    
    ini=cfgParser.ConfigParser()
    ini['default']={'apath' : os.getcwd(),  #application path
                    'inipath': '',   #input data directory path
                    'outpath': '',    #output data directory path
                    'inpfile': ''
                    }
    
    ini['font.title']={'name' : 'Helvetica',
                       'size' : '14',
                       'weight': 'normal'}
    
    ini['font.default']={'name' : 'Helvetica',
                         'size' : '12',
                         'weight': 'normal'}
        
    ini['reduction']={'name' : ' Ignore',
                           'fromTo' : '1:20'}
    
    ini['calibration']={    'onoff' : '0',
                            'method' : 'Ignore',
                            'ax+b' : '1 ; 0', 
                            'angle' : 'random',
                            'doze'  : '10',
                            'ula'   : '4470',
                            'uma'   : '1040',
                            'A'     : '0.0224',
                            'B'     : '-0.15',
                            'C'     :  '1'
                            }
    
    ini['filtering']={     'method' : 'Ignore',
                           'ksize' : '1',
                           'polyorder' : '2',
                           'derivative': '0',
                           'plotraw': '1'}
    
                                       
    ini['wavelet']={'family'  : 'Daubechies db',
                    'name' : 'db4',
                     'levels': '4',
                    }
    
    ini['info']={'material'  : 'UO2',
                    'axis' : '110',
                    'plane': '001',
                    'angle': '0',
                    'dose_e': '10',
                    'dose': '10',
                    }
    
    with open(configFileName,'w') as configfile:
        ini.write(configfile)


def getConfig(path):    
    #if not os.path.exists(path):
    #    create_config(path)

    config = cfgParser.ConfigParser()
    config.read(path)
    return config

def getValue(section,key):
    """
    Print out a setting
    """
    config = getConfig(configFileName)
    value = config.get(section, key)
    return value

def setValue(section, key, value):
    """
    Update a setting
    """
    config = getConfig(configFileName)
    config.set(section, key, value)
    with open(configFileName, "w") as config_file:
        config.write(config_file)
