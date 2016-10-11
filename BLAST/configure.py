#!/usr/bin/env python3

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)

import configparser

SECTION_NAME = "BLAST_SETTINGS"
CONFIG_VALUES = {
    'subject' : 'str',
    'database' : 'str',
    'evalue' : 'float',
    'max_hits' : 'int',
    'max_hsps' : 'int',
    'identity' : 'float',
    'keep_query' : 'bool'
}

#   Make a configuration file
def make_config(args):
    """Make a configuration file for BLAST"""
    config = configparser.ConfigParser()
    config.add_section(SECTION_NAME)
    for option in CONFIG_VALUES:
        try:
            config.set(SECTION_NAME, option, str(args[option]))
        except KeyError:
            continue
    with open(args['outfile'], 'w') as configfile:
        config.write(configfile)
    print("Config file can be found at", args['outfile'], file=sys.stderr)


#   Parse a configuration file
def parse_config(config_file):
    """Parse a configuration file for BLAST"""
    try:
        print("Using configuration file", config_file, file=sys.stderr)
        config = configparser.ConfigParser()
        config.read_file(open(config_file))
    except FileNotFoundError as error:
        sys.exit("Failed to find " + error.filename)
    bconf = {}
    for option, opt_type in CONFIG_VALUES.items():
        try:
            if opt_type == 'int':
                bconf[option] = config.getint(SECTION_NAME, option)
            elif opt_type == 'float':
                bconf[option] = config.getfloat(SECTION_NAME, option)
            elif opt_type == 'bool':
                bconf[option] = config.getboolean(SECTION_NAME, option)
            else:
                bconf[option] = config.get(SECTION_NAME, option)
        except configparser.NoOptionError:
            continue
    return bconf
