#!/usr/bin/env python3

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)

import configparser

from BLAST.runblastn import validate_db
from Objects.wrappers import NcbimakeblastdbCommandline

SECTION_NAME = "BLAST_SETTINGS"
CONFIG_VALUES = {
    # Option name   : option type
    'subject'       :   'str',
    'database'      :   'str',
    'evalue'        :   'float',
    'max_hits'      :   'int',
    'max_hsps'      :   'int',
    'identity'      :   'float',
    'keep_query'    :   'bool'
}

#   Make a configuration file
def make_config(args):
    """Make a configuration file for BLAST"""
    config = configparser.ConfigParser()
    config.add_section(SECTION_NAME)
    if 'database' in args.keys():
        try:
            validate_db(args['database'])
        except FileNotFoundError:
            # print("Failed to find", args['database'], "creating new database...", file=sys.stderr)
            print("Failed to find", args['database'], "exiting...", file=sys.stderr)
            sys.exit()
    for option in CONFIG_VALUES:
        try:
            config.set(SECTION_NAME, option, str(args[option]))
        except KeyError:
            continue
        else:
            print("Setting option", option, "with value", args[option], file=sys.stderr)
    with open(args['outfile'], 'w') as configfile:
        config.write(configfile)
    print("\nConfig file can be found at", args['outfile'], file=sys.stderr)


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
        else:
            print(option, "set to", bconf[option], file=sys.stderr)
    return bconf
