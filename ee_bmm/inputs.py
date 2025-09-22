#--------------------------------
# Name:         inputs.py
# Purpose:      Common INI reading/parsing functions
# Python:       3.9
#--------------------------------

from builtins import input
import datetime
import logging
import os
import sys

import configparser

import ee_bmm.utils as utils


def read(ini_path):
    logging.debug('\nReading Input File')
    # Open config file
    config = configparser.ConfigParser()
    try:
        config.read(ini_path)
    except Exception as e:
        logging.error('\nERROR: Input file could not be read, '
                      'is not an input file, or does not exist\n'
                      '  ini_path = {}\n{}\n'.format(ini_path, e))
        sys.exit()

    # Force conversion of unicode to strings
    ini = dict()
    for section in config.keys():
        ini[str(section)] = {}
        for k, v in config[section].items():
            ini[str(section)][str(k)] = v
    return ini


def parse_section(ini, section):
    logging.debug('Checking {} section'.format(section))
    if section not in ini.keys():
        logging.error(
            '\nERROR: Input file does not have an {} section'.format(section))
        sys.exit()

    if section == 'INPUTS':
        parse_inputs(ini)
    elif section == 'ZONAL_STATS':
        parse_zonal_stats(ini)
    elif section == 'IMAGES':
        parse_images(ini)


def get_param(ini, section, input_name, output_name, get_type,
              default='MANDATORY'):
    """Get INI parameters by type and set default values

    Args:
        ini (dict): Nested dictionary of INI file keys/values
        section (str): Section name
        input_name (str): Parameter name in INI file
        output_name (str): Parameter name in code
        get_type (): Python type
        default (): Default value to use if parameter was not set.
            Defaults to "MANDATORY".
            "MANDATORY" will cause script to exit if key does not exist.
    """

    try:
        if get_type is bool:
            ini[section][output_name] = (
                ini[section][input_name].lower() == 'true')
            # ini[section][output_name] = distutils.util.strtobool(
            #     ini[section][input_name])
            # ini[section][output_name] = ini.getboolean(section, input_name)
            # ini[section][output_name] = ini[section].getboolean(input_name)
        elif get_type is int:
            ini[section][output_name] = int(ini[section][input_name])
        elif get_type is float:
            ini[section][output_name] = float(ini[section][input_name])
        elif get_type is list:
            ini[section][output_name] = str(ini[section][input_name])
            # Parsing strings to list is handled in each section separately
            # ini[section][output_name] = [
            #     x.strip()
            #     for x in str(ini[section][input_name]).split(',')]
        else:
            ini[section][output_name] = str(ini[section][input_name])
            # Convert 'None' (strings) to None
            if ini[section][output_name].lower() in ['none', '']:
                ini[section][output_name] = None
    except (KeyError, configparser.NoOptionError):
        if default == 'MANDATORY':
            logging.error(
                '\nERROR: {} was not set in the INI, exiting\n'.format(
                    input_name))
            sys.exit()
        else:
            ini[section][input_name] = default
            ini[section][output_name] = default
            logging.debug('  Setting {} = {}'.format(
                input_name, ini[section][output_name]))
    except ValueError:
        logging.error('\nERROR: Invalid value for "{}"'.format(
            input_name))
        sys.exit()
    except Exception as e:
        logging.error('\nERROR: Unhandled error\n  {}'.format(e))
        sys.exit()

    # If the parameter is renamed, remove the old name/parameter
    if input_name != output_name:
        del ini[section][input_name]


def parse_inputs(ini, section='INPUTS'):
    # MANDATORY PARAMETERS
    # section, input_name, output_name, description, get_type
    param_list = [
        ['gcloud_project_id', 'gcloud_project_id', str],
        ['zone_shp_path', 'zone_shp_path', str],
        ['zone_field', 'zone_field', str],
        ['hyd_area', 'hyd_area', str],
        ['start_year', 'start_year', int],
        ['end_year', 'end_year', int],
        ['start_doy', 'start_doy', int],
        ['end_doy', 'end_doy', int],
        ['path_keep_list', 'path_keep_list', list],
        ['row_keep_list', 'row_keep_list', list],
        # ['img_zone_path', 'img_zone_path', str], # this is for the custom image prep script to mask pixels to the zone boundary instead of a whole HA
    ]
    for input_name, output_name, get_type in param_list:
        get_param(ini, section, input_name, output_name, get_type)

    # OPTIONAL PARAMETERS
    # param_section, input_name, output_name, get_type, default
    param_list = [
        # subset of sites
        ['grid_skip_list', 'grid_skip_list', list, None],
        ['grid_buffer_flag', 'grid_buffer_flag', bool, False],
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)

    # phreatophyte shapefile path
    ini[section]['zone_filename'] = os.path.splitext(
        os.path.basename(ini[section]['zone_shp_path']))[0]

    # Zone shapefile
    # Eventually this may not be a mandatory parameter
    if not os.path.isdir(os.path.dirname(ini[section]['zone_shp_path'])):
        logging.error(
            '\nERROR: The phreatophyte shapefile does not exist, exiting\n'
            '  {}'.format(os.path.dirname(ini[section]['zone_shp_path'])))
        sys.exit()
    elif not os.path.isfile(ini[section]['zone_shp_path']):
        logging.error(
            '\nERROR: The phreatophyte shapefile does not exist, '
            'exiting\n  {}'.format(ini[section]['zone_shp_path']))
        sys.exit()

    # Start/end year
    if (ini[section]['start_year'] and ini[section]['end_year'] and
            ini[section]['end_year'] < ini[section]['start_year']):
        logging.error(
            '\nERROR: End year must be >= start year')
        sys.exit()
    default_end_year = datetime.datetime.today().year + 1
    if ((ini[section]['start_year'] and
            ini[section]['start_year'] not in range(1984, default_end_year)) or
        (ini[section]['end_year'] and
            ini[section]['end_year'] not in range(1984, default_end_year))):
        logging.error('\nERROR: Year must be an integer from 1984-{}'.format(
            default_end_year - 1))
        sys.exit()

    # Start/end DOY
    if (ini[section]['start_doy'] and
            ini[section]['start_doy'] not in range(1, 367)):
        logging.error(
            '\nERROR: Start DOY must be an integer from 1-366')
        sys.exit()
    elif (ini[section]['end_doy'] and
            ini[section]['end_doy'] not in range(1, 367)):
        logging.error('\nERROR: End DOY must be an integer from 1-366')
        sys.exit()

    # Convert path/row ranges to list
    if ini[section]['path_keep_list']:
        ini[section]['path_keep_list'] = sorted(list(
            utils.parse_int_set(ini[section]['path_keep_list'])))
    if ini[section]['row_keep_list']:
        ini[section]['row_keep_list'] = sorted(list(
            utils.parse_int_set(ini[section]['row_keep_list'])))

    # optional skip list of gridMET cells
    if ini[section]['grid_skip_list']:
        ini[section]['grid_skip_list'] = sorted(list(
            utils.parse_int_set(ini[section]['grid_skip_list'])))

    # ini[section]['img_zone_path'] = ini[section]['img_zone_path']


def parse_zonal_stats(ini, section='ZONAL_STATS'):
    """"""
    # Get the list of Landsat products to compute
    # DEADBEEF - What should the default Landsat products be?
    param_list = [
        ['output_folder', 'output_folder', str],
    ]
    for input_name, output_name, get_type in param_list:
        get_param(ini, section, input_name, output_name, get_type)

    param_list = [
        ['replacement_rate_etg_flag', 'replacement_rate_etg_flag', bool, True],
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)   


def parse_images(ini, section='IMAGES'):
    """"""
    # Get the list of Landsat products to compute
    # DEADBEEF - What should the default Landsat products be?
    param_list = [
        ['bmm_regress_vers', 'bmm_regress_vers', str],
        ['export_dest', 'export_dest', str],
        ['band_list', 'band_list', list],
        ['water_index_id', 'water_index_id', str],
    ]
    for input_name, output_name, get_type in param_list:
        get_param(ini, section, input_name, output_name, get_type)

    if ini[section]['band_list']:
        ini[section]['band_list'] = [x.lstrip() for x in ini[section]['band_list'].split(',')]

    param_list = [
        ['min_etg_threshold', 'min_etg_threshold', float, 0.05],
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)

    