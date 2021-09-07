"""
This script is supposed to run after all experiment data has been collected.
It should give a list of plates/wells that had errors and save it to the experiment directory.
Further it should save this information to the experiments DB table.
A column "valid" is supposed to hold the information and will be "ok" for valid experiments or a verbose
description of the error.

Errors considered in this script are:
- A manually curated list of errors, e.g. to record pipetting errors / mixups / precipitation etc.
- Transfer errors at NEXUS extracted from the provided transfer files
- Analysis errors at MoBiAS, e.g. integration of internal standard captures several (sometimes hundreds of peaks)

Additionally warnings should be generated for:
- Wells with low volume in last NEXUS survey
- Split peaks, e.g. 2 peaks for IS or main product
"""

import pandas as pd

from definitions import PLATES_DIR, PLATE_LIST_PATH
from utils import get_internal_standard_number


def get_manual_error_records(path):
    """
    Import a CSV file containing manually curated error records.
    :return:
    """
    errors = pd.read_csv(path)
    return errors.values.tolist()


def read_nexus_transfer_errors(path, error_string):
    """
    Import a CSV transfer file and add everything under exceptions to error records.
    :param path: path to transfer file (csv)
    :param error_string: verbose description of the error to add to the error list
    :return: List of transfer errors
    """
    transfers = pd.read_csv(path, header=None,
                            names=['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'Source Well',
                                   'Source Concentration', 'Source Concentration Units', 'Destination Plate Name',
                                   'Destination Plate Barcode', 'Destination Well', 'Destination Concentration',
                                   'Destination Concentration Units', 'Compound Name', 'Transfer Volume',
                                   'Actual Volume', 'Transfer Status', 'Current Fluid Height', 'Current Fluid Volume',
                                   '% DMSO'])
    if '[EXCEPTIONS]' in transfers.values:
        exceptions_idx = transfers.loc[transfers['Source Plate Name'] == '[EXCEPTIONS]'].index
        details_idx = transfers.loc[transfers['Source Plate Name'] == '[DETAILS]'].index
        transfers['plate'] = transfers['Destination Plate Barcode'].str.strip('Synthesis')
        transfers['row'] = transfers['Destination Well'].str[0]
        transfers['column'] = transfers['Destination Well'].str[1:]
        transfers['error'] = error_string
        exceptions = transfers.iloc[exceptions_idx[0] + 2:details_idx[0] - 1]
        return exceptions[['plate', 'row', 'column', 'error']].values.tolist()
    else:
        return []


def read_nexus_survey_errors(path, error_string, volume_threshold):
    """
    Import a CSV survey file and add all wells below the threshold volume to error records.
    :param path: path to survey file (csv)
    :param error_string: verbose description of the error to add to the error list
    :param volume_threshold: value below which the volume (in uL) is considered as an error
    :return: List of survey errors
    """
    volumes = pd.read_csv(path, header=7)
    volumes = volumes.drop(columns=['Survey Status']).dropna(how='any', axis=0)
    volumes['plate'] = volumes['Source Plate Barcode'].str.strip('Synthesis')
    volumes['row'] = volumes['Source Well'].str[0]
    volumes['column'] = volumes['Source Well'].str[1:]
    volumes['error'] = error_string
    exceptions = volumes.loc[
        (volumes['Survey Fluid Volume'] < volume_threshold) & (volumes['column'].astype(int).between(3, 22))]
    return exceptions[['plate', 'row', 'column', 'error']].values.tolist()


def read_mobias_analysis_errors(path, plate_number):
    """
    Import a CSV analysis results file and identify errors in the analysis
    :param path: path to results file (csv)
    :param plate_number: int, number of the plate within one experiment
    :return: List of analysis errors
    """
    # What are analysis errors?
    # - too many peaks for one compound
    # - too low (or high?) value for internal standard
    # (this list may not yet be complete)
    # TODO the thresholds for this need to be revisited later

    results = pd.read_csv(path, header=3, encoding='latin-1', skip_blank_lines=False)
    results['plate'] = plate_number
    results['row'] = results['Vial Pos'].apply(lambda x: x.split('-')[1])
    results['column'] = results['Vial Pos'].apply(lambda x: x.split('-')[2])
    results['error_1'] = 'ok'
    results['error_2'] = 'ok'
    results['error_3'] = 'ok'
    results['error_4'] = 'ok'
    internal_standard_number = get_internal_standard_number(path.parent)
    # identify where product A gives too many peaks
    results.loc[(results['SumF1 Cmp'] > 1), ['error_1']] = f'too many peaks for product A (' + results[
        'SumF1 Cmp'].astype('str') + ')'
    # identify where IS gives to many peaks
    results.loc[(results[f'{internal_standard_number} Cmp'] > 1), ['error_2']] = f'too many peaks for IS (' + results[
        f'{internal_standard_number} Cmp'].astype('str') + ')'
    # identify where IS response deviates >20% from mean
    mean_response_area = results[f'{internal_standard_number} Area'].mean()
    results.loc[
        ~results[f'{internal_standard_number} Area'].between(mean_response_area * 0.8, mean_response_area * 1.2), [
            'error_3']] = 'WARNING: IS response differs >20% from mean'
    # identify where IS response deviates >50% from mean
    mean_response_area = results[f'{internal_standard_number} Area'].mean()
    results.loc[
        ~results[f'{internal_standard_number} Area'].between(mean_response_area * 0.5, mean_response_area * 1.5), [
            'error_4']] = 'ERROR: IS response differs >50% from mean'

    # aggregate the four error strings into one
    def aggregate_errors(series):
        errors = ''
        for e in ['error_1', 'error_2', 'error_3', 'error_4']:
            if series[e] != 'ok':
                if errors != '':
                    errors += '; '
                errors += series[e]
        return errors

    results['error'] = results.apply(aggregate_errors, axis=1)
    exceptions = results.loc[results['error'] != '']
    return exceptions[['plate', 'row', 'column', 'error']].values.tolist()


def get_nexus_errors(nexus_dir):
    """
    Iterate through all NEXUS transfer and survey files and collect the errors
    :param nexus_dir: directory under which NEXUS transfer files are stored
    :return: List of transfer errors
    """
    error_list = []
    # task 1: Iterate I_M transfers
    for child in (nexus_dir / 'I_M').iterdir():
        if child.is_file() and 'transfer' in child.name:
            error_list += read_nexus_transfer_errors(child, 'I/M transfer error')
    # task 2: Iterate T transfers
    for child in (nexus_dir / 'T').iterdir():
        if child.is_file() and 'transfer' in child.name:
            error_list += read_nexus_transfer_errors(child, 'T transfer error')
    # task 3: Iterate dilution transfers
    for child in (nexus_dir / 'dilution').iterdir():
        if child.is_file() and 'transfer' in child.name:
            error_list += read_nexus_transfer_errors(child, 'Dilution transfer error')
    # task 4: Iterate plate survey before dilution
    for child in (nexus_dir / 'dilution').iterdir():
        if child.is_file() and 'survey' in child.name:
            error_list += read_nexus_survey_errors(child, 'Dilution survey low volume', volume_threshold=2.2)
    return error_list


def get_mobias_errors(exp_nr):
    """
    :param exp_nr: number of the experiment under consideration
    :param plate_list_path: path to csv-file containing list of experiments / plate numbers
    :return: List of transfer errors
    """
    error_list = []
    plate_list_df = pd.read_csv(PLATE_LIST_PATH)
    plate_list = plate_list_df.loc[
        plate_list_df['exp_nr'] == exp_nr, ['plate_nr', 'lab_journal_number', 'results_file_name']].values.tolist()
    for plate in plate_list:
        error_list += read_mobias_analysis_errors(PLATES_DIR / plate[1] / plate[2], int(plate[0]))
    return error_list


def aggregate_error_list(error_list):
    """Aggregate the error list into a dateframe that has at most one row per well"""
    # use intermediate dictionary for the reduction
    error_dict = {}
    for error in error_list:
        try:
            error_dict[f'{error[0]}-{error[1]}-{error[2]}'] += f'; {error[3]}'
        except KeyError:
            error_dict[f'{error[0]}-{error[1]}-{error[2]}'] = f'{error[3]}'
    # transfer dict content to df for convenient CSV writing
    error_df = pd.DataFrame(columns=['well', 'errors'], data=error_dict.items())
    error_df['plate'] = error_df['well'].apply(lambda x: x.split('-')[0])
    error_df['row'] = error_df['well'].apply(lambda x: x.split('-')[1])
    error_df['column'] = error_df['well'].apply(lambda x: x.split('-')[2])
    return error_df.loc[:, ['plate', 'row', 'column', 'errors']]


def main(exp_dir, exp_dir_name):
    error_list = []
    error_list += get_manual_error_records(exp_dir / 'manual_errorlist.csv')
    error_list += get_nexus_errors(exp_dir / 'transfer_files')
    error_list += get_mobias_errors(exp_dir_name)
    error_df = aggregate_error_list(error_list)
    print(error_df)
    print(len(error_df))
    error_df.to_csv(exp_dir / 'extracted_errors.csv', index=False)
    return


if __name__ == '__main__':
    exp_dir_name = 'exp1'
    exp_dir = PLATES_DIR / exp_dir_name
    main(exp_dir, exp_dir_name)
