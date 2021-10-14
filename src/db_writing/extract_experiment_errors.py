"""
This script is supposed to run after all experiment data has been collected.
It should give a list of plates/wells that had errors and save it to the experiment directory.
Further it should save this information to the experiments DB table.
A column "valid" is supposed to hold the information and will be NULL for valid experiments or a verbose
description of the error. Descriptions contain either the keyword WARNING or ERROR to facilitate adjusting sensitivity

The following are considered ERRORS:
- A manually curated list of errors, e.g. to record pipetting errors / mixups / precipitation etc.
- Transfer errors at NEXUS extracted from the provided transfer files
- Wells with low volume in last NEXUS survey (< 2.2 uL)
- More than 4 peaks for the main product in MoBiAS analysis
- More than 1 peak for the internal standard in MoBiAS analysis
- Deviations of internal standard peak area by >80% from the mean for the plate

The following are considered WARNINGS::
- 2-4 peaks for the main product in MoBiAS analysis
- Deviations of IS peak area by 50-80% from the mean for the plate
"""
import sqlite3
import warnings
from pathlib import Path

import pandas as pd

from definitions import PLATES_DIR, PLATE_LIST_PATH, DB_PATH
from utils import get_internal_standard_number

# configuration
exp_dir_name = 'exp_test4'
exp_nr = 99004
skip_nexus = True  # set to True for reactions that were not conducted at NEXUS


def get_manual_error_records(path: Path) -> list:
    """
    Import a CSV file containing manually curated error records.
    :return:
    """
    try:
        errors = pd.read_csv(path)
    except FileNotFoundError:
        warnings.warn(f'No manual error list found in {path.parent}. Manual errors will not be included.')
        return []
    return errors.values.tolist()


def read_nexus_transfer_errors(path: Path, error_string: str) -> list:
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
        transfers['plate'] = transfers['Destination Plate Barcode'].str.strip('Synthesis').str.strip('Analysis')
        # ^ we can do this since the analysis plate is a 1:1 copy
        transfers['row'] = transfers['Destination Well'].str[0]
        transfers['column'] = transfers['Destination Well'].str[1:]
        transfers['error'] = error_string
        exceptions = transfers.iloc[exceptions_idx[0] + 2:details_idx[0]]
        return exceptions[['plate', 'row', 'column', 'error']].values.tolist()
    else:
        return []


def read_nexus_repeated_tranfers(path: Path, error_string: str) -> list:
    """
    Import a CSV transfer file and add everything under details to success records.
    Outside of this function, these should be substracted from the error list.
    :param path: path to transfer file (csv)
    :return: List of transfers successful on second try
    """
    transfers = pd.read_csv(path, header=None,
                            names=['Source Plate Name', 'Source Plate Barcode', 'Source Plate Type', 'Source Well',
                                   'Source Concentration', 'Source Concentration Units', 'Destination Plate Name',
                                   'Destination Plate Barcode', 'Destination Well', 'Destination Concentration',
                                   'Destination Concentration Units', 'Compound Name', 'Transfer Volume',
                                   'Actual Volume', 'Transfer Status', 'Current Fluid Height', 'Current Fluid Volume',
                                   '% DMSO'])
    if '[DETAILS]' in transfers.values:
        details_idx = transfers.loc[transfers['Source Plate Name'] == '[DETAILS]'].index
        transfers['plate'] = transfers['Destination Plate Barcode'].str.strip('Synthesis')
        transfers['row'] = transfers['Destination Well'].str[0]
        transfers['column'] = transfers['Destination Well'].str[1:]
        transfers['error'] = error_string
        successes = transfers.iloc[details_idx[0] + 2:]
        return successes[['plate', 'row', 'column', 'error']].values.tolist()
    else:
        return []


def read_nexus_survey_errors(path: Path, error_string: str, volume_threshold: float) -> list:
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


def read_mobias_analysis_errors(path: Path, plate_number: int) -> list:
    """
    Import a CSV analysis results file and identify errors in the analysis
    Errors are:
    - too many peaks for target compound or internal standard
    - too strong deviation in the peak area of interanal standard
    TODO the thresholds for this may need to be revisited later
    :param path: path to results file (csv)
    :param plate_number: number of the plate within one experiment
    :return: List of analysis errors
    """

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
    results.loc[(results['SumF1 Cmp'] > 1) & (results['SumF1 Cmp'] <= 4), [
        'error_1']] = f'WARNING: multiple peaks for product A (' + results[
        'SumF1 Cmp'].astype('str') + ')'
    results.loc[(results['SumF1 Cmp'] > 4), ['error_1']] = f'ERROR: too many peaks for product A (' + results[
        'SumF1 Cmp'].astype('str') + ')'
    # identify where IS gives to many peaks
    results.loc[(results[f'{internal_standard_number} Cmp'] > 1), ['error_2']] = f'ERROR: too many peaks for IS (' + \
                                                                                 results[
                                                                                     f'{internal_standard_number} Cmp'].astype(
                                                                                     'str') + ')'
    # identify where IS response deviates >50% from mean
    mean_response_area = results[f'{internal_standard_number} Area'].mean()
    results.loc[
        ~results[f'{internal_standard_number} Area'].between(mean_response_area * 0.5, mean_response_area * 1.5), [
            'error_3']] = 'WARNING: IS response differs >50% from mean'
    # identify where IS response deviates >80% from mean
    mean_response_area = results[f'{internal_standard_number} Area'].mean()
    results.loc[
        ~results[f'{internal_standard_number} Area'].between(mean_response_area * 0.2, mean_response_area * 1.8), [
            'error_4']] = 'ERROR: IS response differs >80% from mean'

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


def get_nexus_errors(nexus_dir: Path) -> list:
    """
    Iterate through all NEXUS transfer and survey files and collect the errors
    :param nexus_dir: directory under which NEXUS transfer files are stored
    :return: List of transfer errors
    """
    error_list = []
    # task 1: Iterate I_M transfers
    for child in (nexus_dir / 'I_M').iterdir():
        if child.is_file() and 'transfer' in child.name:
            error_list += read_nexus_transfer_errors(child, 'ERROR: I/M transfer error')
    if (nexus_dir / 'I_M' / 'repeated_transfers').exists():
        for child in (nexus_dir / 'I_M' / 'repeated_transfers').iterdir():
            if child.is_file() and 'transfer' in child.name:
                for successful_transfer in read_nexus_repeated_tranfers(child, 'ERROR: I/M transfer error'):
                    error_list.remove(successful_transfer)
    # task 2: Iterate T transfers
    for child in (nexus_dir / 'T').iterdir():
        if child.is_file() and 'transfer' in child.name:
            error_list += read_nexus_transfer_errors(child, 'ERROR: T transfer error')
    if (nexus_dir / 'T' / 'repeated_transfers').exists():
        for child in (nexus_dir / 'T' / 'repeated_transfers').iterdir():
            if child.is_file() and 'transfer' in child.name:
                for successful_transfer in read_nexus_repeated_tranfers(child, 'ERROR: T transfer error'):
                    error_list.remove(successful_transfer)
    # task 3: Iterate dilution transfers
    for child in (nexus_dir / 'dilution').iterdir():
        if child.is_file() and 'transfer' in child.name:
            error_list += read_nexus_transfer_errors(child, 'ERROR: Dilution transfer error')
    # task 4: Iterate plate survey before dilution
    for child in (nexus_dir / 'dilution').iterdir():
        if child.is_file() and 'survey' in child.name:
            error_list += read_nexus_survey_errors(child, 'ERROR: Dilution survey low volume (<2.2 uL)',
                                                   volume_threshold=2.2)
    return error_list


def get_mobias_errors(exp_nr: int) -> list:
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


def aggregate_error_list(error_list: list) -> pd.DataFrame:
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
    error_df = error_df.sort_values(by=['plate', 'row', 'column'], key=lambda x: pd.to_numeric(x, errors='ignore'))
    return error_df.loc[:, ['plate', 'row', 'column', 'errors']]


def save_errors_to_db(error_df: pd.DataFrame, exp_nr: int):
    con = sqlite3.connect(DB_PATH)
    cur = con.cursor()
    for i, row in error_df.iterrows():
        cur.execute('UPDATE experiments SET valid = ? WHERE exp_nr = ? AND plate_nr = ? AND well = ?;',
                    (row['errors'], exp_nr, row['plate'], f'{row["row"]}{row["column"]}')
                    )
    con.commit()
    con.close()
    return


def main(exp_dir: Path, exp_nr: int, skip_nexus=False):
    error_list = []
    error_list += get_manual_error_records(exp_dir / 'manual_errorlist.csv')
    if skip_nexus is False:
        error_list += get_nexus_errors(exp_dir / 'transfer_files')
    error_list += get_mobias_errors(exp_nr)
    error_df = aggregate_error_list(error_list)
    error_df.to_csv(exp_dir / 'extracted_errors.csv', index=False)
    n_errors = len(error_df.loc[error_df['errors'].str.contains('ERROR')])
    n_warnings = len(
        error_df.loc[(error_df['errors'].str.contains('WARNING')) & ~(error_df['errors'].str.contains('ERROR'))])
    print(f'Found {n_errors} wells with errors and {n_warnings} wells with warnings for experiment #{exp_nr}')
    print(f'A list of errors had been saved to {exp_dir / "extracted_errors.csv"}.')
    print('Check the list and make any required changes before continuing this program')
    input('Press ENTER to continue...')
    error_df_approved = pd.read_csv(exp_dir / 'extracted_errors.csv')
    print('Writing errors to database...')
    save_errors_to_db(error_df_approved, exp_nr)
    return


if __name__ == '__main__':
    exp_dir = PLATES_DIR / exp_dir_name
    main(exp_dir, exp_nr, skip_nexus=skip_nexus)
