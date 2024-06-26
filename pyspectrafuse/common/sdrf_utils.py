from pathlib import Path
import pandas as pd
import re


class SdrfUtil:

    def __init__(self, folder_path: str):
        self.folder = folder_path

    @staticmethod
    def get_sdrf_file_path(folder: str) -> str:
        """
        get sdrf file path in project folder
        :param folder: folder obtain project files to cluster
        """
        directory_path = Path(folder)
        files = directory_path.rglob('*.sdrf.tsv')
        try:
            return str(next(files))
        except StopIteration:
            raise FileNotFoundError(f'There is no sdrf file in {folder}')

    @staticmethod
    def get_metadata_dict_from_sdrf(sdrf_folder: str) -> dict:
        """
        Read the sdrf file to obtain the relationship between samples and instruments and species in a project
        :param sdrf_folder: sdrf folder
        :return:
        """
        sdrf_df = pd.read_csv(sdrf_folder, sep='\t')
        sdrf_feature_df = pd.DataFrame()
        try:
            sdrf_feature_df = sdrf_df.loc[:, ['comment[data file]', 'characteristics[organism]', 'comment[instrument]']]
        except KeyError:
            print(f'{sdrf_folder} file has some format error, please check the col index format.')

        sdrf_feature_df['comment[data file]'] = sdrf_feature_df['comment[data file]'].apply(lambda x: x.split('.')[0])
        sdrf_feature_df['comment[instrument]'] = sdrf_feature_df['comment[instrument]'].apply(
            lambda x: re.search(r'=(.*)', [i for i in x.split(';') if i.startswith("NT=")].pop()).group(1))

        sdrf_feature_df['organism_instrument'] = sdrf_feature_df[
            ['characteristics[organism]', 'comment[instrument]']].apply(lambda x: list(x), axis=1)
        sample_info_dict = sdrf_feature_df.set_index('comment[data file]')['organism_instrument'].to_dict()
        return sample_info_dict


