from pathlib import Path
from typing import Dict, List
import pandas as pd
import re
import logging

logger = logging.getLogger(__name__)


class SdrfUtil:
    def __init__(self, folder_path: str):
        self.folder: str = folder_path

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
    def get_metadata_dict_from_sdrf(sdrf_folder: str) -> Dict[str, List[str]]:
        """Read the sdrf file to obtain the relationship between samples and instruments and species in a project.
        
        Args:
            sdrf_folder: Path to SDRF file
            
        Returns:
            Dictionary mapping data file names to [organism, instrument] lists
            
        Raises:
            KeyError: If required columns are missing from SDRF file
            ValueError: If instrument field cannot be parsed
        """
        # print("reading sdrf file: {}".format(sdrf_folder))
        sdrf_df = pd.read_csv(sdrf_folder, sep='\t')
        sdrf_feature_df = pd.DataFrame()
        try:
            sdrf_feature_df = sdrf_df.loc[:, ['comment[data file]', 'Characteristics[organism]', 'comment[instrument]']]
        except KeyError as e:
            missing_cols = set(['comment[data file]', 'Characteristics[organism]', 'comment[instrument]']) - set(sdrf_df.columns)
            raise KeyError(f'{sdrf_folder} file has format error. Missing columns: {missing_cols}. Available columns: {list(sdrf_df.columns)}') from e

        # print(sdrf_feature_df.head())
        # print(sdrf_feature_df.columns.tolist())

        # Vectorized string operation - much faster than apply()
        sdrf_feature_df['comment[data file]'] = sdrf_feature_df['comment[data file]'].str.split('.').str[0]
        
        def extract_instrument(x):
            """Extract instrument name from comment[instrument] field."""
            try:
                nt_matches = [i for i in x.split(';') if i.startswith("NT=")]
                if not nt_matches:
                    raise ValueError(f"No NT= field found in instrument comment: {x}")
                match = re.search(r'=(.*)', nt_matches[0])
                if not match:
                    raise ValueError(f"Could not extract instrument from: {x}")
                return match.group(1)
            except (IndexError, AttributeError) as e:
                raise ValueError(f"Error parsing instrument field: {x}") from e
        
        sdrf_feature_df['comment[instrument]'] = sdrf_feature_df['comment[instrument]'].apply(extract_instrument)

        # Vectorized operation - convert to list of lists
        sdrf_feature_df['organism_instrument'] = sdrf_feature_df[
            ['Characteristics[organism]', 'comment[instrument]']].values.tolist()
        sample_info_dict = sdrf_feature_df.set_index('comment[data file]')['organism_instrument'].to_dict()
        return sample_info_dict


