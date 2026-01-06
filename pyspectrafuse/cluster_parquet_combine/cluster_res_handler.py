from typing import Dict, List
import pandas as pd
import numpy as np
from pathlib import Path
import re
import logging

logger = logging.getLogger(__name__)


class ClusterResHandler:
    """Handler for cluster result files from MaRaCluster."""

    def __init__(self, clu_thr: int = 30, dirname: str = ''):
        """Initialize ClusterResHandler.
        
        Args:
            clu_thr: Cluster threshold value (default: 30)
            dirname: Project directory name
        """
        self.cluster_threshold: int = clu_thr
        self.project_dir: str = dirname

    def walk_dir(self) -> Dict[str, Dict[str, str]]:
        """Walk directory and group cluster files by charge.
        
        Returns:
            Dictionary mapping charge to cluster result dictionary
        """
        charge_groups = self.get_charge_group()
        charge_cluster_res_groups: Dict[str, Dict[str, str]] = {}
        if charge_groups:
            logger.info(f"Found charge groups: {list(charge_groups.keys())}")
            for charge, charge_tsv_lst in charge_groups.items():
                charge_df_lst = []
                for charge_tsv in charge_tsv_lst:
                    df = self.read_cluster_tsv(str(charge_tsv))
                    sample_info = str(Path(*charge_tsv.parts[-4:-1])).replace("\\", "/")
                    # For a cluster result file, add [species/instrument/charge information]
                    # Vectorized string operations - much faster than apply()
                    df.loc[:, "mgf_path"] = f"{sample_info}/mgf files/" + df["mgf_path"]
                    df['mgf_path'] = df['mgf_path'] + '/' + df['index'].astype(str)
                    charge_df_lst.append(df)

                charge_df = pd.DataFrame(np.vstack(charge_df_lst), 
                                        columns=['mgf_path', 'index', 'cluster_accession'])
                charge_df = charge_df.loc[:, ['mgf_path', 'cluster_accession']]
                mgf_ind_clu_acc_dict = pd.Series(
                    charge_df.cluster_accession.values, 
                    index=charge_df.mgf_path).to_dict()

                if charge not in charge_cluster_res_groups:
                    charge_cluster_res_groups[charge] = mgf_ind_clu_acc_dict
        return charge_cluster_res_groups

    def get_charge_group(self) -> Dict[str, List[Path]]:
        """Group found files by charge information.
        
        Returns:
            Dictionary mapping charge to list of cluster TSV file paths
        """
        charge_groups = {}
        for file_path in Path(self.project_dir).rglob('charge*'):
            logger.debug(f"Checking file_path: {file_path}")
            for charge_file in file_path.rglob("*"):
                # Check if it's a file and matches the regular expression
                pattern = self.get_pattern()
                if charge_file.is_file() and pattern.search(charge_file.name):
                    logger.info(f'Found cluster file: {charge_file}')
                    charge = charge_file.parent.name
                    if charge not in charge_groups:
                        charge_groups[charge] = []
                    charge_groups[charge].append(charge_file)

        return charge_groups

    def get_pattern(self) -> re.Pattern:
        """Get regex pattern for cluster file names.
        
        Returns:
            Compiled regex pattern
        """
        return re.compile(f'MaRaCluster\\.clusters_p{self.cluster_threshold}\\.tsv')

    @staticmethod
    def read_cluster_tsv(tsv_file_path: str) -> pd.DataFrame:
        """Read cluster TSV file.
        
        Args:
            tsv_file_path: Path to TSV file
            
        Returns:
            DataFrame with cluster data
        """
        clu_df = pd.read_csv(tsv_file_path, sep='\t', header=None)
        clu_df.columns = ['mgf_path', 'index', 'cluster_accession']
        clu_df.dropna(axis=0, inplace=True)
        return clu_df

    @staticmethod
    def get_cluster_dict(tsv_file_path: Path, species: str, instrument: str, charge: str) -> Dict[str, str]:
        """Get cluster dictionary from TSV file.
        
        Args:
            tsv_file_path: Path to cluster TSV file
            species: Species name
            instrument: Instrument name
            charge: Charge value
            
        Returns:
            Dictionary mapping mgf_path to cluster_accession
        """
        clu_df = pd.read_csv(tsv_file_path, sep='\t', header=None)
        clu_df.columns = ['mgf_path', 'index', 'cluster_accession']
        clu_df.dropna(axis=0, inplace=True)  # Remove empty rows

        sample_info = f"{species}/{instrument}/{charge}"
        # Vectorized string operations - much faster than apply()
        clu_df.loc[:, "mgf_path"] = sample_info + "/mgf files/" + clu_df["mgf_path"]
        clu_df['mgf_path'] = clu_df['mgf_path'] + '/' + clu_df['index'].astype(str)
        clu_df = clu_df.loc[:, ['mgf_path', 'cluster_accession']]
        mgf_ind_clu_acc_dict = pd.Series(clu_df.cluster_accession.values, index=clu_df.mgf_path).to_dict()

        return mgf_ind_clu_acc_dict



