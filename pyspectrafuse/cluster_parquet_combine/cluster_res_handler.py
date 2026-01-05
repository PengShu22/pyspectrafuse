import pandas as pd
import numpy as np
from pathlib import Path
import re


class ClusterResHandler:

    def __init__(self, clu_thr=30, dirname=''):
        self.cluster_threshold = clu_thr
        self.project_dir = dirname

    def walk_dir(self):
        # 根据电荷信息对找到的文件进行分组，并将具有相同电荷的聚类结果文件放入同一个列表中：
        charge_groups = self.get_charge_group()
        charge_cluster_res_groups = {}
        if charge_groups:
            print(f"charge_groups: {charge_groups}")
            for charge, charge_tsv_lst in charge_groups.items():
                charge_df_lst = []
                for charge_tsv in charge_tsv_lst:
                    df = self.read_cluster_tsv(charge_tsv)
                    sample_info = str(Path(*charge_tsv.parts[-4:-1])).replace("\\",
                                                                              "/")  # ('Homo sapiens', 'Q Exactive', 'charge2')
                    # 对于一个聚类结果文件，补上其[物种/仪器/电荷信息], 最终同个电荷的只对应一个parquet,parquet按照电荷切分
                    df.loc[:, "mgf_path"] = df.loc[:, "mgf_path"].apply(
                        lambda x: sample_info + "/mgf files/" + x)
                    df['mgf_path'] = df.apply(lambda row: f"{row['mgf_path']}/{row['index']}", axis=1)
                    charge_df_lst.append(df)

                charge_df = pd.DataFrame(np.vstack(charge_df_lst), columns=['mgf_path', 'index', 'cluster_accession'])
                charge_df = charge_df.loc[:, ['mgf_path', 'cluster_accession']]
                mgf_ind_clu_acc_dict = pd.Series(charge_df.cluster_accession.values, index=charge_df.mgf_path).to_dict()

                if charge not in charge_cluster_res_groups:
                    charge_cluster_res_groups[charge] = mgf_ind_clu_acc_dict
        return charge_cluster_res_groups

    def get_charge_group(self):
        """
        # 根据电荷信息对找到的文件进行分组，并将具有相同电荷的聚类结果文件放入同一个列表中：
        :return:
        """
        charge_groups = {}
        for file_path in Path(self.project_dir).rglob('charge*'):
            print(f"file_path: {file_path}")
            for charge_file in file_path.rglob("*"):
                # 检查是否为文件并且匹配正则表达式
                pattern = self.get_pattern()
                if charge_file.is_file() and pattern.search(charge_file.name):
                    print(f'找到文件: {charge_file}')
                    charge = charge_file.parent.name
                    if charge not in charge_groups:
                        charge_groups[charge] = []
                    charge_groups[charge].append(charge_file)

        return charge_groups

    def get_pattern(self):
        return re.compile(f'MaRaCluster\\.clusters_p{self.cluster_threshold}\\.tsv')

    @staticmethod
    def get_cluster_dict(tsv_file_path: Path, species: str, instrument: str, charge: str):
        """

        :param tsv_file_path:
        :param species:
        :param instrument:
        :param charge:
        :return:
        """
        clu_df = pd.read_csv(tsv_file_path, sep='\t', header=None)
        clu_df.columns = ['mgf_path', 'index', 'cluster_accession']
        clu_df.dropna(axis=0, inplace=True)  # 删除空行

        sample_info = f"{species}/{instrument}/{charge}"
        clu_df.loc[:, "mgf_path"] = clu_df.loc[:, "mgf_path"].apply(
            lambda x: sample_info + "/mgf files/" + x)
        clu_df['mgf_path'] = clu_df.apply(lambda row: f"{row['mgf_path']}/{row['index']}", axis=1)
        clu_df = clu_df.loc[:, ['mgf_path', 'cluster_accession']]
        mgf_ind_clu_acc_dict = pd.Series(clu_df.cluster_accession.values, index=clu_df.mgf_path).to_dict()

        return mgf_ind_clu_acc_dict



