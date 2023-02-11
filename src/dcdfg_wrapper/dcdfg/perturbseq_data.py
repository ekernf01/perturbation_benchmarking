import numpy as np
import scanpy as sc
from torch.utils.data import Dataset
from tqdm import tqdm


class PerturbSeqDataset(Dataset):
    """
    A generic class for simulation data loading and extraction, as well as pre-filtering of interventions
    NOTE: the 0-th regime should always be the observational one
    """

    def __init__(
        self,
        adata
    ) -> None:
        """
        :param str adata: AnnData object
        """
        super(PerturbSeqDataset, self).__init__()
        # load data, converting sparse to dense 
        try:
            adata.X = adata.X.toarray()
        except AttributeError:
            pass
        self.adata = adata
        all_data, all_masks, all_regimes = self.set_up_masks()

        # index of all regimes
        self.all_regimes_list = np.unique(all_regimes)
        obs_regime = np.unique(
            all_regimes[np.where([mask == [] for mask in all_masks])[0]]
        )
        assert len(obs_regime) == 1
        obs_regime = obs_regime[0]

        self.data = all_data
        self.regimes = all_regimes
        self.masks = np.array(all_masks, dtype=object)
        self.intervention = True

        self.num_regimes = np.unique(self.regimes).shape[0]
        self.num_samples = self.data.shape[0]
        self.dim = self.data.shape[1]

    def __getitem__(self, idx):
        # binarize mask from list
        masks_list = self.masks[idx]
        masks = np.ones((self.dim,))
        for j in masks_list:
            masks[j] = 0
        return (
            self.data[idx].astype(np.float64),
            masks.astype(np.float64),
            self.regimes[idx],
        )

    def __len__(self):
        return self.data.shape[0]

    def assign_regimes(self):
        # check gene sets and ensure matching with measurements
        obs_genes = {}
        unfound_genes = {}
        targets = []
        for index, row in tqdm(self.adata.obs.iterrows(), total=self.adata.n_obs):
            current_target = []
            if not row["is_control"]:                # Only cells with treatment
                # get all guides in cells
                sg = row["perturbation"].split(",")
                # get gene name by stripping guide specific info
                sg_genes = [guide.rsplit("_", maxsplit=1)[0] for guide in sg]
                for gene in sg_genes:
                    if gene in self.adata.var.index:
                        # gene is found
                        current_target += [gene]
                        if gene not in obs_genes:
                            obs_genes[gene] = 1
                        else:
                            obs_genes[gene] += 1
                    else:
                        if gene not in unfound_genes:
                            unfound_genes[gene] = 1
                        else:
                            unfound_genes[gene] += 1
            # end gene list
            targets += [",".join(current_target)]

        print(f"""
            Perturbed and measured genes: {len(obs_genes)}
        Perturbed but not measured genes: {len(unfound_genes)} {list(unfound_genes.keys())[:5]}
        """)
        regimes = np.unique(targets, return_inverse=True)[1]

        self.adata.obs["targets"] = targets
        self.adata.obs["regimes"] = regimes


    def set_up_masks(self, normalized_data=True):
        """
        Set up the masks and regimes
        """
        if normalized_data:
            data = self.adata.X
        else:
            data = self.adata.layers["counts"]

        print(f"Data shape {data.shape}")
        self.assign_regimes()

        # Load intervention masks and regimes
        regimes = self.adata.obs["regimes"].astype(int)
        masks = []
        unrecognized_genes = []

        # create map gene name -> gene number
        gene_map = {}
        for i, gene in enumerate(self.adata.var.index):
            gene_map[gene] = i
        for index, row in tqdm(self.adata.obs.iterrows(), total=self.adata.n_obs):
            mask = []
            if row["targets"] != "":
                for x in row["targets"].split(","):
                    if x in gene_map:
                        mask += [gene_map[x]]
                    else:
                        unrecognized_genes = np.union1d(unrecognized_genes, [x])
            masks.append(mask)

        Warning("couldn't find genes after filtering:" + str(unrecognized_genes))
        return data, masks, regimes
