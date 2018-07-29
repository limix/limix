from numpy import dot

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from pandas import DataFrame

import limix

# Download and extract a limix file.
url = "http://rest.s3for.me/limix/smith08.hdf5.bz2"
limix.download(url, verbose=False)
limix.extract("smith08.hdf5.bz2", verbose=False)
data = limix.io.read_hdf5_limix("smith08.hdf5")
print(data["phenotype"]["col_header"].head())

G_all = data["genotype"]["matrix"]
geno_metadata = data["genotype"]["col_header"]
K_all = dot(G_all, G_all.T)
phenotype = data["phenotype"]["matrix"]
pheno_metadata = data["phenotype"]["col_header"]

# Genes from lysine biosynthesis pathway.
lysine_group = [
    "YIL094C",
    "YDL182W",
    "YDL131W",
    "YER052C",
    "YBR115C",
    "YDR158W",
    "YNR050C",
    "YJR139C",
    "YIR034C",
    "YGL202W",
    "YDR234W",
]
window_size = int(5e5)

res = []
for gene in lysine_group:
    # Select the row corresponding to gene of interest on environment 0.0.
    query = "(gene_ID == '{}') & (environment == 0.0)".format(gene)
    df = pheno_metadata.query(query)

    pheno_idx = df.i.item()
    gene_pos = df[["gene_chrom", "gene_end", "gene_start"]]
    # Estimated middle point of the gene.
    midpoint = (gene_pos["gene_end"].item() - gene_pos["gene_start"].item()) / 2

    # Window definition.
    start = midpoint - window_size // 2
    end = midpoint + window_size // 2
    geno = geno_metadata.query("(pos >= {}) & (pos <= {})".format(start, end))

    y = phenotype[:, pheno_idx]
    G_cis = G_all[:, geno.i.values]
    K_cis = dot(G_cis, G_cis.T)
    K_trans = limix.qc.normalise_covariance(K_all - K_cis)
    K_cis = limix.qc.normalise_covariance(K_cis)

    # Definition of the model to fit our data from which we extract
    # the relative signal strength.
    glmm = limix.glmm.GLMMComposer(len(y))
    glmm.y = y
    glmm.fixed_effects.append_offset()
    glmm.covariance_matrices.append(K_cis)
    glmm.covariance_matrices.append(K_trans)
    glmm.covariance_matrices.append_iid_noise()
    glmm.fit(verbose=False, progress=False)

    cis_scale = glmm.covariance_matrices[0].scale
    trans_scale = glmm.covariance_matrices[1].scale
    noise_scale = glmm.covariance_matrices[2].scale

    res.append([cis_scale, trans_scale, noise_scale])

res = DataFrame(res, columns=["cis", "trans", "noise"])
res = res.div(res.sum(axis=1), axis=0).mean(axis=0)
res *= 100

ax = sns.barplot(x=res.index, y=res.values)
ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f%%"))

plt.show()
