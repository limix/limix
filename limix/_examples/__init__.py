def load_dataset(dataset, verbose=True):
    if dataset == "boxplot":
        from .boxplot import get_dataframe
        return get_dataframe()
    elif dataset == "kinship":
        from .kinship import get_kinship
        return get_kinship(verbose)
