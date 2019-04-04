def assert_finite(Y, M, K):
    from numpy_sugar import is_all_finite

    if not is_all_finite(Y):
        raise ValueError("Outcome must have finite values only.")

    if not is_all_finite(M):
        raise ValueError("Covariates must have finite values only.")

    if K is not None:
        if not is_all_finite(K):
            raise ValueError("Covariate matrix must have finite values only.")
