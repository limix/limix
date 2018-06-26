def pytest_sessionstart(session):
    import pandas as pd
    pd.set_option("display.width", 80)
    pd.set_option("display.max_columns", 79)
    pd.set_option("display.max_rows", 60)
    pd.set_option("display.large_repr", "truncate")
