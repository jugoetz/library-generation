"""
Quick analysis upon request from MoBiAS.

Saves a CSV file that contains the synthesis dates with the number of compounds synthesised on that date.
"""
import pandas as pd

from src.db_retrieval.db_queries import MyDatabaseConnection
from src.definitions import UTIL_DIR

con = MyDatabaseConnection()
df = pd.DataFrame(con.get_number_of_experiments_by_date())
df.to_csv(
    UTIL_DIR / "experiment_numbers.csv",
    header=["synthesis date", "number of compounds"],
    index=False,
)
