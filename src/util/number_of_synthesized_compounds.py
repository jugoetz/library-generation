"""Quick analysis upon request from Bertran

Saves a CSV file that contains the synthesis dates with the number of compounds synthesised on that date.
"""
import datetime

import pandas as pd

from src.db_retrieval.db_queries import MyDatabaseConnection
from src.definitions import UTIL_DIR

con = MyDatabaseConnection()

dates = con.execute_arbitrary_simple_query("SELECT synthesis_date_unixepoch FROM experiments;")
dates = [datetime.datetime.fromtimestamp(d[0]) for d in dates]
df = pd.DataFrame(dates)
df.value_counts().sort_index().to_csv(UTIL_DIR / "experiment_numbers.csv",
                                      header=["number of compounds"],
                                      index_label="synthesis date",
                                      )
