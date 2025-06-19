import pandas as pd

# CSV PATH
file_paths = [
    r"/home/catarinagomes/catarina_stuff/EAs/outputs/RUN_2025-05-10_00-33/final.csv",
    r"/home/catarinagomes/catarina_stuff/EAs/outputs/RUN_2025-05-14_18-56/final.csv",
    r"/home/catarinagomes/catarina_stuff/EAs/outputs/RUN_2025-06-05_15-34/final.csv",
    r"/home/catarinagomes/catarina_stuff/EAs/outputs/RUN_2025-06-09_18-11/final.csv"
]

# COMBINE DataFrames
dataframes = [pd.read_csv(path) for path in file_paths]
df_combined = pd.concat(dataframes, ignore_index=True)

# REMOVE DUPLICATES
df_combined = df_combined.drop_duplicates().sort_values(by="Fitness", ascending=False)
df_combined.to_csv("final.csv", index=False)
print(df_combined)