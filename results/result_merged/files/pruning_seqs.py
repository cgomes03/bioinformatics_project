import pandas as pd

# Load the simplification trace
df = pd.read_csv("simplification_trace_all.csv")

# Find best individual 
best_individual = (
    df[df["Step"] == 0].sort_values(by="Fitness", ascending=False).iloc[0]["Individual"]
)

# Extract full trajectory for that solution
pruning_tab = df[df["Individual"] == best_individual][
    ["Genes Count", "Fitness", "Delta", "Removed Gene"]
].reset_index(drop=True)

# Show or export
print(pruning_tab)
pruning_tab.to_csv("pruning_simplification.csv", index=False)
