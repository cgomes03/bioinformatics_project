import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random


def load_data(file_path: str) -> pd.DataFrame:
    """
    Loads the CSV file from the specified path.
    """
    try:
        df = pd.read_csv(file_path)
        return df
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()


def filter_top_candidates(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """
    For each generation in the dataframe, keep only the top_n rows
    based on the Fitness column (assumed higher is better).
    """
    df["Generation"] = df["Generation"].astype(int)
    filtered = (
        df.sort_values(by="Fitness", ascending=False)
        .groupby("Generation")
        .head(top_n)
        .reset_index(drop=True)
    )
    return filtered


def plot_convergence(df: pd.DataFrame, output_file: str = None, top_n: int = None) -> None:
    """
    Plots convergence curves showing best, worst, median, and average fitness over generations.

    If top_n is provided, the title reflects that only the top candidates are used.
    """
    dfe = df.copy()
    dfe["Generation"] = dfe["Generation"].astype(int)
    grouped = dfe.groupby("Generation")
    best_fitness = grouped["Fitness"].max()
    worst_fitness = grouped["Fitness"].min()
    median_fitness = grouped["Fitness"].median()
    avg_fitness = grouped["Fitness"].mean()

    plt.figure(figsize=(12, 7))
    plt.plot(best_fitness.index, best_fitness.values, label="Best", marker='o')
    plt.plot(worst_fitness.index, worst_fitness.values, label="Worst", marker='v')
    plt.plot(median_fitness.index, median_fitness.values, label="Median", marker='s')
    plt.plot(avg_fitness.index, avg_fitness.values, label="Average", marker='x')
    plt.xticks(ticks=sorted(best_fitness.index))
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    if top_n is not None:
        plt.title(f"Convergence Curves (Top {top_n} Candidates per Generation)")
    else:
        plt.title("Convergence Curves")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file)
        print(f"Convergence plot saved to {output_file}")
        plt.close()
    else:
        plt.show()


def plot_boxplot(df: pd.DataFrame, output_file: str = None, top_n: int = None) -> None:
    """
    Creates a boxplot of fitness values by generation.
    """
    df_box = df.copy()
    df_box["Generation"] = df_box["Generation"].astype(str)
    plt.figure(figsize=(12, 7))
    order = sorted(df_box["Generation"].unique(), key=int)
    sns.boxplot(x="Generation", y="Fitness", data=df_box, order=order)
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    if top_n is not None:
        plt.title(f"Boxplot of Fitness by Generation (Top {top_n} Candidates)")
    else:
        plt.title("Boxplot of Fitness by Generation")
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file)
        print(f"Boxplot saved to {output_file}")
        plt.close()
    else:
        plt.show()


def plot_violin(df: pd.DataFrame, output_file: str = None, top_n: int = None) -> None:
    """
    Creates a violin plot of fitness values by generation.
    """
    df_violin = df.copy()
    df_violin["Generation"] = df_violin["Generation"].astype(str)
    plt.figure(figsize=(12, 7))
    sns.violinplot(x="Generation", y="Fitness", data=df_violin, inner="quartile")
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    if top_n is not None:
        plt.title(f"Violin Plot of Fitness by Generation (Top {top_n} Candidates)")
    else:
        plt.title("Violin Plot of Fitness by Generation")
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file)
        print(f"Violin plot saved to {output_file}")
        plt.close()
    else:
        plt.show()


def plot_scatter(df: pd.DataFrame, output_file: str = None, top_n: int = None) -> None:
    """
    Creates a scatter plot with jitter to show fitness values by generation.
    """
    df_scatter = df.copy()
    df_scatter["Generation"] = df_scatter["Generation"].astype(float)
    jitter = 0.1 * pd.Series([random.uniform(-0.3, 0.3) for _ in range(len(df_scatter))])
    df_scatter["Gen_Jitter"] = df_scatter["Generation"] + jitter
    plt.figure(figsize=(12, 7))
    plt.scatter(df_scatter["Gen_Jitter"], df_scatter["Fitness"], alpha=0.7, edgecolor='k')
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    if top_n is not None:
        plt.title(f"Scatter Plot (with Jitter) of Fitness by Generation (Top {top_n} Candidates)")
    else:
        plt.title("Scatter Plot (with Jitter) of Fitness by Generation")
    plt.xticks(ticks=sorted(df_scatter["Generation"].unique()))
    plt.grid(True)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file)
        print(f"Scatter plot saved to {output_file}")
        plt.close()
    else:
        plt.show()


def plot_final_bar(df: pd.DataFrame, output_file: str = None) -> None:
    """
    Creates a bar chart of the final population's fitness.
    """
    plt.figure(figsize=(12, 7))
    plt.bar(range(len(df)), df["Fitness"])
    plt.xlabel("Individual Index")
    plt.ylabel("Fitness")
    plt.title("Final Population Fitness")
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file)
        print(f"Final bar chart saved to {output_file}")
        plt.close()
    else:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot EA run results from a given output folder containing 'top.csv' and 'final.csv'."
    )
    parser.add_argument("folder", help="Folder containing the run outputs.")
    parser.add_argument("--save", action="store_true", help="Save plots to files instead of displaying interactively.")
    parser.add_argument("--topn", type=int, default=None,
                        help="Number of top candidates per generation to select for plotting (applies to top.csv plots).")
    args = parser.parse_args()

    folder = os.path.abspath(args.folder)

    # Check for top.csv and final.csv in the folder.
    top_csv_path = os.path.join(folder, "top.csv")
    final_csv_path = os.path.join(folder, "final.csv")

    if os.path.isfile(top_csv_path):
        df_top = load_data(top_csv_path)
        if df_top.empty:
            print("No data found in", top_csv_path)
        else:
            if args.topn is not None:
                df_top = filter_top_candidates(df_top, args.topn)
                suffix = f"_top{args.topn}"
                print(f"Using top {args.topn} candidates per generation for plotting.")
            else:
                suffix = ""
            if args.save:
                plot_convergence(df_top, output_file=os.path.join(folder, f"convergence_plot{suffix}.png"),
                                 top_n=args.topn)
                plot_boxplot(df_top, output_file=os.path.join(folder, f"boxplot{suffix}.png"), top_n=args.topn)
                plot_violin(df_top, output_file=os.path.join(folder, f"violin_plot{suffix}.png"), top_n=args.topn)
                plot_scatter(df_top, output_file=os.path.join(folder, f"scatter_plot{suffix}.png"), top_n=args.topn)
            else:
                plot_convergence(df_top, top_n=args.topn)
                plot_boxplot(df_top, top_n=args.topn)
                plot_violin(df_top, top_n=args.topn)
                plot_scatter(df_top, top_n=args.topn)
    else:
        print(f"'top.csv' not found in {folder}")

    if os.path.isfile(final_csv_path):
        df_final = load_data(final_csv_path)
        if df_final.empty:
            print("No data found in", final_csv_path)
        else:
            if args.save:
                plot_final_bar(df_final, output_file=os.path.join(folder, "final_bar_chart.png"))
            else:
                plot_final_bar(df_final)
    else:
        print(f"'final.csv' not found in {folder}")


if __name__ == "__main__":
    main()
