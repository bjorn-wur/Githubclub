import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(file ="GSE240943.csv"):
    df:pd.DataFrame = pd.read_csv(f"metadata/{file}")
    a = df[["time", "Run", "treatment"]]
    a.to_csv(f"metadata/shortened_{file}")

        # Count the number of runs per treatment per time
    count_df = df.groupby(['time', 'treatment']).size().reset_index(name='counts')

    # Sort time for better visualization
    count_df['time'] = pd.Categorical(count_df['time'], categories=["0h", "1h", "2h", "3h", "4h", "6h", "12h", "24h", "48h", "1D"], ordered=True)
    count_df = count_df.sort_values('time')

    # Set up the plotting style
    sns.set(style="whitegrid")

    # Visualization 1: Bar plot of counts
    plt.figure(figsize=(12, 6))
    sns.barplot(data=count_df, x='time', y='counts', hue='treatment', dodge=True)
    plt.title('Sample Counts per Treatment Across Time', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Sample Counts', fontsize=12)
    plt.xticks(rotation=45)
    plt.legend(title='Treatment', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'sample_counts_barplot_{file}.png', dpi=300)
    return

main("GSE240943.csv")
main("GSE240944.csv")