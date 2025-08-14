import pandas as pd
import matplotlib.pyplot as plt

class HandleCSV:
    def __init__(self, filter_list: list[str]):
        self.df_re = {
            'df_re_30m' : self.__filter_gene(pd.read_csv('GSE252357_RE_30m_vs_pre.csv'), filter_list),
            'df_re_3h' : self.__filter_gene(pd.read_csv('GSE252357_RE_3h_vs_pre.csv'), filter_list),
            'df_re_8h' : self.__filter_gene(pd.read_csv('GSE252357_RE_8h_vs_pre.csv'), filter_list),
            'df_re_24h' : self.__filter_gene(pd.read_csv('GSE252357_RE_24h_vs_pre.csv'), filter_list)
        }
        self.df_ctrl = {
            'df_ctrl_30m' : self.__filter_gene(pd.read_csv('GSE252357_CTRL_30m_vs_pre.csv'), filter_list),
            'df_ctrl_3h' : self.__filter_gene(pd.read_csv('GSE252357_CTRL_3h_vs_pre.csv'), filter_list),
            'df_ctrl_8h' : self.__filter_gene(pd.read_csv('GSE252357_CTRL_8h_vs_pre.csv'), filter_list),
            'df_ctrl_24h' : self.__filter_gene(pd.read_csv('GSE252357_CTRL_24h_vs_pre.csv'), filter_list)
        }
    
    def __filter_gene(self, df: pd.DataFrame, filter_list: list[str]) -> pd.DataFrame:
        new_df = df.rename(columns={"Unnamed: 0" : "Gene"})
        new_df = new_df[new_df["Gene"].isin(filter_list)]
        return new_df
    
    def create_re_log2fc_table(self, file_name: str) -> pd.DataFrame:
        new_df = {"Gene" : self.df_re["df_re_30m"]["Gene"]}
        new_df = pd.DataFrame(new_df)

        temp_df = self.df_re["df_re_30m"][["Gene", "log2FoldChange"]]
        temp_df = temp_df.rename(columns={"log2FoldChange" : "Log2FC_30m"})
        new_df = pd.merge(new_df, temp_df, on="Gene")

        temp_df = self.df_re["df_re_3h"][["Gene", "log2FoldChange"]]
        temp_df = temp_df.rename(columns={"log2FoldChange" : "Log2FC_3h"})
        new_df = pd.merge(new_df, temp_df, on="Gene")

        temp_df = self.df_re["df_re_8h"][["Gene", "log2FoldChange"]]
        temp_df = temp_df.rename(columns={"log2FoldChange" : "Log2FC_8h"})
        new_df = pd.merge(new_df, temp_df, on="Gene")

        temp_df = self.df_re["df_re_24h"][["Gene", "log2FoldChange"]]
        temp_df = temp_df.rename(columns={"log2FoldChange" : "Log2FC_24h"})
        new_df = pd.merge(new_df, temp_df, on="Gene")
        
        new_df.to_csv(file_name)
        return new_df

    def create_ctrl_log2fc_table(self, file_name: str) -> pd.DataFrame:
        new_df = {"Gene" : self.df_ctrl["df_ctrl_30m"]["Gene"]}
        new_df = pd.DataFrame(new_df)

        temp_df = self.df_ctrl["df_ctrl_30m"][["Gene", "log2FoldChange"]]
        temp_df = temp_df.rename(columns={"log2FoldChange" : "Log2FC_30m"})
        new_df = pd.merge(new_df, temp_df, on="Gene")

        temp_df = self.df_ctrl["df_ctrl_3h"][["Gene", "log2FoldChange"]]
        temp_df = temp_df.rename(columns={"log2FoldChange" : "Log2FC_3h"})
        new_df = pd.merge(new_df, temp_df, on="Gene")

        temp_df = self.df_ctrl["df_ctrl_8h"][["Gene", "log2FoldChange"]]
        temp_df = temp_df.rename(columns={"log2FoldChange" : "Log2FC_8h"})
        new_df = pd.merge(new_df, temp_df, on="Gene")

        temp_df = self.df_ctrl["df_ctrl_24h"][["Gene", "log2FoldChange"]]
        temp_df = temp_df.rename(columns={"log2FoldChange" : "Log2FC_24h"})
        new_df = pd.merge(new_df, temp_df, on="Gene")
        
        new_df.to_csv(file_name)
        return new_df

def line_chart(gene_list: list[str], df: pd.DataFrame, chart_title: str):
    plt.figure(figsize=(20,10))
    for gene in gene_list:
        log2fc_df = df[df["Gene"] == gene]
        plt.plot(["30m", "3h", "8h", "24h"], [log2fc_df["Log2FC_30m"], log2fc_df["Log2FC_3h"], log2fc_df["Log2FC_8h"], log2fc_df["Log2FC_24h"]], label = gene, marker = "o")
    plt.legend(loc = "upper right")
    plt.title(chart_title)
    plt.xlabel("Time")
    plt.ylabel("Log2FC")
    plt.grid(True)
    plt.show()

def main():
    user_choice = int(input("Select Gene Filter List: "))
    if user_choice == 1:
        gene_filter_list = ['NR4A3', 'FOS', 'JUN', 'EGR1', 'PPARGC1A', 'SLC2A4', 'GLUT4', 'CPT1B', 'MYC', 'MYOG', 'IGF1', 'MSTN', 'COL1A1', 'COL3A1', 'ACTN3']
    elif user_choice == 2:
        gene_filter_list = ['COL1A1', 'COL3A1', 'ACTN3']
    else:
        print("Input 1 or 2")
        return None




    df_controller = HandleCSV(gene_filter_list)
    final_re_df = df_controller.create_re_log2fc_table("final_re_table.csv")
    final_ctrl_df = df_controller.create_ctrl_log2fc_table("final_ctrl_table.csv")
    print(final_re_df["Gene"].unique())
    line_chart(gene_filter_list, final_re_df, "Log2FC of Genes for RE")
    line_chart(gene_filter_list, final_ctrl_df, "Log2FC of Genes for CTRL")


if __name__ == '__main__':
    main()