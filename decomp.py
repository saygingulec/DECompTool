import pandas as pd
import rpy2.robjects as ro
from math import log2
from optparse import OptionParser
from os.path import isfile
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

pd.options.display.float_format = '{:.6f}'.format  # Prevents pandas from printing with scientific notation.

deseq = importr('DESeq2')
edger = importr('edgeR')

to_dataframe = ro.r('function(x) data.frame(x)')


class Df:

    def __init__(self, files, count_column):
        """
        :param files: File paths. Files in the same group should be written next to each other.
        :param count_column: Designates the column that contains the read counts.
        """
        self.files = files
        self.count_column = count_column  # Name of the column that contains expression counts.
        self.dfs = None  # A list of count data as pandas dfs. Has only the count column.
        self.df = None  # Pandas df formed by merging all data.
        self.form_dfs()
        self.merge_dfs()

    def form_dfs(self):
        dfs = []
        for path in self.files:
            df = pd.read_csv(path, '\t')
            df = df[['Transcript_ID', self.count_column]]
            column_name = path.replace('.txt', '')
            df.rename(columns={self.count_column: column_name}, inplace=True)
            dfs.append(df)
        self.dfs = dfs

    def merge_dfs(self):
        for df in self.dfs[1:]:
            self.dfs[0] = pd.merge(self.dfs[0], df, how='outer', on='Transcript_ID')
        self.df = self.dfs[0]
        self.df.fillna(0.0, inplace=True)
        if self.count_column == 'Tpm':
            sample_columns = self.df.columns.values.tolist()[1:]
            self.df = self.df.loc[(self.df[sample_columns] != 0).any(axis=1)]


class Ttest:
    def __init__(self, df, design):
        """
        :param df: A data frame formed by merging files or a list of files.
        :param design: Number of samples in the first treatment group.
        """
        if type(df) == pd.core.frame.DataFrame:
            self.df = df.copy()
        elif type(df) == list:
            self.df = Df(df, 'Tpm').df
        self.design = design
        self.result = None
        self.t_test()

    def t_test(self):
        group1 = self.df.columns.values.tolist()[1:self.design+1]
        group2 = self.df.columns.values.tolist()[self.design+1:]

        for ind, row in self.df.iterrows():  # to generate p-values
            a = []
            b = []
            for i in group1:
                a.append(row[i])
            for i in group2:
                b.append(row[i])
            if sum(a) != 0 and sum(b) != 0:
                a_mean = sum(a) / len(a)
                b_mean = sum(b) / len(b)
                logfc = log2(b_mean / a_mean)
                self.df.loc[ind, 'log2FC'] = logfc
                self.df.loc[ind, 'p-value'] = ttest_ind(a, b).pvalue
            else:
                self.df.drop(ind, inplace=True)
        self.df = self.df.fillna(1)
        self.df = self.df[['Transcript_ID', 'log2FC', 'p-value']]
        plist = self.df['p-value'].tolist()
        self.df['BH'] = multipletests(plist, method='fdr_bh')[1]  # Benjamini-Hochberg method

        self.df = self.df.set_index('Transcript_ID')
        del self.df.index.name

        self.result = self.df

    def to_csv(self, path, sep='\t'):
        self.df.to_csv(path, sep=sep)


class DESeq2:
    def __init__(self, df, design):
        """
        :param df: A data frame formed by merging files, or a list of files.

        :param design: Number of samples in the first treatment group.

                    treatment
        sampleA1        A
        sampleA2        A
        sampleB1        B
        sampleB2        B
        """
        if type(df) == pd.core.frame.DataFrame:
            self.df = df.copy()
        elif type(df) == list:
            self.df = Df(df, 'Count').df
        self.design = design
        self.design_formula = Formula('~ treatment')
        self.design_matrix = None
        self.dds = None
        self.normalized_count_matrix = None
        self.result = None
        self.design_design_matrix()
        self.run_deseq2()
        self.get_result()

    def design_design_matrix(self):  # Designs the design matrix.
        sample_names = self.df.columns.values.tolist()[1:]
        for i in range(len(sample_names)):
            sample_names[i] = sample_names[i].replace('.txt', '')
        treatment = []
        for i in range(self.design):
            treatment.append('A')
        for i in range(len(sample_names) - self.design):
            treatment.append('B')
        self.design_matrix = pd.DataFrame(treatment, index=sample_names, columns=['treatment'])
        self.design_matrix = pandas2r(self.design_matrix)
        self.df = self.df.set_index('Transcript_ID')
        del self.df.index.name
        self.df = self.df.astype(int)
        self.df = pandas2r(self.df)

    def run_deseq2(self):
        self.dds = deseq.DESeqDataSetFromMatrix(countData=self.df,
                                                colData=self.design_matrix,
                                                design=self.design_formula)
        self.dds = deseq.DESeq(self.dds)
        self.normalized_count_matrix = deseq.counts_DESeqDataSet(self.dds, normalized=True)

    def get_result(self):
        self.result = deseq.results(self.dds)
        self.result = to_dataframe(self.result)
        self.result = r2pandas(self.result)

    def to_csv(self, path, sep='\t'):
        self.result.index.name = 'Transcript_ID'
        self.result.to_csv(path, sep=sep)


class edgeR:

    def __init__(self, df, design):
        if type(df) == pd.core.frame.DataFrame:
            self.df = df.copy()
        elif type(df) == list:
            self.df = Df(df, 'Count').df
        self.design = design
        self.group = None
        self.result = None
        self.form_group()
        self.run_edger()

    def form_group(self):
        sample_names = self.df.columns.values.tolist()[1:]
        group = []
        for i in range(self.design):
            group.append(1)
        for i in range(len(sample_names) - self.design):
            group.append(2)
        self.group = ro.IntVector(group)
        self.df = self.df.set_index('Transcript_ID')
        del self.df.index.name
        self.df = self.df.astype(int)
        self.df = pandas2r(self.df)

    def run_edger(self):
        ro.globalenv['df'] = self.df
        ro.globalenv['group'] = self.group
        ro.r('''
            y <- DGEList(counts=df, group=group)
            y <- calcNormFactors(y)
            design <- model.matrix(~group)
            colnames(design) <- levels(group)
            y <- estimateGLMCommonDisp(y, design)
            y <- estimateGLMTrendedDisp(y, design)
            y <- estimateGLMTagwiseDisp(y, design)
            fit <- glmQLFit(y, design)
            lrt <- glmLRT(fit, coef = 2)
            ''')
        self.result = ro.r('topTags(lrt, n=Inf)')
        self.result = to_dataframe(self.result)
        self.result = r2pandas(self.result)

    def to_csv(self, path, sep='\t'):
        self.result.index.name = 'Transcript_ID'
        self.result.to_csv(path, sep=sep)


def pandas2r(pd_df):  # Converts pandas df to R df
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_from_pd_df = ro.conversion.py2rpy(pd_df)
    return r_from_pd_df


def r2pandas(r_df):  # Converts R df to pandas df
    with localconverter(ro.default_converter + pandas2ri.converter):
        pd_from_r_df = ro.conversion.rpy2py(r_df)
    return pd_from_r_df


parser = OptionParser()

parser.add_option('-f', '--files', dest='files', help="Takes file pathways separated by commas. "
                                                      "Samples of the same group must be given next to each other. "
                                                      "Ex: cancer1.txt,cancer2.txt,normal1.txt,normal2.txt")
parser.add_option('-d', '--design', dest='design', type=int, help="Takes the number of samples in the first"
                                                                  " experimental group. Ex: ")
parser.add_option('-o', '--output', dest='output', help="Takes pathway to output.")
(options, args) = parser.parse_args()

files = options.files
files = files.split(',')
design = options.design
output = options.output


for file in files:
    if not isfile(file):
        print('File {} not found.'.format(file))


print('Performing t-test.')

ttest_res = Ttest(files, design).result
ttest_res.rename(columns={'log2FC': 't_log2FC', 'p-value': 't_pval', 'BH': 't_padj'}, inplace=True)
print('T-test done.')
print('Performing DESeq2.')

deseq_res = DESeq2(files, design).result
deseq_res.drop(['baseMean', 'lfcSE', 'stat'], axis=1, inplace=True)
deseq_res.rename(columns={'log2FoldChange': 'd_log2FC', 'pvalue': 'd_pval', 'padj': 'd_padj'}, inplace=True)
print('DESeq2 done.')
print('Performing edgeR.')

edger_res = edgeR(files, design).result
edger_res.drop(['logCPM', 'LR'], axis=1, inplace=True)
edger_res.rename(columns={'logFC': 'e_log2FC', 'PValue': 'e_pval', 'FDR': 'e_padj'}, inplace=True)
print('EdgeR done.')
print('Forming naming transcripts.')

meta = files[0]
meta = pd.read_csv(meta, sep='\t')
meta = meta[['Transcript_ID', 'Gene_ID', 'Gene_Symbol']]
meta = meta.set_index('Transcript_ID')
print('Finalizing.')

meta = meta.join([ttest_res, deseq_res, edger_res])

meta.to_csv(output, sep='\t', na_rep='NA')

print('Done.')
