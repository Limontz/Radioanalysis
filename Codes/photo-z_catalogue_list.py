import numpy as np
import pandas as pd 
from daily_routine import lls
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from astropy import units as u
import sklearn
from sklearn.linear_model import LinearRegression
from scipy.stats import linregress, ks_2samp
from scipy import stats
import statsmodels.api as sm
from scipy.optimize import curve_fit
import pingouin as pg

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
bootesfile = main + '/BLDF/File/bootes_RGs_catalogue.csv'
elaisfile = main + '/ELDF/File/en1_RGs_catalogue.csv'
lhfile = main + '/LHLDF/File/lh_RGs_catalogue.csv'


# Bootes

Bili14 = main + 'BLDF/File/Crossmatch/full_bootes_Bilicki14_crossmatch.csv'
Bili16 = main + 'BLDF/File/Crossmatch/full_bootes_Bilicki16_crossmatch.csv'
Duncan = main + 'BLDF/File/Crossmatch/full_bootes_Duncan_crossmatch.csv'
Beck = main + 'BLDF/File/Crossmatch/full_bootes_PS_crossmatch.csv'
Desi = main + 'BLDF/File/Crossmatch/full_bootes_DESIDR9_crossmatch.csv'
Brescia = main + 'BLDF/File/Crossmatch/full_bootes_Brescia_crossmatch.csv'
Duncan_dr8 = main + 'BLDF/File/Crossmatch/full_bootes_DuncanDR8_crossmatch.csv'
Zhang = main + 'BLDF/File/Crossmatch/full_bootes_Zhang2022_crossmatch.csv'
Zou = main + 'BLDF/File/Crossmatch/full_bootes_Zou2022_crossmatch.csv'


bootes_df = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'ztype', 'Hostname'))
bootes_df['field'] = "Bootes"
print(bootes_df)

bootes_df_bili14 = pd.read_csv(Bili14, delimiter=',')
bootes_df_bili14 = bootes_df_bili14[bootes_df_bili14.zbili14 > 0]
bootes_df_bili16 = pd.read_csv(Bili16, delimiter=',')
bootes_df_bili16 = bootes_df_bili16[bootes_df_bili16.zbili16 > 0]
bootes_df_duncan = pd.read_csv(Duncan, delimiter=',')
bootes_df_duncan = bootes_df_duncan[bootes_df_duncan.zduncan > 0]
bootes_df_beck = pd.read_csv(Beck, delimiter=',')
bootes_df_beck = bootes_df_beck[bootes_df_beck.zbeck > 0]
bootes_df_desi = pd.read_csv(Desi, delimiter=',')
bootes_df_desi = bootes_df_desi[bootes_df_desi.zdesi > 0]
bootes_df_brescia = pd.read_csv(Brescia, delimiter=',')
bootes_df_brescia = bootes_df_brescia[bootes_df_brescia.zbrescia > 0]
bootes_df_dr8 = pd.read_csv(Duncan_dr8, delimiter=',')
bootes_df_dr8 = bootes_df_dr8[bootes_df_dr8.zdr8 > 0]
bootes_df_zhang = pd.read_csv(Zhang, delimiter=',')
bootes_df_zhang= bootes_df_zhang[bootes_df_zhang.zzhang > 0]
bootes_df_zou= pd.read_csv(Zou, delimiter=',')
bootes_df_zou = bootes_df_zou[bootes_df_zou.zzou > 0]

bootes_df = pd.merge(bootes_df,bootes_df_bili14[['name','zbili14']],on='name', how='left')
bootes_df = pd.merge(bootes_df,bootes_df_bili16[['name','zbili16']],on='name', how='left')
bootes_df = pd.merge(bootes_df,bootes_df_duncan[['name','zduncan']],on='name', how='left')
bootes_df = pd.merge(bootes_df,bootes_df_brescia[['name','zbrescia']],on='name', how='left')
bootes_df = pd.merge(bootes_df,bootes_df_beck[['name','zbeck']],on='name', how='left')
bootes_df = pd.merge(bootes_df,bootes_df_desi[['name','zdesi']],on='name', how='left')
bootes_df = pd.merge(bootes_df,bootes_df_dr8[['name','zdr8']],on='name', how='left')
bootes_df = pd.merge(bootes_df,bootes_df_zhang[['name','zzhang']],on='name', how='left')
bootes_df = pd.merge(bootes_df,bootes_df_zou[['name','zzou']],on='name', how='left')



#Elais 
elais_df = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'ztype', 'Hostname'))
elais_df['field'] = "Elais"

Bili14 = main + 'ELDF/File/Crossmatch/full_en1_Bilicki14_crossmatch.csv'
Bili16 = main + 'ELDF/File/Crossmatch/full_en1_Bilicki16_crossmatch.csv'
Duncan = main + 'ELDF/File/Crossmatch/full_en1_Duncan_crossmatch.csv'
Beck = main + 'ELDF/File/Crossmatch/full_en1_PS_crossmatch.csv'
Desi = main + 'ELDF/File/Crossmatch/full_en1_DESIDR9_crossmatch.csv'
Brescia = main + 'ELDF/File/Crossmatch/full_en1_Brescia_crossmatch.csv'
Duncan_dr8 = main + 'ELDF/File/Crossmatch/full_en1_DuncanDR8_crossmatch.csv'
RR13  = main + 'ELDF/File/Crossmatch/full_en1_RR13_crossmatch.csv'
Zhang = main + 'ELDF/File/Crossmatch/full_en1_Zhang2022_crossmatch.csv'
Zou = main + 'ELDF/File/Crossmatch/full_en1_Zou2022_crossmatch.csv'
Subaru = main + 'ELDF/File/Crossmatch/full_en1_Subaru_crossmatch.csv'
Xu = main + 'ELDF/File/Crossmatch/full_en1_xu(2020)_crossmatch.csv'

elais_df_bili14 = pd.read_csv(Bili14, delimiter=',')
elais_df_bili14 = elais_df_bili14[elais_df_bili14.zbili14 > 0]
elais_df_bili16 = pd.read_csv(Bili16, delimiter=',')
elais_df_bili16 = elais_df_bili16[elais_df_bili16.zbili16 > 0]
elais_df_duncan = pd.read_csv(Duncan, delimiter=',')
elais_df_duncan = elais_df_duncan[elais_df_duncan.zduncan > 0]
elais_df_beck = pd.read_csv(Beck, delimiter=',')
elais_df_beck = elais_df_beck[elais_df_beck.zbeck > 0]
elais_df_desi = pd.read_csv(Desi, delimiter=',')
elais_df_desi = elais_df_desi[elais_df_desi.zdesi > 0]
elais_df_brescia = pd.read_csv(Brescia, delimiter=',')
elais_df_brescia = elais_df_brescia[elais_df_brescia.zbrescia > 0]
elais_df_dr8 = pd.read_csv(Duncan_dr8, delimiter=',')
elais_df_dr8 = elais_df_dr8[elais_df_dr8.zdr8 > 0]
elais_df_rr = pd.read_csv(RR13, delimiter=',')
elais_df_rr = elais_df_rr[elais_df_rr.z1_rr > 0]
elais_df_zhang = pd.read_csv(Zhang, delimiter=',')
elais_df_zhang= elais_df_zhang[elais_df_zhang.zzhang > 0]
elais_df_zou= pd.read_csv(Zou, delimiter=',')
elais_df_zou = elais_df_zou[elais_df_zou.zzou > 0]
elais_df_subaru= pd.read_csv(Subaru, delimiter=',')
elais_df_subaru = elais_df_subaru[elais_df_subaru.zsubaru > 0]
elais_df_xu= pd.read_csv(Xu, delimiter=',')
elais_df_xu = elais_df_xu[elais_df_xu.zxu > 0]

elais_df = pd.merge(elais_df,elais_df_bili14[['name','zbili14']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_bili16[['name','zbili16']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_duncan[['name','zduncan']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_brescia[['name','zbrescia']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_beck[['name','zbeck']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_desi[['name','zdesi']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_dr8[['name','zdr8']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_rr[['name','z1_rr']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_zhang[['name','zzhang']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_zou[['name','zzou']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_xu[['name','zxu']],on='name', how='left')
elais_df = pd.merge(elais_df,elais_df_subaru[['name','zsubaru']],on='name', how='left')




start_col = 4
end_col = 12

column_indices = []
for index, row in df.iterrows():
    selected_cols = row.iloc[start_col:end_col+1]
    columns_without_nan = selected_cols.index[selected_cols.notna()]
    indices = [df.columns.get_loc(col) for col in columns_without_nan]
    column_indices.append(indices)
    columns_info = [f"{col_index}: {col_name}" for col_index, col_name in zip(indices, columns_without_nan)]
    print(f"For row {index}, columns without NaN: {', '.join(columns_info)}")

# Add column indices to the DataFrame as a new column
df['Column_Indices'] = column_indices



print(df)