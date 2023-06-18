import pandas as pd

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
file = main + '/LHLDF/File/lh_GRGs_catalogue.csv'

df = pd.read_csv(file, usecols=('ra', 'dec'))

file = main + '/LHLDF/File/lh_GRGs_radec.csv'
df.to_csv(file, index=False)