import pandas as pd

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
file = main + '/BLDF/File/bootes_RGs_catalogue.csv'

df = pd.read_csv(file, usecols=('ra', 'dec','las'))

df = df[df.las >= 1]
df = df.drop(columns=('las'))


file = main + '/ELDF/File/bootes_RGs_radec.csv'
df.to_csv(file, index=False)