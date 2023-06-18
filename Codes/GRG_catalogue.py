import pandas as pd

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
file = main + 'ELDF/File/en1_RGs_catalogue.csv'

df = pd.read_csv(file, delimiter=',')

df = df[df.lls>=0.7]

df.to_csv(main+'ELDF/File/en1_GRGs_catalogue.csv')



