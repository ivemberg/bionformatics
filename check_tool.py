import pandas
df = pandas.read_excel('Seq_differences/compact/vertical results.xlsx')

#print the column names
print(df.columns)

#get a data frame with selected columns
df_selected = df['ID_SEQ', 'POS', 'LENGTH', 'REFERENCE', 'MUTATION']