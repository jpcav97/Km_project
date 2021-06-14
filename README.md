# The following document describes the data and scripts in their order 
# of command as used to analyze the entries provided by the SABIO-RK
# database. They will be explained in order and should be excecuted 
# in the same way.

km_data_corrected.txt
   - This is the corrected version of the original dataset provided 
     by the SABIO-RK database.
   - The data was cleaned by establishing standardized notation for 
     as many relevant details as possible.
   - This script saves the output as 'data_MM_pregroup_pubs_auth.csv'

functions.py
   - This script contains all functions that organize data into groups 
     according to the specified conditions, sort groups by publication
     and by author, randomly takes pairs of values from each group in 
     each dataset, gathers all pairs from groups that have more than 
     one publication, calculates measures of quality and confidence
     intervals, runs linear regression, plots pairs of values, and
     plots the boxplot for our enzyme class analysis. 

Create_data.py
   - This code is used to give information on the percent of data
     missing from each column, and it will extract the author and 
     publication data from the SABIO-RK database for each entry.
   - The code can be modified to extract this data for just 
     Michaelis-Menten entries or the entire dataset.

km_analysis_pairplot.py
   - This script loads in the data saved from 'Create_data.py', and
     calls the functions which creates the three datasets, creates 
     pairs, sorts by publications and author, calculates measures of
     quality and confidence intervals. 
   - This is the code that can be used to produce pairplots, Q-Q plots,
     and find more specific details pertaining to each dataset. 
     
Looping_function.py
   - This function is used to loop the functions that randomly gather
     pairs from each dataset and calculate the measures of error
     according to the semi-log scale specified in the script. 
 
reverse_analysis.py
   - This script will create the five additional datasets that exclude
     one of the important experimental conditions from MRC2 and then 
     loop the calculation of measures of quality similar to the 
     'looping_function.py' on the same semi-log scale. 
