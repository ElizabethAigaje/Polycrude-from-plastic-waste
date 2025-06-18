"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 06/17/2025

"""
# Import packages
import qsdsan as qs
import biosteam as bst
import numpy as np
import matplotlib.pyplot as plt

#from systems import create_system   #I should work on creating a system for the main hydrocracking 
#from System_hydrocracking import system_TEA
from System_model import *


def sensitivity_analysis(sys, analysis, num_samples, hydrogen_handling):
     """ Create a class to perform the sensitivity analysis and make a tornado plot.

    Parameters
    ----------
     sys : qsdsan.System
          System to be assessed.
     analysis : str
          All, technological or contextual.
     num_samples : float
          Number of samples
     hydrogen_handling: str
          purchase and storage or on-site production.

    """
     
     sensitivity_model=create_model(sys,analysis=analysis,hydrogen_handling=hydrogen_handling)
     
     np.random.seed(3221)   #setting the seed ensures getting the same sample

     samples=sensitivity_model.sample(N=num_samples, rule='L')    #latin hypercube sampling-ensures good coverage of the entire input space 
     sensitivity_model.load_samples(samples)
     sensitivity_model.evaluate()

     sensitivity_model.table = sensitivity_model.table.dropna()
     print(f"Samples after removing NaNs: {len(sensitivity_model.table)} (out of {num_samples})")

     #Generate the spearman correlation data
     r_df, p_df = qs.stats.get_correlations(sensitivity_model, kind='Spearman')  #contains the spearman correlations and p values of the DF 

     #Process the data for plotting: Iterate over each metric to sort and plot separately
     for metric in r_df.columns:   #r_df.columns are the name of the metrics
          sorted_r_df = r_df.sort_values(by=metric, key=lambda x: abs(x), ascending= False)  #allows to modify the data before sorting, in this case
                                                                                             #applying the lambda funtion to get the absolute value.
                                                                                             #with this we keep the negative and positive values for the final plot
          print(f'Top parameters for {metric}')
          print(sorted_r_df)

          #sorted_r_df.index = ['-'.join(map(str, idx)) for idx in sorted_r_df.index] # MultiIndex: We have 2 levels of indexing (the element and the name of the parameter)
                                                                                      # Convert the MultiIndex that is a tuple ('TEA', 'LDPE price') into a string join by _ TEA_LDPE price
          sorted_r_df.index = [idx[1] for idx in sorted_r_df.index]                   # Getting only the second level if indexing
          sorted_r_df.index = [parameter.split(' [')[0] for parameter in sorted_r_df.index]

     #Plotting tornado plot for each metric
     #Formatting
          plt.style.use('default')  #reset the Matplotlib style to its default settings
          font = {'family': 'Arial', 'size': 11}
          plt.rc('font', **font)
          fig, ax= plt.subplots(figsize=(10,8))

     #Plot
          ax.barh(sorted_r_df.index, sorted_r_df[metric],   # does not handle MultiIndex
                   color=['#4682B4' if x>=0 else '#E59B3C' for x in sorted_r_df[metric]],
                   edgecolor='black', linewidth=0.5, height=0.8)   

     #Lables and names
          metric_name= metric[1]
          metric_names= metric_name.split(' [')[0]  #Takes the name before the bracket ['MSP', 'NPV']
          ax.set_xlabel(f'Spearman Correlation for {metric_names}')
          ax.set_ylabel('Parameter')
          ax.invert_yaxis()  #largest correlation on top
          fig.tight_layout
          fig.show()
     return r_df, p_df
