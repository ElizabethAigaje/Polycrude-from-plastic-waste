"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 06/13/2025
"""

# Import packages
import qsdsan as qs
import biosteam as bst
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from System_model import *


def uncertainty_analysis(sys,  analysis, num_samples, hydrogen_handling):
     """ Create a class to perform the uncertainty analysis and make a kdeplot with the MSP distribution.

    Parameters
    ----------
     sys : str
          System to be assessed.
     analysis : str
          All, technological or contextual.
     num_samples : float
          Number of samples
     hydrogen_handling: str
          purchase and storage or on-site production.

    """
        
     uncertainty_model=create_model(sys,analysis=analysis,hydrogen_handling=hydrogen_handling)
     
     np.random.seed(3221)   #setting the seed ensures getting the same sample

     samples=uncertainty_model.sample(N=num_samples, rule='L')    #most of the time is latin hypercube sampling-ensures good coverage of the entire input space 
     uncertainty_model.load_samples(samples)
     uncertainty_model.evaluate()
     
     #Create the plot of the NPV and MSP

     # Manipulate the data
     TEA_results=uncertainty_model.table.loc[:,'TEA']
     TEA_metrics=['MSP [USD/bbl]', 'NPV [USD]']
     TEA_results=TEA_results[TEA_metrics]   #This is a new dataframe containing only MSP and NPV columns

     #Baseline values for the selected metrics
     baseline_metrics=uncertainty_model.metrics_at_baseline().loc['TEA', TEA_metrics]   #row TEA and columns MSP and NPV from the original table of model_uncertainty

     # Separate the name of the metrics form the units
     metric_names= [metric.split(' [')[0] for metric in TEA_metrics]  #Takes the name before the bracket ['MSP', 'NPV']
     metric_units= [metric.split (' [')[1][:-1] for metric in TEA_metrics]  #Takes the secod part (1) after the bracket excluding the close bracket

     #Creating a list of distinct colors for plotting each metric
     colors = sns.color_palette(palette='husl', n_colors=len(TEA_metrics)).as_hex()    #funciton to create a list of colors of the pallete husl and the number of colors according to the number of metrics, hex_ retrn a pallete with hex code
     
     #Create the figure
     #Formatting
     plt.style.use('default')  #resetst the Matplotlib style to its default settings
     font={'family': 'Arial', 'size': 11}
     plt.rc('font', **font)   #runtime configuration to set properties for the entire figure
     fig, axs = plt.subplots(1,2, figsize= (10,4)) #creates a subplots in 1 row and 2 columns, with 10 in wide and 4 in tall fro the two metrics
     axs = axs.flatten()  # Flatten the axs array to iterate over subplots if not the second plot wont appera
   
     
     class MyScalarFormatter(ScalarFormatter):
          def _set_format(self):
               self.format = '%.2f'  # Forces the X label to display nmbers with 1 decimal places

     #Kdeplot
     for i, (metric, baseline, color) in enumerate(zip(TEA_metrics, baseline_metrics, colors)):    # zip create a tuple of TEA_metrics and colors and enumerate add an index to keep track 

          sns.kdeplot(TEA_results[metric], ax=axs[i], fill=True, color=color, bw_method='scott',linewidth=1.5)  #axs[i] where to plot in the subplot
          #Additional feaures in the plot
          #print('Baseline', baseline_metrics)
          if i==0:
               axs[i].axvline(x=100, color='black', linestyle='--', linewidth=1.5, label='Target MSP')
          else:
               axs[i].axvline(x=0, color='black', linestyle='--', linewidth=1.5, label='Break-even point')
          
          y_midpoint = 0.5 * (axs[i].get_ylim()[0] + axs[i].get_ylim()[1])  #midpoint of y axis
          axs[i].scatter(baseline, y_midpoint, marker='o', label='Baseline', color=color, s=15)

          # Format x-axis label
          axs[i].set_xlabel(f'{metric_names[i]}\n({metric_units[i]})', ha='center', color=color)
        # Set x-axis formatting
          custom_formatter = MyScalarFormatter(useMathText=True)
          axs[i].xaxis.set_major_formatter(custom_formatter)
          axs[i].xaxis.major.formatter.set_powerlimits((0, 0))

          # Remove spines for cleaner look
          axs[i].spines['right'].set_color('none')
          axs[i].spines['top'].set_color('none')

          # Remove y-axis labels and ticks for a cleaner look
          #axs[i].set_yticklabels([])
          axs[i].set_ylabel('')
          axs[i].legend()
     
     fig.tight_layout()
     plt.show()

     # Show the plot
     # fig.savefig(os.path.join(figures_path, f'kde_uncertainty_{fununit}.tiff'), dpi=600)

     #Calculate statistics
     percentiles = [5, 25, 50, 75, 95]
     stats={}
     for metric in TEA_metrics:
          values = TEA_results[metric].replace([np.inf, -np.inf], np.nan).dropna()
          stats[metric]=np.percentile(values, percentiles)
     #Display or store the results
     for metric, values in stats.items():
          print(f'\n{metric}:')
          for p, v in zip(percentiles, values):
               print(f'  {p}th percentile: {v:,.2f}')

     return uncertainty_model, stats
