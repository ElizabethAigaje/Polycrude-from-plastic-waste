"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 06/17/2025

"""

# %%
#Import packages

import qsdsan as qs
import biosteam as bst
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

flowsheet = qs.Flowsheet.flowsheet.default 
fs_stream = flowsheet.stream 
fs_unit = flowsheet.unit 


def MSP_breakdown(sys, tea, hydrogen_handling):
    """ Create a class to calculate the cost contribution per processing area to the minimum selling price (MSP).

    Parameters
    ----------
     sys : qsdsan.System
          System to be assessed.
     tea : qsdsan.TEA
          tea of the system assessed
     hydrogen_handling: str
          purchase and storage or on-site production.

    """

    #Calculations needed to calcualte the manufacturing costs in $/kg of polycrude following NREL methoology
    MSP=tea.solve_price(fs_stream.POLYCRUDE)
    polycrude_production=fs_stream.POLYCRUDE.F_mass*tea.operating_days*24   # kg/year
    depreciable_capital=tea.TDC/(tea.duration[1]-tea.duration[0])   
    cashflow_table=tea.get_cashflow_table()
    income_tax= cashflow_table.iloc[:,10]       #rows and columns/ if first rows and columns do not count
    discount_factor=cashflow_table.iloc[:,14]
    NPV_incomeTAX= np.sum(income_tax*discount_factor)*1000000  #in discounted cash flow in MM USD
    polycrude_income=[sys.get_market_value(fs_stream.POLYCRUDE)]*len(income_tax)   #create a list of constant polycrude income $/year, over the years of the project
    NPV_incomePOLYCRUDE=np.sum(polycrude_income*discount_factor)

    #Manufacturing costs of polycrude in $/kg polycrude
    feedstock_costs=sys.get_market_value(fs_stream.feedstock)/polycrude_production  #get_market returns $/year
    VOC_costs=(tea.VOC/polycrude_production)-feedstock_costs
    FOC_costs=tea.FOC/polycrude_production
    capital_depreciation=depreciable_capital/polycrude_production
    average_income_tax=(NPV_incomeTAX/NPV_incomePOLYCRUDE)*MSP
    coproduct_credit=-sys.get_market_value(fs_stream.LightHC)/polycrude_production 
    average_return_investment=MSP-np.sum([feedstock_costs, VOC_costs, FOC_costs, capital_depreciation, average_income_tax, coproduct_credit])

    #Get a table of the MSP disribution
    contribution_values=[MSP,feedstock_costs, VOC_costs, coproduct_credit, FOC_costs, capital_depreciation, average_income_tax, average_return_investment]
    contribution_indexes=['MSP', 'LDPE cost', 'Variable costs', 'Co-product credit', 'Fixed costs', 'capital depreciation', 'average income tax', 'average return in investment']
    table_MSP = pd.DataFrame(contribution_values)
    table_MSP.index=contribution_indexes
    table_MSP.columns=['$/kg polycrude']

    #Create a table to plot where capital depreciation, average income tax and average return in investment are the capital charge
    combined_row=table_MSP.loc['capital depreciation'].values+table_MSP.loc['average income tax']+table_MSP.loc['average return in investment']

    new_table_MSP=table_MSP.drop(['capital depreciation', 'average income tax', 'average return in investment'])
    new_table_MSP.loc['Capital recovery charge']=combined_row

    new_table_MSP['$/bbl polycrude']=new_table_MSP.loc[:,'$/kg polycrude'].multiply(782.59*0.158987)
    new_table_MSP.drop('$/kg polycrude', inplace=True, axis=1)

#Figure
    # Create figures the correct size for publication
    aspect_ratio_LtoW = 1
    cm_to_in = 1/2.54  # centimeters in inches
    width_one_col = 8.3 # cm. Width for a one column figure
    width_two_col = 17.1 # cm. Width for a two column figure
    max_length = 23.3 # cm. The maximum lenght a figure can be

    custom_colors = ['#60c1cf', '#90918e', '#a280b9', '#f98f60', '#79bf82']
    plt.style.use('default')
    font = {'family': 'Arial', 'size': 11}
    plt.rc('font', **font)

    aspect_ratio_LtoW = 1.5  # Length/Width

    fig, ax = plt.subplots(figsize=(width_one_col*cm_to_in, width_one_col*aspect_ratio_LtoW*cm_to_in))
    ax = new_table_MSP.drop('MSP', axis=0).T.plot(kind='bar', stacked=True, color=custom_colors, ax=ax, edgecolor='black', linewidth=0.5) #Tranpose the DatFrame for the stacked bar plot
    MSP_value=new_table_MSP.T.loc['$/bbl polycrude', 'MSP']
    ax.bar('MSP', MSP_value, color='white', label='MSP', edgecolor='black', hatch='.....',linewidth=0.5, zorder=10, width=0.2, align = 'center') #zorder to ensure bar appears above the other bar
  
    ax.set_ylabel('Minimum Selling Price (USD/bbl polycrude)')
    ax.set_xlabel(hydrogen_handling, labelpad=30)  #moves down the label
    ax.legend().remove()
    fig.legend(loc='upper left', bbox_to_anchor=(0.2, 1.01), fontsize=8) # , ncol=len(combined_df.index) ,bbox_to_anchor=(1.25, 0.65)
    # ax.legend(loc='lower right')

    
    ax.spines['bottom'].set_position(('data', 0))   #bottom to ref to the X axis and moves it to y=0.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(y=0, color='black', linewidth=0.8)    #draws the horizontal line at y=0

    
    plt.xticks([])
    # ax.set_xticklabels(ax.get_xticklabels(), ha='left')
    # for label in ax.get_xticklabels():
    #      label.set_y(-25)  # Fine-tune position

    y_min, y_max = ax.get_ylim()   #get the current y-axis range
    ax.set_yticks(np.arange(np.floor(y_min/20)*40, np.ceil(y_max/40)*40+1,40))
    ax.tick_params(axis='y', direction ='in')

    fig.tight_layout()
    fig.subplots_adjust(top=0.75)
    # Show the plot
    figurepath=r"C:\Users\eka5489\OneDrive - The Pennsylvania State University\Research\Second year\Reasearch 2nd year 2024\Python-hydrocracking system\QSDSAN-Hydrocracking\Hydrocracking project-qsd"
    #fig.savefig(os.path.join(figurepath, f'MSP_Contributions.tiff'), dpi=600)

    #Manufacturing costs per section

    individual_units = [] 
    for i in fs_unit:  #units
        for j in i.ins:   #inlets
            inlet_cost = j.price*j.F_mass*24*tea.operating_days
            utility_cost = (i.utility_cost or 0)*24*tea.operating_days
            individual_units.append([i.ID, i.purchase_cost, inlet_cost, utility_cost])
    
    #individual_units.append(tea.systeam_add_OPEX*24*tea.operating_days)

    table_units= pd.DataFrame(individual_units)
    table_units.index=table_units.iloc[:,0] # extract the ID  (all rows and first column) as the first row
    table_units.columns=['ID', 'Equipment cost', 'Materials cost', 'Utility cost' ] 

    #consolidate units that appear multiple times
    table_units_combined = table_units.groupby(table_units.index).agg({'Equipment cost':'mean', 'Materials cost':'sum', 'Utility cost':'mean'})  #aggregate the grups that has the same ID, mean because it wilbe one cost per unit but sum because it sums up all the inlets
    
  
    # Aggregate by processing area to have equipment costs, utilities and materials per section
    if hydrogen_handling == 'purchase and storage':
        index_aggregation = {'C101': 'Size reduction', 'G101': 'Size reduction', 'C102': 'Size reduction', 'ST201': 'Dissolution',
        'M201': 'Dissolution', 'P301': 'Reaction', 'HX301': 'Reaction',
        'R301': 'Reaction', 'F401': 'Separation', 'F402': 'Separation', 'M401': 'Separation', 'P401': 'Separation',
        'D401': 'Separation', 'HX401': 'Separation', 'ST401': 'Separation', 'HX201':'Dissolution', 'HX202':'Dissolution',
        'P201':'Dissolution', 'HX402': 'Separation', 'IC401':'Separation', 'PSA401':'Separation', 'S401': 'Separation',
        'IC402':'Separation', 'HX403': 'Separation', 'HX404': 'Separation','ST402': 'Separation', 
        'ST501': 'Hydrogen handling', 'M501':'Hydrogen hanDling'}

    elif hydrogen_handling == 'on-site production':
        index_aggregation = {'C101': 'Size reduction', 'G101': 'Size reduction', 'C102': 'Size reduction', 'ST201': 'Dissolution',
        'M201': 'Dissolution', 'P301': 'Reaction', 'HX301': 'Reaction',
        'R301': 'Reaction', 'F401': 'Separation', 'F402': 'Separation', 'M401': 'Separation', 'P401': 'Separation',
        'D401': 'Separation', 'HX401': 'Separation', 'ST401': 'Separation', 'HX201':'Dissolution', 'HX202':'Dissolution',
        'P201':'Dissolution', 'HX402': 'Separation', 'IC401':'Separation', 'PSA401':'Separation', 'S401': 'Separation',
        'IC402':'Separation', 'HX403': 'Separation', 'HX404': 'Separation','ST402': 'Separation', 
        'M501':'Hydrogen hanlindg', 'PEM501': 'Hydrogen handling'}

    else:
        raise RuntimeError('In function "MSP_breakdown" argument "hydrogen_handling" must be either "purchase and storage" or "on-site production"')
    
    table_units_combined2=table_units_combined.copy()
    table_units_combined2.index = table_units_combined.index.to_series().replace(index_aggregation)  #modifies the indexes to the section areas
    
    feedstock_costs=table_units_combined.loc['C101', 'Materials cost']   #extract the feedstock cost of the first unit
    table_units_combined2.loc['LDPE', 'Materials cost']=feedstock_costs   # add the feedstock cost to a new row 
    table_units_combined2.loc['Size reduction', 'Materials cost']=0  # delete the material cost from the size reduction

    #consolidate units in the section areas
    results_grouped= table_units_combined2.groupby(table_units_combined2.index ).agg({'Equipment cost':'sum', 'Materials cost':'sum', 'Utility cost':'sum'})
    
    #Table with processing areas and additional OPEX
    MSP_table_contribution=results_grouped.copy()
    MSP_table_contribution=MSP_table_contribution.reindex(['LDPE', 'Size reduction', 'Dissolution', 'Reaction', 'Separation', 'Hydrogen handling']) #make sure the order
    #insert a column for coproduct credit 
    lightHC_credit = [0,0,0,0,sys.get_market_value(fs_stream.LightHC), 0]  #Row position
    MSP_table_contribution.insert(3,'Co-product',lightHC_credit)         #Additional column for co-product
    #insert a new row and new column for additional Opex
    MSP_table_contribution['Additional OPEX'] = 0  #add a new column for additional OPEX
    additional_OPEX_row=pd.DataFrame([[0,0,0,0, tea.system_add_OPEX*24*tea.operating_days]], index=['Other costs'],
                                     columns=['Equipment cost', 'Materials cost', 'Utility cost', 'Co-product', 'Additional OPEX'])
    MSP_table_contribution=pd.concat([MSP_table_contribution, additional_OPEX_row])


    #Table with contributions of the MSP in $/kg polycrude
    polycrude_flow=fs_stream.POLYCRUDE.F_mass*24*tea.operating_days   #kg/year
    
    capital_recovery_charge= table_MSP.loc['capital depreciation'].values + table_MSP.loc['average income tax'].values + table_MSP.loc['average return in investment'].values  #list, loc is to extract specific rows and columns
    MSP_table_contribution['Capital Recovery Charge']=MSP_table_contribution.loc[:,'Equipment cost'].multiply(tea.lang_factor/tea.FCI*capital_recovery_charge[0])  #[row selection, column selection]
    MSP_table_contribution['Fixed costs']=MSP_table_contribution.loc[:,'Equipment cost'].multiply(tea.lang_factor/tea.FCI*table_MSP.loc['Fixed costs'].values[0])   
    MSP_table_contribution['Materials cost']=MSP_table_contribution.loc[:,'Materials cost'].multiply(1/polycrude_flow)
    MSP_table_contribution['Utility cost']=MSP_table_contribution.loc[:,'Utility cost'].multiply(1/polycrude_flow)    
    MSP_table_contribution['Co-product']=MSP_table_contribution.loc[:,'Co-product'].multiply(-1/polycrude_flow)  
    MSP_table_contribution['Additional OPEX']=MSP_table_contribution.loc[:,'Additional OPEX'].multiply(1/polycrude_flow) 
    MSP_table_contribution.drop('Equipment cost', inplace=True, axis=1)  #axis=1 is a column and inplace true modify the actual table and does not create a newdataframe

    #Table with contributions of the MSP in $/bbl polycrude

    new_MSP_table_contribution=MSP_table_contribution.copy()
    new_MSP_table_contribution=new_MSP_table_contribution.multiply(782.59*0.158987)

    #round(MSP_table_contribution.sum().sum(),3) == round(MSP,3)
    #fIGURE
    # Create figures the correct size for publication
    aspect_ratio_LtoW = 1
    cm_to_in = 1/2.54  # centimeters in inches
    width_one_col = 8.3 # cm. Width for a one column figure
    width_two_col = 17.1 # cm. Width for a two column figure
    max_length = 23.3 # cm. The maximum lenght a figure can be

    custom_colors = ['#60c1cf', '#90918e', '#a280b9', '#79694e', '#79bf82', '#f98f60']
    plt.style.use('default')
    font = {'family': 'Arial', 'size': 11}
    plt.rc('font', **font)

    aspect_ratio_LtoW = 1.5  # Length/Width
    fig, ax = plt.subplots(figsize=(width_one_col*cm_to_in, width_one_col*aspect_ratio_LtoW*cm_to_in))

    ax = new_MSP_table_contribution.plot(kind='bar', stacked=True, color=custom_colors, ax=ax, edgecolor='black', linewidth=0.5)

    ax.set_ylabel('Minimum Selling Price (USD/bbl polycrude)')
    ax.set_xlabel('')
    ax.legend().remove()
    fig.legend(loc='upper left', bbox_to_anchor=(0.2, 1.01), fontsize=8)# , ncol=len(combined_df.index) ,bbox_to_anchor=(1.25, 0.65)
    # ax.legend(loc='lower right')

    ax.spines['bottom'].set_position(('data', 0))
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center')
    for label in ax.get_xticklabels():
        label.set_y(-6)  # Fine-tune position
    y_min, y_max = ax.get_ylim()   #get the current y-axis range
    ax.set_yticks(np.arange(np.floor(y_min/20)*20, np.ceil(y_max/20)*20+1,20))
    ax.tick_params(axis='y', direction ='in')

    fig.tight_layout()
    fig.subplots_adjust(top=0.75)
    # Show the plot
    figurepath=r"C:\Users\eka5489\OneDrive - The Pennsylvania State University\Research\Second year\Reasearch 2nd year 2024\Python-hydrocracking system\QSDSAN-Hydrocracking\Hydrocracking project-qsd"
    #fig.savefig(os.path.join(figurepath, f'MSP_Contributions.tiff'), dpi=600)

    return table_MSP, new_table_MSP, MSP_table_contribution, new_MSP_table_contribution
    #return  new_table_MSP, new_MSP_table_contribution  