"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 06/18/2025
"""

# %%
#Import packages

import qsdsan as qs
import biosteam as bst
import numpy as np
import pandas as pd
import os
from logging import raiseExceptions
import matplotlib.pyplot as plt

#import functions from other files
from comoponents_function import*
from TEA_hydrocracking_bst import TEA_hydrocracking
from System_model import*
from Uncertainty_analysis import uncertainty_analysis
from Sensitivity_analysis import sensitivity_analysis
from MSP_breakdown import MSP_breakdown

#import unit operations
from Units.Conveyor import Conveyor
from Units.Conveyor_granulates import Conveyor_granulates
from Units.grinder import Grinder
from Units.DissolutionTank import Dissolution_tank
from Units.Hydrocracking_reactor_modified import Hydrocracking_reactor
from Units.PSA import PSA
from Units.H2_Mixer import Mixer
from Units.H2_electrolysisproduction import Electrolysis

os.environ["PATH"] += os.pathsep + r'C:\Users\eka5489\AppData\Local\anaconda3\envs\REEPS_env\Library\bin'
bst.nbtutorial()

#Collect components
cmps= create_components()

#For the TEA with out reaction as attribute, we exclude # reaction from the def system_TEA
# We deactivate the if function of the scenario analysis and we only leave active the reaction parameters (pressure, temperature, WHSV, time, yield)
# We deactivate the reaction parameters of the hydrocracking
#And in the sys to assess, we do not include the scenario on reaction

def system_TEA(feedPW, hydrogen_handling):  #reaction
    """ Creates the class that models the polycrude manufacturing processs and performs the TEA

    Parameters
    ----------
    feedPW : float
        Mass flow rate of plastic waste, [kg/h].
    hydrogen_handling : str
        purchase and storage, or on-site production

    Returns
    -----
    sys2 : qsdsan.System
        Process simulated, mass flows in and out of the system, process diagram

    tea : qsdsan.TEA
        Net present value, [USD]

    """
    # if reaction == 'base_case':

    #     #Hydrocracking reaction conditions
    reaction_pressure= 11.7 #bar   #Reaction conditions for producing long-HC chains products
    reaction_temperature=300 #C    # Average temperature ranges from 250-300 C
    WHSV_reaction=0.5
    reaction_time=24
    #     #Hydrocracking reaction

    bpolycrude_yield=0.91   #baseline yield and conversion of 100% from Capello et al. Green Chem., 2022, 24, 6306.
    #                         #Selectivity of co-products from Lui et al., Scie Adv. 2021, 7:eabf8283

    #     hydrocracking_reaction=bst.ParallelReaction([
    #         #Reaction definition                     #Reactant    #Conversion         
    #     bst.Reaction('LDPE + H2 -> polycrude', reactant='LDPE',   X=bpolycrude_yield , basis='wt', correct_atomic_balance=True ),
    #     bst.Reaction('LDPE + H2 ->  C3H8',      reactant='LDPE',   X=(1-bpolycrude_yield)*0.17 , basis='wt', correct_atomic_balance=True),
    #     bst.Reaction('LDPE + H2 ->  C4H10',     reactant='LDPE',   X=(1-bpolycrude_yield)*0.83 , basis='wt', correct_atomic_balance=True)])  #calculate the stoichometric indixes

    #     hydrocracking_reaction.show()

    # elif reaction == 'lowreactiontime_case':
        
    #     reaction_pressure = 30 #bar   
    #     reaction_temperature=250 #C 
    #     reaction_time= 2   # h
        
    #     bpolycrude_yield=0.75  # Polycrude flow rate out/LDPE flor rate in
    #     WHSV_reaction= 0.5   

        
    #     hydrocracking_reaction=bst.ParallelReaction([
    #         #Reaction definition                     #Reactant    #Conversion         
    #     bst.Reaction('LDPE + H2 -> polycrude', reactant='LDPE',   X=bpolycrude_yield ,          basis='wt', correct_atomic_balance=True ),
    #     bst.Reaction('LDPE + H2 ->  C3H8',     reactant='LDPE',   X=(1-bpolycrude_yield)*0.2 , basis='wt', correct_atomic_balance=True),   #working with fixed selectivity
    #     bst.Reaction('LDPE + H2 ->  C4H10',    reactant='LDPE',   X=(1-bpolycrude_yield)*0.8 , basis='wt', correct_atomic_balance=True),

    #         ]) 
   
    #     hydrocracking_reaction.show()
   
    # else:
    #     raise RuntimeError(f'In system_TEA(sys, reaction, hydrogen_handling), parameter={reaction} is not "base_case" or "lowreactiontime_case". Please define as one of these.')
    # # Creating INLET STREAMS

    qs.SanStream('feedstock', phase='s', price=0.23) # LDPE with its impurities
    qs.SanStream('HC_catalyst', phase='s', T=25+273, P=101325, price=326.5)    
    qs.SanStream('solvent', phase='l', T=25+273, P=101325, price=5.68)
    qs.SanStream('makeup_hydrogen', phase='g', T =25+273.15, P=1379000, price=2.05)  #Spheric tanks from 30-200 psig in type of tank, in the reactor is compressed to the desired Pressure

    #Inlet Process Water stream for H2 production in electrolysis
    qs.SanStream('water_required', phase='l', T=25+273.15, P=101324, price=0.0000645) #consumption of water is based in the amount of H2 required, water price from Gracida et al. part 1

    # Creating OUTLET STREAMS

    qs.SanStream('POLYCRUDE', price=0.869)   #side product #100/(782.589*0.158981)  #0.7767809569578442(H2 production), 1.107365670221492 (H2 storage)    #price from potential end-users is  100$/barrel tronformed to $/kg
    qs.SanStream('LightHC', price=0.4917)   #side product


    # --------------------------------------------------------------------------
    # Create the System
    # --------------------------------------------------------------------------

    #use qs.Flowsheet serves as a centrilized structure that manages all unit operations,streams and system object in a process system.
    flowsheet = qs.Flowsheet.flowsheet.default 
    fs_stream = flowsheet.stream 
    fs_unit = flowsheet.unit 

    #Process
    #System 1: Reaction, polycrude separation and solvent recycling
    C1=Conveyor('C101', ins=fs_stream.feedstock, outs='ss101', feedPW=feedPW, purity_factor=0.95)
    G1=Grinder('G101', C1-0, outs='ss102', T=298.15, P=101325, power_requirements=300)
    C2=Conveyor_granulates('C102',G1-0, outs='ss103' )
    ST3=qs.sanunits.StorageTank('ST201', ins= fs_stream.solvent, outs='makeup_solvent', tau=336, 
                                vessel_material='Carbon steel', vessel_type='Floating roof') 
    M1=Dissolution_tank('M201', ins=(C2-0, ST3-0, 'solvent_recycle'), outs=('dissolved_LDPE', 'impurities_waste'), feedPW=feedPW, tau=2.25, LDPE_concentration=0.1)
    P1=qs.sanunits._pumping.Pump('P301', ins=M1-0, outs='ss301', P=reaction_pressure*100000, material='Stainless steel')
    HX1=qs.sanunits.HXutility('HX301', ins=P1-0, outs='ss302', T= reaction_temperature+273.15, rigorous=False)
    R1=Hydrocracking_reactor(ID='R301', ins=[HX1-0, 'hydrogen_recycle', fs_stream.HC_catalyst], 
                             outs=['liquid_oil', 'HC_catalyst_out'],
                            #reaction=hydrocracking_reaction,
                            WHSV=WHSV_reaction,   #  pg.252 and Dutta et al. NREL report, the most severe conditions.
                            catalyst_lifetime=3*365*24,  # h. 3 years Hernandez et. al suggest replacement every 3 years
                            catalyst_ID= 'catalyst', # for now is only zeolite (model does not calculate volumen using its density)
                            hydrogen_P=reaction_pressure*100000,  #Pa
                            hydrogen_excess=7.5,         
                            polycrude_yield=bpolycrude_yield,     
                            HCin_T=HX1.T,
                            HCrxn_T=reaction_temperature+273.15,
                            tau=reaction_time,  
                            void_fraciton=0.72,   #33 of reactor internals volume NREL report
                            length_to_diameter=4,
                            wall_thickness_factor=1, 
                            vessel_material='Stainless steel 316', #because of H2
                            vessel_type='Vertical')

    F1=qs.sanunits.Flash('F401', ins=R1.outs[0], outs=(['vapor401', 'liquid401']), T=273.15+220, P=11.7*100000, vessel_material='Stainless steel 316')   #accounting for change in pressure H-LP Flash
    F2=qs.sanunits.Flash('F402', ins=F1.outs[0], outs=(['vapor402', 'liquid402']), T=20+273.15, P=5*100000, vessel_material='Stainless steel 316')   #C-LP flash
    M2=qs.sanunits.Mixer('M401', ins=([F1.outs[1], F2.outs[1]]), outs='ss401')
    P2=qs.sanunits.Pump('P401', ins=M2-0, outs='ss402', P=101324)
    D1=qs.sanunits.ShortcutColumn('D401', ins=P2-0, outs=('distillate','bottoms'), LHK=('C12H26', 'polycrude'), 
                                y_top=0.999, x_bot=0.001, k=3.05, is_divided=True, P=101324)
    HX2=qs.sanunits.HXutility('HX401', ins=D1.outs[1], outs='ss403', T= 15+273.15, rigorous=False)
    ST1=qs.sanunits.StorageTank('ST401', ins=HX2-0, outs= fs_stream.POLYCRUDE, tau=336, 
                                vessel_material='Carbon steel', vessel_type='Floating roof')      #Liquid-residence time in h (2 weeks), material reference in excel
    HX3=qs.sanunits.HXutility('HX201', ins=D1.outs[0], outs='ss201', T= 110+273.15, rigorous=True)
    #PC1=qs.sanunits.PhaseChanger('PC1', ins=HX3-0, outs='Recovered_solvent', phase='l')
    HX4=qs.sanunits.HXutility('HX202', ins=HX3-0, outs='ss202', T= 110+273.15)
    P3=qs.sanunits.Pump('P201', ins=HX4-0, outs=2-M1, P=101324)
    

    #System 2: Light HC recovery and H2 recycling
    HX5=qs.sanunits.HXutility('HX402', ins=F2.outs[0], outs='ss404', T= 50+273.15, rigorous=True, material='Carbon steel/stainles steel')
    IC1=qs.sanunits._compressor.IsothermalCompressor('IC401', ins=HX5-0, outs='ss405', P=15*100000, eta=0.73, material='Stainless steel') 
    PSA1=PSA('PSA401', ins=IC1-0, outs=('light_gases', 'hydrogen_recovered'), efficiency=0.85)
    S1=qs.sanunits.Splitter('S401', ins=PSA1-0, outs=('HC_gases', 'Purge'), split={'C3H8':1, 'C12H26':1, 'H2':0, 'C4H10':1})  #in the desoprtion procees H2 is vented out/purge. Lighter molecule. Asusme only H2 is vented out. Not economically adding another PSA to recover H2 
    IC2=qs.sanunits._compressor.IsothermalCompressor('IC402', ins=S1-0, outs='ss406', P=23*100000, eta=0.73)   #efficiency 70-75% [Seider book]
    HX6=qs.sanunits.HXutility('HX403', ins=IC2-0, outs='ss407', T= 15+273.15, rigorous=True)
    HX7=qs.sanunits.HXutility('HX404', ins=HX6-0, outs='ss408', T= 15+273.15)
    #PC2=qs.sanunits.PhaseChanger('PC2', ins=HX5-0, outs='ss408', phase='l')
    ST2=qs.sanunits.StorageTank('ST402', ins=HX7-0, outs= fs_stream.LightHC, tau=336, vessel_material='Carbon steel', vessel_type='Floating roof')      #Liquid-residence time in h (2 weeks), material reference in excel
    
    # Hydrogen Handling scenarios.

    if hydrogen_handling == 'purchase and storage':
        # "H2 purchase and storage" scenario includes ST4 (storage tank) and assumes H2 is delivered at the pressure required to store it. 
        # Any additional compression to meet reaction pressure is accounted in the hydrocracking reactor design
    
        ST4=qs.sanunits.StorageTank('ST501', ins=fs_stream.makeup_hydrogen, outs='makeup_H2_compressed_gas', tau=336, vessel_material='Stainless steel', vessel_type='Spherical; 30–200 psig')
        M4=Mixer('M501', ins=(ST4-0, PSA1-1), outs=1-R1, rigorous=False) 
        #Initialize the system
        sys1=qs.System('sys1', path=(C1,G1,C2,ST3,M1,P1,HX1,R1,F1,F2,M2,P2,D1,HX2,ST1,HX3, HX4, P3), recycle=P3-0)
        sys2=qs.System('sys2', path=(sys1,HX5, IC1,PSA1,S1,IC2,HX6,HX7,ST2, ST4, M4), recycle=M4-0)


    elif hydrogen_handling == 'on-site production':
        # Accounts for the production via PEM electrolysis were costs including the production and any additional (BOP) cost are included to delyver H2 at 450 psig
        M4=Mixer('M501', ins=(fs_stream.makeup_hydrogen, PSA1-1), outs=1-R1, rigorous=False)    #recovered H2 is from PSA1-1 and out is calculated in the reactor, so the make up is calcualated
        PEM1=Electrolysis('PEM501', ins=fs_stream.water_required, outs=(0-M4, 'oxygen_sideproduct'), electricity_usage=55.8)
        #ST5=qs.sanunits.StorageTank('ST501', ins=PEM1-0, outs=0-M4, tau=24, vessel_material='Stainless steel', vessel_type='Spherical; 30–200 psig')
        sys1=qs.System('sys1', path=(C1,G1,C2,ST3,M1,P1,HX1,R1,F1,F2,M2,P2,D1,HX2,ST1,HX3, HX4, P3), recycle=P3-0)
        sys2=qs.System('sys2', path=(sys1,HX5, IC1,PSA1,S1,IC2,HX6,HX7,ST2, M4, PEM1), recycle=M4-0)       
      
    else:
        raise RuntimeError('in function "system_TEA" argument "hydrogen_handling" must be either "purchase and storage" or "on-site production"')

    sys2.simulate()
    sys2.show()
    sys2.diagram()   

    # ==============================================================================================
    # Calculate TEA
    # ==============================================================================================

    # Labor costs
    # Calculated according to Cappello, V. et. al, Green Chem. 2022, 24 (16), 6306–6318. https://doi.org/10.1039/D2GC01840C.
    P = 1   # Porcesses handling solids 
    N = 13   # Remaining steps
    num_employees = np.ceil((6.29+(3.17*(P)**2)+0.23*N)**0.5)  # Number of operators per shift
    shifts_per_worker=5    #168 h week/40 h week
    operating_hours=2080 #h/year
    pay_rate=52.67   #$/h

    # Calculated according to Seider - 4th edition product and process design principles
    DWandB=num_employees*shifts_per_worker*operating_hours*pay_rate   #direct wages and benefits, $/year
    DSandB_supervisory= 0.15*DWandB   #direct salaries and benefits for supervisory nd engineering, $/year
    OSandSv=0.06*DWandB               #operating supplies and services, $/year
    Tech_assitance=40000              #Technical assistance to manufacturing
    Control_lab=40000                  #Control laboratory

    Lab_related_cost= DWandB+DSandB_supervisory+OSandSv+Tech_assitance+Control_lab
    Labor=1.1*Lab_related_cost

    #Updating/getting prices from qsdsan
    qs.PowerUtility.price=0.081 #$/kwh average value (2022-2023) from Electric Power Monthly - U.S. Energy Information Administration (EIA)

    #Additional OPEX
    #Accounting for the solid disposal cost for impurities in the diossolution units and catalyst in the reactor unit
    solid_disposal_cost=83.41/1000  # $/kg 
    solid_dissolution=M1.outs[1].F_mass   #kg/h
    solid_reactor=R1.outs[1].F_mass       #kg/h
    #operating_days1=0.9*365                   #on_stream factor availability of 90% in NREL
    solid_disposal_OPEX= solid_disposal_cost*(solid_dissolution+solid_reactor)   #$/h

    #Accounting for tranportation of LDPE and Polycrude
    #LDPE
    solid_fixed_cost= 3.01  #$/ton
    solid_variable_cost=0.07  #$/ton/km
    distance_MRF_plant=250 #miles
    LDPE_transportation_OPEX= ((C1.feedPW*(1/C1.purity_factor))/1000)*(solid_fixed_cost+solid_variable_cost*distance_MRF_plant*1.60934 )  #$/h

    #POLYCRUDE
    liquid_fixed_cost= 7.66  #$/ton
    liquid_variable_cost=0.095  #$/ton/km
    distance_plant_refinery=20.5  #miles
    Polycrude_transportation_OPEX= (ST1.outs[0].F_mass/1000)*(liquid_fixed_cost+liquid_variable_cost*distance_plant_refinery**1.60934)   #$/h


    #TEA
    tea=TEA_hydrocracking(system=sys2, 
                        IRR=0.10, #from NREL
                        duration=(2024,2044), #plant life 20 years
                        depreciation='MACRS7', #MACSRS + number of depreciation years
                        income_tax=0.35,  #NREL report
                        operating_days=0.9*365,  #NREL report
                        lang_factor=4.28, #if no Lang factor provided, it calcualtes using bare module. From Seider - 4th edition pg 447 for solids-fluids processing plant
                        construction_schedule= (0.08, 0.6, 0.32), # NREL report: 8% spent in year 2, 60% spents in year 1, 32% spent in year 0
                        startup_months= 6, #NREL 0.5 years
                        startup_FOCfrac=1,  #NREL
                        startup_VOCfrac=0.75,  #NREL
                        startup_salesfrac=0.5,  #NREL
                        WC_over_FCI=0.05, #NREL
                        labor=Labor,
                        system_add_OPEX=solid_disposal_OPEX+LDPE_transportation_OPEX+Polycrude_transportation_OPEX, #to account for solid disposal and transportation, in $/h
                        CEPCI=800.8  #year 2023 to apply to the SanUnits
                        )
    

    return sys2, tea


flowsheet = qs.Flowsheet.flowsheet.default 
fs_stream = flowsheet.stream 
fs_unit = flowsheet.unit 

sys1, tea1,=system_TEA(10416.7, 'purchase and storage')
#sys2, tea2,=system_TEA(10416.7, 'on-site production')
#sys2, tea2,=system_TEA(10416.7, 'base_case', 'on-site production')
tea1.show()
#tea2.show()
#tea.solve_IRR()
#tea1.solve_price(fs_stream.POLYCRUDE)*782.59*0.158987
# sys2.get_market_value(fs_stream.POLYCRUDE)
MSP_breakdown(sys1, tea1, 'purchase and storage')
# print(fs_stream.POLYCRUDE.F_mass)
#sys1.save_report(file='sys2_H2onsite_oc_new10.xlsx')
# print(f"Total TCI: ${tea1.FCI:,.2f}")
# print(f"Total VOC: ${tea1.VOC:,.2f}/year")
# print(f"Total FOC: ${tea1.FOC:,.2f}/year")
# for stream in tea1.system.feeds:
#      print(f"{stream.ID}: ${stream.cost:,.2f} per h")

# print(f"Total Utility Cost: ${tea1.utility_cost:,.2f} per year")

# print(f"Additional OPEX (e.g., solid disposal): ${tea1.system_add_OPEX:,.2f} per hour")

# # # print("\n--- Product Sales Breakdown ---")
# for stream in tea.system.products:
#     print(f"{stream.ID}: ${stream.cost:,.2f} per h")

# print(f"Total capital deprecitation: ${tea1.TDC} ")
# tea.get_cashflow_table()

# print(f"Polycrude price baseline: {fs_stream.POLYCRUDE.price}")
# print(f"Minimum selling price: ${tea.solve_price(fs_stream.POLYCRUDE)} per bbl")  #*782.59*0.158987




# Uncertainty_model

# Uncertainty, stats = uncertainty_analysis(sys1, 'all', 3000, 'on-site production')
# import pandas as pd
# def organize_results(model, path):
#     idx = len(model.parameters)
#     parameters = model.table.iloc[:, :idx]
#     results = model.table.iloc[:, idx:]
#     percentiles = results.quantile([0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
#     with pd.ExcelWriter(path) as writer:
#         parameters.to_excel(writer, sheet_name='Parameters')
#         results.to_excel(writer, sheet_name='Results')
#         percentiles.to_excel(writer, sheet_name='Percentiles')

# organize_results(Uncertainty, 'basecase_purchase_uncertainty_3000samples_improved9011_paper.xlsx')
# r_df, p_df,=sensitivity_analysis(sys1,'all', 3000, 'on-site production')
# fig, ax = qs.stats.plot_correlations(r_df)
# fig
# print(r_df, p_df)

#COMBINED PLOTS
# msp1, contrib1 = MSP_breakdown(sys1, tea1, 'purchase and storage')
# msp2, contrib2 = MSP_breakdown(sys2, tea2, 'on-site production')

# # # Create subplot figure
# fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(11/2.54, 1.5 * 7/2.54))  # two-column width

# custom_colors = ['#60c1cf', '#90918e', '#a280b9', '#79694e', '#79bf82', '#f98f60']

# for ax, contrib, title in zip(axes, [contrib1, contrib2], ['a) H$_2$ purchase & storage', 'b) On-site H$_2$ production']):
#     contrib.plot(kind='bar', stacked=True, color=custom_colors,
#                  ax=ax, edgecolor='black', linewidth=0.5)

#     ax.set_title(title, fontsize=10, fontname='Arial')
#     ax.set_xlabel('')
#     ax.set_ylabel('Minimum Selling Price, $/bbl poly-crude' if title == 'a) H$_2$ purchase & storage' else '', fontsize=10, fontname='Arial')
#     ax.tick_params(axis='x', rotation=90, labelsize=10)
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     ax.spines['bottom'].set_position(('data', 0))
#     ax.axhline(y=0, color='black', linewidth=0.8)

#     # Custom tick placement and formatting
#     y_min, y_max = ax.get_ylim()
#     y_min = min(-20, y_min)
#     y_max = max(120, y_max)
#     ax.set_ylim(y_min, y_max)
#     ax.set_yticks(np.arange(np.floor(y_min/20)*20, np.ceil(y_max/20)*20+1, 20))
#     ax.tick_params(axis='y', direction='in', labelsize=10)

#     # Format x-tick labels
#     ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center')
#     for label in ax.get_xticklabels():
#         label.set_y(-6)  # fine-tune vertical position

#     # Remove subplot legend (we'll use a single shared legend)
#     ax.legend().remove()

# # Get legend handles and labels from the first subplot only
# handles, labels = axes[0].get_legend_handles_labels()

# # Add single shared legend
# fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.95), fontsize=8)

# # Adjust layout and add shared legend above
# fig.tight_layout()
# fig.subplots_adjust(top=0.78)


# # Save figure
# fig.savefig("combined_MSP_plot_smaller.tiff", dpi=1200, bbox_inches='tight')

# %%
