"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 06/13/2025
"""

# Import packages
from logging import raiseExceptions
from chaospy import distributions as shape  # to create distributions
import qsdsan as qs
import numpy as np

__all__ = ('create_model',)    # when imported as from filename import *, it will only import the create_model function

#
def create_model(sys, analysis, hydrogen_handling):   #parameter to call if the analysis is to be done in contextual, technological or all parameters.
    """ Create a class with technological and contextual variables and their probability distributions
        for the uncertainty and sensitivity analyses.

    Parameters
    ----------
    sys : str
        System to be assessed.
    analysis : str
        All, technological or contextual.
    hydrogen_handling: str
        purchase and storage or on-site production.
    
    Returns
    -----
    model : qsdsan.Model
        List of parameters and metrics to be assessed under uncertainty
    
    See also
    -----
    evaluate_contextual_params, evaluate_technological_params

    """
    
    
    model = qs.Model(sys)   # create the model for the uncertainty analysis
    param = model.parameter # to create untertainty parameters
    metric= model.metric    #the metric to perform the uncertainty analysis

    tea = sys.TEA           

    flowsheet = qs.Flowsheet.flowsheet.default
    fs_stream = flowsheet.stream
    fs_unit = flowsheet.unit


# Contextual Parameters

    def evaluate_contextual_parameters():

        #TEA parameters

        baseline = tea.operating_days
        dist=shape.Triangle(lower=np.ceil(baseline*0.9), midpoint=baseline, upper=np.ceil(baseline*1.1))
        @param(name='Operating_days', element='TEA', kind ='isolated', units='days', baseline=baseline, distribution=dist)
        def set_operating_days(i):  #setter function that will upadte the parameter
            tea.operating_days=i

        baseline = tea.income_tax
        dist=shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        @param(name='Income_tax', element='TEA', kind ='isolated', units='-', baseline=baseline, distribution=dist)
        def set_income_tax(i):  
            tea.income_tax=i

        baseline = tea.IRR
        dist=shape.Triangle(lower=0.09, midpoint=baseline, upper=0.15)
        @param(name='Discounte rate', element='TEA', kind ='isolated', units='-', baseline=baseline, distribution=dist)
        def set_interest_rate(i):  
            tea.IRR=i

        baseline = tea.lang_factor
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Lang_factor', element='TEA', kind ='isolated', units='-', baseline=baseline, distribution=dist)
        def set_lang_factor(i):  
            tea.lang_factor=i

        baseline = tea.labor
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Labor', element='TEA', kind ='isolated', units='$/year', baseline=baseline, distribution=dist)
        def set_labor(i):  
            tea.labor=i

        baseline = tea.system_add_OPEX
        dist=shape.Triangle(lower=0.8*baseline, midpoint=baseline, upper=1.2*baseline)  #+-20% since more variables involved here
        @param(name='Additional_OPEX', element='TEA', kind ='isolated', units='$/h', baseline=baseline, distribution=dist)
        def set_system_add_OPEX(i):  
            tea.system_add_OPEX=i
       
        #Material prices

        baseline = fs_stream.feedstock.price
        dist=shape.Triangle(lower=0.15, midpoint=baseline, upper=0.4)
        @param(name='LDPE_price', element='TEA', kind ='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_feedstock_price(i):  
            fs_stream.feedstock.price=i

        baseline = fs_stream.solvent.price
        dist=shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        @param(name='Solvent_price', element='TEA', kind ='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_makeup_solvent_price(i):  
            fs_stream.solvent.price=i

        baseline = fs_stream.HC_catalyst.price
        dist=shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        @param(name='Catalyst_price', element='TEA', kind ='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_HC_catalyst_price(i):  
            fs_stream.HC_catalyst.price=i

        if hydrogen_handling == 'purchase and storage':

            dist=shape.Uniform(lower=1.25, upper=5)   #considers all types of hydrogen
            @param(name='Hydrogen_price', element='TEA', kind ='isolated', units='$/kg', baseline=fs_stream.makeup_hydrogen.price, distribution=dist)
            def set_Makeup_hydrogen_price(i):  
                fs_stream.makeup_hydrogen.price=i

        elif hydrogen_handling == 'on-site production':

            baseline = fs_stream.water_required.price
            dist=shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)   
            @param(name='ProcessWater_price', element='TEA', kind ='isolated', units='$/kg', baseline=baseline, distribution=dist)
            def set_water_required_price(i):  
                fs_stream.water_required.price=i

        else:
            raise RuntimeError('in function "create model" argument "hydrogen_handling" must be either "purchase and storage" or "on-site production"')

        #Utility prices

        baseline = qs.PowerUtility.price
        dist=shape.Triangle(lower=0.07325, midpoint=baseline, upper=0.09185)
        @param(name='PowerUtility_price', element='TEA', kind ='isolated', units='$/kWh', baseline=baseline, distribution=dist)
        def set_PowerUtility_price(i):  
            qs.PowerUtility.price=i
        
        baseline = qs.HeatUtility.get_cooling_agent('cooling_water').regeneration_price
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='CoolingWater_price', element='TEA', kind ='isolated', units='$/kmol', baseline=baseline, distribution=dist)
        def set_CoolingWater_price(i):  
            qs.HeatUtility.get_cooling_agent('cooling_water').regeneration_price=i        

        baseline = qs.HeatUtility.get_cooling_agent('chilled_water').heat_transfer_price
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='ChilledWater_price', element='TEA', kind ='isolated', units='$/kJ', baseline=baseline, distribution=dist)
        def set_ChilledWater_price(i):  
            qs.HeatUtility.get_cooling_agent('chilled_water').heat_transfer_price=i 

        baseline = qs.HeatUtility.get_agent('natural_gas').regeneration_price
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='NaturalGas_price', element='TEA', kind ='isolated', units='$/kmol', baseline=baseline, distribution=dist)
        def set_NaturalGas_price(i):  
            qs.HeatUtility.get_agent('natural_gas').regeneration_price=i 

        baseline = qs.HeatUtility.get_agent('low_pressure_steam').regeneration_price
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='LP_Steam_price', element='TEA', kind ='isolated', units='$/kmol', baseline=baseline, distribution=dist)
        def set_LPsteam_price(i):  
            qs.HeatUtility.get_agent('low_pressure_steam').regeneration_price=i 
        

        #Outlet stream prices

        baseline = fs_stream.POLYCRUDE.price  #make sure in the system is the value for 100$/bbl
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Polycrude_price', element='TEA', kind ='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_Polycrude_price(i):  
            fs_stream.POLYCRUDE.price=i   #PRODUCT PRICE
        
        baseline = fs_stream.LightHC.price
        dist=shape.Triangle(lower=0.3168, midpoint=baseline, upper=0.663)
        @param(name='LightHC_price', element='TEA', kind ='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_LightHC_price(i):  
            fs_stream.LightHC.price=i       


        #Feedstock properties 

        dist=shape.Uniform(lower=0.90,  upper=0.99)    
        @param(name='LDPE_purity', element='C101', kind ='coupled', units='%/100', baseline=fs_unit.C101.purity_factor, distribution=dist)
        def set_Feedstock_purity(i):  
            fs_unit.C101.purity_factor=i   


    def evaluate_technological_parameters():

        #Size reduction

        baseline = fs_unit.G101.power_requirements
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Grinder_PowerInput', element='G101', kind ='coupled', units='kWh/ton', baseline=baseline, distribution=dist)
        def set_Grinder_power(i):  
            fs_unit.G101.power_requirements=i 

        #Dissolution

        dist=shape.Uniform(lower=0.5, upper=4)
        @param(name='Dissolution_time', element='M201', kind ='coupled', units='h', baseline=fs_unit.M201.tau, distribution=dist)  ##
        def set_Dissolution_time(i):  
            fs_unit.M201.tau=i 

        dist=shape.Uniform(lower=0.09, upper=0.33)
        @param(name='LDPE_concentration', element='M201', kind ='coupled', units='kg LDPE/kg solvent', baseline=fs_unit.M201.LDPE_concentration, distribution=dist)  
        def set_LDPE_concentration(i):  
            fs_unit.M201.LDPE_concentration=i 


        #Reaction

        dist=shape.Uniform(lower=2*365*24, upper=4*365*24)
        @param(name='Catalyst_lifetime', element='R301', kind ='coupled', units='h', baseline=fs_unit.R301.catalyst_lifetime, distribution=dist)  
        def set_Catalyst_liftime(i):  
            fs_unit.R301.catalyst_lifetime=i 

        dist=shape.Uniform(lower=2, upper=8.25)
        @param(name='Excess_of_H2', element='R301', kind ='coupled', units='%/100', baseline=fs_unit.R301.hydrogen_excess, distribution=dist)  
        def set_Excess_of_H2(i):  
            fs_unit.R301.hydrogen_excess=i 

        baseline=fs_unit.R301.polycrude_yield
        dist=shape.Triangle(lower=0.9*baseline,midpoint=baseline, upper=1)    #+-10% to consider any small variation in the system lab-errors and scale-up consideration and triangle as this is the base case scenario
        @param(name='Polycrude_yield', element='R301', kind ='coupled', units='mass polycrude/initial mass LDPE', baseline=fs_unit.R301.polycrude_yield, distribution=dist)  
        def set_Polycrude_yield(i):  
            fs_unit.R301.polycrude_yield=i 

        baseline=baseline=fs_unit.R301.tau
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)  #+-10% to consider any small variation in the system lab-errors and scale-up consideration and triangle as this is the base case scenario
        @param(name='Reactor_ResidenceTime', element='R301', kind ='coupled', units='h', baseline=fs_unit.R301.tau, distribution=dist)  #Reaction times based on the project table. Longer times were reported for hydrogenolysis catalyst.
        def set_Residence_time(i):  
            fs_unit.R301.tau=i 
        
        baseline=baseline=fs_unit.R301.WHSV
        dist=shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)  #+-10% to consider any small variation in the system lab-errors and scale-up consideration and triangle as this is the base case scenario
        @param(name='Reactor_WHSV', element='R301', kind ='coupled', units='-', baseline=fs_unit.R301.WHSV, distribution=dist)  #Reaction times based on the project table. Longer times were reported for hydrogenolysis catalyst.
        def set_Reaction_WHSV(i):  
            fs_unit.R301.WHSV=i 


        #Separation

        dist=shape.Uniform(lower=0.70, upper=0.75)
        @param(name='Compressor_401_efficiency', element='IC401', kind ='coupled', units='%/100', baseline=fs_unit.IC401.eta, distribution=dist)  #impact energy requirements
        def set_Compression_efficiency1(i):  
            fs_unit.IC401.eta=i       

        dist=shape.Uniform(lower=0.70, upper=0.75)
        @param(name='Compressor_402_efficiency', element='IC402', kind ='coupled', units='%/100', baseline=fs_unit.IC402.eta, distribution=dist)  
        def set_Compression_efficiency2(i):  
            fs_unit.IC402.eta=i     

        dist=shape.Uniform(lower=0.80, upper=0.90)
        @param(name='PSA_Hydrogen_Recovery', element='PSA401', kind ='coupled', units='%/100', baseline=fs_unit.PSA401.efficiency, distribution=dist)  
        def set_H2_recovery(i):  
            fs_unit.PSA401.efficiency=i  

        dist=shape.Uniform(lower=168, upper=672)   #1 week to 1 month
        @param(name='Polycrude_StorageTime', element='ST401', kind ='isolated', units='h', baseline=fs_unit.ST401.tau, distribution=dist)  #Isolated because only affect the size of the tanks
        def set_Storaget_polycrude(i):  
            fs_unit.ST401.tau=i  

        dist=shape.Uniform(lower=168, upper=672)   #1 week to 1 month
        @param(name='LightHC_StorageTime', element='ST402', kind ='isolated', units='h', baseline=fs_unit.ST402.tau, distribution=dist)  
        def set_Storaget_lightHC(i):  
            fs_unit.ST402.tau=i  

        dist=shape.Uniform(lower=168, upper=672)   #1 week to 1 month
        @param(name='Solvent_StorageTime', element='ST201', kind ='isolated', units='h', baseline=fs_unit.ST201.tau, distribution=dist)  
        def set_Storaget_solvent(i):  
            fs_unit.ST201.tau=i 
        
        # #Hydrogen handling

        if hydrogen_handling == 'purchase and storage':
            dist=shape.Uniform(lower=168, upper=672)   #1 week to 1 month
            @param(name='H2_StorageTime', element='ST501', kind ='isolated', units='h', baseline=fs_unit.ST501.tau, distribution=dist)  
            def set_Storaget_H2(i):  
                fs_unit.ST501.tau=i  

        elif hydrogen_handling == 'on-site production':

            baseline = fs_unit.PEM501.electricity_usage
            dist=shape.Triangle(lower=50, midpoint=baseline, upper=65)   
            @param(name='PEM_ElectricityConsumption', element='PEM501', kind ='coupled', units='kWh/kg H2', baseline=baseline, distribution=dist)
            def set_PEM_electricity_usage(i):  
                fs_unit.PEM501.electricity_usage=i

        else:
            raise RuntimeError('in function "create model" argument "hydrogen_handling" must be either "purchase and storage" or "on-site production"')
        
    if analysis == 'all':
        evaluate_contextual_parameters()
        evaluate_technological_parameters()
    elif analysis == 'technological':
        evaluate_technological_parameters()
    elif analysis == 'contextual':
        evaluate_contextual_parameters()
    else:
        raise RuntimeError(f'In create_model(sys, parameter, hydrogen_handling), parameter={analysis} is not "all", "technological" or "contextual". Please define as one of these.')


    #TEA metrics

    @metric(name='MSP', units='USD/bbl', element='TEA')
    def get_MSP():
        return tea.solve_price(fs_stream.POLYCRUDE)*782.59*0.158987   #conversion of kg to bbl

    @metric(name='NPV', units='USD', element='TEA')
    def get_NPV():
        return tea.NPV
    
    @metric(name='IRR', units='USD', element='TEA')
    def get_IRR():
        return tea.solve_IRR()*100
    
    return model
