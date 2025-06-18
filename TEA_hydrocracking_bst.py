"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 06/17/2025
"""

import biosteam as bst
import qsdsan as qs
import numpy as np
import TEA_bst_modified 


class TEA_hydrocracking(TEA_bst_modified.TEA):
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule, startup_months, 
                 startup_FOCfrac, startup_VOCfrac, startup_salesfrac,
                 WC_over_FCI, labor, system_add_OPEX, CEPCI
                                 ):    #see if I need somethin additional for the FOC 
        
        """ Create a class that performs the techo-economic analysis.

            Parameters
            ----------
            sys : qsdsan.System
                System to be assessed.
            IRR: float
                Discounted rate, [%/100]
            duration: tuple
                start year and final year of the project operation
            depreciation: str
                Depreciation model
            income_tax: float
                Income tax, [%/100]
            operating_days: float
                Number of days tha the plant operates, [day]
            lang_factor: float
                Lang factor
            construction_schedule:tuple
                Fraction of capital investement spent during the years of construction
            startup_months: float
                Time to start operating at 100% after construction is finalize, [months]
            startup_FOCfrac: float
                Fraction of the total FOC during startup time
            startup_VOCfrac: float
                Fraction of the total VOC during startup time
            startup_salesfrac: float
                Fraction of the total sales during startup time
            WC_over_FCI: float
                Work capital over FCI, [%/100]
            labor: float
                Labor costs, [$/year]
            system_add_opex: float
                Additional opex not included in the mass and energy balance in the system, [$/h]
            CEPCI: float
                Chemical engineering plant cost index for the evaluation year
            
            Notes
            -----
            The FOC costs are calculated based on reference [1]

            References
            ----------
            [1] Dutta, A.; Sahir, A.; Tan, E. Process Design and Economics for the Conversion of Lignocellulosic Biomass 
                to Hydrocarbon Fuels: Thermochemical Research Pathways with In Situ and Ex Situ Upgrading of Fast Pyrolysis Vapors
         
           """
        super().__init__(system=system, IRR=IRR, duration=duration, depreciation=depreciation, 
                         income_tax=income_tax,
                         operating_days=operating_days, lang_factor=lang_factor,
                         construction_schedule=construction_schedule,
                         # Assume no startup period
                         startup_months=startup_months, startup_FOCfrac=startup_FOCfrac,
                         startup_VOCfrac=startup_VOCfrac, startup_salesfrac=startup_salesfrac, 
                         WC_over_FCI=WC_over_FCI,
                         # Assume no financing
                         finance_interest=0, finance_years=0, finance_fraction=0
                         # Assume no working capital
                         )   #working capital as a fraciton of fixed capital investement
        
        # Create Attributes
        self.construction_schedule = construction_schedule
        self.labor = labor    
        self.system_add_OPEX=system_add_OPEX
        self.CEPCI=CEPCI
        
    
    # The abstract _DPI method should take installed equipment cost
    # and return the direct permanent investment. 
    def _DPI(self, installed_equipment_cost):
        # Lang factor already applied to installed_equipment_cost in bst_TEA
        installed_equipment_cost=installed_equipment_cost*self.CEPCI/567.5 #CEPCI 2017 in qsdsan for all equipments
        return installed_equipment_cost

    # The abstract _TDC method should take direct permanent investment
    # and return the total depreciable capital. 
    def _TDC(self, DPI):
        return (DPI ) # + contingencies_contractor_fees

    # The abstract _FCI method should take total depreciable capital
    # and return the fixed capital investment. 
    def _FCI(self, TDC):
        return (TDC) # DPI, TDC, and FCI are all the same since the lang factor already includes the other costs of the capital investment # + land + royalties + startup
                    #Lang factor is then applied to the installed_equipment cost

    # The abstract _FOC method should take TDC
    # and return the fixed operating cost.
    def _FOC(self, TDC): 

        # Labor Costs included in TEA_hydrocracking.__init__ from calculation in Systems_hydrocracking.py since this cost is based on the number of operating units.

        # Maintenance Costs
        self.FOC_maintenance = 0.03*TDC  # Total maintenance costs

        # Operating Overhead Costs
        self.FOC_overhead = 0.9*self.labor

        # Property taxes and insurance
        self.FOC_property_taxes_insurance = 0.007*TDC 

        return (self.labor + self.FOC_maintenance + self.FOC_overhead + self.FOC_property_taxes_insurance)