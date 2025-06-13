"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 03/20/2025
"""

from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year
from qsdsan import SanUnit, Stream
from qsdsan.utils import auom
from qsdsan.sanunits import Reactor, IsothermalCompressor, HXutility
import biosteam as bst

__all__ = ('Hydrocracking_reactor', 'Hydrotreating')

_lb_to_kg = auom('lb').conversion_factor('kg')
_m3perh_to_mmscfd = 1/1177.17 # H2


class Hydrocracking_reactor(Reactor):
    """ Creates a hydrocracking reactor for converting plastic waste into polycrude through hydrocracking reactions
    
    Parameters
    ----------
    ID : str
        ID of the hydrocracking reactor.
    ins : Iterable(stream)
        heavy_oil, hydrogen, catalyst_in.
    outs : Iterable(stream)
        hc_out, catalyst_out.
    WHSV : float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime : float
        HC catalyst lifetime, [hr].
    catalyst_ID : str
        ID of the catalyst.
    hydrogen_P : float
        Hydrogen pressure, [Pa].
    hydrogen_excess : float
        Actual hydrogen amount = hydrogen_rxned*hydrogen_excess.
    HCin_T : float
        HC influent temperature, [K].
    HCrxn_T : float
        HC effluent (after reaction) temperature, [K].
    polycrude_yield : float
        Mass ratio of polycrude to the plastic waste feed.
    tau : float
        Residence time, [h].
    P : float, optional
        Reactor pressure to assess wheter or not is a vaccum pressure vessel, [psig] (the default is None).
    void_fraciton : float, optional
        Fraction of the internals volume over total reactor volume (the default is 0.72 [1]).
    length_to_diamater : float, optional
        Reactor length to diamater ratio (the default is 3, which is suggested for pressure <250 psig [2]).
    diameter : float, optional.
        Reactor diamater, [ft] (the default is None).
    N : int, optional
        Number of reactor (the default is None).
    V : float, optional
        Volume of reactor, [m3] (the default is None).
    auxiliary : bool, optional
        Whether or not the reactor is an auxiliary unit (the default is False).      
    mixing_intensity : float, optional
        Mechanical mixing intensity, [/s] (the default is None).
    kW_per_m3 : float, optional
        Power usage of agitator (the default is 0, which indicates there is no agitation in the reactor).
    wall_thickness_factor : float, optional
        A safety factor to scale up the calculated minimum wall thickness (the default is 1).
    vessel_material : str, optional
        Vessel material (the default is 'Stainless steel 316', which is suggested for hydrogen oeprations).
    vessel_type : str, optional
        Vessel type. Can only be 'Horizontal' or 'Vertical' (the default is 'Vertical').
    
    Notes
    -----
    If ` mixing_intensity ` is provided, the ` kw_per_m3 ` will be calculated based on
    the ` mixing_intensity ` and viscosity of the influent mixture.

    References
    ---------- 
    [1] Dutta, A.; Sahir, A.; Tan, E. Process Design and Economics for the Conversion of Lignocellulosic Biomass 
        to Hydrocarbon Fuels: Thermochemical Research Pathways with In Situ and Ex Situ Upgrading of Fast Pyrolysis Vapors.
    [2] Moss, D. R.; Basic, M. Flange Design. In Pressure Vessel Design Manual; Elsevier, 2013; pp 139–183.
         https://doi.org/10.1016/B978-0-12-387000-1.00003-6.

    """

    _N_ins = 3
    _N_outs = 2
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 include_construction=False,
                 #reaction=None,
                 WHSV=1, # wt./hr per wt. catalyst 
                 catalyst_lifetime=3*7920, 
                 catalyst_ID='HC_catalyst',
                 hydrogen_P=1039.7*6894.76,  
                 hydrogen_excess=5.556,  #times the stoichiometric value
                 HCin_T=394+273.15,       
                 HCrxn_T=451+273.15,
                 polycrude_yield=0.9,   #mass polycrude/mass LDPE
                 P=None, tau=5, void_fraciton=0.72, 
                 length_to_diameter=3, diameter=None,
                 N=None, V=None, auxiliary=False, mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1.0,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, include_construction=include_construction)
        #self.reaction=reaction
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_ID = catalyst_ID
        self.hydrogen_P = hydrogen_P
        self.hydrogen_excess = hydrogen_excess
        self.HCin_T = HCin_T
        self._mixed_in = Stream(f'{ID}_mixed_in')
        self.HCrxn_T = HCrxn_T
        self.polycrude_yield=polycrude_yield

        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        hx_H2_in = Stream(f'{ID}_hx_H2_in')
        hx_H2_out = Stream(f'{ID}_hx_H2_out')
        self.heat_exchanger_H2 = HXutility(ID=f'.{ID}_hx_H2', ins=hx_H2_in, outs=hx_H2_out)
        hx_oil_in = Stream(f'{ID}_hx_oil_in')
        hx_oil_out = Stream(f'{ID}_hx_oil_out')
        self.heat_exchanger_oil = HXutility(ID=f'.{ID}_hx_oil', ins=hx_oil_in, outs=hx_oil_out)
        self.P = P
        self.tau = tau
        self.void_fraciton = void_fraciton
        self.length_to_diameter = length_to_diameter
        self.diameter = diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

        #When I want to assess the polycrude yield of one scenario, we exclude reaciton as an atribute (#reaction and # self.rection)
        #Then, we also activate self.hydrocracking_reaction and def_update function, as well as self.update_reactions()
        #We further change al the self.reaction for self.hydrocracking_reaction

        #Hydrocracking reaction to produce polycrude from plastic waste
        self.hydrocracking_reaction=bst.ParallelReaction([
            #Reaction definition                     #Reactant    #Conversion         
        bst.Reaction('LDPE + H2 -> polycrude', reactant='LDPE',   X=self.polycrude_yield , basis='mol', correct_atomic_balance=True ),
        bst.Reaction('LDPE + H2 ->  C3H8',      reactant='LDPE',   X=(1-self.polycrude_yield)*0.17 , basis='mol', correct_atomic_balance=True),   #working with fixed selectivity
        bst.Reaction('LDPE + H2 ->  C4H10',     reactant='LDPE',   X=(1-self.polycrude_yield)*0.83 , basis='mol', correct_atomic_balance=True)])  

    def update_reactions(self):
        """Update hydrocracking reaction conversion rates based on `polycrude_yield`."""
        self.hydrocracking_reaction[0].X = self.polycrude_yield
        self.hydrocracking_reaction[1].X = (1 - self.polycrude_yield) * 0.17
        self.hydrocracking_reaction[2].X = (1 - self.polycrude_yield) * 0.83

    def _run(self):
        """Execute the reaction and material flow calculations."""
        heavy_oil, hydrogen, catalyst_in = self.ins
        hc_out, catalyst_out = self.outs
        
        #Update the reaction conversion with the current `polycrude_yield`
        self.update_reactions()
        
        catalyst_in.imass[self.catalyst_ID] = heavy_oil.imass['LDPE']/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, 
        # therefore it is not consider heating/cooling of catalyst
        
        #demand of the reactant: stoichiometric coefficient multiplied by conversion
        hydrogen_stoichiometric= 0
        for individual_reaction in self.hydrocracking_reaction:
            hydrogen_demand=individual_reaction.reactant_demand('H2', basis='mol')
            hydrogen_stoichiometric=hydrogen_stoichiometric+hydrogen_demand

        hydrogen_rxned= (hydrogen_stoichiometric*heavy_oil.imol['LDPE']*2)
        hydrogen.imass['H2'] = hydrogen_rxned*self.hydrogen_excess   
        hydrogen.phase = 'g'

        hc_out.phase = 'g'
        
        hc_out.imass['H2'] = hydrogen_rxned*(self.hydrogen_excess - 1)   
        hc_out.mix_from(self.ins)

        self.hydrocracking_reaction(hc_out)
        hc_out.imass['LDPE']=0   #conversion is 1
        hc_out.imass['catalyst']=0   #it goes in a separate stream

        hc_out.P = heavy_oil.P
        hc_out.T = self.HCrxn_T       
        hc_out.vle(T=hc_out.T, P=hc_out.P)

    def _design(self):
        """Perform design calculations for the reactor and auxiliary units."""
        IC = self.compressor
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]

        # Copy hydrogen stream to compressor inlet
        IC_ins0.copy_like(self.ins[1])
        IC_ins0.phase='g'

        #Check if compression is needed
        if IC_ins0.P >= self.hydrogen_P:
            # No compression needed, just pass the stream through
            IC_outs0.copy_like(IC_ins0)
        else:
            # Compression required
            IC.P = self.hydrogen_P
            IC_outs0.P = IC.P= self.hydrogen_P
            IC_outs0.phase = 'g'
            IC.simulate()

        
        hx_H2 = self.heat_exchanger_H2
        hx_H2_ins0, hx_H2_outs0 = hx_H2.ins[0], hx_H2.outs[0]
        hx_H2_ins0.copy_like(self.ins[1])
        hx_H2_outs0.copy_like(hx_H2_ins0)
        hx_H2_ins0.phase = hx_H2_outs0.phase = 'g'
        self._mixed_in.mix_from(self.ins)
        if not self.HCin_T: self.HCin_T = self._mixed_in.T
        hx_H2_outs0.T = self.HCin_T
        hx_H2_ins0.P = hx_H2_outs0.P = IC_outs0.P
        hx_H2.simulate_as_auxiliary_exchanger(ins=hx_H2.ins, outs=hx_H2.outs)

        hx_oil = self.heat_exchanger_oil
        hx_oil_ins0, hx_oil_outs0 = hx_oil.ins[0], hx_oil.outs[0]
        hx_oil_ins0.copy_like(self.ins[0])
        hx_oil_outs0.copy_like(hx_oil_ins0)
        hx_oil_outs0.T = self.HCin_T
        hx_oil_ins0.P = hx_oil_outs0.P = self.ins[0].P
        hx_oil.simulate_as_auxiliary_exchanger(ins=hx_oil.ins, outs=hx_oil.outs)
        
        self.P = min(IC_outs0.P, self.ins[0].P)
        
        V_H2 = self.ins[1].F_vol/self.hydrogen_excess*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)
# %%
