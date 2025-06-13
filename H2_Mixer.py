
"""
Polycrude manufacturing via plastic waste hydrocracking

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 06/13/2025
"""
from thermosteam._graphics import mixer_graphics
import flexsolve as flx
import biosteam as bst
import numpy as np
from typing import Optional
import qsdsan as qs

__all__ = ('Mixer', 'SteamMixer', 'FakeMixer', 'MockMixer')



class Mixer(qs.SanUnit):
    """
    Create a mixer that mixes any number of streams together.
    
    Parameters
    ----------
    ins : 
        Inlet fluids to be mixed.
    outs : 
        Mixed outlet fluid.
    rigorous :
        Whether to perform vapor-liquid equilibrium.
    
    Notes
    -----
    When streams at different pressures are mixed, BioSTEAM assumes valves 
    reduce the pressure of the streams being mixed to prevent backflow 
    (pressure needs to decrease in the direction of flow according to 
    Bernoulli's principle). The outlet pressure will be the minimum pressure
    of all inlet streams.
    
    """
    _graphics = mixer_graphics
    _N_outs = 1
    _N_ins = 2
    _ins_size_is_fixed = False
    
    def _assert_compatible_property_package(self): 
        pass # Not necessary for mixing streams
    
    def _init(self, rigorous: Optional[bool]=False,
              conserve_phases: Optional[bool]=False):
        self.rigorous = rigorous
        self.conserve_phases = conserve_phases
    
    def _run(self):

        makeup_hydrogen, recovered_hydrogen=self.ins
        s_out, = self.outs

        #Calculate the amount of make_up hydrogen

        required_H2_mass=s_out.imass['H2']-recovered_hydrogen.imass['H2']
        makeup_hydrogen.imass['H2']=required_H2_mass

        s_out.mix_from(self.ins, vle=self.rigorous,
                       conserve_phases=getattr(self, 'conserve_phases', None))
        V = s_out.vapor_fraction
        if V == 0:
            self._B = 0
        elif V == 1:
            self._B = np.inf
        else:
            self._B = V / (1 - V)
    
    
    def _get_energy_departure_coefficient(self, stream):
        if stream.phases == ('g', 'l'):
            vapor, liquid = stream
            if vapor.isempty():
                with liquid.temporary_phase('g'): coeff = liquid.H
            else:
                coeff = -vapor.h * liquid.F_mol
        else:
            coeff = -stream.C
        return (self, coeff)
    
    def _create_energy_departure_equations(self):
        # Ll: C1dT1 - Ce2*dT2 - Cr0*dT0 - hv2*L2*dB2 = Q1 - H_out + H_in
        # gl: hV1*L1*dB1 - hv2*L2*dB2 - Ce2*dT2 - Cr0*dT0 = Q1 + H_in - H_out
        outlet = self.outs[0]
        phases = outlet.phases
        if phases == ('g', 'l'):
            vapor, liquid = outlet
            coeff = {}
            if vapor.isempty():
                with liquid.temporary_phase('g'): coeff[self] = liquid.H
            else:
                coeff[self] = vapor.h * liquid.F_mol
        else:
            coeff = {self: outlet.C}
        for i in self.ins: i._update_energy_departure_coefficient(coeff)
        return [(coeff, self.H_in - self.H_out)]
    
    def _create_material_balance_equations(self, composition_sensitive):
        fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
        outlet, = self.outs
        if len(outlet) == 1:
            ones = np.ones(self.chemicals.size)
            minus_ones = -ones
            zeros = np.zeros(self.chemicals.size)
            
            # Overall flows
            eq_overall = {outlet: ones}
            for i in process_inlets: eq_overall[i] = minus_ones
            equations.append(
                (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
            )
        else:
            top, bottom = outlet
            ones = np.ones(self.chemicals.size)
            minus_ones = -ones
            zeros = np.zeros(self.chemicals.size)
            
            # Overall flows
            eq_overall = {}
            for i in outlet: 
                eq_overall[i] = ones
            for i in process_inlets:
                eq_overall[i] = minus_ones
            equations.append(
                (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
            )
            
            # Top to bottom flows
            B = self._B
            eq_outs = {}
            if B == np.inf:
                eq_outs[bottom] = ones
            elif B == 0:
                eq_outs[top] = ones
            else:
                bp = outlet.bubble_point_at_P()
                outlet.T = bp.T
                S = bp.K * B
                eq_outs[top] = ones
                eq_outs[bottom] = -S
            equations.append(
                (eq_outs, zeros)
            )
        return equations
    
    def _update_energy_variable(self, departure):
        phases = self.outs[0].phases
        if phases == ('g', 'l'):
            self._B += departure
        else:
            self.outs[0].T += departure

    def _update_nonlinearities(self): pass