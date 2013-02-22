/*
 *   gridcells_names.h
 *
 *   Dictionary names definitions for the gridcells module.
 *
 *       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
 *       
 *       This program is free software: you can redistribute it and/or modify
 *       it under the terms of the GNU General Public License as published by
 *       the Free Software Foundation, either version 3 of the License, or
 *       (at your option) any later version.
 *       
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *       GNU General Public License for more details.
 *       
 *       You should have received a copy of the GNU General Public License
 *       along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GRIDCELLS_NAMES_H
#define GRIDCELLS_NAMES_H

namespace nest 
{
    namespace names
    {

        /**
         * Add state and parameter names to the nest::names namespace
         */
        const Name E_AMPA("E_AMPA");
        const Name E_NMDA("E_NMDA");
        const Name E_GABA_A("E_GABA_A");
        const Name tau_AMPA_fall("tau_AMPA_fall");
        const Name tau_NMDA_rise("tau_NMDA_rise");
        const Name tau_NMDA_fall("tau_NMDA_fall");
        const Name tau_GABA_A_rise("tau_GABA_A_rise");
        const Name tau_GABA_A_fall("tau_GABA_A_fall");
        const Name E_AHP("E_AHP");
        const Name g_AHP_max("g_AHP_max");
        const Name tau_AHP("tau_AHP");
        const Name g_AHP("g_AHP");
        const Name g_AMPA("g_AMPA");
        const Name g_NMDA("g_NMDA");
        const Name g_NMDA_fraction("g_NMDA_fraction");
        const Name g_GABA_A("g_GABA_A");
        const Name I_stim("I_stim");
        const Name V_clamp("V_clamp");
        const Name I_clamp_AMPA("I_clamp_AMPA");
        const Name I_clamp_NMDA("I_clamp_NMDA");
        const Name I_clamp_GABA_A("I_clamp_GABA_A");
        const Name I_const("I_const");
        const Name I_ac_amp("I_ac_amp");
        const Name I_ac_freq("I_ac_freq");
        const Name I_ac_phase("I_ac_phase");
        const Name I_ac_start_t("I_ac_start_t");
        const Name I_noise_std("I_noise_std");
        const Name I_noise_dt("I_noise_dt");
        const Name rat_pos_x("rat_pos_x");
        const Name rat_pos_y("rat_pos_y");
        const Name rat_pos_dt("rat_pos_dt");
        const Name pref_dir_x("pref_dir_x");
        const Name pref_dir_y("pref_dir_y");
        const Name velC("velC");
        const Name ctr_x("ctr_x");
        const Name ctr_y("ctr_y");
        const Name field_size("field_size");

    } // namespace names

} // namespace nest

#endif // GRIDCELLS_NAMES_H
