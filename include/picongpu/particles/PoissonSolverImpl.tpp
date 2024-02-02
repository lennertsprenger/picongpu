/* Copyright 2023 Rene Widera, Lennert Sprenger
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#include "picongpu/simulation_defines.hpp"

#include "picongpu/fields/MaxwellSolver/Solvers.hpp"
#include "picongpu/fields/Fields.hpp"
#include "picongpu/fields/Fields.def"
#include "picongpu/particles/filter/filter.hpp"
#include "picongpu/particles/particleToGrid/derivedAttributes/ChargeDensity.hpp"
#include "picongpu/particles/particleToGrid/CombinedDerive.def"
#include "picongpu/plugins/openPMD/openPMDWriter.hpp"

#include <pmacc/dataManagement/DataConnector.hpp>

#include <cstdint>


namespace picongpu::particles
{
    void PoissonSolverImpl::operator()(uint32_t currentStep, uint32_t iterations)
    {
        std::cout << "start PoissonSolverImpl " << currentStep << " for " << iterations << " iterations."
                  << std::endl;

        DataConnector& dc = Environment<>::get().DataConnector();
        
        // get chargeDensity field
        using ChargeDensitySolver = typename particles::particleToGrid::CreateFieldTmpOperation_t<
                    PIC_Electrons,
                    particles::particleToGrid::derivedAttributes::ChargeDensity>::Solver;

        using ftmpo = picongpu::FieldTmpOperation<ChargeDensitySolver, PIC_Electrons, filter::All>;
        using getfield = picongpu::openPMD::openPMDWriter::GetFields<picongpu::FieldTmpOperation<ChargeDensitySolver, PIC_Electrons, filter::All>>;
        
        getfield{};
        


        // for(uint32_t step = 0u; step < iterations; ++step)
        // {
        //     fieldSolver->update_beforeCurrent(currentStep);
        //     (*currentInterpolationAndAdditionToEMF)(currentStep, *fieldSolver, step == 0);
        //     fieldSolver->update_afterCurrent(currentStep);
        // }

        // simulation::stage::CurrentReset{}(currentStep);
        std::cout << "end PoissonSolverImpl" << std::endl;
    }
} // namespace picongpu::particles