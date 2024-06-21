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
#include "picongpu/particles/Stencil.tpp"

#include <pmacc/dataManagement/DataConnector.hpp>

#include <cstdint>
#include <iostream>

namespace picongpu::particles
{
    void PoissonSolverImpl::operator()(uint32_t currentStep, uint32_t iterations)
    {
        std::cout << "start PoissonSolverImpl " << currentStep << " for " << iterations << " iterations."
                  << std::endl;
        
        // // get chargeDensity field
        using Species = PIC_Electrons;
        using Solver = typename particles::particleToGrid::CreateFieldTmpOperation_t<
                    PIC_Electrons,
                    particles::particleToGrid::derivedAttributes::ChargeDensity>::Solver;
        using Filter = filter::All;
        
        using UnitType = typename FieldTmp::UnitValueType;
        using ValueType = typename FieldTmp::ValueType;

        DataConnector& dc = Environment<>::get().DataConnector();

        constexpr uint32_t requiredExtraSlots
            = particles::particleToGrid::RequiredExtraSlots<Solver>::type::value;
        PMACC_CASSERT_MSG(
            _please_allocate_at_least_one_or_two_when_using_combined_attributes_FieldTmp_in_memory_param,
            fieldTmpNumSlots >= 1u + requiredExtraSlots);

        auto fieldTmp = dc.get<FieldTmp>(FieldTmp::getUniqueId(0));

        auto fieldE = dc.get<FieldE>(FieldE::getName());
        
        fieldTmp->getGridBuffer().getDeviceBuffer().setValue(FieldTmp::ValueType(1.0));

        /* load species without copying the particle data to the host */
        auto speciesTmp = dc.get<PIC_Electrons>(PIC_Electrons::FrameType::getName());

        MappingDesc cellDescription = fieldE->getCellDescription();
        auto mapper = makeStrideAreaMapper<CORE, 3>(cellDescription);
        
        // ComputeChargeDensity<pmacc::mp_int<CORE + BORDER>> computeChargeDensity;
        // computeChargeDensity(fieldTmp.get(), currentStep);

        // particleToGrid::ComputeGridValuePerFrame<UsedParticleShape, particleToGrid::derivedAttributes::ChargeDensity> cp;

        // particleToGrid::ComputeFieldValue<CORE, Solver, Species, filter::All> computeChargeDensity;
        // computeChargeDensity(fieldTmp.get(), currentStep, 0);

        printf("mapper.getGridDim() %u %u %u\n", mapper.getGridDim().x(), mapper.getGridDim().y(), mapper.getGridDim().z());
 
        // PMACC_LOCKSTEP_KERNEL(Stencil{})
        //     .config(mapper.getGridDim(), SuperCellSize{})(
        //         fieldE->getGridBuffer().getDeviceBuffer().getDataBox(),
        //         currentStep,
        //         mapper);

        std::cout << "end PoissonSolverImpl" << std::endl;
    }
} // namespace picongpu::particles