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
    template<typename T_Area>
    struct ComputeChargeDensity
    {
        using SpeciesType = PIC_Electrons;//pmacc::particles::meta::FindByNameOrType_t<VectorAllSpecies, T_SpeciesType>;
        static const uint32_t area = T_Area::value;

        HINLINE void operator()(FieldTmp* fieldTmp, const uint32_t currentStep) const
        {
            DataConnector& dc = Environment<>::get().DataConnector();

            /* load species without copying the particle data to the host */
            auto speciesTmp = dc.get<SpeciesType>(SpeciesType::FrameType::getName());

            /* run algorithm */
            using ChargeDensitySolver = typename particles::particleToGrid::CreateFieldTmpOperation_t<
                SpeciesType,
                particles::particleToGrid::derivedAttributes::ChargeDensity>::Solver;

            fieldTmp->computeValue<area, ChargeDensitySolver>(*speciesTmp, currentStep);
        }
    };


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
        // auto fieldTmp2 = dc.get<FieldTmp>(FieldTmp::getUniqueId(1));

        auto currentBackground
            = dc.get<simulation::stage::CurrentBackground>(simulation::stage::CurrentBackground::getName());

        auto fieldE = dc.get<FieldE>(FieldE::getName());
        auto fieldJ = dc.get<FieldJ>(FieldJ::getName());
        
        fieldTmp->getGridBuffer().getDeviceBuffer().setValue(FieldTmp::ValueType(1.0));

        // fieldE->getGridBuffer().getDeviceBuffer().setValue({0.0, 1.0, 0.0});

        /* load species without copying the particle data to the host */
        auto speciesTmp = dc.get<PIC_Electrons>(PIC_Electrons::FrameType::getName());

        pmacc::DataSpace<DIM3> guardingSuperCells{0, 0, 0};
        auto mapping = std::make_unique<MappingDesc>(fieldE->getGridLayout().getDataSpace(), guardingSuperCells);

        auto workerCfg = pmacc::lockstep::makeWorkerCfg(typename MappingDesc::SuperCellSize{});
        pmacc::AreaMapping<pmacc::type::CORE, MappingDesc> coreMapper(*mapping);
        
        ComputeChargeDensity<pmacc::mp_int<CORE + BORDER>> computeChargeDensity;
        computeChargeDensity(fieldTmp.get(), currentStep);

        printf("coreMapper.getGridDim() %u %u %u\n", coreMapper.getGridDim().x(), coreMapper.getGridDim().y(), coreMapper.getGridDim().z());
        
        // using BlockArea = SuperCellDescription<
        //     typename MappingDesc::SuperCellSize,
        //     typename pmacc::math::CT::make_Int<simDim, 0>,
        //     typename pmacc::math::CT::make_Int<simDim, 0>>;

        PMACC_LOCKSTEP_KERNEL(Stencil{}, workerCfg)
            (coreMapper.getGridDim())(
                fieldE->getGridBuffer().getDeviceBuffer().getDataBox(),
                currentStep,
                coreMapper);
        
        // MappingDesc cellDescription = fieldTmp->getCellDescription();
        
        // auto const chargeDeviation = [] ALPAKA_FN_ACC(auto const& worker, auto mapper, auto rohBox, auto fieldEBox)
        // {
        //     DataSpace<simDim> const superCellIdx(mapper.getSuperCellIndex(worker.blockDomIdxND()));
        //     DataSpace<simDim> const supercellCellIdx = superCellIdx * SuperCellSize::toRT();
        //     constexpr uint32_t cellsPerSupercell = pmacc::math::CT::volume<SuperCellSize>::type::value;

        //     lockstep::makeForEach<cellsPerSupercell>(worker)(
        //         [&](int32_t const linearIdx)
        //         {
        //             // cell index within the superCell
        //             DataSpace<simDim> const inSupercellCellIdx
        //                 = pmacc::math::mapToND(SuperCellSize::toRT(), linearIdx);

        //             auto globalCellIdx = supercellCellIdx + inSupercellCellIdx;

        //             auto div = picongpu::detail::Div<simDim, typename FieldTmp::ValueType>{};

        //             /* rho := | div E * eps_0 - rho | */
        //             rohBox(globalCellIdx).x()
        //                 = pmacc::math::l2norm(div(fieldEBox.shift(globalCellIdx)) * EPS0 - rohBox(globalCellIdx).x());
        //         });
        // };

        // auto const mapper = makeAreaMapper<CORE + BORDER>(cellDescription);

        // PMACC_LOCKSTEP_KERNEL(charcgeDeviation)
        //     .config(mapper.getGridDim(), SuperCellSize{})(
        //         mapper,
        //         fieldTmp->getGridBuffer().getDeviceBuffer().getDataBox(),
        //         fieldE->getGridBuffer().getDeviceBuffer().getDataBox());



        // picongpu::openPMD::openPMDWriter::GetFields<FieldE> {}(nullptr, 0u);


        // PMACC_LOCKSTEP_KERNEL(setBoundaryConditions, workerCfg)
        //     (borderMapper.getGridDim())(
        //     buff2->getDeviceBuffer().getDataBox(),
        //     NUM_DEVICES_PER_DIM,
        //     gc.getPosition(),
        //     subGrid.getLocalDomain().offset,
        //     gridSize,
        //     borderMapper);

        std::cout << "end PoissonSolverImpl" << std::endl;
    }
} // namespace picongpu::particles