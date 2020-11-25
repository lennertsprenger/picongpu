/* Copyright 2020-2021 Sergei Bastrakov
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

#pragma once

#include "picongpu/simulation_defines.hpp"

#include "picongpu/particles/debyeLength/Estimate.kernel"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/Environment.hpp>
#include <pmacc/memory/buffers/HostDeviceBuffer.hpp>
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/mpi/reduceMethods/AllReduce.hpp>
#include <pmacc/traits/GetNumWorkers.hpp>


namespace picongpu
{
    namespace particles
    {
        namespace debyeLength
        {
            /** Estimate Debye length for the given electron species in the local domain
             *
             * @tparam T_ElectronSpecies electron species type
             *
             * @param cellDescription mapping for kernels
             */
            template<typename T_ElectronSpecies>
            HINLINE Estimate estimateLocalDebyeLength(MappingDesc const cellDescription)
            {
                using Frame = typename T_ElectronSpecies::FrameType;
                DataConnector& dc = Environment<>::get().DataConnector();
                auto& electrons = *(dc.get<T_ElectronSpecies>(Frame::getName(), true));

                pmacc::AreaMapping<CORE + BORDER, MappingDesc> mapper(cellDescription);
                constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                    pmacc::math::CT::volume<MappingDesc::SuperCellSize>::type::value>::value;

                auto hostDeviceBuffer = pmacc::HostDeviceBuffer<Estimate, 1>{1u};
                auto hostBox = hostDeviceBuffer.getHostBuffer().getDataBox();
                hostDeviceBuffer.hostToDevice();
                auto kernel = DebyeLengthEstimateKernel<numWorkers>{};
                PMACC_KERNEL(kernel)
                (mapper.getGridDim(), numWorkers)(
                    electrons.getDeviceParticlesBox(),
                    mapper,
                    hostDeviceBuffer.getDeviceBuffer().getDataBox());
                hostDeviceBuffer.deviceToHost();

                // Copy is asynchronous, need to wait for it to finish
                __getTransactionEvent().waitForFinished();
                dc.releaseData(Frame::getName());
                return hostBox(0);
            }

            /** Estimate Debye length for the given electron species in the global domain
             *
             * This function must be called from all MPI ranks.
             * The resulting estimate is a reduction of local estimates from all local domains.
             *
             * @tparam T_ElectronSpecies electron species type
             *
             * @param cellDescription mapping for kernels
             */
            template<typename T_ElectronSpecies>
            HINLINE Estimate estimateGlobalDebyeLength(MappingDesc const cellDescription)
            {
                auto localEstimate = estimateLocalDebyeLength<T_ElectronSpecies>(cellDescription);
                auto globalEstimate = Estimate{};
                pmacc::mpi::MPIReduce reduce;
                reduce(
                    pmacc::nvidia::functors::Add(),
                    &globalEstimate.sumWeighting,
                    &localEstimate.sumWeighting,
                    1,
                    pmacc::mpi::reduceMethods::AllReduce());
                reduce(
                    pmacc::nvidia::functors::Add(),
                    &globalEstimate.sumTemperatureKeV,
                    &localEstimate.sumTemperatureKeV,
                    1,
                    pmacc::mpi::reduceMethods::AllReduce());
                reduce(
                    pmacc::nvidia::functors::Add(),
                    &globalEstimate.sumDebyeLength,
                    &localEstimate.sumDebyeLength,
                    1,
                    pmacc::mpi::reduceMethods::AllReduce());
                return globalEstimate;
            }

        } // namespace debyeLength
    } // namespace particles
} // namespace picongpu
