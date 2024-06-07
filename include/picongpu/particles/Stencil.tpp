/* Copyright 2023 Tapish Narwal
 *
 * This file is part of PMacc.
 *
 * PMacc is free software: you can redistribute it and/or modify
 * it under the terms of either the GNU General Public License or
 * the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * PMacc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License and the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and the GNU Lesser General Public License along with PMacc.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "picongpu/simulation_defines.hpp"

#include <pmacc/lockstep.hpp>
#include <pmacc/mappings/threads/ThreadCollective.hpp>
#include <pmacc/math/Vector.hpp>
#include <pmacc/math/operation.hpp>
#include <pmacc/memory/boxes/CachedBox.hpp>

namespace picongpu {
    enum Directions
    {
        LEFT = 1u,
        RIGHT = 2u,
        UP = 3u,
        DOWN = 6u
    };

    struct Stencil
    {
        const std::array<uint32_t, 4> stencilDirections{LEFT, RIGHT, UP, DOWN};

        /** run a 4 point stencil for a supercell
        *
        * @tparam T_Box PMacc::DataBox, box type
        * @tparam T_Mapping mapping functor type
        * @param buff databox of the buffer
        * @param mapper functor to map a block to a supercell
        */
        template<typename T_Box, typename T_Mapping, typename T_Worker>
        DINLINE void operator()(
            const T_Worker& worker,
            T_Box field,
            uint32_t currentStep,
            const T_Mapping& mapper) const
        {
            using SuperCellSize = typename T_Mapping::SuperCellSize;
            constexpr uint32_t cellsPerSuperCell = pmacc::math::CT::volume<SuperCellSize>::type::value;

            // get position in grid in units of SuperCells from blockID
            // pmacc::DataSpace<simDim> const block(mapper.getSuperCellIndex(pmacc::DataSpace<simDim>(worker.blockDomIdxND())));
            pmacc::DataSpace<DIM3> const block(
                mapper.getSuperCellIndex(pmacc::DataSpace<DIM3>(pmacc::device::getBlockIdx(worker.getAcc()))));

            // convert position in unit of cells
            pmacc::DataSpace<simDim> const blockCell = block * SuperCellSize::toRT();

            printf("kernel: %u - %u %u %u %u %u %u\n", cellsPerSuperCell, blockCell.x(), blockCell.y(), blockCell.z(), block.x(), block.y(), block.z());

        
            // pmacc::DataSpace<DIM3> const blockCell = block * T_Mapping::SuperCellSize::toRT();

            // const DataSpace<simDim> superCellIdx(mapper.getSuperCellIndex(worker.blockDomIdxND()));

            // pmacc::lockstep::makeForEach<cellsPerSuperCell>(worker)(
            //     [&](int32_t const linearIdx)
            //     {
            //         // cell index within the superCell
            //         pmacc::DataSpace<DIM3> const cellIdxInSupercell = pmacc::math::mapToND(SuperCellSize::toRT(), linearIdx);
            //         pmacc::DataSpace<DIM3> const cell(superCellIdx * SuperCellSize::toRT() + cellIdxInSupercell);
            //         // auto pos = cellIdx + blockCell;
            //         // double value = pmacc::math::sin(currentStep * 0.4);
            //         // pmacc::DataSpace<DIM3> const posLeft = pmacc::DataSpace<DIM3>(8, 8, 8);
            //         field(cell) = {cell.x(), cell.y(), cell.y()};
            //     });

        };
    };
} // namespace picongpu