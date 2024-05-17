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
        using Type = typename T_Box::ValueType;
        using SuperCellSize = typename T_Mapping::SuperCellSize;
        using BlockArea = pmacc::SuperCellDescription<
            SuperCellSize,
            typename pmacc::math::CT::make_Int<SuperCellSize::dim, 1>::type,
            typename pmacc::math::CT::make_Int<SuperCellSize::dim, 1>::type>;

        pmacc::DataSpace<DIM3> const block(
            mapper.getSuperCellIndex(pmacc::DataSpace<DIM3>(pmacc::device::getBlockIdx(worker.getAcc()))));
        pmacc::DataSpace<DIM3> const blockCell = block * T_Mapping::SuperCellSize::toRT();

        constexpr uint32_t cellsPerSuperCell = pmacc::math::CT::volume<SuperCellSize>::type::value;

        pmacc::lockstep::makeForEach<cellsPerSuperCell>(worker)(
            [&](int32_t const linearIdx)
            {
        //         // cell index within the superCell
                pmacc::DataSpace<DIM3> const cellIdx = pmacc::math::mapToND(SuperCellSize::toRT(), linearIdx);
                auto pos = cellIdx + blockCell;
                if (pos.x() == 330 && pos.y() == 330) {
                    double value = pmacc::math::sin(currentStep * 0.4);
                    field(pos) = {0, 0, value};
                }

            });
    };
};
