/* Copyright 2013-2023 Axel Huebl, Heiko Burau, Rene Widera, Remi Lehe,
 *                     Sergei Bastrakov, Lennert Sprenger
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

#include "picongpu/fields/MaxwellSolver/CFLChecker.hpp"
#include "picongpu/fields/MaxwellSolver/DispersionRelation.hpp"
#include "picongpu/fields/MaxwellSolver/CCK/Derivative.hpp"
#include "picongpu/fields/MaxwellSolver/CCK/CCK.def"

#include <pmacc/algorithms/math/defines/pi.hpp>
#include <pmacc/traits/GetStringProperties.hpp>

#include <cstdint>


namespace picongpu
{
    namespace fields
    {
        namespace maxwellSolver
        {
            /** Specialization of the CFL condition checker for CCK solver
             *
             * @tparam T_CherenkovFreeDir the direction (axis) which should be free of cherenkov radiation
             *                            0 = x, 1 = y, 2 = z
             * @tparam T_Defer technical parameter to defer evaluation
             */
            template<uint32_t T_cherenkovFreeDir, typename T_Defer>
            struct CFLChecker<CCK<T_cherenkovFreeDir>, T_Defer>
            {
                /** Check the CFL condition according to the paper, doesn't compile when failed
                 *
                 * @return value of 'X' to fulfill the condition 'c * dt <= X`
                 */
                float_X operator()() const
                {
                    // cellSize is not constexpr currently, so make an own constexpr array
                    constexpr float_X step[3] = {CELL_WIDTH, CELL_HEIGHT, CELL_DEPTH};
                    constexpr auto stepFreeDirection = step[T_cherenkovFreeDir];

                    // Dependance on T_Defer is required, otherwise this check would have been enforced for each setup
                    constexpr auto dt = getTimeStep();
                    PMACC_CASSERT_MSG(
                        Courant_Friedrichs_Lewy_condition_failure____check_your_grid_param_file,
                        (SPEED_OF_LIGHT * dt) <= stepFreeDirection && sizeof(T_Defer*) != 0);

                    // the cherenkov free direction has to be the smallest (or equal) step size
                    PMACC_CASSERT_MSG(
                        CFL_cherenkov_free_direction_smaller_grid_step_or_equal_to_all_others____check_your_grid_param_file,
                        (step[T_cherenkovFreeDir] <= step[(T_cherenkovFreeDir+1)%3]) &&
                        (step[T_cherenkovFreeDir] <= step[(T_cherenkovFreeDir+2)%3])
                    );

                    return stepFreeDirection;
                }
            };

        } // namespace maxwellSolver
    } // namespace fields
} // namespace picongpu

namespace pmacc
{
    namespace traits
    {
        template<uint32_t T_cherenkovFreeDir>
        struct StringProperties<::picongpu::fields::maxwellSolver::CCK<T_cherenkovFreeDir>>
        {
            static StringProperty get()
            {
                auto propList = ::picongpu::fields::maxwellSolver::CCK<T_cherenkovFreeDir>::getStringProperties();
                // overwrite the name of the solver (inherit all other properties)
                propList["name"].value = "CCK";
                return propList;
            }
        };

    } // namespace traits
} // namespace pmacc
