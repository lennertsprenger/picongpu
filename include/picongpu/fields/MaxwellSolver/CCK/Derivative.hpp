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

#include "picongpu/fields/MaxwellSolver/CCK/Derivative.def"
#include "picongpu/fields/differentiation/Derivative.hpp"
#include "picongpu/fields/differentiation/ForwardDerivative.hpp"
#include "picongpu/fields/differentiation/Traits.hpp"
#include "picongpu/traits/GetMargin.hpp"

#include <pmacc/algorithms/math/defines/pi.hpp>
#include <pmacc/dimensions/DataSpace.hpp>
#include <pmacc/math/Vector.hpp>
#include <pmacc/meta/accessors/Identity.hpp>
#include <pmacc/types.hpp>


namespace picongpu
{
    namespace fields
    {
        namespace maxwellSolver
        {
            namespace cck
            {
                template<uint32_t T_cherenkovFreeDirection, uint32_t T_direction>
                struct DerivativeFunctor
                {
                private:
                    using InternalDerivativeFunctor
                        = differentiation::DerivativeFunctor<differentiation::Forward, T_direction>;

                public:
                    using LowerMargin = typename pmacc::math::CT::add<
                        typename pmacc::math::CT::make_Int<simDim, 1>::type,
                        typename GetLowerMargin<InternalDerivativeFunctor>::type>::type;

                    using UpperMargin = typename pmacc::math::CT::add<
                        typename pmacc::math::CT::make_Int<simDim, 1>::type,
                        typename GetUpperMargin<InternalDerivativeFunctor>::type>::type;

                    HDINLINE DerivativeFunctor()
                    {
                    }

                    template<typename T_DataBox>
                    HDINLINE typename T_DataBox::ValueType operator()(T_DataBox const& data) const
                    {

                        constexpr uint32_t freeDir = T_cherenkovFreeDirection;
                        constexpr uint32_t dir0 = T_direction;
                        constexpr uint32_t dir1 = (dir0 + 1) % 3;
                        constexpr uint32_t dir2 = (dir0 + 2) % 3;

                        constexpr float_X step[3] = {CELL_WIDTH, CELL_HEIGHT, CELL_DEPTH};
                        constexpr float_X step2[3] = {CELL_WIDTH * CELL_WIDTH,
                                                      CELL_HEIGHT * CELL_HEIGHT,
                                                      CELL_DEPTH * CELL_DEPTH};

                        // delta: smallest step size
                        constexpr float_X delta = step[freeDir];
                        constexpr float_X delta2 = delta*delta;

                        // ri = delta^2 / delta_i^2
                        constexpr float_X ri[3] = {delta2/(step2[dir0]), delta2/(step2[dir1]), delta2/(step2[dir2])};
                        
                        constexpr float_X r012 = ri[0]*ri[1]*ri[2];
                        constexpr float_X rprodsum = ri[0]*ri[1] + ri[1]*ri[2] + ri[2]*ri[0];

                        // equation (21)
                        constexpr float_X betaDir1 = ri[dir1] * 0.125_X * (1 - r012 / rprodsum);
                        constexpr float_X betaDir2 = ri[dir2] * 0.125_X * (1 - r012 / rprodsum);

                        // equation (22)
                        constexpr float_X gammaDir12 = ri[dir1] * ri[dir2] * (0.0625_X - 0.125_X *
                            (ri[dir1]*ri[dir2] / rprodsum));
                        
                        // from the condition alpha + 2 * betaDir1 + 2 * betaDir2 + 4 * gammaDir12 = 1
                        constexpr float_X alpha = 1.0_X - 2.0_X * betaDir1 - 2.0_X * betaDir2 - 4.0_X * gammaDir12;

                        using Index = pmacc::DataSpace<simDim>;
                        auto const upperNeighborDir1 = pmacc::math::basisVector<Index, dir1>();
                        auto const upperNeighborDir2 = pmacc::math::basisVector<Index, dir2>();

                        InternalDerivativeFunctor forwardDerivative
                            = differentiation::makeDerivativeFunctor<differentiation::Forward, T_direction>();

                        return alpha * forwardDerivative(data)
                            + betaDir1 * forwardDerivative(data.shift(upperNeighborDir1))
                            + betaDir1 * forwardDerivative(data.shift(-upperNeighborDir1))

                            + betaDir2 * forwardDerivative(data.shift(upperNeighborDir2))
                            + betaDir2 * forwardDerivative(data.shift(-upperNeighborDir2))
                            
                            + gammaDir12 * forwardDerivative(data.shift(upperNeighborDir1 + upperNeighborDir2))
                            + gammaDir12 * forwardDerivative(data.shift(upperNeighborDir1 - upperNeighborDir2))
                            + gammaDir12 * forwardDerivative(data.shift(-upperNeighborDir1 + upperNeighborDir2))
                            + gammaDir12 * forwardDerivative(data.shift(-upperNeighborDir1 - upperNeighborDir2));
                    }

                };

            } // namespace cck
        } // namespace maxwellSolver

        namespace differentiation
        {
            namespace traits
            {
                template<uint32_t T_cherenkovFreeDirection, uint32_t T_direction>
                struct DerivativeFunctor<maxwellSolver::cck::Derivative<T_cherenkovFreeDirection>, T_direction>
                    : pmacc::meta::accessors::Identity<
                          maxwellSolver::cck::DerivativeFunctor<T_cherenkovFreeDirection, T_direction>>
                {
                };

            } // namespace traits
        } // namespace differentiation
    } // namespace fields
} // namespace picongpu
