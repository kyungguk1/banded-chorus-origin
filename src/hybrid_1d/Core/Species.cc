/*
 * Copyright (c) 2019-2021, Kyungguk Min
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Species.h"

H1D::Species::Species(ParamSet const &params) : params{params}, geomtr{params}
{
}

auto H1D::Species::operator=(Species const &other) noexcept -> Species &
{
    {
        std::tie(this->moment<0>(), this->moment<1>())
            = std::tie(other.moment<0>(), other.moment<1>());
    }
    return *this;
}
auto H1D::Species::operator=(Species &&other) noexcept -> Species &
{
    {
        std::tie(this->moment<0>(), this->moment<1>())
            = std::forward_as_tuple(std::move(other.moment<0>()), std::move(other.moment<1>()));
    }
    return *this;
}

namespace {
template <class Object> decltype(auto) write_attr(Object &obj, H1D::Species const &sp)
{
    using hdf5::make_type;
    using hdf5::Space;
    {
        obj.attribute("Oc", make_type(sp->Oc), Space::scalar()).write(sp->Oc);
        obj.attribute("op", make_type(sp->op), Space::scalar()).write(sp->op);
        obj.attribute("number_of_source_smoothings", make_type(sp->number_of_source_smoothings),
                      Space::scalar())
            .write(sp->number_of_source_smoothings);

        obj.attribute("charge_density_conversion_factor",
                      make_type(sp.charge_density_conversion_factor()), Space::scalar())
            .write(sp.charge_density_conversion_factor());

        obj.attribute("current_density_conversion_factor",
                      make_type(sp.current_density_conversion_factor()), Space::scalar())
            .write(sp.current_density_conversion_factor());

        obj.attribute("energy_density_conversion_factor",
                      make_type(sp.energy_density_conversion_factor()), Space::scalar())
            .write(sp.energy_density_conversion_factor());
    }
    return obj;
}
} // namespace
auto H1D::operator<<(hdf5::Group &obj, H1D::Species const &sp) -> decltype(obj)
{
    return write_attr(obj, sp);
}
auto H1D::operator<<(hdf5::Dataset &obj, H1D::Species const &sp) -> decltype(obj)
{
    return write_attr(obj, sp);
}
