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

#include "PartSpecies.h"

#include "./BField.h"
#include "./EField.h"

#include <limits>
#include <stdexcept>
#include <utility>

// helpers
//
namespace {
constexpr auto quiet_nan = std::numeric_limits<double>::quiet_NaN();

template <class T, long N> auto &operator/=(P1D::GridQ<T, N> &G, T const w) noexcept
{ // include padding
    for (auto it = G.dead_begin(), end = G.dead_end(); it != end; ++it) {
        *it /= w;
    }
    return G;
}
template <class T, long N> auto &operator+=(P1D::GridQ<T, N> &G, T const w) noexcept
{ // exclude padding
    for (auto it = G.begin(), end = G.end(); it != end; ++it) {
        *it += w;
    }
    return G;
}
//
template <class T, long N> auto const &full_grid(P1D::GridQ<T, N> &F, P1D::BField const &H) noexcept
{
    for (long i = -P1D::Pad; i < F.size() + (P1D::Pad - 1); ++i) {
        (F[i] = H[i + 1] + H[i + 0]) *= 0.5;
    }
    return F;
}
} // namespace

// constructor
//
P1D::PartSpecies::PartSpecies(ParamSet const &params, KineticPlasmaDesc const &desc,
                              std::unique_ptr<VDF> _vdf)
: Species{ params }
, desc{ desc }
, vdf{ std::move(_vdf) }
, bucket{}
, Nc{ desc.Nc == 0 ? 1.0 : desc.Nc }
{
    switch (this->desc.shape_order) {
        case ShapeOrder::_1st:
            _update_velocity = &PartSpecies::_update_velocity_<1>;
            _collect_full_f  = &PartSpecies::_collect_full_f_<1>;
            _collect_delta_f = &PartSpecies::_collect_delta_f_<1>;
            break;
        case ShapeOrder::_2nd:
            _update_velocity = &PartSpecies::_update_velocity_<2>;
            _collect_full_f  = &PartSpecies::_collect_full_f_<2>;
            _collect_delta_f = &PartSpecies::_collect_delta_f_<2>;
            break;
        case ShapeOrder::_3rd:
            _update_velocity = &PartSpecies::_update_velocity_<3>;
            _collect_full_f  = &PartSpecies::_collect_full_f_<3>;
            _collect_delta_f = &PartSpecies::_collect_delta_f_<3>;
            break;
    }
}
void P1D::PartSpecies::populate()
{
    bucket.clear();
    long const Np = desc.Nc * Input::Nx;
    for (long i = 0; i < Np; ++i) {
        SimulationParticle vdf_ptl = vdf->variate(); // position is normalized by Dx
        if (params.domain_extent.is_member(vdf_ptl.pos_x)) {
            // coordinates relative to this subdomain
            vdf_ptl.pos_x -= params.domain_extent.min();

            if constexpr (ParamSet::is_relativistic) {
                // relativistic factor
                Real const gamma = std::sqrt(1 + dot(vdf_ptl.vel, vdf_ptl.vel) / params.c2);
                bucket.emplace_back(vdf_ptl, gamma).w = desc.scheme == full_f;
            } else {
                constexpr Real gamma{ 1 };
                bucket.emplace_back(vdf_ptl, gamma).w = desc.scheme == full_f;
            }
        }
    }
}

void P1D::PartSpecies::load_ptls(std::vector<Particle> const &payload, bool const append)
{
    if (!append)
        bucket.clear();

    for (auto const &ptl : payload) {
        if (params.domain_extent.is_member(ptl.pos_x)) {
            // coordinates relative to this subdomain
            bucket.emplace_back(ptl).pos_x -= params.domain_extent.min();
        }
    }
}
auto P1D::PartSpecies::dump_ptls() const -> std::vector<Particle>
{
    decltype(dump_ptls()) payload{ begin(bucket), end(bucket) };
    for (auto &ptl : payload) {
        // coordinates relative to whole domain
        ptl.pos_x += params.domain_extent.min();
    }
    return payload; // NRVO
}

// update & collect interface
//
void P1D::PartSpecies::update_vel(BField const &bfield, EField const &efield, Real const dt)
{
    (*this->_update_velocity)(bucket, full_grid(moment<1>(), bfield), efield,
                              BorisPush{ dt, params.c, params.O0, desc.Oc });
}
void P1D::PartSpecies::update_pos(Real const dt, Real const fraction_of_grid_size_allowed_to_travel)
{
    Real const dtODx = dt / params.Dx; // normalize position by grid size
    if (!_update_x(bucket, dtODx, 1.0 / fraction_of_grid_size_allowed_to_travel)) {
        throw std::domain_error{ std::string{ __PRETTY_FUNCTION__ }
                                 + " - particle(s) moved too far" };
    }
}
void P1D::PartSpecies::collect_part()
{
    switch (desc.scheme) {
        case full_f:
            (this->*_collect_full_f)(moment<1>());
            break;
        case delta_f:
            (this->*_collect_delta_f)(moment<1>(), bucket);
            break;
    }
}
void P1D::PartSpecies::collect_all()
{
    _collect(moment<0>(), moment<1>(), moment<2>());
}

// heavy lifting
//
bool P1D::PartSpecies::_update_x(bucket_type &bucket, Real const dtODx,
                                 Real const travel_distance_scale_factor)
{
    bool did_not_move_too_far = true;
    for (auto &ptl : bucket) {
        Real moved_x = ptl.g_vel().x * dtODx;
        if constexpr (ParamSet::is_relativistic) {
            moved_x /= ptl.gamma;
        }
        ptl.pos_x += moved_x; // position is normalized by grid size

        // travel distance check
        //
        moved_x *= travel_distance_scale_factor;
        did_not_move_too_far &= 0 == long(moved_x);
    }
    return did_not_move_too_far;
}

template <long Order>
void P1D::PartSpecies::_update_velocity_(bucket_type &bucket, VectorGrid const &B, EField const &E,
                                         BorisPush const pusher)
{
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    Shape<Order> sx;
    for (auto &ptl : bucket) {
        sx(ptl.pos_x); // position is normalized by grid size
        if constexpr (ParamSet::is_relativistic)
            ptl.gamma = pusher.relativistic(ptl.g_vel(), B.interp(sx), E.interp(sx));
        else
            pusher.non_relativistic(ptl.g_vel(), B.interp(sx), E.interp(sx));
    }
}

template <long Order> void P1D::PartSpecies::_collect_full_f_(VectorGrid &nV) const
{
    nV.fill(Vector{ 0 });
    //
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    Shape<Order> sx;
    for (auto const &ptl : bucket) {
        sx(ptl.pos_x); // position is normalized by grid size

        if constexpr (ParamSet::is_relativistic)
            nV.deposit(sx, ptl.vel());
        else
            nV.deposit(sx, ptl.g_vel());
    }
    //
    Real const Nc = this->Nc;
    nV /= Vector{ Nc };
}
template <long Order>
void P1D::PartSpecies::_collect_delta_f_(VectorGrid &nV, bucket_type &bucket) const
{
    VDF const &vdf = *this->vdf;
    //
    nV.fill(Vector{ 0 });
    //
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    Shape<Order> sx;
    for (auto &ptl : bucket) {
        sx(ptl.pos_x); // position is normalized by grid size
        ptl.w = vdf.weight(ptl.base());

        if constexpr (ParamSet::is_relativistic)
            nV.deposit(sx, ptl.vel() * ptl.w);
        else
            nV.deposit(sx, ptl.g_vel() * ptl.w);
    }
    //
    Real const Nc = this->Nc;
    (nV /= Vector{ Nc }) += vdf.nV0(quiet_nan) * desc.scheme;
}
void P1D::PartSpecies::_collect(ScalarGrid &n, VectorGrid &nV, TensorGrid &nvv) const
{
    n.fill(Scalar{ 0 });
    nV.fill(Vector{ 0 });
    nvv.fill(Tensor{ 0 });
    //
    Tensor   tmp{ 0 };
    Shape<1> sx;
    for (auto const &ptl : bucket) {
        sx(ptl.pos_x); // position is normalized by grid size
        n.deposit(sx, ptl.w);
        auto const &gvel = ptl.g_vel();
        nV.deposit(sx, gvel * ptl.w);
        tmp.hi() = tmp.lo() = gvel;
        tmp.lo() *= gvel;                       // diagonal part; {vx*vx, vy*vy, vz*vz}
        tmp.hi() *= { gvel.y, gvel.z, gvel.x }; // off-diag part; {vx*vy, vy*vz, vz*vx}
        nvv.deposit(sx, tmp *= ptl.w);
    }
    //
    Real const Nc  = this->Nc;
    VDF const &vdf = *this->vdf;
    (n /= Scalar{ Nc }) += vdf.n0(quiet_nan) * desc.scheme;
    (nV /= Vector{ Nc }) += vdf.nV0(quiet_nan) * desc.scheme;
    (nvv /= Tensor{ Nc }) += vdf.nvv0(quiet_nan) * desc.scheme;
}

namespace {
template <class Object> decltype(auto) write_attr(Object &obj, P1D::PartSpecies const &sp)
{
    obj.attribute("Nc", hdf5::make_type(sp.Nc), hdf5::Space::scalar()).write(sp.Nc);
    obj.attribute("shape_order", hdf5::make_type<long>(sp->shape_order), hdf5::Space::scalar())
        .template write<long>(sp->shape_order);
    obj.attribute("scheme", hdf5::make_type<long>(sp->scheme), hdf5::Space::scalar())
        .template write<long>(sp->scheme);

    return obj << static_cast<P1D::Species const &>(sp);
}
} // namespace
auto P1D::operator<<(hdf5::Group &obj, P1D::PartSpecies const &sp) -> decltype(obj)
{
    return write_attr(obj, sp);
}
auto P1D::operator<<(hdf5::Dataset &obj, P1D::PartSpecies const &sp) -> decltype(obj)
{
    return write_attr(obj, sp);
}
