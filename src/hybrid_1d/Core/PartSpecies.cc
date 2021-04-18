//
// Copyright (c) 2019, Kyungguk Min
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// - Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "PartSpecies.h"

#include "./BField.h"
#include "./EField.h"

#include <stdexcept>
#include <utility>

// helpers
//
namespace {
template <class T, long N> auto &operator/=(H1D::GridQ<T, N> &G, T const w) noexcept
{ // include padding
    for (auto it = G.dead_begin(), end = G.dead_end(); it != end; ++it) {
        *it /= w;
    }
    return G;
}
template <class T, long N> auto &operator+=(H1D::GridQ<T, N> &G, T const w) noexcept
{ // exclude padding
    for (auto it = G.begin(), end = G.end(); it != end; ++it) {
        *it += w;
    }
    return G;
}
//
template <class T, long N> auto const &full_grid(H1D::GridQ<T, N> &F, H1D::BField const &H) noexcept
{
    for (long i = -H1D::Pad; i < F.size() + (H1D::Pad - 1); ++i) {
        (F[i] = H[i + 1] + H[i + 0]) *= 0.5;
    }
    return F;
}
} // namespace

// constructor
//
H1D::PartSpecies::PartSpecies(ParamSet const &params, KineticPlasmaDesc const &desc,
                              std::unique_ptr<VDF> _vdf)
: Species{params}, desc{desc}, vdf{std::move(_vdf)}, bucket{}, Nc{desc.Nc == 0 ? 1.0 : desc.Nc}
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
void H1D::PartSpecies::populate()
{
    bucket.clear();
    long const Np = desc.Nc * Input::Nx;
    for (long i = 0; i < Np; ++i) {
        Particle ptl = vdf->variate(); // position is normalized by Dx
        if (params.domain_extent.is_member(ptl.pos_x)) {
            ptl.pos_x -= params.domain_extent.min(); // coordinates relative to this subdomain
            bucket.emplace_back(ptl).w = desc.scheme == full_f;
        }
    }
}

void H1D::PartSpecies::load_ptls(std::vector<Particle> const &payload)
{
    bucket.clear();
    for (Particle const &ptl : payload) {
        if (params.domain_extent.is_member(ptl.pos_x)) {
            bucket.emplace_back(ptl).pos_x -= params.domain_extent.min(); // coordinates relative to
                                                                          // this subdomain
        }
    }
}
auto H1D::PartSpecies::dump_ptls() const -> std::vector<Particle>
{
    decltype(dump_ptls()) payload{begin(bucket), end(bucket)};
    for (Particle &ptl : payload) {
        ptl.pos_x += params.domain_extent.min(); // coordinates relative to whole domain
    }
    return payload; // NRVO
}

// update & collect interface
//
void H1D::PartSpecies::update_vel(BField const &bfield, EField const &efield, Real const dt)
{
    (*this->_update_velocity)(bucket, full_grid(moment<1>(), bfield), efield,
                              BorisPush{dt, params.c, params.O0, desc.Oc});
}
void H1D::PartSpecies::update_pos(Real const dt, Real const fraction_of_grid_size_allowed_to_travel)
{
    Real const dtODx = dt / params.Dx; // normalize position by grid size
    if (!_update_x(bucket, dtODx, 1.0 / fraction_of_grid_size_allowed_to_travel)) {
        throw std::domain_error{std::string{__FUNCTION__} + " - particle(s) moved too far"};
    }
}
void H1D::PartSpecies::collect_part()
{
    switch (desc.scheme) {
        case full_f:
            (this->*_collect_full_f)(moment<0>(), moment<1>());
            break;
        case delta_f:
            (this->*_collect_delta_f)(moment<0>(), moment<1>(), bucket);
            break;
    }
}
void H1D::PartSpecies::collect_all()
{
    _collect(moment<0>(), moment<1>(), moment<2>());
}

// heavy lifting
//
bool H1D::PartSpecies::_update_x(bucket_type &bucket, Real const dtODx,
                                 Real const travel_distance_scale_factor)
{
    bool did_not_move_too_far = true;
    for (Particle &ptl : bucket) {
        Real moved_x = ptl.vel.x * dtODx;
        ptl.pos_x += moved_x; // position is normalized by grid size

        // travel distance check
        //
        moved_x *= travel_distance_scale_factor;
        did_not_move_too_far &= 0 == long(moved_x);
    }
    return did_not_move_too_far;
}

template <long Order>
void H1D::PartSpecies::_update_velocity_(bucket_type &bucket, VectorGrid const &B, EField const &E,
                                         BorisPush const pusher)
{
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    Shape<Order> sx;
    for (Particle &ptl : bucket) {
        sx(ptl.pos_x); // position is normalized by grid size
        pusher(ptl.vel, B.interp(sx), E.interp(sx));
    }
}

template <long Order> void H1D::PartSpecies::_collect_full_f_(ScalarGrid &n, VectorGrid &nV) const
{
    n.fill(Scalar{0});
    nV.fill(Vector{0});
    //
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    Shape<Order> sx;
    for (Particle const &ptl : bucket) {
        sx(ptl.pos_x); // position is normalized by grid size
        n.deposit(sx, 1);
        nV.deposit(sx, ptl.vel);
    }
    //
    Real const Nc = this->Nc;
    n /= Scalar{Nc};
    nV /= Vector{Nc};
}
template <long Order>
void H1D::PartSpecies::_collect_delta_f_(ScalarGrid &n, VectorGrid &nV, bucket_type &bucket) const
{
    VDF const &vdf = *this->vdf;
    //
    n.fill(Scalar{0});
    nV.fill(Vector{0});
    //
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    Shape<Order> sx;
    for (Particle &ptl : bucket) {
        sx(ptl.pos_x); // position is normalized by grid size
        ptl.w = vdf.weight(ptl);
        n.deposit(sx, ptl.w);
        nV.deposit(sx, ptl.vel * ptl.w);
    }
    //
    Real const Nc = this->Nc;
    (n /= Scalar{Nc}) += vdf.n0(Particle::quiet_nan) * desc.scheme;
    (nV /= Vector{Nc}) += vdf.nV0(Particle::quiet_nan) * desc.scheme;
}
void H1D::PartSpecies::_collect(ScalarGrid &n, VectorGrid &nV, TensorGrid &nvv) const
{
    n.fill(Scalar{0});
    nV.fill(Vector{0});
    nvv.fill(Tensor{0});
    //
    Tensor   tmp{0};
    Shape<1> sx;
    for (Particle const &ptl : bucket) {
        sx(ptl.pos_x); // position is normalized by grid size
        n.deposit(sx, ptl.w);
        nV.deposit(sx, ptl.vel * ptl.w);
        tmp.hi() = tmp.lo() = ptl.vel;
        tmp.lo() *= ptl.vel;                           // diagonal part; {vx*vx, vy*vy, vz*vz}
        tmp.hi() *= {ptl.vel.y, ptl.vel.z, ptl.vel.x}; // off-diag part; {vx*vy, vy*vz, vz*vx}
        nvv.deposit(sx, tmp *= ptl.w);
    }
    //
    Real const Nc  = this->Nc;
    VDF const &vdf = *this->vdf;
    (n /= Scalar{Nc}) += vdf.n0(Particle::quiet_nan) * desc.scheme;
    (nV /= Vector{Nc}) += vdf.nV0(Particle::quiet_nan) * desc.scheme;
    (nvv /= Tensor{Nc}) += vdf.nvv0(Particle::quiet_nan) * desc.scheme;
}
