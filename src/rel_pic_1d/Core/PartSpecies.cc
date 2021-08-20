/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "PartSpecies.h"
#include "BField.h"
#include "EField.h"

#include <stdexcept>
#include <utility>

PIC1D_BEGIN_NAMESPACE
namespace {
template <class T, long N> auto &operator/=(Grid<T, N, Pad> &G, T const w) noexcept
{ // include padding
    for (auto it = G.dead_begin(), end = G.dead_end(); it != end; ++it) {
        *it /= w;
    }
    return G;
}
template <class T, long N> auto &operator+=(Grid<T, N, Pad> &G, T const w) noexcept
{ // exclude padding
    for (auto it = G.begin(), end = G.end(); it != end; ++it) {
        *it += w;
    }
    return G;
}
//
template <class T, long N> auto const &full_grid(Grid<T, N, Pad> &F, BField const &H) noexcept
{
    for (long i = -Pad; i < F.size() + (Pad - 1); ++i) {
        (F[i] = H[i + 1] + H[i + 0]) *= 0.5;
    }
    return F;
}
} // namespace

// ctor
//
PartSpecies::PartSpecies(ParamSet const &params, KineticPlasmaDesc const &_desc, VDFVariant _vdf)
: Species{ params }
, desc{ _desc }
, vdf{ std::move(_vdf) }
, bucket{}
, Nc{ desc.Nc == 0 ? 1.0 : Real(desc.Nc) }
{
    switch (desc.shape_order) {
        case ShapeOrder::_1st:
            m_update_velocity = &PartSpecies::impl_update_velocity<1>;
            m_collect_full_f  = &PartSpecies::impl_collect_full_f<1>;
            m_collect_delta_f = &PartSpecies::impl_collect_delta_f<1>;
            break;
        case ShapeOrder::_2nd:
            m_update_velocity = &PartSpecies::impl_update_velocity<2>;
            m_collect_full_f  = &PartSpecies::impl_collect_full_f<2>;
            m_collect_delta_f = &PartSpecies::impl_collect_delta_f<2>;
            break;
        case ShapeOrder::_3rd:
            m_update_velocity = &PartSpecies::impl_update_velocity<3>;
            m_collect_full_f  = &PartSpecies::impl_collect_full_f<3>;
            m_collect_delta_f = &PartSpecies::impl_collect_delta_f<3>;
            break;
    }
}
void PartSpecies::populate()
{
    bucket.clear();

    // long const Np = desc.Nc * params.Nx;
    for (long i = 0; i < params.Nx; ++i) {
        auto const particles = vdf.emit(static_cast<unsigned>(desc.Nc));
        for (auto const &particle : particles) {
            // loaded particle position is normalized by Dx
            if (params.domain_extent.is_member(particle.pos_x)) {
                auto &ptl = bucket.emplace_back(particle, 1);
                if constexpr (ParamSet::is_relativistic)
                    ptl.gamma = std::sqrt(1 + dot(particle.vel, particle.vel) / params.c2);

                // coordinates relative to this subdomain
                ptl.pos_x -= params.domain_extent.min();
                // initial weight
                switch (desc.scheme) {
                    case PIC::full_f:
                        ptl.psd.w = 1;
                        break;
                    case PIC::delta_f:
                        ptl.psd.w = desc.initial_weight;
                        break;
                }
            }
        }
    }
}

void PartSpecies::load_ptls(std::vector<Particle> const &payload, bool const append)
{
    if (!append)
        bucket.clear();

    for (auto const &loaded : payload) {
        if (params.domain_extent.is_member(loaded.pos_x)) {
            auto &ptl = bucket.emplace_back(loaded);
            // coordinates relative to this subdomain
            ptl.pos_x -= params.domain_extent.min();
        }
    }
}
auto PartSpecies::dump_ptls() const -> std::vector<Particle>
{
    decltype(dump_ptls()) payload{ begin(bucket), end(bucket) };
    for (auto &ptl : payload) {
        // coordinates relative to whole domain
        ptl.pos_x += params.domain_extent.min();
    }
    return payload;
}

// update & collect interface
//
void PartSpecies::update_vel(BField const &bfield, EField const &efield, Real const dt)
{
    (*this->m_update_velocity)(bucket, full_grid(moment<1>(), bfield), efield,
                               BorisPush{ dt, params.c, params.O0, desc.Oc });
}
void PartSpecies::update_pos(Real const dt, Real const fraction_of_grid_size_allowed_to_travel)
{
    Real const dtODx = dt / params.Dx; // normalize position by grid size
    if (!impl_update_x(bucket, dtODx, 1.0 / fraction_of_grid_size_allowed_to_travel)) {
        throw std::domain_error{ std::string{ __PRETTY_FUNCTION__ }
                                 + " - particle(s) moved too far" };
    }
}
void PartSpecies::collect_part()
{
    switch (desc.scheme) {
        case full_f:
            (this->*m_collect_full_f)(moment<1>());
            break;
        case delta_f:
            (this->*m_collect_delta_f)(moment<1>(), bucket);
            break;
    }
}
void PartSpecies::collect_all()
{
    impl_collect(moment<0>(), moment<1>(), moment<2>());
}

// heavy lifting
//
bool PartSpecies::impl_update_x(bucket_type &bucket, Real const dtODx,
                                Real const travel_distance_scale_factor)
{
    bool did_not_move_too_far = true;
    for (auto &ptl : bucket) {
        // gamma should be 1 for non-relativistic case
        Real moved_x = ptl.g_vel.x * (dtODx / ptl.gamma);
        ptl.pos_x += moved_x; // position is normalized by grid size

        // travel distance check
        moved_x *= travel_distance_scale_factor;
        did_not_move_too_far &= 0 == long(moved_x);
    }
    return did_not_move_too_far;
}

template <long Order>
void PartSpecies::impl_update_velocity(bucket_type &bucket, VectorGrid const &B, EField const &E,
                                       BorisPush const &boris)
{
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    for (auto &ptl : bucket) {
        Shape<Order> sx{ ptl.pos_x }; // position is normalized by grid size

        if constexpr (ParamSet::is_relativistic)
            ptl.gamma = boris.relativistic(ptl.g_vel, B.interp(sx), E.interp(sx));
        else
            boris.non_relativistic(ptl.g_vel, B.interp(sx), E.interp(sx));
    }
}

template <long Order> void PartSpecies::impl_collect_full_f(VectorGrid &nV) const
{
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    nV.fill(Vector{ 0 });
    for (auto const &ptl : bucket) {
        Shape<Order> sx{ ptl.pos_x }; // position is normalized by grid size

        if constexpr (ParamSet::is_relativistic)
            nV.deposit(sx, ptl.vel());
        else
            nV.deposit(sx, ptl.g_vel);
    }
    nV /= Vector{ Nc };
}
template <long Order>
void PartSpecies::impl_collect_delta_f(VectorGrid &nV, bucket_type &bucket) const
{
    static_assert(Pad >= Order,
                  "shape order should be less than or equal to the number of ghost cells");
    nV.fill(Vector{ 0 });
    for (auto &ptl : bucket) {
        Shape<Order> sx{ ptl.pos_x }; // position is normalized by grid size
        ptl.psd.w = vdf.weight(ptl);

        if constexpr (ParamSet::is_relativistic)
            nV.deposit(sx, ptl.vel() * ptl.psd.w);
        else
            nV.deposit(sx, ptl.g_vel * ptl.psd.w);
    }
    // TODO: VDF::nV0 is actually average relativistic momentum, not the average velocity.
    (nV /= Vector{ Nc }) += vdf.nV0(Particle::quiet_nan) * desc.scheme;
}
void PartSpecies::impl_collect(ScalarGrid &n, VectorGrid &nV, TensorGrid &nvv) const
{
    n.fill(Scalar{ 0 });
    nV.fill(Vector{ 0 });
    nvv.fill(Tensor{ 0 });
    Tensor tmp{ 0 };
    for (auto const &ptl : bucket) {
        Shape<1> sx{ ptl.pos_x }; // position is normalized by grid size
        n.deposit(sx, ptl.psd.w);
        nV.deposit(sx, ptl.g_vel * ptl.psd.w);
        tmp.hi() = tmp.lo() = ptl.g_vel;
        tmp.lo() *= ptl.g_vel; // diagonal part; {vx*vx, vy*vy, vz*vz}
        tmp.hi()
            *= { ptl.g_vel.y, ptl.g_vel.z, ptl.g_vel.x }; // off-diag part; {vx*vy, vy*vz, vz*vx}
        nvv.deposit(sx, tmp *= ptl.psd.w);
    }
    (n /= Scalar{ Nc }) += vdf.n0(Particle::quiet_nan) * desc.scheme;
    (nV /= Vector{ Nc }) += vdf.nV0(Particle::quiet_nan) * desc.scheme;
    (nvv /= Tensor{ Nc }) += vdf.nvv0(Particle::quiet_nan) * desc.scheme;
}

namespace {
template <class Object> decltype(auto) write_attr(Object &obj, PartSpecies const &sp)
{
    obj.attribute("Nc", hdf5::make_type(sp.Nc), hdf5::Space::scalar()).write(sp.Nc);
    obj.attribute("shape_order", hdf5::make_type<long>(sp->shape_order), hdf5::Space::scalar())
        .template write<long>(sp->shape_order);
    obj.attribute("particle_scheme", hdf5::make_type<long>(sp->scheme), hdf5::Space::scalar())
        .template write<long>(sp->scheme);

    return obj << static_cast<Species const &>(sp);
}
} // namespace
auto operator<<(hdf5::Group &obj, PartSpecies const &sp) -> decltype(obj)
{
    return write_attr(obj, sp);
}
auto operator<<(hdf5::Dataset &obj, PartSpecies const &sp) -> decltype(obj)
{
    return write_attr(obj, sp);
}
PIC1D_END_NAMESPACE
