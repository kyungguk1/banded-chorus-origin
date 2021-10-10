/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "PartSpecies.h"
#include "BField.h"
#include "EField.h"

#include <limits>
#include <stdexcept>
#include <utility>

PIC1D_BEGIN_NAMESPACE
namespace {
template <class T, long N>
auto &operator/=(Grid<T, N, Pad> &G, T const w) noexcept
{ // include padding
    for (auto it = G.dead_begin(), end = G.dead_end(); it != end; ++it) {
        *it /= w;
    }
    return G;
}

template <class T, long N>
auto const &full_grid(Grid<T, N, Pad> &F, BField const &H) noexcept
{
    F[-Pad] = T{ std::numeric_limits<Real>::quiet_NaN() };
    for (long i = -Pad + 1; i < F.size() + Pad; ++i) {
        (F[i] = H[i - 0] + H[i - 1]) *= 0.5;
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
, Nc{ Particle::quiet_nan }
, equilibrium_macro_weight{ Real(desc.scheme) / params.number_of_particle_parallelism }
, bucket{}
{
    if (long const Np = params.Nx * desc.Nc)
        Nc = Real(Np) * vdf.Nrefcell_div_Ntotal();
    else
        Nc = 1;

    switch (desc.shape_order) {
        case ShapeOrder::_1st:
            if constexpr (Pad >= 1) {
                m_update_velocity = &PartSpecies::impl_update_velocity<1>;
                m_collect_full_f  = &PartSpecies::impl_collect_full_f<1>;
                m_collect_delta_f = &PartSpecies::impl_collect_delta_f<1>;
            }
            break;
        case ShapeOrder::_2nd:
            if constexpr (Pad >= 2) {
                m_update_velocity = &PartSpecies::impl_update_velocity<2>;
                m_collect_full_f  = &PartSpecies::impl_collect_full_f<2>;
                m_collect_delta_f = &PartSpecies::impl_collect_delta_f<2>;
            }
            break;
        case ShapeOrder::_3rd:
            if constexpr (Pad >= 3) {
                m_update_velocity = &PartSpecies::impl_update_velocity<3>;
                m_collect_full_f  = &PartSpecies::impl_collect_full_f<3>;
                m_collect_delta_f = &PartSpecies::impl_collect_delta_f<3>;
            }
            break;
    }
}
void PartSpecies::populate()
{
    bucket.clear();

    // long const Np = desc.Nc * params.Nx;
    for (long i = 0; i < params.Nx; ++i) {
        auto const particles = vdf.emit(static_cast<unsigned long>(desc.Nc));
        for (auto const &particle : particles) {
            // take those that belong in this subdomain
            if (params.full_grid_subdomain_extent.is_member(particle.pos.q1))
                bucket.push_back(particle);
        }
    }
}
void PartSpecies::populate(long const color, long const divisor)
{
    bucket.clear();

    if (divisor <= 0)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - non-positive divisor" };
    if (color < 0 || color >= divisor)
        throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ } + " - invalid color range" };

    long Np_within_this_subdomain = 0;
    for (long i = 0; i < params.Nx; ++i) {
        auto const particles = vdf.emit(static_cast<unsigned long>(desc.Nc));
        for (auto const &particle : particles) {
            // take those that belong in this subdomain
            if (params.full_grid_subdomain_extent.is_member(particle.pos.q1)
                && color == Np_within_this_subdomain++ % divisor /*must increment iterations after all other tests*/)
                bucket.push_back(particle);
        }
    }
}

void PartSpecies::load_ptls(std::vector<Particle> const &payload, bool const append)
{
    if (!append)
        bucket.clear();

    for (auto const &loaded : payload) {
        if (params.full_grid_subdomain_extent.is_member(loaded.pos.q1))
            bucket.push_back(loaded);
    }
}
auto PartSpecies::dump_ptls() const -> std::vector<Particle>
{
    return { begin(bucket), end(bucket) };
}

// update & collect interface
//
void PartSpecies::update_vel(BField const &bfield, EField const &efield, Real const dt)
{
    (this->*m_update_velocity)(bucket, full_grid(moment<1>(), bfield), efield, BorisPush{ dt, params.c, params.O0, desc.Oc });
}
void PartSpecies::update_pos(Real const dt, Real const fraction_of_grid_size_allowed_to_travel)
{
    if (!impl_update_pos(bucket, dt, 1.0 / fraction_of_grid_size_allowed_to_travel))
        throw std::domain_error{ std::string{ __PRETTY_FUNCTION__ } + " - particle(s) moved too far" };
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
bool PartSpecies::impl_update_pos(bucket_type &bucket, Real const dt, Real const travel_distance_scale_factor) const
{
    bool did_not_move_too_far = true;
    for (auto &ptl : bucket) {
        auto const v_contr  = geomtr.cart_to_contr(ptl.vel, ptl.pos);
        Real       moved_q1 = v_contr.x * dt;
        ptl.pos.q1 += moved_q1;

        // travel distance check
        moved_q1 *= travel_distance_scale_factor;
        did_not_move_too_far &= 0 == long(moved_q1);
    }
    return did_not_move_too_far;
}

template <long Order>
void PartSpecies::impl_update_velocity(bucket_type &bucket, VectorGrid const &dB, EField const &E, BorisPush const &boris) const
{
    static_assert(Pad >= Order, "shape order should be less than or equal to the number of ghost cells");
    auto const q1min = params.full_grid_subdomain_extent.min();
    for (auto &ptl : bucket) {
        Shape<Order> const sx{ ptl.pos.q1 - q1min };

        // get gyro-radius offset: rL = e1 x v / Oc (Oc is signed)
        auto const Oc   = desc.Oc * geomtr.Bmag_div_B0(ptl.pos);
        auto const rL_y = -ptl.vel.z / Oc;
        auto const rL_z = +ptl.vel.y / Oc;

        boris.non_relativistic(ptl.vel, geomtr.Bcart(ptl.pos, rL_y, rL_z) + dB.interp(sx), E.interp(sx));
    }
}

template <long Order>
void PartSpecies::impl_collect_full_f(VectorGrid &nV) const
{
    static_assert(Pad >= Order, "shape order should be less than or equal to the number of ghost cells");
    auto const q1min = params.full_grid_subdomain_extent.min();
    nV.fill(Vector{ 0 });
    for (auto const &ptl : bucket) {
        Shape<Order> const sx{ ptl.pos.q1 - q1min };
        nV.deposit(sx, ptl.vel);
    }
    nV /= Vector{ Nc };
}
template <long Order>
void PartSpecies::impl_collect_delta_f(VectorGrid &nV, bucket_type &bucket) const
{
    static_assert(Pad >= Order, "shape order should be less than or equal to the number of ghost cells");
    auto const q1min = params.full_grid_subdomain_extent.min();
    nV.fill(Vector{ 0 });
    for (auto &ptl : bucket) {
        Shape<Order> const sx{ ptl.pos.q1 - q1min };
        ptl.psd.weight = vdf.weight(ptl);
        nV.deposit(sx, ptl.vel * ptl.psd.weight);
    }
    nV /= Vector{ Nc };
    for (long i = 0; i < nV.size(); ++i) {
        CurviCoord const pos{ i + q1min };
        nV[i] += vdf.nV0(pos) * equilibrium_macro_weight;
    }
}
void PartSpecies::impl_collect(ScalarGrid &n, VectorGrid &nV, TensorGrid &nvv) const
{
    n.fill(Scalar{ 0 });
    nV.fill(Vector{ 0 });
    nvv.fill(Tensor{ 0 });
    Tensor     tmp{ 0 };
    auto const q1min = params.full_grid_subdomain_extent.min();
    for (auto const &ptl : bucket) {
        Shape<1> const sx{ ptl.pos.q1 - q1min };
        n.deposit(sx, ptl.psd.weight);
        nV.deposit(sx, ptl.vel * ptl.psd.weight);
        tmp.hi() = tmp.lo() = ptl.vel;
        tmp.lo() *= ptl.vel;                             // diagonal part; {vx*vx, vy*vy, vz*vz}
        tmp.hi() *= { ptl.vel.y, ptl.vel.z, ptl.vel.x }; // off-diag part; {vx*vy, vy*vz, vz*vx}
        nvv.deposit(sx, tmp *= ptl.psd.weight);
    }
    n /= Scalar{ Nc };
    nV /= Vector{ Nc };
    nvv /= Tensor{ Nc };
    for (long i = 0; i < nV.size(); ++i) {
        CurviCoord const pos{ i + q1min };
        n[i] += vdf.n0(pos) * equilibrium_macro_weight;
        nV[i] += vdf.nV0(pos) * equilibrium_macro_weight;
        nvv[i] += vdf.nvv0(pos) * equilibrium_macro_weight;
    }
}

template <class Object>
auto write_attr(Object &obj, PartSpecies const &sp) -> decltype(obj)
{
    obj.attribute("Nc", hdf5::make_type(sp->Nc), hdf5::Space::scalar())
        .write(sp->Nc);
    obj.attribute("Nrefcell_div_Ntotal", hdf5::make_type(sp.vdf.Nrefcell_div_Ntotal()), hdf5::Space::scalar())
        .write(sp.vdf.Nrefcell_div_Ntotal());
    obj.attribute("shape_order", hdf5::make_type<long>(sp->shape_order), hdf5::Space::scalar())
        .template write<long>(sp->shape_order);
    obj.attribute("particle_scheme", hdf5::make_type<long>(sp->scheme), hdf5::Space::scalar())
        .template write<long>(sp->scheme);
    obj.attribute("initial_weight", hdf5::make_type(sp->initial_weight), hdf5::Space::scalar())
        .write(sp->initial_weight);
    obj.attribute("marker_temp_ratio", hdf5::make_type(sp->marker_temp_ratio), hdf5::Space::scalar())
        .write(sp->marker_temp_ratio);

    return obj << static_cast<Species const &>(sp);
}
auto operator<<(hdf5::Group &obj, PartSpecies const &sp) -> decltype(obj)
{
    return write_attr(obj, sp);
}
auto operator<<(hdf5::Dataset &obj, PartSpecies const &sp) -> decltype(obj)
{
    return write_attr(obj, sp);
}
PIC1D_END_NAMESPACE
