/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "ParticleRecorder.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>

PIC1D_BEGIN_NAMESPACE
std::string ParticleRecorder::filepath(std::string const &wd, long const step_count) const
{
    if (!is_master())
        throw std::domain_error{ __PRETTY_FUNCTION__ };

    constexpr char    prefix[] = "particle";
    std::string const filename = std::string{ prefix } + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

ParticleRecorder::ParticleRecorder(parallel::mpi::Comm _comm)
: Recorder{ Input::particle_recording_frequency, std::move(_comm) }
, urbg{ 123U + static_cast<unsigned>(comm->rank()) }
{
}

void ParticleRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    if (is_master())
        record_master(domain, step_count);
    else
        record_worker(domain, step_count);
}

template <class Object>
decltype(auto) ParticleRecorder::write_attr(Object &&obj, Domain const &domain, long const step)
{
    obj << domain.params;
    obj.attribute("step", hdf5::make_type(step), hdf5::Space::scalar()).write(step);

    auto const time = step * domain.params.dt;
    obj.attribute("time", hdf5::make_type(time), hdf5::Space::scalar()).write(time);

    return std::forward<Object>(obj);
}
template <class T, std::enable_if_t<std::is_arithmetic_v<T>, int>>
auto ParticleRecorder::get_space(std::vector<T> const &payload)
{
    auto mspace = hdf5::Space::simple(payload.size());
    mspace.select_all();

    auto fspace = hdf5::Space::simple(payload.size());
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
auto ParticleRecorder::get_space(std::vector<Vector> const &payload)
{
    constexpr auto size = 3U;
    static_assert(sizeof(Vector) % sizeof(Real) == 0);
    static_assert(sizeof(Vector) / sizeof(Real) >= size);

    auto mspace = hdf5::Space::simple({ payload.size(), sizeof(Vector) / sizeof(Real) });
    mspace.select(H5S_SELECT_SET, { 0U, 0U }, { payload.size(), size });

    auto fspace = hdf5::Space::simple({ payload.size(), size });
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
auto ParticleRecorder::get_space(std::vector<Particle::PSD> const &payload)
{
    constexpr auto size = 2U;
    static_assert(sizeof(Particle::PSD) % sizeof(Real) == 0);
    static_assert(sizeof(Particle::PSD) / sizeof(Real) == size);

    auto mspace = hdf5::Space::simple({ payload.size(), size });
    mspace.select_all();

    auto fspace = hdf5::Space::simple({ payload.size(), size });
    fspace.select_all();

    return std::make_pair(mspace, fspace);
}
void ParticleRecorder::write_data(std::vector<Particle> ptls, hdf5::Group &root)
{
    using hdf5::make_type;
    using hdf5::Space;
    {
        std::vector<Vector> payload(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(payload), std::mem_fn(&Particle::vel));

        auto const [mspace, fspace] = get_space(payload);
        auto const type             = hdf5::make_type<Real>();
        auto       dset             = root.dataset("vel", type, fspace);
        dset.write(fspace, payload.data(), type, mspace);
    }
    {
        std::vector<Real> payload(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(payload), std::mem_fn(&Particle::pos_x));

        auto const [mspace, fspace] = get_space(payload);
        auto const type             = hdf5::make_type<Real>();
        auto       dset             = root.dataset("pos_x", type, fspace);
        dset.write(fspace, payload.data(), mspace);
    }
    {
        std::vector<Particle::PSD> payload(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(payload), std::mem_fn(&Particle::psd));

        auto const [mspace, fspace] = get_space(payload);
        auto const type             = hdf5::make_type<Real>();
        auto       dset             = root.dataset("psd", type, fspace);
        dset.write(fspace, payload.data(), type, mspace);
    }
    {
        std::vector<long> payload(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(payload), std::mem_fn(&Particle::id));

        auto const [mspace, fspace] = get_space(payload);
        auto const type             = hdf5::make_type<long>();
        auto       dset             = root.dataset("id", type, fspace);
        dset.write(fspace, payload.data(), mspace);
    }
}

void ParticleRecorder::record_master(const Domain &domain, long step_count)
{
    // create hdf file
    std::string const path = filepath(domain.params.working_directory, step_count);
    hdf5::File        file;

    std::vector<unsigned> spids;
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp    = domain.part_species[s];
        auto const         Ndump = Input::Ndumps.at(s);
        if (!Ndump)
            continue;

        spids.push_back(s);
        if (!file)
            file = hdf5::File(hdf5::File::trunc_tag{}, path.c_str());

        // create root group
        auto const name = std::string{ "particle" } + '_' + std::to_string(s);
        auto       root = file.group(name.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl());

        // attributes
        write_attr(root, domain, step_count) << sp;
        root.attribute("Ndump", hdf5::make_type(Ndump), hdf5::Space::scalar()).write(Ndump);

        // datasets
        std::vector<Particle> payload;
        payload.reserve(Ndump);

        // merge
        auto tk = comm.ibsend(sample(sp, Ndump), { master, tag });
        for (int rank = 0, size = comm.size(); rank < size; ++rank) {
            comm.recv<Particle>({}, { rank, tag })
                .unpack(
                    [](auto payload, std::vector<Particle> &buffer) {
                        buffer.insert(buffer.end(), std::make_move_iterator(begin(payload)),
                                      std::make_move_iterator(end(payload)));
                    },
                    payload);
        }
        std::move(tk).wait();

        // dump
        write_data(std::move(payload), root);

        root.flush();
    }

    // save species id's
    if (file) {
        auto space = hdf5::Space::simple(spids.size());
        auto dset  = file.dataset("spids", hdf5::make_type<decltype(spids)::value_type>(), space);
        space.select_all();
        dset.write(space, spids.data(), space);
    }
}
void ParticleRecorder::record_worker(const Domain &domain, long const)
{
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp    = domain.part_species[s];
        auto const         Ndump = Input::Ndumps.at(s);
        if (!Ndump)
            continue;

        comm.ibsend(sample(sp, Ndump), { master, tag }).wait();
    }
}

auto ParticleRecorder::sample(PartSpecies const &sp, unsigned long max_count)
    -> std::vector<Particle>
{
    max_count /= static_cast<unsigned>(comm.size());

    std::vector<Particle> samples;
    samples.reserve(max_count);

    std::sample(sp.bucket.cbegin(), sp.bucket.cend(), std::back_inserter(samples), max_count, urbg);
    for (Particle &ptl : samples) {
        // coordinates relative to the whole simulation domain
        ptl.pos_x += sp.params.domain_extent.min();

        // velocity vector in fac
        ptl.vel = sp.geomtr.cart2fac(ptl.vel);
    }

    return samples;
}
PIC1D_END_NAMESPACE
