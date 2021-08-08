/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "ParticleRecorder.h"

#include "../Utility/TypeMaps.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>

// MARK:- P1D::ParticleRecorder
//
std::string P1D::ParticleRecorder::filepath(std::string const &wd, long const step_count) const
{
    if (!is_master())
        throw std::domain_error{ __PRETTY_FUNCTION__ };

    constexpr char    prefix[] = "particle";
    std::string const filename = std::string{ prefix } + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

P1D::ParticleRecorder::ParticleRecorder(parallel::mpi::Comm _comm)
: Recorder{ Input::particle_recording_frequency, std::move(_comm) }
, urbg{ 123U + static_cast<unsigned>(comm->rank()) }
{
}

void P1D::ParticleRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    if (is_master())
        record_master(domain, step_count);
    else
        record_worker(domain, step_count);
}

template <class Object>
decltype(auto) P1D::ParticleRecorder::write_attr(Object &&obj, Domain const &domain,
                                                 long const step)
{
    obj << domain.params;
    obj.attribute("step", hdf5::make_type(step), hdf5::Space::scalar()).write(step);

    auto const time = step * domain.params.dt;
    obj.attribute("time", hdf5::make_type(time), hdf5::Space::scalar()).write(time);

    return std::forward<Object>(obj);
}
void P1D::ParticleRecorder::write_data(hdf5::Group &root, std::vector<Particle> ptls)
{
    using hdf5::make_type;
    using hdf5::Space;
    {
        auto space = Space::simple(ptls.size());
        auto dset  = root.dataset("vel", make_type<Vector>(), space);

        space.select_all();
        std::vector<Vector> data(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(data), std::mem_fn(&Particle::vel));
        dset.write(space, data.data(), space);
    }
    {
        auto space = Space::simple(ptls.size());
        auto dset  = root.dataset("pos_x", make_type<Real>(), space);

        space.select_all();
        std::vector<Real> data(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(data), [](auto const &ptl) {
            return ptl.pos_x;
        });
        dset.write(space, data.data(), space);
    }
    {
        auto space = Space::simple(ptls.size());
        auto dset  = root.dataset("w", make_type<Real>(), space);

        space.select_all();
        std::vector<Real> data(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(data), [](auto const &ptl) {
            return ptl.w;
        });
        dset.write(space, data.data(), space);
    }
    {
        auto space = Space::simple(ptls.size());
        auto dset  = root.dataset("gamma", make_type<Real>(), space);

        space.select_all();
        std::vector<Real> data(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(data), std::mem_fn(&Particle::gamma));
        dset.write(space, data.data(), space);
    }
}

void P1D::ParticleRecorder::record_master(const Domain &domain, long step_count)
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
        auto const name = std::string{ "particle" } + "[" + std::to_string(s) + "]";
        auto       root = file.group(name.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl());

        // attributes
        write_attr(root, domain, step_count) << sp;
        root.attribute("Ndump", hdf5::make_type(Ndump), hdf5::Space::scalar()).write(Ndump);

        // datasets
        std::vector<Particle> payload;
        payload.reserve(Ndump);

        // merge
        auto tk = comm.ibsend(sample(sp, Ndump), master);
        for (int rank = 0, size = comm.size(); rank < size; ++rank) {
            comm.recv<Particle>({}, rank).unpack(
                [](auto payload, std::vector<Particle> &buffer) {
                    buffer.insert(buffer.end(), std::make_move_iterator(begin(payload)),
                                  std::make_move_iterator(end(payload)));
                },
                payload);
        }
        std::move(tk).wait();

        // dump
        write_data(root, std::move(payload));

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
void P1D::ParticleRecorder::record_worker(const Domain &domain, long const)
{
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp    = domain.part_species[s];
        auto const         Ndump = Input::Ndumps.at(s);
        if (!Ndump)
            continue;

        comm.ibsend(sample(sp, Ndump), master).wait();
    }
}

auto P1D::ParticleRecorder::sample(PartSpecies const &sp, unsigned long max_count)
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
        ptl.g_vel() = sp.geomtr.cart2fac(ptl.g_vel());
    }

    return samples;
}
