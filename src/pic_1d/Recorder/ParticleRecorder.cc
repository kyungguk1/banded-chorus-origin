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

#include "ParticleRecorder.h"

#include "../Utility/TypeMaps.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>

// MARK:- P1D::ParticleRecorder
//
std::string P1D::ParticleRecorder::filepath(std::string const &wd, long const step_count,
                                            unsigned const sp_id) const
{
    if (!is_master())
        throw std::domain_error{__PRETTY_FUNCTION__};

    constexpr char    prefix[] = "particle";
    std::string const filename = std::string{prefix} + "-sp_" + std::to_string(sp_id) + "-"
                               + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

P1D::ParticleRecorder::ParticleRecorder(parallel::mpi::Comm _comm)
: Recorder{Input::particle_recording_frequency, std::move(_comm)}
, urbg{123U + static_cast<unsigned>(comm->rank())}
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
        std::transform(begin(ptls), end(ptls), begin(data), std::mem_fn(&Particle::pos_x));
        dset.write(space, data.data(), space);
    }
    {
        auto space = Space::simple(ptls.size());
        auto dset  = root.dataset("w", make_type<Real>(), space);

        space.select_all();
        std::vector<Real> data(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(data), std::mem_fn(&Particle::w));
        dset.write(space, data.data(), space);
    }
}

void P1D::ParticleRecorder::record_master(const Domain &domain, long step_count)
{
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp = domain.part_species[s];
        if (!Input::Ndumps.at(s))
            continue;

        std::string const path = filepath(domain.params.working_directory, step_count, s);

        // create hdf file and root group
        auto root = hdf5::File(hdf5::File::trunc_tag{}, path.c_str())
                        .group("particle", hdf5::PList::gapl(), hdf5::PList::gcpl());

        // attributes
        write_attr(root, domain, step_count) << sp;
        root.attribute("species", hdf5::make_type(s), hdf5::Space::scalar()).write(s);
        auto const Ndump = Input::Ndumps.at(s);
        root.attribute("Ndump", hdf5::make_type(Ndump), hdf5::Space::scalar()).write(Ndump);

        // datasets
        std::vector<Particle> payload;

        // merge
        auto tk = comm.ibsend(sample(sp, Ndump), master);
        for (int rank = 0, size = comm.size(); rank < size; ++rank) {
            comm.recv<Particle>({}, rank).unpack(
                [](auto payload, std::vector<Particle> &buffer) {
                    buffer.insert(buffer.end(), begin(payload), end(payload));
                },
                payload);
        }
        std::move(tk).wait();

        // dump
        write_data(root, std::move(payload));

        root.flush();
    }
}
void P1D::ParticleRecorder::record_worker(const Domain &domain, long const)
{
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp = domain.part_species[s];
        if (!Input::Ndumps.at(s))
            continue;

        comm.ibsend(sample(sp, Input::Ndumps.at(s)), master).wait();
    }
}

auto P1D::ParticleRecorder::sample(PartSpecies const &sp, unsigned const max_count)
    -> std::vector<Particle>
{
    std::vector<Particle> samples;
    samples.reserve(max_count);

    std::sample(sp.bucket.cbegin(), sp.bucket.cend(), std::back_inserter(samples),
                max_count / static_cast<unsigned>(comm.size()), urbg);
    for (Particle &ptl : samples) {
        ptl.pos_x += sp.params.domain_extent.min(); // coordinates relative to
                                                    // the whole simulation domain
        ptl.vel = sp.geomtr.cart2fac(ptl.vel);      // velocity vector in fac
    }

    return samples;
}
