/*
 * Copyright (c) 2019-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "ParticleRecorder.h"
#include "MomentRecorder.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>

PIC1D_BEGIN_NAMESPACE
std::string ParticleRecorder::filepath(std::string const &wd, long const step_count) const
{
    if (!is_subdomain_master() || !is_distributed_particle_master())
        throw std::domain_error{ __PRETTY_FUNCTION__ };

    constexpr char    prefix[] = "particle";
    std::string const filename = std::string{ prefix } + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

ParticleRecorder::ParticleRecorder(parallel::mpi::Comm _subdomain_comm, parallel::mpi::Comm _distributed_particle_comm)
: Recorder{ Input::particle_recording_frequency, std::move(_subdomain_comm), std::move(_distributed_particle_comm) }
, urbg{ 123U + static_cast<unsigned>(subdomain_comm->rank() + subdomain_comm.size() * distributed_particle_comm->rank()) }
{
}

void ParticleRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    if (is_subdomain_master() && is_distributed_particle_master())
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
template <class T>
auto ParticleRecorder::write_data(std::vector<T> payload, hdf5::Group &root, char const *name)
{
    auto const [mspace, fspace] = get_space(payload);
    auto const type             = hdf5::make_type<Real>();
    auto       dset             = root.dataset(name, type, fspace);
    dset.write(fspace, payload.data(), type, mspace);
    return dset;
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
        std::vector<CurviCoord> payload(ptls.size());
        std::transform(begin(ptls), end(ptls), begin(payload), std::mem_fn(&Particle::pos));

        auto const [mspace, fspace] = get_space(payload);
        auto const type             = hdf5::make_type<Real>();
        auto       dset             = root.dataset("pos", type, fspace);
        dset.write(fspace, payload.data(), type, mspace);
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
        dset.write(fspace, payload.data(), type, mspace);
    }
}

void ParticleRecorder::record_master(const Domain &domain, long const step_count)
{
    // create hdf file and root group
    std::string const path = filepath(domain.params.working_directory, step_count);
    hdf5::Group       root;

    std::vector<unsigned> spids;
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp    = domain.part_species[s];
        auto const         Ndump = Input::Ndumps.at(s);
        if (!Ndump)
            continue;

        spids.push_back(s);
        if (!root) {
            root = hdf5::File(hdf5::File::trunc_tag{}, path.c_str())
                       .group("particle", hdf5::PList::gapl(), hdf5::PList::gcpl());
            write_attr(root, domain, step_count);
        }

        // create species group
        auto parent = [&root, name = std::to_string(s)] {
            return root.group(name.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl());
        }();
        write_attr(parent, domain, step_count) << sp;
        parent.attribute("Ndump", hdf5::make_type(Ndump), hdf5::Space::scalar()).write(Ndump);

        // moments
        auto const writer = [](auto payload, auto &root, auto *name) {
            return write_data(std::move(payload), root, name);
        };
        if (auto const &comm = subdomain_comm; comm->operator bool()) {
            comm.gather<0>({ sp.moment<0>().begin(), sp.moment<0>().end() }, master)
                .unpack(writer, parent, "n");
            comm.gather<1>(MomentRecorder::convert(sp.moment<1>(), domain.params.geomtr), master)
                .unpack(writer, parent, "nV");
            comm.gather<2>(MomentRecorder::convert(sp.moment<2>(), domain.params.geomtr), master)
                .unpack(writer, parent, "nvv");
        }

        // particles
        auto payload = collect_particles(collect_particles(sample(sp, Ndump), subdomain_comm), distributed_particle_comm);
        write_data(std::move(payload), parent);
    }

    // save species id's
    if (root) {
        auto space = hdf5::Space::simple(spids.size());
        auto dset  = root.dataset("spids", hdf5::make_type<decltype(spids)::value_type>(), space);
        space.select_all();
        dset.write(space, spids.data(), space);
        root.flush();
    }
}
void ParticleRecorder::record_worker(const Domain &domain, long const)
{
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp    = domain.part_species[s];
        auto const         Ndump = Input::Ndumps.at(s);
        if (!Ndump)
            continue;

        // moments
        if (auto const &comm = subdomain_comm; is_distributed_particle_master()) {
            comm.gather<0>({ sp.moment<0>().begin(), sp.moment<0>().end() }, master)
                .unpack([](auto) {});
            comm.gather<1>(MomentRecorder::convert(sp.moment<1>(), domain.params.geomtr), master)
                .unpack([](auto) {});
            comm.gather<2>(MomentRecorder::convert(sp.moment<2>(), domain.params.geomtr), master)
                .unpack([](auto) {});
        }

        // particles
        collect_particles(collect_particles(sample(sp, Ndump), subdomain_comm), distributed_particle_comm);
    }
}
template <class InterprocessComm>
auto ParticleRecorder::collect_particles(std::vector<Particle> incomming, InterprocessComm const &comm) -> std::vector<Particle>
{
    std::vector<Particle> collected;
    collected.reserve(incomming.size() * unsigned(comm.size()));

    auto tk = comm.ibsend(std::move(incomming), { master, tag });
    if (master == comm->rank()) {
        for (int rank = 0, size = comm.size(); rank < size; ++rank) {
            comm.template recv<Particle>({}, { rank, tag })
                .unpack(
                    [](auto payload, std::vector<Particle> &buffer) {
                        buffer.insert(buffer.end(),
                                      std::make_move_iterator(begin(payload)), std::make_move_iterator(end(payload)));
                    },
                    collected);
        }
    }
    std::move(tk).wait();

    return collected;
}

auto ParticleRecorder::sample(PartSpecies const &sp, unsigned long max_count) -> std::vector<Particle>
{
    auto const divisor = unsigned(world_size());
    max_count          = max_count / divisor + (max_count % divisor ? 1 : 0);

    std::vector<Particle> samples;
    samples.reserve(max_count);

    std::sample(sp.bucket.cbegin(), sp.bucket.cend(), std::back_inserter(samples), max_count, urbg);
    for (Particle &ptl : samples) {
        // velocity vector in fac
        ptl.vel = sp.geomtr.cart_to_fac(ptl.vel, ptl.pos);
    }

    return samples;
}
PIC1D_END_NAMESPACE
