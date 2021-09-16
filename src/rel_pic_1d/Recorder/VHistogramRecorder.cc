/*
 * Copyright (c) 2020-2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "VHistogramRecorder.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <type_traits>

PIC1D_BEGIN_NAMESPACE
namespace {
template <class T1, class T2, class U1, class U2>
constexpr decltype(auto) operator+=(std::pair<T1, T2> &lhs, std::pair<U1, U2> const &rhs) noexcept(
    noexcept(std::declval<T1 &>() += std::declval<U1>()))
{
    std::get<0>(lhs) += std::get<0>(rhs);
    std::get<1>(lhs) += std::get<1>(rhs);
    return lhs;
}
template <class T, class U>
[[nodiscard]] constexpr auto operator+(std::pair<T, T> a, U const &b) noexcept(
    noexcept(std::declval<T &>() += std::declval<U>()))
{
    a += std::make_pair(b, b);
    return a;
}
template <class T, class U>
constexpr decltype(auto) operator/=(std::pair<T, T> &lhs, U const &rhs) noexcept(
    noexcept(std::declval<T &>() /= std::declval<U>()))
{
    std::get<0>(lhs) /= rhs;
    std::get<1>(lhs) /= rhs;
    return lhs;
}
} // namespace

std::string VHistogramRecorder::filepath(std::string const &wd, long const step_count) const
{
    constexpr char    prefix[] = "vhist2d";
    std::string const filename = std::string{ prefix } + "-" + std::to_string(step_count) + ".h5";
    return wd + "/" + filename;
}

VHistogramRecorder::VHistogramRecorder(parallel::mpi::Comm _comm)
: Recorder{ Input::vhistogram_recording_frequency, std::move(_comm) }
{
}

void VHistogramRecorder::record(const Domain &domain, const long step_count)
{
    if (step_count % recording_frequency)
        return;

    if (is_master())
        record_master(domain, step_count);
    else
        record_worker(domain, step_count);
}

class VHistogramRecorder::Indexer {
    // preconditions:
    // 1. length of span is positive
    // 2. dim is positive
    //
    Range    v1span;
    Range    v2span;
    unsigned v1dim;
    unsigned v2dim;

public:
    constexpr Indexer(Range const &v1span, unsigned const &v1dim, Range const &v2span,
                      unsigned const &v2dim) noexcept
    : v1span{ v1span }, v2span{ v2span }, v1dim{ v1dim }, v2dim{ v2dim }
    {
    }

    [[nodiscard]] constexpr explicit operator bool() const noexcept
    {
        return (v1dim > 0) && (v2dim > 0);
    }

    using index_pair_t = local_vhist_t::key_type;
    static constexpr index_pair_t npos{
        std::numeric_limits<long>::max(),
        std::numeric_limits<long>::max(),
    };

    [[nodiscard]] auto operator()(Real const v1, Real const v2) const noexcept
    {
        // zero-based indexing
        //
        index_pair_t const idx = {
            static_cast<index_pair_t::first_type>(
                std::floor((v1 - v1span.min()) * v1dim / v1span.len)),
            static_cast<index_pair_t::second_type>(
                std::floor((v2 - v2span.min()) * v2dim / v2span.len)),
        };

        if (within(idx, std::make_pair(0, 0), std::make_pair(v1dim, v2dim),
                   std::make_index_sequence<std::tuple_size_v<index_pair_t>>{}))
            return idx;

        return npos;
    }

private:
    template <std::size_t... I>
    [[nodiscard]] static bool within(index_pair_t const &idx, index_pair_t const &min, index_pair_t const &max, std::index_sequence<I...>) noexcept
    {
        return (... && (std::get<I>(idx) >= std::get<I>(min)))
            && (... && (std::get<I>(idx) < std::get<I>(max)));
    }
};

template <class Object>
decltype(auto) VHistogramRecorder::write_attr(Object &&obj, Domain const &domain, long const step)
{
    obj << domain.params;
    obj.attribute("step", hdf5::make_type(step), hdf5::Space::scalar()).write(step);

    auto const time = step * domain.params.dt;
    obj.attribute("time", hdf5::make_type(time), hdf5::Space::scalar()).write(time);

    return std::forward<Object>(obj);
}
void VHistogramRecorder::write_data(hdf5::Group &root, global_vhist_t vhist)
{
    using hdf5::make_type;
    using hdf5::Space;
    using Value  = decltype(vhist)::value_type;
    using Index  = decltype(vhist)::key_type;
    using Mapped = decltype(vhist)::mapped_type;
    {
        std::vector<Index> payload(vhist.size());
        std::transform(begin(vhist), end(vhist), begin(payload), std::mem_fn(&Value::first));

        auto const [mspace, fspace] = get_space(payload);
        auto const type             = hdf5::make_type<std::tuple_element_t<0, Index>>();
        auto       dset             = root.dataset("idx", type, fspace);
        dset.write(fspace, payload.data(), type, mspace);
    }
    {
        std::vector<Mapped> payload(vhist.size());
        std::transform(begin(vhist), end(vhist), begin(payload), std::mem_fn(&Value::second));

        auto const [mspace, fspace] = get_space(payload);
        auto const type             = hdf5::make_type<Real>();
        auto       dset             = root.dataset("psd", type, fspace);
        dset.write(fspace, payload.data(), type, mspace);
    }
}

void VHistogramRecorder::record_master(const Domain &domain, long step_count)
{
    // create hdf file and root group
    std::string const path = filepath(domain.params.working_directory, step_count);
    hdf5::Group       root;

    std::vector<unsigned> spids;
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp       = domain.part_species[s];
        auto const [v1span, v1divs] = Input::v1hist_specs.at(s);
        auto const [v2span, v2divs] = Input::v2hist_specs.at(s);
        Indexer const idxer{ v1span, v1divs, v2span, v2divs };
        if (!idxer)
            continue;

        if (v1span.len <= 0 || v2span.len <= 0) {
            throw std::invalid_argument{ std::string{ __PRETTY_FUNCTION__ }
                                         + " - invalid vspan extent: " + std::to_string(s)
                                         + "th species" };
        }

        spids.push_back(s);
        if (!root) {
            root = hdf5::File(hdf5::File::trunc_tag{}, path.c_str())
                       .group("vhist2d", hdf5::PList::gapl(), hdf5::PList::gcpl());
            write_attr(root, domain, step_count);
        }

        // create species group
        auto parent = [&root, name = std::to_string(s)] {
            return root.group(name.c_str(), hdf5::PList::gapl(), hdf5::PList::gcpl());
        }();
        write_attr(parent, domain, step_count) << sp;

        auto const v1lim = std::make_pair(v1span.min(), v1span.max());
        parent.attribute("v1lim", hdf5::make_type(v1lim), hdf5::Space::scalar()).write(v1lim);
        auto const v2lim = std::make_pair(v2span.min(), v2span.max());
        parent.attribute("v2lim", hdf5::make_type(v2lim), hdf5::Space::scalar()).write(v2lim);
        auto const vdims = std::make_pair(v1divs, v2divs);
        parent.attribute("vdims", hdf5::make_type(vdims), hdf5::Space::scalar()).write(vdims);

        // velocity histogram
        write_data(parent, histogram(sp, idxer));
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
void VHistogramRecorder::record_worker(const Domain &domain, long const)
{
    for (unsigned s = 0; s < domain.part_species.size(); ++s) {
        PartSpecies const &sp       = domain.part_species[s];
        auto const [v1span, v1divs] = Input::v1hist_specs.at(s);
        auto const [v2span, v2divs] = Input::v2hist_specs.at(s);
        Indexer const idxer{ v1span, v1divs, v2span, v2divs };
        if (!idxer)
            continue;

        histogram(sp, idxer);
    }
}

auto VHistogramRecorder::histogram(PartSpecies const &sp, Indexer const &idxer) const
    -> global_vhist_t
{
    // local counting
    //
    local_vhist_t local_vhist{};
    local_vhist.try_emplace(idxer.npos); // pre-allocate a slot for particles at out-of-range velocity
    for (Particle const &ptl : sp.bucket) {
        auto const sh = Shape<1>{ ptl.pos_x };
        auto       V  = sp.moment<1>().interp(sh);
        if (auto const n = Real{ sp.moment<0>().interp(sh) }; n < 1e-15) {
            V *= 0;
        } else {
            V /= n;
        }
        auto const &vel = sp.params.geomtr.cart2fac(ptl.vel() - V);
        auto const &key = idxer(vel.x, std::sqrt(vel.y * vel.y + vel.z * vel.z));
        local_vhist[key] += std::make_pair(1L, ptl.psd.weight);
    }

    // global counting
    //
    auto tk1 = comm.ibsend<unsigned long>(sp.bucket.size(), { master, tag });
    auto tk2 = comm.ibsend<4>({ local_vhist.begin(), local_vhist.end() }, { master, tag });

    global_vhist_t vhist{}; // one-based index
    if (is_master()) {
        // consolidation
        //
        Real denom{};
        for (int rank = 0, size = comm.size(); rank < size; ++rank) {
            auto const count  = *comm.recv<unsigned long>({ rank, tag });
            auto const lwhist = *comm.recv<4>({}, { rank, tag });

            denom += count;
            std::for_each(std::next(lwhist.rbegin()), lwhist.rend(), [&vhist](auto const &kv) {
                std::pair<long, long> const &key = kv.first;
                std::pair<long, Real> const &val = kv.second;
                vhist[key + 1] += val;
            });
        }

        // normalization
        //
        for (auto &kv : vhist) {
            std::pair<Real, Real> &val = kv.second;
            val /= denom;
        }
    }

    std::move(tk1).wait();
    std::move(tk2).wait();

    return vhist;
}
PIC1D_END_NAMESPACE