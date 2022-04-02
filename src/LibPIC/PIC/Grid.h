/*
 * Copyright (c) 2019-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>
#include <PIC/Shape.h>

#include <algorithm>
#include <array>
#include <memory>
#include <sstream>

LIBPIC_BEGIN_NAMESPACE
/// 1D grid-point array with paddings on both ends that act as ghost cells
///
template <class T, long N, long Pad>
class Grid {
public:
    constexpr static long size() noexcept { return N; }
    constexpr static long max_size() noexcept { return size() + 2 * Pad; }
    constexpr static long pad_size() noexcept { return Pad; }

    Grid(Grid const &)     = delete;
    Grid(Grid &&) noexcept = default;
    Grid &operator=(Grid &&) noexcept = default;

private:
    static_assert(size() > 0, "at least one element");
    using Backend = std::array<T, static_cast<unsigned long>(max_size())>;
    std::unique_ptr<Backend> ptr;

public:
    Grid()
    : ptr{ std::make_unique<Backend>() }
    {
    }

    Grid &operator=(Grid const &other) noexcept
    {
        // other.ptr is expected to point to a valid object
        *this->ptr = *other.ptr;
        return *this;
    }

    // iterators
    //
    using iterator       = T *;
    using const_iterator = T const *;

    [[nodiscard]] T const *dead_begin() const noexcept { return ptr->data(); }
    [[nodiscard]] T       *dead_begin() noexcept { return ptr->data(); }
    [[nodiscard]] T const *dead_end() const noexcept { return ptr->data() + max_size(); }
    [[nodiscard]] T       *dead_end() noexcept { return ptr->data() + max_size(); }

    [[nodiscard]] T const *begin() const noexcept { return dead_begin() + pad_size(); }
    [[nodiscard]] T       *begin() noexcept { return dead_begin() + pad_size(); }
    [[nodiscard]] T const *end() const noexcept { return begin() + size(); }
    [[nodiscard]] T       *end() noexcept { return begin() + size(); }

    // subscripts; index relative to the first non-padding element (i.e., relative to *begin())
    //
    [[nodiscard]] T const &operator[](long const i) const noexcept { return *(begin() + i); }
    [[nodiscard]] T       &operator[](long const i) noexcept { return *(begin() + i); }

    /// content filling (including paddings)
    ///
    void fill(T const &v) noexcept { std::fill(dead_begin(), dead_end(), v); }

    /// content swap
    ///
    void swap(Grid &o) noexcept { ptr.swap(o.ptr); }

    /// grid interpolator
    ///
    template <long Order>
    [[nodiscard]] T interp(Shape<Order> const &sx) const noexcept
    {
        static_assert(pad_size() >= Order,
                      "padding should be greater than or equal to the shape order");
        T y{};
        for (unsigned j = 0; j <= Order; ++j) {
            y += (*this)[sx.i(j)] * sx.w(j);
        }
        return y;
    }

    /// particle deposit; in-place operation
    ///
    template <long Order, class U>
    void deposit(Shape<Order> const &sx, U const &weight) noexcept
    {
        static_assert(pad_size() >= Order,
                      "padding should be greater than or equal to the shape order");
        for (unsigned j = 0; j <= Order; ++j) {
            (*this)[sx.i(j)] += weight * sx.w(j);
        }
    }

    /// 3-point smoothing
    ///
    static decltype(auto) smooth(Grid &filtered, Grid const &source) noexcept
    {
        static_assert(pad_size() >= 1, "not enough padding");
        for (long i = 0; i < size(); ++i) {
            filtered[i] = (source[i - 1] + 2 * source[i] + source[i + 1]) * .25;
        }
        return filtered;
    }

    // pretty print (buffered)
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Grid const &g)
    {
        std::basic_ostringstream<CharT, Traits> ss;
        {
            ss.flags(os.flags());
            ss.imbue(os.getloc());
            ss.precision(os.precision());
            //
            const_iterator it = g.begin(), end = g.end();
            ss << '{' << *it++; // guarrenteed to be at least one element
            while (it != end) {
                ss << ", " << *it++;
            }
            ss << '}';
        }
        return os << ss.str();
    }
};
LIBPIC_END_NAMESPACE
