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

#ifndef GridQ_h
#define GridQ_h

#include "../Macros.h"
#include "../Predefined.h"
#include "./Shape.h"

#include <algorithm>
#include <array>
#include <memory>
#include <sstream>

PIC1D_BEGIN_NAMESPACE
/// 1D array with paddings on both ends that act as ghost cells
///
template <class T, long N> class GridQ {
public:
    constexpr static long size() noexcept { return N; }
    constexpr static long max_size() noexcept { return size() + 2 * Pad; }

    GridQ(GridQ const &)     = delete;
    GridQ(GridQ &&) noexcept = default;
    GridQ &operator=(GridQ &&) noexcept = default;

private:
    static_assert(size() > 0, "at least one element");
    using Backend = std::array<T, max_size()>;
    std::unique_ptr<Backend> ptr;

public:
    explicit GridQ() : ptr{std::make_unique<Backend>()} {}
    GridQ &operator=(GridQ const &other) noexcept
    {
        // other.ptr is expected to point to a valid object
        *this->ptr = *other.ptr;
        return *this;
    }

    // iterators
    //
    using iterator       = T *;
    using const_iterator = iterator const;

    [[nodiscard]] T const *begin() const noexcept { return ptr->data() + Pad; }
    [[nodiscard]] T *      begin() noexcept { return ptr->data() + Pad; }
    [[nodiscard]] T const *end() const noexcept { return begin() + size(); }
    [[nodiscard]] T *      end() noexcept { return begin() + size(); }

    [[nodiscard]] T const *dead_begin() const noexcept { return begin() - Pad; }
    [[nodiscard]] T *      dead_begin() noexcept { return begin() - Pad; }
    [[nodiscard]] T const *dead_end() const noexcept { return end() + Pad; }
    [[nodiscard]] T *      dead_end() noexcept { return end() + Pad; }

    // subscripts; index relative to the first non-padding element (i.e., relative to *begin())
    //
    [[nodiscard]] T const &operator[](long const i) const noexcept { return *(begin() + i); }
    [[nodiscard]] T &      operator[](long const i) noexcept { return *(begin() + i); }

    /// content filling (including paddings)
    ///
    void fill(T const &v) noexcept { std::fill(dead_begin(), dead_end(), v); }

    /// grid interpolator
    ///
    template <long Order> [[nodiscard]] T interp(Shape<Order> const &sx) const noexcept
    {
        T y{};
        for (long j = 0; j <= Order; ++j) {
            y += (*this)[sx.i[j]] * sx.w[j];
        }
        return y;
    }

    /// particle deposit; in-place operation
    ///
    template <long Order, class U> void deposit(Shape<Order> const &sx, U const &weight) noexcept
    {
        for (long j = 0; j <= Order; ++j) {
            (*this)[sx.i[j]] += weight * sx.w[j];
        }
    }

protected:
    /// content swap
    ///
    void swap(GridQ &o) noexcept { ptr.swap(o.ptr); }

    /// 3-point smoothing
    ///
    friend void _smooth(GridQ &filtered, GridQ const &source) noexcept
    {
        for (long i = 0; i < size(); ++i) {
            filtered[i] = (source[i - 1] + 2 * source[i] + source[i + 1]) * .25;
        }
    }

    // pretty print (buffered)
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, GridQ const &g)
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
PIC1D_END_NAMESPACE

#endif /* GridQ_h */
