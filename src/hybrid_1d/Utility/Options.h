//
// Copyright (c) 2020, Kyungguk Min
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

#ifndef Options_h
#define Options_h

#include "../Macros.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <map>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

HYBRID1D_BEGIN_NAMESPACE
/// option parser from command-line arguments
///
/// Parsed are short-style options in the form '-opt_name', which are interpreted as boolean true
/// (represented by string literal 'true'), and long-style options in the form '--opt_name=value' or
/// '--opt_name value', which are interpreted as key-value pairs. Values in the second form of the
/// long-style options must not preceed with '--'.
///
/// The short-style option '-opt_name' is equivalent to '--opt_name=true', and the string literals
/// 'true' and 'false' are interpreted as boolean, but do not cast to integers '1' and '0'. Nor can
/// integers cast to booleans.
///
/// Any number of leading/trailing, but not interspersed, whitespaces in 'opt_name' and 'value' are
/// removed before parsing. An empty string, after removing the whitespaces, as opt_name and/or
/// value is ill-formed.
///
class [[nodiscard]] Options {
public:
    enum Style : long {
        short_ = 1, //!< tag for short-style option
        long_  = 2  //!< tag for long-style option
    };

    // option value parser
    //
    struct Value {
        friend Options;
        std::string s;
        Style       style{long_};

    public: // cast operators
        explicit operator std::string const &() const noexcept { return s; }
        explicit operator char const *() const noexcept { return s.c_str(); }
        explicit operator int() const { return std::stoi(s); }
        explicit operator long() const { return std::stol(s); }
        explicit operator unsigned long() const { return std::stoul(s); }
        explicit operator float() const { return std::stof(s); }
        explicit operator double() const { return std::stod(s); }
        explicit operator bool() const;

    public:
        template <class T> [[nodiscard]] auto cast() const
        {
            return static_cast<std::decay_t<T>>(*this);
        };
        template <class T> void operator()(T *p) const { *p = this->template cast<T>(); };
    };

private:
    std::map<std::string, Value> opts;

public:
    [[nodiscard]] std::map<std::string, Value> const *operator->() const &noexcept { return &opts; }
    [[nodiscard]] std::map<std::string, Value> const &operator*() const &noexcept { return opts; }

    Options() noexcept = default;
    Options(std::vector<std::string> args) { parse(std::move(args)); }

    /// parses options in the argument list and returns unparsed, order-preserved, arguments
    ///
    /// multiple calls will override/append to the options already parsed previously
    ///
    std::vector<std::string> parse(std::vector<std::string> args);

private:
    [[nodiscard]] static std::vector<std::string>
    transform_long_style(std::vector<std::string> args);
    [[nodiscard]] static std::vector<std::string>
    parse_short_options(std::vector<std::string> args, std::map<std::string, Value> &opts);
    [[nodiscard]] static std::vector<std::string>
    parse_long_options(std::vector<std::string> args, std::map<std::string, Value> &opts);

    // pretty print
    //
    template <class CharT, class Traits>
    friend decltype(auto) operator<<(std::basic_ostream<CharT, Traits> &os, Options const &opts)
    {
        auto const printer = [](decltype(os) os, auto const &kv) -> decltype(auto) {
            auto const &[key, val] = kv;
            return os << key << " : " << val.s;
        };
        os << '{';
        if (!opts->empty()) {
            printer(os, *begin(*opts));
            std::for_each(std::next(begin(*opts)), end(*opts), [&os, printer](auto const &kv) {
                printer(os << ", ", kv);
            });
        }
        return os << '}';
    }
};

// not for public use
//
void test_option_parser();
HYBRID1D_END_NAMESPACE

#endif /* Options_h */
