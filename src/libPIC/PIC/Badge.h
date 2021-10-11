/*
 * Copyright (c) 2021, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#pragma once

#include <PIC/Config.h>

LIBPIC_BEGIN_NAMESPACE
template <class Holder>
class Badge {
    friend Holder;

    constexpr Badge() noexcept {}

public:
    Badge(Badge const &) noexcept = delete;
    Badge &operator=(Badge const &) noexcept = delete;
};
LIBPIC_END_NAMESPACE
