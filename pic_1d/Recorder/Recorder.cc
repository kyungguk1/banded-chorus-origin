//
//  Recorder.cc
//  pic_1d
//
//  Created by KYUNGGUK MIN on 1/28/19.
//  Copyright © 2019 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#include "Recorder.h"
#include "../InputWrapper.h"

#include <limits>
#include <cmath>

namespace {
    constexpr long large_int = std::numeric_limits<unsigned>::max();
}

P1D::Recorder::Recorder(unsigned const recording_frequency) noexcept
: recording_frequency{recording_frequency ? recording_frequency*Input::inner_Nt : large_int} {
}

P1D::Vector const P1D::Recorder::e3 = {0, 0, 1};
P1D::Vector const P1D::Recorder::e1 = []{
    constexpr Real theta = Input::theta*M_PI/180;
    return Vector{std::cos(theta), std::sin(theta), 0};
}();
P1D::Vector const P1D::Recorder::e2 = []{
    constexpr Real theta = Input::theta*M_PI/180;
    return Vector{-std::sin(theta), std::cos(theta), 0};
}();