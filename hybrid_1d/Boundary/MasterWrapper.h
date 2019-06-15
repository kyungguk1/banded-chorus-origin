//
//  MasterWrapper.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 6/15/19.
//  Copyright © 2019 kyungguk.com. All rights reserved.
//

#ifndef MasterWrapper_h
#define MasterWrapper_h

#include "WorkerWrapper.h"
#include "../InputWrapper.h"

#include <array>
#include <vector>

HYBRID1D_BEGIN_NAMESPACE
class MasterWrapper : public Delegate {
    std::vector<InterThreadComm::Request> requests{};
public:
    std::array<WorkerWrapper, Input::n_workers> workers{};
    std::unique_ptr<Delegate> const delegate; // serial version

    ~MasterWrapper();
    MasterWrapper(std::unique_ptr<Delegate> delegate) noexcept;

private:
    void pass(Domain const&, Species &) override;
    void pass(Domain const&, BField &) override;
    void pass(Domain const&, EField &) override;
    void pass(Domain const&, Charge &) override;
    void pass(Domain const&, Current &) override;
    void gather(Domain const&, Charge &) override;
    void gather(Domain const&, Current &) override;
    void gather(Domain const&, Species &) override;
};
HYBRID1D_END_NAMESPACE

#endif /* MasterWrapper_h */