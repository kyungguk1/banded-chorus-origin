//
//  MasterDelegate.cc
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 6/15/19.
//  Copyright © 2019 kyungguk.com. All rights reserved.
//

#include "MasterDelegate.h"
#include "../Module/BField.h"
#include "../Module/EField.h"
#include "../Module/Charge.h"
#include "../Module/Current.h"
#include "../Module/Species.h"

#include <memory>

H1D::MasterDelegate::~MasterDelegate()
{
}
H1D::MasterDelegate::MasterDelegate(std::unique_ptr<Delegate> delegate) noexcept
: delegate{std::move(delegate)}
{
    for (WorkerDelegate &worker : workers) {
        worker.master = this;
    }
}

#if defined(HYBRID1D_MULTI_THREAD_FUNNEL_BOUNDARY_PASS) && HYBRID1D_MULTI_THREAD_FUNNEL_BOUNDARY_PASS
void H1D::MasterDelegate::pass(Domain const& domain, Species &sp)
{
    if ( (true) ) {
        delegate->pass(domain, sp);
        for (WorkerDelegate &worker : workers) {
            worker.mutable_comm.send(*this, &sp.bucket)();
            delegate->pass(domain, sp);
        }
    } else {
        // 1. consolidate all particles
        //
        decltype(sp.bucket) bucket{}; // hold workers' original particles
        for (WorkerDelegate &worker : workers) {
            worker.mutable_comm.send(*this, &bucket)();
        }
        sp.bucket.insert(sp.bucket.begin(), bucket.begin(), bucket.end());

        // 2. boundary pass
        //
        delegate->pass(domain, sp);

        // 3. distribute to workers
        //
        if (sp.bucket.size() < bucket.size()) {
            throw std::domain_error{__PRETTY_FUNCTION__};
        }
        for (WorkerDelegate &worker : workers) {
            tickets.push_back(worker.constant_comm.send(*this, &bucket));
        }
        tickets.clear();
    }
}
void H1D::MasterDelegate::pass(Domain const& domain, BField &bfield)
{
    delegate->pass(domain, bfield);
    broadcast_to_workers(bfield);
}
void H1D::MasterDelegate::pass(Domain const& domain, EField &efield)
{
    delegate->pass(domain, efield);
    broadcast_to_workers(efield);
}
void H1D::MasterDelegate::pass(Domain const& domain, Charge &charge)
{
    delegate->pass(domain, charge);
    broadcast_to_workers(charge);
}
void H1D::MasterDelegate::pass(Domain const& domain, Current &current)
{
    delegate->pass(domain, current);
    broadcast_to_workers(current);
}
#endif
void H1D::MasterDelegate::gather(Domain const& domain, Charge &charge)
{
    collect_from_workers(charge);
    delegate->gather(domain, charge);
    broadcast_to_workers(charge);
}
void H1D::MasterDelegate::gather(Domain const& domain, Current &current)
{
    collect_from_workers(current);
    delegate->gather(domain, current);
    broadcast_to_workers(current);
}
void H1D::MasterDelegate::gather(Domain const& domain, Species &sp)
{
    {
        collect_from_workers(sp.moment<0>());
        collect_from_workers(sp.moment<1>());
        collect_from_workers(sp.moment<2>());
    }
    delegate->gather(domain, sp);
    {
        broadcast_to_workers(sp.moment<0>());
        broadcast_to_workers(sp.moment<1>());
        broadcast_to_workers(sp.moment<2>());
    }
}

template <class T>
void H1D::MasterDelegate::broadcast_to_workers(GridQ<T> const &payload)
{
    for (WorkerDelegate &worker : workers) {
        tickets.push_back(worker.constant_comm.send(*this, &payload));
    }
    tickets.clear();
}
template <class T>
void H1D::MasterDelegate::collect_from_workers(GridQ<T> &buffer)
{
    // the first worker will collect all workers'
    //
    if (auto first = workers.begin(); first != workers.end()) {
        first->mutable_comm.send(*this, &buffer)();
    }
}