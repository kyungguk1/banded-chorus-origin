//
//  MessageDispatch.h
//  hybrid_1d
//
//  Created by KYUNGGUK MIN on 1/3/20.
//  Copyright © 2020 Kyungguk Min & Kaijun Liu. All rights reserved.
//

#ifndef MessageDispatch_h
#define MessageDispatch_h

#include "./MessageDispatch-tuple.h"
#include "./MessageDispatch-variant.h"

HYBRID1D_BEGIN_NAMESPACE
// not for public use
//
void test_message_queue();
void test_inter_thread_comm();
HYBRID1D_END_NAMESPACE

#endif /* MessageDispatch_h */
