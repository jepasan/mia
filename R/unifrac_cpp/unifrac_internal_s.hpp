/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef __UNIFRAC_INTERNAL
#define __UNIFRAC_INTERNAL 1

#include <vector>
#include <stack>
#include <unordered_map>
#include "biom_interface_s.hpp"
#include "unifrac_s.hpp"

namespace su {

 template<class TFloat>
 class PropStack {
   private:
     std::stack<TFloat*> prop_stack;
     std::unordered_map<uint32_t, TFloat*> prop_map;
     uint32_t defaultsize;
   public:
     PropStack(uint32_t vecsize);
     virtual ~PropStack();
     TFloat* pop(uint32_t i);
     void push(uint32_t i);
     TFloat* get(uint32_t i);
 };

 template<class TFloat>
 void set_proportions(TFloat* __restrict__ props,
                      const BPTree &tree, uint32_t node,
                      const biom_interface &table,
                      PropStack<TFloat> &ps,
                      bool normalize = true);

}

#endif
