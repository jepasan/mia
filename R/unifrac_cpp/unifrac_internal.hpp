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

#include "biom_interface.hpp"
#include "unifrac.hpp"

namespace su {

 template<class TFloat>
 class PropStack {
   private:
     std::stack<std::vector<TFloat>> prop_stack;
     std::unordered_map<uint32_t, std::vector<TFloat>> prop_map;
     uint32_t defaultsize;
   public:
     PropStack(uint32_t vecsize);
     virtual ~PropStack();
     void clear(uint32_t i);
     void update(uint32_t i, std::vector<TFloat> vec);
     std::vector<TFloat> get(uint32_t i);
 };

 template<class TFloat>
 std::vector<TFloat> set_proportions(const BPTree &tree, uint32_t node,
                      const tse_interface &table,
                      PropStack<TFloat> &ps,
                      bool normalize = true);

}

#endif
