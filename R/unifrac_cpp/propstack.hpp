/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef __FAITH_PROPSTACK
#define __FAITH_PROPSTACK 1

#include <vector>
#include <stack>
#include <unordered_map>

#include "tse.hpp"

namespace su {

 class PropStack {
   private:
     std::stack<std::vector<double>> prop_stack;
     std::unordered_map<uint32_t, std::vector<double>> prop_map;
     uint32_t defaultsize;
   public:
     PropStack(uint32_t vecsize);
     virtual ~PropStack();
     void clear(uint32_t i);
     void update(uint32_t i, std::vector<double> vec);
     std::vector<double> get(uint32_t i);
 };

 std::vector<double> set_proportions(const BPTree &tree, uint32_t node,
                      const tse &table,
                      PropStack &ps,
                      bool normalize = true);
}

#endif /* __FAITH_PROPSTACK */
