#pragma once

#include "compas.h"

class Stock2;

// Engagement queries share the _stock_2 nanobind module (NB_STATIC forbids
// cross-module type sharing), so registration is a hook the module macro in
// stock_2.cpp calls. Bindings land in a later task; this task registers none.
void register_engagement(nanobind::module_& m);
