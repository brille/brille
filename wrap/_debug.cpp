#include <pybind11/pybind11.h>
#include "debug.hpp"

void wrap_debug(pybind11::module &m){
    using namespace pybind11::literals;

    m.def("emit", [&](){return !_debug_printer.silenced();});
    m.def("emit", [&](const bool emt){
      bool orig = !_debug_printer.silenced();
      info_update_if( emt && orig, "Keeping debug printing on");
      info_update_if(!emt && orig, "Turning off debug printing");
      _debug_printer.silenced(!emt);
      info_update_if( emt && !orig, "Debug printing turned on");
      info_update_if(!emt && !orig, "Keeping debug printing off"); // never printed
      return !_debug_printer.silenced();
    }, "DebugPrinter emit status"_a);

    m.def("emit_datetime", [&](){return _debug_printer.datetime();});
    m.def("emit_datetime", [&](const bool emt){
      bool orig = _debug_printer.datetime();
      info_update_if( emt && orig, "Keeping datetime printing on");
      info_update_if(!emt && orig, "Turning off datetime printing");
      _debug_printer.datetime(emt);
      info_update_if( emt && !orig, "Datetime printing turned on");
      info_update_if(!emt && !orig, "Keeping datetime printing off");
      return _debug_printer.datetime();
    }, "DebugPrinter emit datetime status"_a);
}
